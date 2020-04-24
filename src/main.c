#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/ioctl.h>

#include "htslib/ketopt.h" // command-line argument parser
#include "htslib/kseq.h" // FASTA/Q parser
#include "htslib/khashl.h" // hash table
#include "htslib/cgranges.h" // read bed
#include "htslib/faidx.h"
#include "zlib.h"

KSEQ_INIT(gzFile, gzread)
KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

#define max(a, b) (((a) < (b)) ? (b) : (a))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

const uint8_t kmertochar[5] = { 'A', 'C', 'G', 'T', 'N' };

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char sig3[][4] = {
	"ACA", "ACC", "ACG", "ACT",
	"CCA", "CCC", "CCG", "CCT",
	"GCA", "GCC", "GCG", "GCT",
	"TCA", "TCC", "TCG", "TCT",
	"ATA", "ATC", "ATG", "ATT",
	"CTA", "CTC", "CTG", "CTT",
	"GTA", "GTC", "GTG", "GTT",
	"TTA", "TTC", "TTG", "TTT"
};

char key3[][4] = {
	"ACA", "ACC", "ACG", "ACT",
	"CCA", "CCC", "CCG", "AGG",
	"GCA", "GCC", "CGC", "AGC",
	"TCA", "GGA", "CGA", "AGA",
	"ATA", "ATC", "ATG", "AAT",
	"CTA", "CTC", "CAG", "AAG",
	"GTA", "GAC", "CAC", "AAC",
	"TAA", "GAA", "CAA", "AAA"
};

/*
uint64_t s64u(const char* s)
{
	size_t k = 0;
	uint64_t x_i, x = 0;
	while (*s && k++ < 4 * sizeof(uint64_t)) {
		x_i = seq_nt4_table[(uint8_t) *s];
		if (x_i > 3) break;
		x = (x << 2) | x_i;
	}
	return x;
}
*/
uint64_t s64u(const char* s)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<6) - 1, shift = 4;
	for (i = 0, x[0] = x[1] = 0; i < 3; ++i)
	{
		int c = seq_nt4_table[(uint8_t)s[i]];
		x[0] = (x[0] << 2 | c) & mask;                  // forward strand
		x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
	}
	return x[0] < x[1]? x[0] : x[1];
}


static void count_seq(kc_c1_t *h, int k, int len, char *seq) // insert k-mers in $seq to hash table $h
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				khint_t itr;
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				itr = kc_c1_put(h, y, &absent); // only add one strand!
				if (absent) kh_val(h, itr) = 0;
				++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static kc_c1_t *count_file(const char *fn, int k)
{
	gzFile fp;
	kseq_t *ks;
	kc_c1_t *h;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = kseq_init(fp);
	h = kc_c1_init();
	while (kseq_read(ks) >= 0)
		count_seq(h, k, ks->seq.l, ks->seq.s);
	kseq_destroy(ks);
	gzclose(fp);
	return h;
}

void u64s(uint64_t x, char* s)
{
	size_t i = 0;
	for (i = 0; i < 3; ++i)
	{
		s[2 - i] = "ACGT"[(uint8_t)x & 0x3];
		x >>= 2;
	}
	s[i] = '\0';
}

void u64rc(uint64_t x, char* s)
{
	size_t i = 0;
	for (i = 2; i > 0; --i)
	{
		s[2 - i] = "ACGT"[3 - (uint8_t)x & 0x3];
		x >>= 2;
	}
	s[2] = "ACGT"[3 - (uint8_t)x & 0x3];
	s[3] = '\0';
}

static void get_count(const kc_c1_t *h, const char *s, int *c)
{
	khint_t k;
	for (k = 0; k < kh_end(h); ++k)
	{
		if (kh_exist(h, k))
		{
			if (s64u(s) == kh_key(h, k))
			{
				*c = kh_val(h, k);
				break;
			}
		}
	}
}

/*
 * https://github.com/parklab/SigMA/blob/master/R/get_trinuc_norm.R
 *
 * The get_trinuc_norm() can be used to determine the
 * normalization to be used to match the frequency in the
 * sequencing platform to the frequency in the whole genomes
 *
 * @param bed_file the path to the bed file that defines the
 * library of the sequencing platform
 * @param do_MC_sampling if set to TRUE, to speed up
 * the calculation the regions are sampled randomly rather
 * than providing a full count
 */
static void get_trinuc_norm(const kc_c1_t *h, const kc_c1_t *h_wg, FILE *fp)
{
	int i, j;
	/*
	 * counts_expanded <- c(rep(counts[1:16], 3), rep(counts[17:32], 3))
	 * norm <- counts_expanded/counts_trinuc_genome
	 * return(norm=100*norm/sum(norm))
	 */
	double norm[32] = {0}, sum = 0;
	int c[32] = {0}, c_wg[32] = {0};
	for (i = 0; i < 32; ++i)
	{
		get_count(h, key3[i], c + i);
		get_count(h_wg, key3[i], c_wg + i);
		norm[i] = (double)c[i] / (c_wg[i] ? c_wg[i] : 1);
		sum += norm[i] * 3;
	}
	sum = sum ? sum : 1;
	for (j = 0; j < 3; ++j)
		for (i = 0; i < 16; ++i)
			fprintf(fp, "%s\t%f\n", sig3[i], norm[i] / sum * 100);
	for (j = 0; j < 3; ++j)
		for (i = 16; i < 32; ++i)
			fprintf(fp, "%s\t%f\n", sig3[i], norm[i] / sum * 100);
}

static void print_sig3(const kc_c1_t *h, FILE *fp)
{
	int i;
	for (i = 0; i < 32; ++i)
	{
		int c = 0;
		get_count(h, key3[i], &c);
		fprintf(fp, "%s\t%d\n", sig3[i], c);
	}
}

// bed functions
typedef struct {
	int32_t l;
	char *s;
} bed_rest1_t;

typedef struct {
	int64_t n, m;
	bed_rest1_t *a;
} bed_rest_t;

static char *parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	if (r) *r = 0;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) {
				en = atol(p);
				if (r && c != 0) *r = q, *q = c;
			}
			++i, p = q + 1;
			if (i == 3 || c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

char *parse_bed3(char *s, int32_t *st_, int32_t *en_)
{
	return parse_bed3b(s, st_, en_, 0);
}

static cgranges_t *read_bed(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	cgranges_t *cr;
	kstring_t str = {0,0,0};
	int64_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = ks_init(fp);
	cr = cr_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
	{
		char *ctg;
		int32_t st, en;
		ctg = parse_bed3(str.s, &st, &en);
		if (ctg) cr_add(cr, ctg, st, en, k++);
	}
	if (k > INT32_MAX)
		fprintf(stderr, "WARNING: more than %d records; some functionality may not work!\n", INT32_MAX);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return cr;
}

void context(const char *ref, const char *bed, const char *out)
{
	size_t i;
	char *chr;
	int beg, end, ret;
	cgranges_t *cr = read_bed(bed); assert(cr);
	if (!cr_is_sorted(cr)) cr_sort(cr);
	//cr_merge_pre_index(cr);
	FILE *fp = fopen(out, "w");
	//fputs("#kmer\trevcom\tcount\n", fp);
	faidx_t *fai = fai_load(ref); assert(fai);
	kc_c1_t *h = kc_c1_init();
	// iterate over the merged regions
	for (i = 0; i < cr->n_r; ++i)
	{
		chr = cr->ctg[cr->r[i].x >> 32].name;
		const int len = faidx_seq_len(fai, chr);
		beg = max(0, (int32_t)cr->r[i].x - 2);
		end = min(len, cr->r[i].y + 1);
		char *seq = faidx_fetch_seq(fai, chr, beg, end - 1, &ret);
		count_seq(h, 3, end - beg, seq);
		free(seq);
	}
	kc_c1_t *h_wg = kc_c1_init();
	int n = faidx_nseq(fai);
	for (i = 0; i < n; ++i)
	{
		const char *chr = faidx_iseq(fai, i);
		const int len = faidx_seq_len(fai, chr);
		char *seq = faidx_fetch_seq(fai, chr, 0, len - 1, &ret);
		count_seq(h_wg, 3, len, seq);
		free(seq);
	}
	fai_destroy(fai);
	cr_destroy(cr);
	// output
	get_trinuc_norm(h, h_wg, fp);
	kc_c1_destroy(h);
	kc_c1_destroy(h_wg);
	fclose(fp);
}

int main(int argc, char *argv[])
{
	kc_c1_t *h;
	int c;
	char *b = 0, *o = 0, *r = 0;
	ketopt_t opt = KETOPT_INIT;
	while ((c = ketopt(&opt, argc, argv, 1, "b:o:", 0)) >= 0)
	{
		if (c == 'b') b = opt.arg;
		if (c == 'o') o = opt.arg;
	}
	if (argc - opt.ind < 1)
	{
		fputc('\n', stderr);
		fprintf(stderr, "usage: \e[1;31m%s\e[0;0m [options] <in.fa>\n", basename(argv[0]));
		fputc('\n', stderr);
		fputs("  -b gene panel bed file\n", stderr);
		fputs("  -o normalization relative to genome\n", stderr);
		fputs("\nnotes: fai index is required for fa\n", stderr);
		fputs("\nreference: https://github.com/lh3/kmer-cnt\n", stderr);
		fputs("           https://github.com/lh3/cgranges\n", stderr);
		fputs("           https://github.com/parklab/SigMA\n", stderr);
		fputs("\ncontact liujc@geneplus.org.cn for bug reports\n", stderr);
		fputc('\n', stderr);
		return 1;
	}
	if (!b || (access(b, R_OK) == -1))
	{
		fputs("required region file no specified or can't be found\n", stderr);
		exit(1);
	}
	r = argv[opt.ind];
	if ((access(r, R_OK) == -1))
	{
		fputs("reference fa specified can't be accessed\n", stderr);
		exit(1);
	}
	context(r, b, o);
	return 0;
}
