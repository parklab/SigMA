/*  faidx.c -- FASTA random access.

    Copyright (C) 2008, 2009, 2013-2015 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/hfile.h"
#include "htslib/khash.h"

typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline void fai_insert_index(faidx_t *idx, const char *name, int len, int line_len, int line_blen, uint64_t offset)
{
    khint_t k;
    int ret;
    faidx1_t t;
    if (idx->n == idx->m) {
        idx->m = idx->m? idx->m<<1 : 16;
        idx->name = (char**)realloc(idx->name, sizeof(char*) * idx->m);
    }
    idx->name[idx->n] = strdup(name);
    k = kh_put(s, idx->hash, idx->name[idx->n], &ret);
    t.len = len; t.line_len = line_len; t.line_blen = line_blen; t.offset = offset;
    kh_value(idx->hash, k) = t;
    ++idx->n;
}

faidx_t *fai_read(FILE *fp)
{
    faidx_t *fai;
    char *buf, *p;
    int len, line_len, line_blen;
#ifdef _WIN32
    long offset;
#else
    long long offset;
#endif
    fai = (faidx_t*)calloc(1, sizeof(faidx_t));
    fai->hash = kh_init(s);
    buf = (char*)calloc(0x10000, 1);
    while (!feof(fp) && fgets(buf, 0x10000, fp)) {
        for (p = buf; *p && isgraph(*p); ++p);
        *p = 0; ++p;
#ifdef _WIN32
        sscanf(p, "%d%ld%d%d", &len, &offset, &line_blen, &line_len);
#else
        sscanf(p, "%d%lld%d%d", &len, &offset, &line_blen, &line_len);
#endif
        fai_insert_index(fai, buf, len, line_len, line_blen, offset);
    }
    free(buf);
    return fai;
}

void fai_destroy(faidx_t *fai)
{
    int i;
    for (i = 0; i < fai->n; ++i) free(fai->name[i]);
    free(fai->name);
    kh_destroy(s, fai->hash);
    if (fai->bgzf) bgzf_close(fai->bgzf);
    free(fai);
}

faidx_t *fai_load(const char *fn)
{
    char *str;
    FILE *fp;
    faidx_t *fai;
    str = (char*)calloc(strlen(fn) + 5, 1);
    sprintf(str, "%s.fai", fn);

    fp = fopen(str, "rb");

        if (fp == 0) {
            fprintf(stderr, "[fai_load] fail to open FASTA index.\n");
            free(str);
            return 0;
        }

    fai = fai_read(fp);
    fclose(fp);

    fai->bgzf = bgzf_open(fn, "rb");
    free(str);
    if (fai->bgzf == 0) {
        fprintf(stderr, "[fai_load] fail to open FASTA file.\n");
        return 0;
    }
    if ( fai->bgzf->is_compressed==1 )
    {
        if ( bgzf_index_load(fai->bgzf, fn, ".gzi") < 0 )
        {
            fprintf(stderr, "[fai_load] failed to load .gzi index: %s[.gzi]\n", fn);
            fai_destroy(fai);
            return NULL;
        }
    }
    return fai;
}

char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
    char *s;
    int c, i, l, k, name_end;
    khiter_t iter;
    faidx1_t val;
    khash_t(s) *h;
    int beg, end;

    beg = end = -1;
    h = fai->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    // remove space
    for (i = k = 0; i < l; ++i)
        if (!isspace(str[i])) s[k++] = str[i];
    s[k] = 0; l = k;
    // determine the sequence name
    for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
    if (i >= 0) name_end = i;
    if (name_end < l) { // check if this is really the end
        int n_hyphen = 0;
        for (i = name_end + 1; i < l; ++i) {
            if (s[i] == '-') ++n_hyphen;
            else if (!isdigit(s[i]) && s[i] != ',') break;
        }
        if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
        s[name_end] = 0;
        iter = kh_get(s, h, s);
        if (iter == kh_end(h)) { // cannot find the sequence name
            iter = kh_get(s, h, str); // try str as the name
            if (iter == kh_end(h)) {
                *len = 0;
            free(s); return 0;
            } else s[name_end] = ':', name_end = l;
        }
    } else iter = kh_get(s, h, str);
    if(iter == kh_end(h)) {
        fprintf(stderr, "[fai_fetch] Warning - Reference %s not found in FASTA file, returning empty sequence\n", str);
        free(s);
        *len = -2;
        return 0;
    };
    val = kh_value(h, iter);
    // parse the interval
    if (name_end < l) {
        for (i = k = name_end + 1; i < l; ++i)
            if (s[i] != ',') s[k++] = s[i];
        s[k] = 0;
        beg = atoi(s + name_end + 1);
        for (i = name_end + 1; i != k; ++i) if (s[i] == '-') break;
        end = i < k? atoi(s + i + 1) : val.len;
        if (beg > 0) --beg;
    } else beg = 0, end = val.len;
    if (beg >= val.len) beg = val.len;
    if (end >= val.len) end = val.len;
    if (beg > end) beg = end;
    free(s);

    // now retrieve the sequence
    int ret = bgzf_useek(fai->bgzf, val.offset + beg / val.line_blen * val.line_len + beg % val.line_blen, SEEK_SET);
    if ( ret<0 )
    {
        *len = -1;
        fprintf(stderr, "[fai_fetch] Error: fai_fetch failed. (Seeking in a compressed, .gzi unindexed, file?)\n");
        return NULL;
    }
    l = 0;
    s = (char*)malloc(end - beg + 2);
    while ( (c=bgzf_getc(fai->bgzf))>=0 && l < end - beg )
        if (isgraph(c)) s[l++] = c;
    s[l] = '\0';
    *len = l;
    return s;
}

int faidx_fetch_nseq(const faidx_t *fai)
{
    return fai->n;
}

int faidx_nseq(const faidx_t *fai)
{
    return fai->n;
}

const char *faidx_iseq(const faidx_t *fai, int i)
{
    return fai->name[i];
}

int faidx_seq_len(const faidx_t *fai, const char *seq)
{
    khint_t k = kh_get(s, fai->hash, seq);
    if ( k == kh_end(fai->hash) ) return -1;
    return kh_val(fai->hash, k).len;
}

char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    int l, c;
    khiter_t iter;
    faidx1_t val;
    char *seq=NULL;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);
    if (iter == kh_end(fai->hash))
    {
        *len = -2;
        fprintf(stderr, "[fai_fetch_seq] The sequence \"%s\" not found\n", c_name);
        return NULL;
    }
    val = kh_value(fai->hash, iter);
    if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence
    int ret = bgzf_useek(fai->bgzf, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
    if ( ret<0 )
    {
        *len = -1;
        fprintf(stderr, "[fai_fetch_seq] Error: fai_fetch failed. (Seeking in a compressed, .gzi unindexed, file?)\n");
        return NULL;
    }
    l = 0;
    seq = (char*)malloc(p_end_i - p_beg_i + 2);
    while ( (c=bgzf_getc(fai->bgzf))>=0 && l < p_end_i - p_beg_i + 1)
        if (isgraph(c)) seq[l++] = c;
    seq[l] = '\0';
    *len = l;
    return seq;
}

int faidx_has_seq(const faidx_t *fai, const char *seq)
{
    khiter_t iter = kh_get(s, fai->hash, seq);
    if (iter == kh_end(fai->hash)) return 0;
    return 1;
}

