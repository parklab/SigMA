/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013, 2014 Genome Research Ltd

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <sys/types.h>
#include <inttypes.h>

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"

#define BGZF_CACHE
#define BGZF_MT

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


/* BGZF/GZIP header (speciallized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

#ifdef BGZF_CACHE
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
#include "htslib/khash.h"
KHASH_MAP_INIT_INT64(cache, cache_t)
#endif

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

struct __bgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};

void bgzf_index_destroy(BGZF *fp);
int bgzf_index_add_block(BGZF *fp);

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static BGZF *bgzf_read_init(hFILE *hfpr)
{
    BGZF *fp;
    uint8_t magic[18];
    ssize_t n = hpeek(hfpr, magic, 18);
    if (n < 0) return NULL;

    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) return NULL;

    fp->is_write = 0;
    fp->is_compressed = (n==2 && magic[0]==0x1f && magic[1]==0x8b);
    fp->uncompressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->compressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
    fp->is_compressed = (n==18 && magic[0]==0x1f && magic[1]==0x8b) ? 1 : 0;
    fp->is_gzip = ( !fp->is_compressed || ((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0) ) ? 0 : 1;
#ifdef BGZF_CACHE
    fp->cache = kh_init(cache);
#endif
    return fp;
}

// get the compress level from the mode string: compress_level==-1 for the default level, -2 plain uncompressed
static int mode2level(const char *__restrict mode)
{
    int i, compress_level = -1;
    for (i = 0; mode[i]; ++i)
        if (mode[i] >= '0' && mode[i] <= '9') break;
    if (mode[i]) compress_level = (int)mode[i] - '0';
    if (strchr(mode, 'u')) compress_level = -2;
    return compress_level;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
    BGZF *fp = 0;
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hopen(path, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; }
        fp->fp = fpr;
    }
	else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_dopen(int fd, const char *mode)
{
    BGZF *fp = 0;
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hdopen(fd, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; } // FIXME this closes fd
        fp->fp = fpr;
    }
	else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_hopen(hFILE *hfp, const char *mode)
{
    BGZF *fp = NULL;
    if (strchr(mode, 'r')) {
        fp = bgzf_read_init(hfp);
        if (fp == NULL) return NULL;
    }
	else { errno = EINVAL; return 0; }

    fp->fp = hfp;
    fp->is_be = ed_is_big();
    return fp;
}

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = (Bytef*)fp->compressed_block + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = (Bytef*)fp->uncompressed_block;
    zs.avail_out = BGZF_MAX_BLOCK_SIZE;

    if (inflateInit2(&zs, -15) != Z_OK) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
        inflateEnd(&zs);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    if (inflateEnd(&zs) != Z_OK) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return zs.total_out;
}

static int inflate_gzip_block(BGZF *fp, int cached)
{
    int ret = Z_OK;
    do
    {
        if ( !cached && fp->gz_stream->avail_out!=0 )
        {
            fp->gz_stream->avail_in = hread(fp->fp, fp->compressed_block, BGZF_BLOCK_SIZE);
            if ( fp->gz_stream->avail_in<=0 ) return fp->gz_stream->avail_in;
            if ( fp->gz_stream->avail_in==0 ) break;
            fp->gz_stream->next_in = fp->compressed_block;
        }
        else cached = 0;
        do
        {
            fp->gz_stream->next_out = (Bytef*)fp->uncompressed_block + fp->block_offset;
            fp->gz_stream->avail_out = BGZF_MAX_BLOCK_SIZE - fp->block_offset;
            ret = inflate(fp->gz_stream, Z_NO_FLUSH);
            if ( ret==Z_BUF_ERROR ) continue;   // non-critical error
            if ( ret<0 ) return -1;
            unsigned int have = BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
            if ( have ) return have;
        }
        while ( fp->gz_stream->avail_out == 0 );
    }
    while (ret != Z_STREAM_END);
    return BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

#ifdef BGZF_CACHE
static void free_cache(BGZF *fp)
{
    khint_t k;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (fp->is_write) return;
    for (k = kh_begin(h); k < kh_end(h); ++k)
        if (kh_exist(h, k)) free(kh_val(h, k).block);
    kh_destroy(cache, h);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, BGZF_MAX_BLOCK_SIZE);
    if ( hseek(fp->fp, p->end_offset, SEEK_SET) < 0 )
    {
        // todo: move the error up
        fprintf(stderr,"Could not hseek to %"PRId64"\n", p->end_offset);
        exit(1);
    }
    return p->size;
}

static void cache_block(BGZF *fp, int size)
{
    int ret;
    khint_t k;
    cache_t *p;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (BGZF_MAX_BLOCK_SIZE >= fp->cache_size) return;
    if ((kh_size(h) + 1) * BGZF_MAX_BLOCK_SIZE > (uint32_t)fp->cache_size) {
        /* A better way would be to remove the oldest block in the
         * cache, but here we remove a random one for simplicity. This
         * should not have a big impact on performance. */
        for (k = kh_begin(h); k < kh_end(h); ++k)
            if (kh_exist(h, k)) break;
        if (k < kh_end(h)) {
            free(kh_val(h, k).block);
            kh_del(cache, h, k);
        }
    }
    k = kh_put(cache, h, fp->block_address, &ret);
    if (ret == 0) return; // if this happens, a bug!
    p = &kh_val(h, k);
    p->size = fp->block_length;
    p->end_offset = fp->block_address + size;
    p->block = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    memcpy(kh_val(h, k).block, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
}
#else
static void free_cache(BGZF *fp) {}
static int load_block_from_cache(BGZF *fp, int64_t block_address) {return 0;}
static void cache_block(BGZF *fp, int size) {}
#endif

int bgzf_read_block(BGZF *fp)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // Reading an uncompressed file
    if ( !fp->is_compressed )
    {
        count = hread(fp->fp, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
        if ( count==0 )
        {
            fp->block_length = 0;
            return 0;
        }
        if (fp->block_length != 0) fp->block_offset = 0;
        fp->block_address += count;
        fp->block_length = count;
        return 0;
    }

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);
    if ( fp->is_gzip && fp->gz_stream ) // is this is a initialized gzip stream?
    {
        count = inflate_gzip_block(fp, 0);
        if ( count<0 )
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        return 0;
    }
    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hread(fp->fp, header, sizeof(header));
    if (count == 0) { // no data read
        fp->block_length = 0;
        return 0;
    }
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        fp->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if ( ret==-1 )
    {
        // GZIP, not BGZF
        uint8_t *cblock = (uint8_t*)fp->compressed_block;
        memcpy(cblock, header, sizeof(header));
        count = hread(fp->fp, cblock+sizeof(header), BGZF_BLOCK_SIZE - sizeof(header)) + sizeof(header);
        int nskip = 10;

        // Check optional fields to skip: FLG.FNAME,FLG.FCOMMENT,FLG.FHCRC,FLG.FEXTRA
        // Note: Some of these fields are untested, I did not have appropriate data available
        if ( header[3] & 0x4 ) // FLG.FEXTRA
        {
            nskip += unpackInt16(&cblock[nskip]) + 2;
        }
        if ( header[3] & 0x8 ) // FLG.FNAME
        {
            while ( nskip<BGZF_BLOCK_SIZE && cblock[nskip] ) nskip++;
            if ( nskip==BGZF_BLOCK_SIZE )
            {
                fp->errcode |= BGZF_ERR_HEADER;
                return -1;
            }
            nskip++;
        }
        if ( header[3] & 0x10 ) // FLG.FCOMMENT
        {
            while ( nskip<BGZF_BLOCK_SIZE && cblock[nskip] ) nskip++;
            if ( nskip==BGZF_BLOCK_SIZE )
            {
                fp->errcode |= BGZF_ERR_HEADER;
                return -1;
            }
            nskip++;
        }
        if ( header[3] & 0x2 ) nskip += 2;  //  FLG.FHCRC

        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*) calloc(1,sizeof(z_stream));
        int ret = inflateInit2(fp->gz_stream, -15);
        if (ret != Z_OK)
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->gz_stream->avail_in = count - nskip;
        fp->gz_stream->next_in  = cblock + nskip;
        count = inflate_gzip_block(fp, 1);
        if ( count<0 )
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        if ( fp->idx_build_otf ) return -1; // cannot build index for gzip
        return 0;
    }
    size = count;
    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    compressed_block = (uint8_t*)fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    if ((count = inflate_block(fp, block_length)) < 0) return -1;
    if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
    fp->block_address = block_address;
    fp->block_length = count;
    if ( fp->idx_build_otf )
    {
        bgzf_index_add_block(fp);
        fp->idx->ublock_addr += count;
    }
    cache_block(fp, size);
    return 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
{
    ssize_t bytes_read = 0;
    uint8_t *output = (uint8_t*)data;
    if (length <= 0) return 0;
    assert(fp->is_write == 0);
    while (bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
        uint8_t *buffer;
        if (available <= 0) {
            if (bgzf_read_block(fp) != 0) return -1;
            available = fp->block_length - fp->block_offset;
            if (available <= 0) break;
        }
        copy_length = length - bytes_read < available? length - bytes_read : available;
        buffer = (uint8_t*)fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
        fp->block_address = htell(fp->fp);
        fp->block_offset = fp->block_length = 0;
    }
    fp->uncompressed_address += bytes_read;
    return bytes_read;
}

ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)
{
    return hread(fp->fp, data, length);
}

int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if ( fp->is_gzip )
    {
        (void)inflateEnd(fp->gz_stream);
        free(fp->gz_stream);
    }
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    free_cache(fp);
    free(fp);
    return 0;
}

void bgzf_set_cache_size(BGZF *fp, int cache_size)
{
    if (fp) fp->cache_size = cache_size;
}

int bgzf_check_EOF(BGZF *fp)
{
    uint8_t buf[28];
    off_t offset = htell(fp->fp);
    if (hseek(fp->fp, -28, SEEK_END) < 0) {
        if (errno == ESPIPE) { hclearerr(fp->fp); return 2; }
        else return -1;
    }
    if ( hread(fp->fp, buf, 28) != 28 ) return -1;
    if ( hseek(fp->fp, offset, SEEK_SET) < 0 ) return -1;
    return (memcmp("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", buf, 28) == 0)? 1 : 0;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
    int block_offset;
    int64_t block_address;

    if (fp->is_write || where != SEEK_SET) {
        fp->errcode |= BGZF_ERR_MISUSE;
        return -1;
    }
    block_offset = pos & 0xFFFF;
    block_address = pos >> 16;
    if (hseek(fp->fp, block_address, SEEK_SET) < 0) {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    fp->block_length = 0;  // indicates current block has not been loaded
    fp->block_address = block_address;
    fp->block_offset = block_offset;
    return 0;
}

int bgzf_is_bgzf(const char *fn)
{
    uint8_t buf[16];
    int n;
    hFILE *fp;
    if ((fp = hopen(fn, "r")) == 0) return 0;
    n = hread(fp, buf, 16);
    if ( hclose(fp) < 0 ) return -1;
    if (n != 16) return 0;
    return memcmp(g_magic, buf, 16) == 0? 1 : 0;
}

int bgzf_getc(BGZF *fp)
{
    int c;
    if (fp->block_offset >= fp->block_length) {
        if (bgzf_read_block(fp) != 0) return -2; /* error */
        if (fp->block_length == 0) return -1; /* end-of-file */
    }
    c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
        fp->block_address = htell(fp->fp);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    fp->uncompressed_address++;
    return c;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

int bgzf_getline(BGZF *fp, int delim, kstring_t *str)
{
    int l, state = 0;
    unsigned char *buf = (unsigned char*)fp->uncompressed_block;
    str->l = 0;
    do {
        if (fp->block_offset >= fp->block_length) {
            if (bgzf_read_block(fp) != 0) { state = -2; break; }
            if (fp->block_length == 0) { state = -1; break; }
        }
        for (l = fp->block_offset; l < fp->block_length && buf[l] != delim; ++l);
        if (l < fp->block_length) state = 1;
        l -= fp->block_offset;
        if (str->l + l + 1 >= str->m) {
            str->m = str->l + l + 2;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
        }
        memcpy(str->s + str->l, buf + fp->block_offset, l);
        str->l += l;
        fp->block_offset += l + 1;
        if (fp->block_offset >= fp->block_length) {
            fp->block_address = htell(fp->fp);
            fp->block_offset = 0;
            fp->block_length = 0;
        }
    } while (state == 0);
    if (str->l == 0 && state < 0) return state;
    fp->uncompressed_address += str->l;
    if ( delim=='\n' && str->l>0 && str->s[str->l-1]=='\r' ) str->l--;
    str->s[str->l] = 0;
    return str->l;
}

void bgzf_index_destroy(BGZF *fp)
{
    if ( !fp->idx ) return;
    free(fp->idx->offs);
    free(fp->idx);
    fp->idx = NULL;
    fp->idx_build_otf = 0;
}

int bgzf_index_build_init(BGZF *fp)
{
    bgzf_index_destroy(fp);
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if ( !fp->idx ) return -1;
    fp->idx_build_otf = 1;  // build index on the fly
    return 0;
}

int bgzf_index_add_block(BGZF *fp)
{
    fp->idx->noffs++;
    if ( fp->idx->noffs > fp->idx->moffs )
    {
        fp->idx->moffs = fp->idx->noffs;
        kroundup32(fp->idx->moffs);
        fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
        if ( !fp->idx->offs ) return -1;
    }
    fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->ublock_addr;
    fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->block_address;
    return 0;
}

int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)
{
    char *tmp = NULL;
    if ( suffix )
    {
        int blen = strlen(bname);
        int slen = strlen(suffix);
        tmp = (char*) malloc(blen + slen + 1);
        if ( !tmp ) return -1;
        memcpy(tmp,bname,blen);
        memcpy(tmp+blen,suffix,slen+1);
    }

    FILE *idx = fopen(tmp?tmp:bname,"rb");
    if ( tmp ) free(tmp);
    if ( !idx ) return -1;

    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    uint64_t x;
    if ( fread(&x, 1, sizeof(x), idx) != sizeof(x) ) return -1;

    fp->idx->noffs = fp->idx->moffs = 1 + (fp->is_be ? ed_swap_8(x) : x);
    fp->idx->offs  = (bgzidx1_t*) malloc(fp->idx->moffs*sizeof(bgzidx1_t));
    fp->idx->offs[0].caddr = fp->idx->offs[0].uaddr = 0;

    int i;
    if ( fp->is_be )
    {
        int ret = 0;
        for (i=1; i<fp->idx->noffs; i++)
        {
            ret += fread(&x, 1, sizeof(x), idx); fp->idx->offs[i].caddr = ed_swap_8(x);
            ret += fread(&x, 1, sizeof(x), idx); fp->idx->offs[i].uaddr = ed_swap_8(x);
        }
        if ( ret != sizeof(x)*2*(fp->idx->noffs-1) ) return -1;
    }
    else
    {
        int ret = 0;
        for (i=1; i<fp->idx->noffs; i++)
        {
            ret += fread(&x, 1, sizeof(x), idx); fp->idx->offs[i].caddr = x;
            ret += fread(&x, 1, sizeof(x), idx); fp->idx->offs[i].uaddr = x;
        }
        if ( ret != sizeof(x)*2*(fp->idx->noffs-1) ) return -1;
    }
    fclose(idx);
    return 0;

}

int bgzf_useek(BGZF *fp, long uoffset, int where)
{
    if ( !fp->is_compressed )
    {
        if (hseek(fp->fp, uoffset, SEEK_SET) < 0)
        {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = uoffset;
        fp->block_offset = 0;
        bgzf_read_block(fp);
        fp->uncompressed_address = uoffset;
        return 0;
    }

    if ( !fp->idx )
    {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }

    // binary search
    int ilo = 0, ihi = fp->idx->noffs - 1;
    while ( ilo<=ihi )
    {
        int i = (ilo+ihi)*0.5;
        if ( uoffset < fp->idx->offs[i].uaddr ) ihi = i - 1;
        else if ( uoffset >= fp->idx->offs[i].uaddr ) ilo = i + 1;
        else break;
    }
    int i = ilo-1;
    if (hseek(fp->fp, fp->idx->offs[i].caddr, SEEK_SET) < 0)
    {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    fp->block_length = 0;  // indicates current block has not been loaded
    fp->block_address = fp->idx->offs[i].caddr;
    fp->block_offset = 0;
    if ( bgzf_read_block(fp) < 0 ) return -1;
    if ( uoffset - fp->idx->offs[i].uaddr > 0 )
    {
        fp->block_offset = uoffset - fp->idx->offs[i].uaddr;
        assert( fp->block_offset <= fp->block_length );     // todo: skipped, unindexed, blocks
    }
    fp->uncompressed_address = uoffset;
    return 0;
}

long bgzf_utell(BGZF *fp)
{
    return fp->uncompressed_address;    // currently maintained only when reading
}

