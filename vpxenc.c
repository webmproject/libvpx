/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/* This is a simple program that encodes YV12 files and generates ivf
 * files using the new interface.
 */
#if defined(_WIN32) || !CONFIG_OS_SUPPORT
#define USE_POSIX_MMAP 0
#else
#define USE_POSIX_MMAP 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include "vpx/vpx_encoder.h"
#if USE_POSIX_MMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#endif
#include "vpx_version.h"
#include "vpx/vp8cx.h"
#include "vpx_ports/mem_ops.h"
#include "vpx_ports/vpx_timer.h"
#include "tools_common.h"
#include "y4minput.h"
#include "libmkv/EbmlWriter.h"
#include "libmkv/EbmlIDs.h"

/* Need special handling of these functions on Windows */
#if defined(_MSC_VER)
/* MSVS doesn't define off_t, and uses _f{seek,tell}i64 */
typedef __int64 off_t;
#define fseeko _fseeki64
#define ftello _ftelli64
#elif defined(_WIN32)
/* MinGW defines off_t, and uses f{seek,tell}o64 */
#define fseeko fseeko64
#define ftello ftello64
#endif

#if defined(_MSC_VER)
#define LITERALU64(n) n
#else
#define LITERALU64(n) n##LLU
#endif

/* We should use 32-bit file operations in WebM file format
 * when building ARM executable file (.axf) with RVCT */
#if !CONFIG_OS_SUPPORT
typedef long off_t;
#define fseeko fseek
#define ftello ftell
#endif

static const char *exec_name;

static const struct codec_item
{
    char const              *name;
    const vpx_codec_iface_t *iface;
    unsigned int             fourcc;
} codecs[] =
{
#if CONFIG_VP8_ENCODER
    {"vp8",  &vpx_codec_vp8_cx_algo, 0x30385056},
#endif
};

static void usage_exit();

void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    usage_exit();
}

static void ctx_exit_on_error(vpx_codec_ctx_t *ctx, const char *s)
{
    if (ctx->err)
    {
        const char *detail = vpx_codec_error_detail(ctx);

        fprintf(stderr, "%s: %s\n", s, vpx_codec_error(ctx));

        if (detail)
            fprintf(stderr, "    %s\n", detail);

        exit(EXIT_FAILURE);
    }
}

/* This structure is used to abstract the different ways of handling
 * first pass statistics.
 */
typedef struct
{
    vpx_fixed_buf_t buf;
    int             pass;
    FILE           *file;
    char           *buf_ptr;
    size_t          buf_alloc_sz;
} stats_io_t;

int stats_open_file(stats_io_t *stats, const char *fpf, int pass)
{
    int res;

    stats->pass = pass;

    if (pass == 0)
    {
        stats->file = fopen(fpf, "wb");
        stats->buf.sz = 0;
        stats->buf.buf = NULL,
                   res = (stats->file != NULL);
    }
    else
    {
#if 0
#elif USE_POSIX_MMAP
        struct stat stat_buf;
        int fd;

        fd = open(fpf, O_RDONLY);
        stats->file = fdopen(fd, "rb");
        fstat(fd, &stat_buf);
        stats->buf.sz = stat_buf.st_size;
        stats->buf.buf = mmap(NULL, stats->buf.sz, PROT_READ, MAP_PRIVATE,
                              fd, 0);
        res = (stats->buf.buf != NULL);
#else
        size_t nbytes;

        stats->file = fopen(fpf, "rb");

        if (fseek(stats->file, 0, SEEK_END))
        {
            fprintf(stderr, "First-pass stats file must be seekable!\n");
            exit(EXIT_FAILURE);
        }

        stats->buf.sz = stats->buf_alloc_sz = ftell(stats->file);
        rewind(stats->file);

        stats->buf.buf = malloc(stats->buf_alloc_sz);

        if (!stats->buf.buf)
        {
            fprintf(stderr, "Failed to allocate first-pass stats buffer (%lu bytes)\n",
                    (unsigned long)stats->buf_alloc_sz);
            exit(EXIT_FAILURE);
        }

        nbytes = fread(stats->buf.buf, 1, stats->buf.sz, stats->file);
        res = (nbytes == stats->buf.sz);
#endif
    }

    return res;
}

int stats_open_mem(stats_io_t *stats, int pass)
{
    int res;
    stats->pass = pass;

    if (!pass)
    {
        stats->buf.sz = 0;
        stats->buf_alloc_sz = 64 * 1024;
        stats->buf.buf = malloc(stats->buf_alloc_sz);
    }

    stats->buf_ptr = stats->buf.buf;
    res = (stats->buf.buf != NULL);
    return res;
}


void stats_close(stats_io_t *stats, int last_pass)
{
    if (stats->file)
    {
        if (stats->pass == last_pass)
        {
#if 0
#elif USE_POSIX_MMAP
            munmap(stats->buf.buf, stats->buf.sz);
#else
            free(stats->buf.buf);
#endif
        }

        fclose(stats->file);
        stats->file = NULL;
    }
    else
    {
        if (stats->pass == last_pass)
            free(stats->buf.buf);
    }
}

void stats_write(stats_io_t *stats, const void *pkt, size_t len)
{
    if (stats->file)
    {
        if(fwrite(pkt, 1, len, stats->file));
    }
    else
    {
        if (stats->buf.sz + len > stats->buf_alloc_sz)
        {
            size_t  new_sz = stats->buf_alloc_sz + 64 * 1024;
            char   *new_ptr = realloc(stats->buf.buf, new_sz);

            if (new_ptr)
            {
                stats->buf_ptr = new_ptr + (stats->buf_ptr - (char *)stats->buf.buf);
                stats->buf.buf = new_ptr;
                stats->buf_alloc_sz = new_sz;
            }
            else
            {
                fprintf(stderr,
                        "\nFailed to realloc firstpass stats buffer.\n");
                exit(EXIT_FAILURE);
            }
        }

        memcpy(stats->buf_ptr, pkt, len);
        stats->buf.sz += len;
        stats->buf_ptr += len;
    }
}

vpx_fixed_buf_t stats_get(stats_io_t *stats)
{
    return stats->buf;
}

enum video_file_type
{
    FILE_TYPE_RAW,
    FILE_TYPE_IVF,
    FILE_TYPE_Y4M
};

struct detect_buffer {
    char buf[4];
    size_t buf_read;
    size_t position;
};


#define IVF_FRAME_HDR_SZ (4+8) /* 4 byte size + 8 byte timestamp */
static int read_frame(FILE *f, vpx_image_t *img, unsigned int file_type,
                      y4m_input *y4m, struct detect_buffer *detect)
{
    int plane = 0;
    int shortread = 0;

    if (file_type == FILE_TYPE_Y4M)
    {
        if (y4m_input_fetch_frame(y4m, f, img) < 1)
           return 0;
    }
    else
    {
        if (file_type == FILE_TYPE_IVF)
        {
            char junk[IVF_FRAME_HDR_SZ];

            /* Skip the frame header. We know how big the frame should be. See
             * write_ivf_frame_header() for documentation on the frame header
             * layout.
             */
            if(fread(junk, 1, IVF_FRAME_HDR_SZ, f));
        }

        for (plane = 0; plane < 3; plane++)
        {
            unsigned char *ptr;
            int w = (plane ? (1 + img->d_w) / 2 : img->d_w);
            int h = (plane ? (1 + img->d_h) / 2 : img->d_h);
            int r;

            /* Determine the correct plane based on the image format. The for-loop
             * always counts in Y,U,V order, but this may not match the order of
             * the data on disk.
             */
            switch (plane)
            {
            case 1:
                ptr = img->planes[img->fmt==VPX_IMG_FMT_YV12? VPX_PLANE_V : VPX_PLANE_U];
                break;
            case 2:
                ptr = img->planes[img->fmt==VPX_IMG_FMT_YV12?VPX_PLANE_U : VPX_PLANE_V];
                break;
            default:
                ptr = img->planes[plane];
            }

            for (r = 0; r < h; r++)
            {
                size_t needed = w;
                size_t buf_position = 0;
                const size_t left = detect->buf_read - detect->position;
                if (left > 0)
                {
                    const size_t more = (left < needed) ? left : needed;
                    memcpy(ptr, detect->buf + detect->position, more);
                    buf_position = more;
                    needed -= more;
                    detect->position += more;
                }
                if (needed > 0)
                {
                    shortread |= (fread(ptr + buf_position, 1, needed, f) < needed);
                }

                ptr += img->stride[plane];
            }
        }
    }

    return !shortread;
}


unsigned int file_is_y4m(FILE      *infile,
                         y4m_input *y4m,
                         char       detect[4])
{
    if(memcmp(detect, "YUV4", 4) == 0)
    {
        return 1;
    }
    return 0;
}

#define IVF_FILE_HDR_SZ (32)
unsigned int file_is_ivf(FILE *infile,
                         unsigned int *fourcc,
                         unsigned int *width,
                         unsigned int *height,
                         struct detect_buffer *detect)
{
    char raw_hdr[IVF_FILE_HDR_SZ];
    int is_ivf = 0;

    if(memcmp(detect->buf, "DKIF", 4) != 0)
        return 0;

    /* See write_ivf_file_header() for more documentation on the file header
     * layout.
     */
    if (fread(raw_hdr + 4, 1, IVF_FILE_HDR_SZ - 4, infile)
        == IVF_FILE_HDR_SZ - 4)
    {
        {
            is_ivf = 1;

            if (mem_get_le16(raw_hdr + 4) != 0)
                fprintf(stderr, "Error: Unrecognized IVF version! This file may not"
                        " decode properly.");

            *fourcc = mem_get_le32(raw_hdr + 8);
        }
    }

    if (is_ivf)
    {
        *width = mem_get_le16(raw_hdr + 12);
        *height = mem_get_le16(raw_hdr + 14);
        detect->position = 4;
    }

    return is_ivf;
}


static void write_ivf_file_header(FILE *outfile,
                                  const vpx_codec_enc_cfg_t *cfg,
                                  unsigned int fourcc,
                                  int frame_cnt)
{
    char header[32];

    if (cfg->g_pass != VPX_RC_ONE_PASS && cfg->g_pass != VPX_RC_LAST_PASS)
        return;

    header[0] = 'D';
    header[1] = 'K';
    header[2] = 'I';
    header[3] = 'F';
    mem_put_le16(header + 4,  0);                 /* version */
    mem_put_le16(header + 6,  32);                /* headersize */
    mem_put_le32(header + 8,  fourcc);            /* headersize */
    mem_put_le16(header + 12, cfg->g_w);          /* width */
    mem_put_le16(header + 14, cfg->g_h);          /* height */
    mem_put_le32(header + 16, cfg->g_timebase.den); /* rate */
    mem_put_le32(header + 20, cfg->g_timebase.num); /* scale */
    mem_put_le32(header + 24, frame_cnt);         /* length */
    mem_put_le32(header + 28, 0);                 /* unused */

    if(fwrite(header, 1, 32, outfile));
}


static void write_ivf_frame_header(FILE *outfile,
                                   const vpx_codec_cx_pkt_t *pkt)
{
    char             header[12];
    vpx_codec_pts_t  pts;

    if (pkt->kind != VPX_CODEC_CX_FRAME_PKT)
        return;

    pts = pkt->data.frame.pts;
    mem_put_le32(header, pkt->data.frame.sz);
    mem_put_le32(header + 4, pts & 0xFFFFFFFF);
    mem_put_le32(header + 8, pts >> 32);

    if(fwrite(header, 1, 12, outfile));
}


typedef off_t EbmlLoc;


struct cue_entry
{
    unsigned int time;
    uint64_t     loc;
};


struct EbmlGlobal
{
    int debug;

    FILE    *stream;
    int64_t last_pts_ms;
    vpx_rational_t  framerate;

    /* These pointers are to the start of an element */
    off_t    position_reference;
    off_t    seek_info_pos;
    off_t    segment_info_pos;
    off_t    track_pos;
    off_t    cue_pos;
    off_t    cluster_pos;

    /* This pointer is to a specific element to be serialized */
    off_t    track_id_pos;

    /* These pointers are to the size field of the element */
    EbmlLoc  startSegment;
    EbmlLoc  startCluster;

    uint32_t cluster_timecode;
    int      cluster_open;

    struct cue_entry *cue_list;
    unsigned int      cues;

};


void Ebml_Write(EbmlGlobal *glob, const void *buffer_in, unsigned long len)
{
    if(fwrite(buffer_in, 1, len, glob->stream));
}


void Ebml_Serialize(EbmlGlobal *glob, const void *buffer_in, unsigned long len)
{
    const unsigned char *q = (const unsigned char *)buffer_in + len - 1;

    for(; len; len--)
        Ebml_Write(glob, q--, 1);
}


/* Need a fixed size serializer for the track ID. libmkv provdes a 64 bit
 * one, but not a 32 bit one.
 */
static void Ebml_SerializeUnsigned32(EbmlGlobal *glob, unsigned long class_id, uint64_t ui)
{
    unsigned char sizeSerialized = 4 | 0x80;
    Ebml_WriteID(glob, class_id);
    Ebml_Serialize(glob, &sizeSerialized, 1);
    Ebml_Serialize(glob, &ui, 4);
}


static void
Ebml_StartSubElement(EbmlGlobal *glob, EbmlLoc *ebmlLoc,
                          unsigned long class_id)
{
    //todo this is always taking 8 bytes, this may need later optimization
    //this is a key that says lenght unknown
    unsigned long long unknownLen =  LITERALU64(0x01FFFFFFFFFFFFFF);

    Ebml_WriteID(glob, class_id);
    *ebmlLoc = ftello(glob->stream);
    Ebml_Serialize(glob, &unknownLen, 8);
}

static void
Ebml_EndSubElement(EbmlGlobal *glob, EbmlLoc *ebmlLoc)
{
    off_t pos;
    uint64_t size;

    /* Save the current stream pointer */
    pos = ftello(glob->stream);

    /* Calculate the size of this element */
    size = pos - *ebmlLoc - 8;
    size |=  LITERALU64(0x0100000000000000);

    /* Seek back to the beginning of the element and write the new size */
    fseeko(glob->stream, *ebmlLoc, SEEK_SET);
    Ebml_Serialize(glob, &size, 8);

    /* Reset the stream pointer */
    fseeko(glob->stream, pos, SEEK_SET);
}


static void
write_webm_seek_element(EbmlGlobal *ebml, unsigned long id, off_t pos)
{
    uint64_t offset = pos - ebml->position_reference;
    EbmlLoc start;
    Ebml_StartSubElement(ebml, &start, Seek);
    Ebml_SerializeBinary(ebml, SeekID, id);
    Ebml_SerializeUnsigned64(ebml, SeekPosition, offset);
    Ebml_EndSubElement(ebml, &start);
}


static void
write_webm_seek_info(EbmlGlobal *ebml)
{

    off_t pos;

    /* Save the current stream pointer */
    pos = ftello(ebml->stream);

    if(ebml->seek_info_pos)
        fseeko(ebml->stream, ebml->seek_info_pos, SEEK_SET);
    else
        ebml->seek_info_pos = pos;

    {
        EbmlLoc start;

        Ebml_StartSubElement(ebml, &start, SeekHead);
        write_webm_seek_element(ebml, Tracks, ebml->track_pos);
        write_webm_seek_element(ebml, Cues,   ebml->cue_pos);
        write_webm_seek_element(ebml, Info,   ebml->segment_info_pos);
        Ebml_EndSubElement(ebml, &start);
    }
    {
        //segment info
        EbmlLoc startInfo;
        uint64_t frame_time;

        frame_time = (uint64_t)1000 * ebml->framerate.den
                     / ebml->framerate.num;
        ebml->segment_info_pos = ftello(ebml->stream);
        Ebml_StartSubElement(ebml, &startInfo, Info);
        Ebml_SerializeUnsigned(ebml, TimecodeScale, 1000000);
        Ebml_SerializeFloat(ebml, Segment_Duration,
                            ebml->last_pts_ms + frame_time);
        Ebml_SerializeString(ebml, 0x4D80,
            ebml->debug ? "vpxenc" : "vpxenc" VERSION_STRING);
        Ebml_SerializeString(ebml, 0x5741,
            ebml->debug ? "vpxenc" : "vpxenc" VERSION_STRING);
        Ebml_EndSubElement(ebml, &startInfo);
    }
}


static void
write_webm_file_header(EbmlGlobal                *glob,
                       const vpx_codec_enc_cfg_t *cfg,
                       const struct vpx_rational *fps)
{
    {
        EbmlLoc start;
        Ebml_StartSubElement(glob, &start, EBML);
        Ebml_SerializeUnsigned(glob, EBMLVersion, 1);
        Ebml_SerializeUnsigned(glob, EBMLReadVersion, 1); //EBML Read Version
        Ebml_SerializeUnsigned(glob, EBMLMaxIDLength, 4); //EBML Max ID Length
        Ebml_SerializeUnsigned(glob, EBMLMaxSizeLength, 8); //EBML Max Size Length
        Ebml_SerializeString(glob, DocType, "webm"); //Doc Type
        Ebml_SerializeUnsigned(glob, DocTypeVersion, 2); //Doc Type Version
        Ebml_SerializeUnsigned(glob, DocTypeReadVersion, 2); //Doc Type Read Version
        Ebml_EndSubElement(glob, &start);
    }
    {
        Ebml_StartSubElement(glob, &glob->startSegment, Segment); //segment
        glob->position_reference = ftello(glob->stream);
        glob->framerate = *fps;
        write_webm_seek_info(glob);

        {
            EbmlLoc trackStart;
            glob->track_pos = ftello(glob->stream);
            Ebml_StartSubElement(glob, &trackStart, Tracks);
            {
                unsigned int trackNumber = 1;
                uint64_t     trackID = 0;

                EbmlLoc start;
                Ebml_StartSubElement(glob, &start, TrackEntry);
                Ebml_SerializeUnsigned(glob, TrackNumber, trackNumber);
                glob->track_id_pos = ftello(glob->stream);
                Ebml_SerializeUnsigned32(glob, TrackUID, trackID);
                Ebml_SerializeUnsigned(glob, TrackType, 1); //video is always 1
                Ebml_SerializeString(glob, CodecID, "V_VP8");
                {
                    unsigned int pixelWidth = cfg->g_w;
                    unsigned int pixelHeight = cfg->g_h;
                    float        frameRate   = (float)fps->num/(float)fps->den;

                    EbmlLoc videoStart;
                    Ebml_StartSubElement(glob, &videoStart, Video);
                    Ebml_SerializeUnsigned(glob, PixelWidth, pixelWidth);
                    Ebml_SerializeUnsigned(glob, PixelHeight, pixelHeight);
                    Ebml_SerializeFloat(glob, FrameRate, frameRate);
                    Ebml_EndSubElement(glob, &videoStart); //Video
                }
                Ebml_EndSubElement(glob, &start); //Track Entry
            }
            Ebml_EndSubElement(glob, &trackStart);
        }
        // segment element is open
    }
}


static void
write_webm_block(EbmlGlobal                *glob,
                 const vpx_codec_enc_cfg_t *cfg,
                 const vpx_codec_cx_pkt_t  *pkt)
{
    unsigned long  block_length;
    unsigned char  track_number;
    unsigned short block_timecode = 0;
    unsigned char  flags;
    int64_t        pts_ms;
    int            start_cluster = 0, is_keyframe;

    /* Calculate the PTS of this frame in milliseconds */
    pts_ms = pkt->data.frame.pts * 1000
             * (uint64_t)cfg->g_timebase.num / (uint64_t)cfg->g_timebase.den;
    if(pts_ms <= glob->last_pts_ms)
        pts_ms = glob->last_pts_ms + 1;
    glob->last_pts_ms = pts_ms;

    /* Calculate the relative time of this block */
    if(pts_ms - glob->cluster_timecode > SHRT_MAX)
        start_cluster = 1;
    else
        block_timecode = pts_ms - glob->cluster_timecode;

    is_keyframe = (pkt->data.frame.flags & VPX_FRAME_IS_KEY);
    if(start_cluster || is_keyframe)
    {
        if(glob->cluster_open)
            Ebml_EndSubElement(glob, &glob->startCluster);

        /* Open the new cluster */
        block_timecode = 0;
        glob->cluster_open = 1;
        glob->cluster_timecode = pts_ms;
        glob->cluster_pos = ftello(glob->stream);
        Ebml_StartSubElement(glob, &glob->startCluster, Cluster); //cluster
        Ebml_SerializeUnsigned(glob, Timecode, glob->cluster_timecode);

        /* Save a cue point if this is a keyframe. */
        if(is_keyframe)
        {
            struct cue_entry *cue, *new_cue_list;

            new_cue_list = realloc(glob->cue_list,
                                   (glob->cues+1) * sizeof(struct cue_entry));
            if(new_cue_list)
                glob->cue_list = new_cue_list;
            else
            {
                fprintf(stderr, "\nFailed to realloc cue list.\n");
                exit(EXIT_FAILURE);
            }

            cue = &glob->cue_list[glob->cues];
            cue->time = glob->cluster_timecode;
            cue->loc = glob->cluster_pos;
            glob->cues++;
        }
    }

    /* Write the Simple Block */
    Ebml_WriteID(glob, SimpleBlock);

    block_length = pkt->data.frame.sz + 4;
    block_length |= 0x10000000;
    Ebml_Serialize(glob, &block_length, 4);

    track_number = 1;
    track_number |= 0x80;
    Ebml_Write(glob, &track_number, 1);

    Ebml_Serialize(glob, &block_timecode, 2);

    flags = 0;
    if(is_keyframe)
        flags |= 0x80;
    if(pkt->data.frame.flags & VPX_FRAME_IS_INVISIBLE)
        flags |= 0x08;
    Ebml_Write(glob, &flags, 1);

    Ebml_Write(glob, pkt->data.frame.buf, pkt->data.frame.sz);
}


static void
write_webm_file_footer(EbmlGlobal *glob, long hash)
{

    if(glob->cluster_open)
        Ebml_EndSubElement(glob, &glob->startCluster);

    {
        EbmlLoc start;
        int i;

        glob->cue_pos = ftello(glob->stream);
        Ebml_StartSubElement(glob, &start, Cues);
        for(i=0; i<glob->cues; i++)
        {
            struct cue_entry *cue = &glob->cue_list[i];
            EbmlLoc start;

            Ebml_StartSubElement(glob, &start, CuePoint);
            {
                EbmlLoc start;

                Ebml_SerializeUnsigned(glob, CueTime, cue->time);

                Ebml_StartSubElement(glob, &start, CueTrackPositions);
                Ebml_SerializeUnsigned(glob, CueTrack, 1);
                Ebml_SerializeUnsigned64(glob, CueClusterPosition,
                                         cue->loc - glob->position_reference);
                //Ebml_SerializeUnsigned(glob, CueBlockNumber, cue->blockNumber);
                Ebml_EndSubElement(glob, &start);
            }
            Ebml_EndSubElement(glob, &start);
        }
        Ebml_EndSubElement(glob, &start);
    }

    Ebml_EndSubElement(glob, &glob->startSegment);

    /* Patch up the seek info block */
    write_webm_seek_info(glob);

    /* Patch up the track id */
    fseeko(glob->stream, glob->track_id_pos, SEEK_SET);
    Ebml_SerializeUnsigned32(glob, TrackUID, glob->debug ? 0xDEADBEEF : hash);

    fseeko(glob->stream, 0, SEEK_END);
}


/* Murmur hash derived from public domain reference implementation at
 *   http://sites.google.com/site/murmurhash/
 */
static unsigned int murmur ( const void * key, int len, unsigned int seed )
{
    const unsigned int m = 0x5bd1e995;
    const int r = 24;

    unsigned int h = seed ^ len;

    const unsigned char * data = (const unsigned char *)key;

    while(len >= 4)
    {
        unsigned int k;

        k  = data[0];
        k |= data[1] << 8;
        k |= data[2] << 16;
        k |= data[3] << 24;

        k *= m;
        k ^= k >> r;
        k *= m;

        h *= m;
        h ^= k;

        data += 4;
        len -= 4;
    }

    switch(len)
    {
    case 3: h ^= data[2] << 16;
    case 2: h ^= data[1] << 8;
    case 1: h ^= data[0];
            h *= m;
    };

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
}

#include "math.h"

static double vp8_mse2psnr(double Samples, double Peak, double Mse)
{
    double psnr;

    if ((double)Mse > 0.0)
        psnr = 10.0 * log10(Peak * Peak * Samples / Mse);
    else
        psnr = 60;      // Limit to prevent / 0

    if (psnr > 60)
        psnr = 60;

    return psnr;
}


#include "args.h"

static const arg_def_t debugmode = ARG_DEF("D", "debug", 0,
        "Debug mode (makes output deterministic)");
static const arg_def_t outputfile = ARG_DEF("o", "output", 1,
        "Output filename");
static const arg_def_t use_yv12 = ARG_DEF(NULL, "yv12", 0,
                                  "Input file is YV12 ");
static const arg_def_t use_i420 = ARG_DEF(NULL, "i420", 0,
                                  "Input file is I420 (default)");
static const arg_def_t codecarg = ARG_DEF(NULL, "codec", 1,
                                  "Codec to use");
static const arg_def_t passes           = ARG_DEF("p", "passes", 1,
        "Number of passes (1/2)");
static const arg_def_t pass_arg         = ARG_DEF(NULL, "pass", 1,
        "Pass to execute (1/2)");
static const arg_def_t fpf_name         = ARG_DEF(NULL, "fpf", 1,
        "First pass statistics file name");
static const arg_def_t limit = ARG_DEF(NULL, "limit", 1,
                                       "Stop encoding after n input frames");
static const arg_def_t deadline         = ARG_DEF("d", "deadline", 1,
        "Deadline per frame (usec)");
static const arg_def_t best_dl          = ARG_DEF(NULL, "best", 0,
        "Use Best Quality Deadline");
static const arg_def_t good_dl          = ARG_DEF(NULL, "good", 0,
        "Use Good Quality Deadline");
static const arg_def_t rt_dl            = ARG_DEF(NULL, "rt", 0,
        "Use Realtime Quality Deadline");
static const arg_def_t verbosearg       = ARG_DEF("v", "verbose", 0,
        "Show encoder parameters");
static const arg_def_t psnrarg          = ARG_DEF(NULL, "psnr", 0,
        "Show PSNR in status line");
static const arg_def_t framerate        = ARG_DEF(NULL, "fps", 1,
        "Stream frame rate (rate/scale)");
static const arg_def_t use_ivf          = ARG_DEF(NULL, "ivf", 0,
        "Output IVF (default is WebM)");
static const arg_def_t *main_args[] =
{
    &debugmode,
    &outputfile, &codecarg, &passes, &pass_arg, &fpf_name, &limit, &deadline,
    &best_dl, &good_dl, &rt_dl,
    &verbosearg, &psnrarg, &use_ivf, &framerate,
    NULL
};

static const arg_def_t usage            = ARG_DEF("u", "usage", 1,
        "Usage profile number to use");
static const arg_def_t threads          = ARG_DEF("t", "threads", 1,
        "Max number of threads to use");
static const arg_def_t profile          = ARG_DEF(NULL, "profile", 1,
        "Bitstream profile number to use");
static const arg_def_t width            = ARG_DEF("w", "width", 1,
        "Frame width");
static const arg_def_t height           = ARG_DEF("h", "height", 1,
        "Frame height");
static const arg_def_t timebase         = ARG_DEF(NULL, "timebase", 1,
        "Stream timebase (frame duration)");
static const arg_def_t error_resilient  = ARG_DEF(NULL, "error-resilient", 1,
        "Enable error resiliency features");
static const arg_def_t lag_in_frames    = ARG_DEF(NULL, "lag-in-frames", 1,
        "Max number of frames to lag");

static const arg_def_t *global_args[] =
{
    &use_yv12, &use_i420, &usage, &threads, &profile,
    &width, &height, &timebase, &framerate, &error_resilient,
    &lag_in_frames, NULL
};

static const arg_def_t dropframe_thresh   = ARG_DEF(NULL, "drop-frame", 1,
        "Temporal resampling threshold (buf %)");
static const arg_def_t resize_allowed     = ARG_DEF(NULL, "resize-allowed", 1,
        "Spatial resampling enabled (bool)");
static const arg_def_t resize_up_thresh   = ARG_DEF(NULL, "resize-up", 1,
        "Upscale threshold (buf %)");
static const arg_def_t resize_down_thresh = ARG_DEF(NULL, "resize-down", 1,
        "Downscale threshold (buf %)");
static const struct arg_enum_list end_usage_enum[] = {
    {"vbr", VPX_VBR},
    {"cbr", VPX_CBR},
    {"cq",  VPX_CQ},
    {NULL, 0}
};
static const arg_def_t end_usage          = ARG_DEF_ENUM(NULL, "end-usage", 1,
        "Rate control mode", end_usage_enum);
static const arg_def_t target_bitrate     = ARG_DEF(NULL, "target-bitrate", 1,
        "Bitrate (kbps)");
static const arg_def_t min_quantizer      = ARG_DEF(NULL, "min-q", 1,
        "Minimum (best) quantizer");
static const arg_def_t max_quantizer      = ARG_DEF(NULL, "max-q", 1,
        "Maximum (worst) quantizer");
static const arg_def_t undershoot_pct     = ARG_DEF(NULL, "undershoot-pct", 1,
        "Datarate undershoot (min) target (%)");
static const arg_def_t overshoot_pct      = ARG_DEF(NULL, "overshoot-pct", 1,
        "Datarate overshoot (max) target (%)");
static const arg_def_t buf_sz             = ARG_DEF(NULL, "buf-sz", 1,
        "Client buffer size (ms)");
static const arg_def_t buf_initial_sz     = ARG_DEF(NULL, "buf-initial-sz", 1,
        "Client initial buffer size (ms)");
static const arg_def_t buf_optimal_sz     = ARG_DEF(NULL, "buf-optimal-sz", 1,
        "Client optimal buffer size (ms)");
static const arg_def_t *rc_args[] =
{
    &dropframe_thresh, &resize_allowed, &resize_up_thresh, &resize_down_thresh,
    &end_usage, &target_bitrate, &min_quantizer, &max_quantizer,
    &undershoot_pct, &overshoot_pct, &buf_sz, &buf_initial_sz, &buf_optimal_sz,
    NULL
};


static const arg_def_t bias_pct = ARG_DEF(NULL, "bias-pct", 1,
                                  "CBR/VBR bias (0=CBR, 100=VBR)");
static const arg_def_t minsection_pct = ARG_DEF(NULL, "minsection-pct", 1,
                                        "GOP min bitrate (% of target)");
static const arg_def_t maxsection_pct = ARG_DEF(NULL, "maxsection-pct", 1,
                                        "GOP max bitrate (% of target)");
static const arg_def_t *rc_twopass_args[] =
{
    &bias_pct, &minsection_pct, &maxsection_pct, NULL
};


static const arg_def_t kf_min_dist = ARG_DEF(NULL, "kf-min-dist", 1,
                                     "Minimum keyframe interval (frames)");
static const arg_def_t kf_max_dist = ARG_DEF(NULL, "kf-max-dist", 1,
                                     "Maximum keyframe interval (frames)");
static const arg_def_t kf_disabled = ARG_DEF(NULL, "disable-kf", 0,
                                     "Disable keyframe placement");
static const arg_def_t *kf_args[] =
{
    &kf_min_dist, &kf_max_dist, &kf_disabled, NULL
};


#if CONFIG_VP8_ENCODER
static const arg_def_t noise_sens = ARG_DEF(NULL, "noise-sensitivity", 1,
                                    "Noise sensitivity (frames to blur)");
static const arg_def_t sharpness = ARG_DEF(NULL, "sharpness", 1,
                                   "Filter sharpness (0-7)");
static const arg_def_t static_thresh = ARG_DEF(NULL, "static-thresh", 1,
                                       "Motion detection threshold");
#endif

#if CONFIG_VP8_ENCODER
static const arg_def_t cpu_used = ARG_DEF(NULL, "cpu-used", 1,
                                  "CPU Used (-16..16)");
#endif


#if CONFIG_VP8_ENCODER
static const arg_def_t token_parts = ARG_DEF(NULL, "token-parts", 1,
                                     "Number of token partitions to use, log2");
static const arg_def_t auto_altref = ARG_DEF(NULL, "auto-alt-ref", 1,
                                     "Enable automatic alt reference frames");
static const arg_def_t arnr_maxframes = ARG_DEF(NULL, "arnr-maxframes", 1,
                                        "AltRef Max Frames");
static const arg_def_t arnr_strength = ARG_DEF(NULL, "arnr-strength", 1,
                                       "AltRef Strength");
static const arg_def_t arnr_type = ARG_DEF(NULL, "arnr-type", 1,
                                   "AltRef Type");
static const struct arg_enum_list tuning_enum[] = {
    {"psnr", VP8_TUNE_PSNR},
    {"ssim", VP8_TUNE_SSIM},
    {NULL, 0}
};
static const arg_def_t tune_ssim = ARG_DEF_ENUM(NULL, "tune", 1,
                                   "Material to favor", tuning_enum);
static const arg_def_t cq_level = ARG_DEF(NULL, "cq-level", 1,
                                   "Constrained Quality Level");

static const arg_def_t *vp8_args[] =
{
    &cpu_used, &auto_altref, &noise_sens, &sharpness, &static_thresh,
    &token_parts, &arnr_maxframes, &arnr_strength, &arnr_type,
    &tune_ssim, &cq_level, NULL
};
static const int vp8_arg_ctrl_map[] =
{
    VP8E_SET_CPUUSED, VP8E_SET_ENABLEAUTOALTREF,
    VP8E_SET_NOISE_SENSITIVITY, VP8E_SET_SHARPNESS, VP8E_SET_STATIC_THRESHOLD,
    VP8E_SET_TOKEN_PARTITIONS,
    VP8E_SET_ARNR_MAXFRAMES, VP8E_SET_ARNR_STRENGTH , VP8E_SET_ARNR_TYPE,
    VP8E_SET_TUNING, VP8E_SET_CQ_LEVEL, 0
};
#endif

static const arg_def_t *no_args[] = { NULL };

static void usage_exit()
{
    int i;

    fprintf(stderr, "Usage: %s <options> -o dst_filename src_filename \n",
            exec_name);

    fprintf(stderr, "\nOptions:\n");
    arg_show_usage(stdout, main_args);
    fprintf(stderr, "\nEncoder Global Options:\n");
    arg_show_usage(stdout, global_args);
    fprintf(stderr, "\nRate Control Options:\n");
    arg_show_usage(stdout, rc_args);
    fprintf(stderr, "\nTwopass Rate Control Options:\n");
    arg_show_usage(stdout, rc_twopass_args);
    fprintf(stderr, "\nKeyframe Placement Options:\n");
    arg_show_usage(stdout, kf_args);
#if CONFIG_VP8_ENCODER
    fprintf(stderr, "\nVP8 Specific Options:\n");
    arg_show_usage(stdout, vp8_args);
#endif
    fprintf(stderr, "\n"
           "Included encoders:\n"
           "\n");

    for (i = 0; i < sizeof(codecs) / sizeof(codecs[0]); i++)
        fprintf(stderr, "    %-6s - %s\n",
               codecs[i].name,
               vpx_codec_iface_name(codecs[i].iface));

    exit(EXIT_FAILURE);
}

#define ARG_CTRL_CNT_MAX 10


int main(int argc, const char **argv_)
{
    vpx_codec_ctx_t        encoder;
    const char                  *in_fn = NULL, *out_fn = NULL, *stats_fn = NULL;
    int                    i;
    FILE                  *infile, *outfile;
    vpx_codec_enc_cfg_t    cfg;
    vpx_codec_err_t        res;
    int                    pass, one_pass_only = 0;
    stats_io_t             stats;
    vpx_image_t            raw;
    const struct codec_item  *codec = codecs;
    int                    frame_avail, got_data;

    struct arg               arg;
    char                   **argv, **argi, **argj;
    int                      arg_usage = 0, arg_passes = 1, arg_deadline = 0;
    int                      arg_ctrls[ARG_CTRL_CNT_MAX][2], arg_ctrl_cnt = 0;
    int                      arg_limit = 0;
    static const arg_def_t **ctrl_args = no_args;
    static const int        *ctrl_args_map = NULL;
    int                      verbose = 0, show_psnr = 0;
    int                      arg_use_i420 = 1;
    unsigned long            cx_time = 0;
    unsigned int             file_type, fourcc;
    y4m_input                y4m;
    struct vpx_rational      arg_framerate = {30, 1};
    int                      arg_have_framerate = 0;
    int                      write_webm = 1;
    EbmlGlobal               ebml = {0};
    uint32_t                 hash = 0;
    uint64_t                 psnr_sse_total = 0;
    uint64_t                 psnr_samples_total = 0;
    double                   psnr_totals[4] = {0, 0, 0, 0};
    int                      psnr_count = 0;

    exec_name = argv_[0];
    ebml.last_pts_ms = -1;

    if (argc < 3)
        usage_exit();


    /* First parse the codec and usage values, because we want to apply other
     * parameters on top of the default configuration provided by the codec.
     */
    argv = argv_dup(argc - 1, argv_ + 1);

    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        arg.argv_step = 1;

        if (arg_match(&arg, &codecarg, argi))
        {
            int j, k = -1;

            for (j = 0; j < sizeof(codecs) / sizeof(codecs[0]); j++)
                if (!strcmp(codecs[j].name, arg.val))
                    k = j;

            if (k >= 0)
                codec = codecs + k;
            else
                die("Error: Unrecognized argument (%s) to --codec\n",
                    arg.val);

        }
        else if (arg_match(&arg, &passes, argi))
        {
            arg_passes = arg_parse_uint(&arg);

            if (arg_passes < 1 || arg_passes > 2)
                die("Error: Invalid number of passes (%d)\n", arg_passes);
        }
        else if (arg_match(&arg, &pass_arg, argi))
        {
            one_pass_only = arg_parse_uint(&arg);

            if (one_pass_only < 1 || one_pass_only > 2)
                die("Error: Invalid pass selected (%d)\n", one_pass_only);
        }
        else if (arg_match(&arg, &fpf_name, argi))
            stats_fn = arg.val;
        else if (arg_match(&arg, &usage, argi))
            arg_usage = arg_parse_uint(&arg);
        else if (arg_match(&arg, &deadline, argi))
            arg_deadline = arg_parse_uint(&arg);
        else if (arg_match(&arg, &best_dl, argi))
            arg_deadline = VPX_DL_BEST_QUALITY;
        else if (arg_match(&arg, &good_dl, argi))
            arg_deadline = VPX_DL_GOOD_QUALITY;
        else if (arg_match(&arg, &rt_dl, argi))
            arg_deadline = VPX_DL_REALTIME;
        else if (arg_match(&arg, &use_yv12, argi))
        {
            arg_use_i420 = 0;
        }
        else if (arg_match(&arg, &use_i420, argi))
        {
            arg_use_i420 = 1;
        }
        else if (arg_match(&arg, &verbosearg, argi))
            verbose = 1;
        else if (arg_match(&arg, &limit, argi))
            arg_limit = arg_parse_uint(&arg);
        else if (arg_match(&arg, &psnrarg, argi))
            show_psnr = 1;
        else if (arg_match(&arg, &framerate, argi))
        {
            arg_framerate = arg_parse_rational(&arg);
            arg_have_framerate = 1;
        }
        else if (arg_match(&arg, &use_ivf, argi))
            write_webm = 0;
        else if (arg_match(&arg, &outputfile, argi))
            out_fn = arg.val;
        else if (arg_match(&arg, &debugmode, argi))
            ebml.debug = 1;
        else
            argj++;
    }

    /* Ensure that --passes and --pass are consistent. If --pass is set and --passes=2,
     * ensure --fpf was set.
     */
    if (one_pass_only)
    {
        /* DWIM: Assume the user meant passes=2 if pass=2 is specified */
        if (one_pass_only > arg_passes)
        {
            fprintf(stderr, "Warning: Assuming --pass=%d implies --passes=%d\n",
                   one_pass_only, one_pass_only);
            arg_passes = one_pass_only;
        }

        if (arg_passes == 2 && !stats_fn)
            die("Must specify --fpf when --pass=%d and --passes=2\n", one_pass_only);
    }

    /* Populate encoder configuration */
    res = vpx_codec_enc_config_default(codec->iface, &cfg, arg_usage);

    if (res)
    {
        fprintf(stderr, "Failed to get config: %s\n",
                vpx_codec_err_to_string(res));
        return EXIT_FAILURE;
    }

    /* Change the default timebase to a high enough value so that the encoder
     * will always create strictly increasing timestamps.
     */
    cfg.g_timebase.den = 1000;

    /* Never use the library's default resolution, require it be parsed
     * from the file or set on the command line.
     */
    cfg.g_w = 0;
    cfg.g_h = 0;

    /* Now parse the remainder of the parameters. */
    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        arg.argv_step = 1;

        if (0);
        else if (arg_match(&arg, &threads, argi))
            cfg.g_threads = arg_parse_uint(&arg);
        else if (arg_match(&arg, &profile, argi))
            cfg.g_profile = arg_parse_uint(&arg);
        else if (arg_match(&arg, &width, argi))
            cfg.g_w = arg_parse_uint(&arg);
        else if (arg_match(&arg, &height, argi))
            cfg.g_h = arg_parse_uint(&arg);
        else if (arg_match(&arg, &timebase, argi))
            cfg.g_timebase = arg_parse_rational(&arg);
        else if (arg_match(&arg, &error_resilient, argi))
            cfg.g_error_resilient = arg_parse_uint(&arg);
        else if (arg_match(&arg, &lag_in_frames, argi))
            cfg.g_lag_in_frames = arg_parse_uint(&arg);
        else if (arg_match(&arg, &dropframe_thresh, argi))
            cfg.rc_dropframe_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_allowed, argi))
            cfg.rc_resize_allowed = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_up_thresh, argi))
            cfg.rc_resize_up_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_down_thresh, argi))
            cfg.rc_resize_down_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &resize_down_thresh, argi))
            cfg.rc_resize_down_thresh = arg_parse_uint(&arg);
        else if (arg_match(&arg, &end_usage, argi))
            cfg.rc_end_usage = arg_parse_enum_or_int(&arg);
        else if (arg_match(&arg, &target_bitrate, argi))
            cfg.rc_target_bitrate = arg_parse_uint(&arg);
        else if (arg_match(&arg, &min_quantizer, argi))
            cfg.rc_min_quantizer = arg_parse_uint(&arg);
        else if (arg_match(&arg, &max_quantizer, argi))
            cfg.rc_max_quantizer = arg_parse_uint(&arg);
        else if (arg_match(&arg, &undershoot_pct, argi))
            cfg.rc_undershoot_pct = arg_parse_uint(&arg);
        else if (arg_match(&arg, &overshoot_pct, argi))
            cfg.rc_overshoot_pct = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_sz, argi))
            cfg.rc_buf_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_initial_sz, argi))
            cfg.rc_buf_initial_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &buf_optimal_sz, argi))
            cfg.rc_buf_optimal_sz = arg_parse_uint(&arg);
        else if (arg_match(&arg, &bias_pct, argi))
        {
            cfg.rc_2pass_vbr_bias_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                fprintf(stderr,
                        "Warning: option %s ignored in one-pass mode.\n",
                        arg.name);
        }
        else if (arg_match(&arg, &minsection_pct, argi))
        {
            cfg.rc_2pass_vbr_minsection_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                fprintf(stderr,
                        "Warning: option %s ignored in one-pass mode.\n",
                        arg.name);
        }
        else if (arg_match(&arg, &maxsection_pct, argi))
        {
            cfg.rc_2pass_vbr_maxsection_pct = arg_parse_uint(&arg);

            if (arg_passes < 2)
                fprintf(stderr,
                        "Warning: option %s ignored in one-pass mode.\n",
                        arg.name);
        }
        else if (arg_match(&arg, &kf_min_dist, argi))
            cfg.kf_min_dist = arg_parse_uint(&arg);
        else if (arg_match(&arg, &kf_max_dist, argi))
            cfg.kf_max_dist = arg_parse_uint(&arg);
        else if (arg_match(&arg, &kf_disabled, argi))
            cfg.kf_mode = VPX_KF_DISABLED;
        else
            argj++;
    }

    /* Handle codec specific options */
#if CONFIG_VP8_ENCODER

    if (codec->iface == &vpx_codec_vp8_cx_algo)
    {
        ctrl_args = vp8_args;
        ctrl_args_map = vp8_arg_ctrl_map;
    }

#endif

    for (argi = argj = argv; (*argj = *argi); argi += arg.argv_step)
    {
        int match = 0;

        arg.argv_step = 1;

        for (i = 0; ctrl_args[i]; i++)
        {
            if (arg_match(&arg, ctrl_args[i], argi))
            {
                match = 1;

                if (arg_ctrl_cnt < ARG_CTRL_CNT_MAX)
                {
                    arg_ctrls[arg_ctrl_cnt][0] = ctrl_args_map[i];
                    arg_ctrls[arg_ctrl_cnt][1] = arg_parse_enum_or_int(&arg);
                    arg_ctrl_cnt++;
                }
            }
        }

        if (!match)
            argj++;
    }

    /* Check for unrecognized options */
    for (argi = argv; *argi; argi++)
        if (argi[0][0] == '-' && argi[0][1])
            die("Error: Unrecognized option %s\n", *argi);

    /* Handle non-option arguments */
    in_fn = argv[0];

    if (!in_fn)
        usage_exit();

    if(!out_fn)
        die("Error: Output file is required (specify with -o)\n");

    memset(&stats, 0, sizeof(stats));

    for (pass = one_pass_only ? one_pass_only - 1 : 0; pass < arg_passes; pass++)
    {
        int frames_in = 0, frames_out = 0;
        unsigned long nbytes = 0;
        struct detect_buffer detect;

        /* Parse certain options from the input file, if possible */
        infile = strcmp(in_fn, "-") ? fopen(in_fn, "rb")
                                    : set_binary_mode(stdin);

        if (!infile)
        {
            fprintf(stderr, "Failed to open input file\n");
            return EXIT_FAILURE;
        }

        /* For RAW input sources, these bytes will applied on the first frame
         *  in read_frame().
         */
        detect.buf_read = fread(detect.buf, 1, 4, infile);
        detect.position = 0;

        if (detect.buf_read == 4 && file_is_y4m(infile, &y4m, detect.buf))
        {
            if (y4m_input_open(&y4m, infile, detect.buf, 4) >= 0)
            {
                file_type = FILE_TYPE_Y4M;
                cfg.g_w = y4m.pic_w;
                cfg.g_h = y4m.pic_h;

                /* Use the frame rate from the file only if none was specified
                 * on the command-line.
                 */
                if (!arg_have_framerate)
                {
                    arg_framerate.num = y4m.fps_n;
                    arg_framerate.den = y4m.fps_d;
                }

                arg_use_i420 = 0;
            }
            else
            {
                fprintf(stderr, "Unsupported Y4M stream.\n");
                return EXIT_FAILURE;
            }
        }
        else if (detect.buf_read == 4 &&
                 file_is_ivf(infile, &fourcc, &cfg.g_w, &cfg.g_h, &detect))
        {
            file_type = FILE_TYPE_IVF;
            switch (fourcc)
            {
            case 0x32315659:
                arg_use_i420 = 0;
                break;
            case 0x30323449:
                arg_use_i420 = 1;
                break;
            default:
                fprintf(stderr, "Unsupported fourcc (%08x) in IVF\n", fourcc);
                return EXIT_FAILURE;
            }
        }
        else
        {
            file_type = FILE_TYPE_RAW;
        }

        if(!cfg.g_w || !cfg.g_h)
        {
            fprintf(stderr, "Specify stream dimensions with --width (-w) "
                            " and --height (-h).\n");
            return EXIT_FAILURE;
        }

#define SHOW(field) fprintf(stderr, "    %-28s = %d\n", #field, cfg.field)

        if (verbose && pass == 0)
        {
            fprintf(stderr, "Codec: %s\n", vpx_codec_iface_name(codec->iface));
            fprintf(stderr, "Source file: %s Format: %s\n", in_fn,
                    arg_use_i420 ? "I420" : "YV12");
            fprintf(stderr, "Destination file: %s\n", out_fn);
            fprintf(stderr, "Encoder parameters:\n");

            SHOW(g_usage);
            SHOW(g_threads);
            SHOW(g_profile);
            SHOW(g_w);
            SHOW(g_h);
            SHOW(g_timebase.num);
            SHOW(g_timebase.den);
            SHOW(g_error_resilient);
            SHOW(g_pass);
            SHOW(g_lag_in_frames);
            SHOW(rc_dropframe_thresh);
            SHOW(rc_resize_allowed);
            SHOW(rc_resize_up_thresh);
            SHOW(rc_resize_down_thresh);
            SHOW(rc_end_usage);
            SHOW(rc_target_bitrate);
            SHOW(rc_min_quantizer);
            SHOW(rc_max_quantizer);
            SHOW(rc_undershoot_pct);
            SHOW(rc_overshoot_pct);
            SHOW(rc_buf_sz);
            SHOW(rc_buf_initial_sz);
            SHOW(rc_buf_optimal_sz);
            SHOW(rc_2pass_vbr_bias_pct);
            SHOW(rc_2pass_vbr_minsection_pct);
            SHOW(rc_2pass_vbr_maxsection_pct);
            SHOW(kf_mode);
            SHOW(kf_min_dist);
            SHOW(kf_max_dist);
        }

        if(pass == (one_pass_only ? one_pass_only - 1 : 0)) {
            if (file_type == FILE_TYPE_Y4M)
                /*The Y4M reader does its own allocation.
                  Just initialize this here to avoid problems if we never read any
                   frames.*/
                memset(&raw, 0, sizeof(raw));
            else
                vpx_img_alloc(&raw, arg_use_i420 ? VPX_IMG_FMT_I420 : VPX_IMG_FMT_YV12,
                              cfg.g_w, cfg.g_h, 1);
        }

        outfile = strcmp(out_fn, "-") ? fopen(out_fn, "wb")
                                      : set_binary_mode(stdout);

        if (!outfile)
        {
            fprintf(stderr, "Failed to open output file\n");
            return EXIT_FAILURE;
        }

        if(write_webm && fseek(outfile, 0, SEEK_CUR))
        {
            fprintf(stderr, "WebM output to pipes not supported.\n");
            return EXIT_FAILURE;
        }

        if (stats_fn)
        {
            if (!stats_open_file(&stats, stats_fn, pass))
            {
                fprintf(stderr, "Failed to open statistics store\n");
                return EXIT_FAILURE;
            }
        }
        else
        {
            if (!stats_open_mem(&stats, pass))
            {
                fprintf(stderr, "Failed to open statistics store\n");
                return EXIT_FAILURE;
            }
        }

        cfg.g_pass = arg_passes == 2
                     ? pass ? VPX_RC_LAST_PASS : VPX_RC_FIRST_PASS
                 : VPX_RC_ONE_PASS;
#if VPX_ENCODER_ABI_VERSION > (1 + VPX_CODEC_ABI_VERSION)

        if (pass)
        {
            cfg.rc_twopass_stats_in = stats_get(&stats);
        }

#endif

        if(write_webm)
        {
            ebml.stream = outfile;
            write_webm_file_header(&ebml, &cfg, &arg_framerate);
        }
        else
            write_ivf_file_header(outfile, &cfg, codec->fourcc, 0);


        /* Construct Encoder Context */
        vpx_codec_enc_init(&encoder, codec->iface, &cfg,
                           show_psnr ? VPX_CODEC_USE_PSNR : 0);
        ctx_exit_on_error(&encoder, "Failed to initialize encoder");

        /* Note that we bypass the vpx_codec_control wrapper macro because
         * we're being clever to store the control IDs in an array. Real
         * applications will want to make use of the enumerations directly
         */
        for (i = 0; i < arg_ctrl_cnt; i++)
        {
            if (vpx_codec_control_(&encoder, arg_ctrls[i][0], arg_ctrls[i][1]))
                fprintf(stderr, "Error: Tried to set control %d = %d\n",
                        arg_ctrls[i][0], arg_ctrls[i][1]);

            ctx_exit_on_error(&encoder, "Failed to control codec");
        }

        frame_avail = 1;
        got_data = 0;

        while (frame_avail || got_data)
        {
            vpx_codec_iter_t iter = NULL;
            const vpx_codec_cx_pkt_t *pkt;
            struct vpx_usec_timer timer;
            int64_t frame_start, next_frame_start;

            if (!arg_limit || frames_in < arg_limit)
            {
                frame_avail = read_frame(infile, &raw, file_type, &y4m,
                                         &detect);

                if (frame_avail)
                    frames_in++;

                fprintf(stderr,
                        "\rPass %d/%d frame %4d/%-4d %7ldB \033[K", pass + 1,
                        arg_passes, frames_in, frames_out, nbytes);
            }
            else
                frame_avail = 0;

            vpx_usec_timer_start(&timer);

            frame_start = (cfg.g_timebase.den * (int64_t)(frames_in - 1)
                          * arg_framerate.den) / cfg.g_timebase.num / arg_framerate.num;
            next_frame_start = (cfg.g_timebase.den * (int64_t)(frames_in)
                                * arg_framerate.den)
                                / cfg.g_timebase.num / arg_framerate.num;
            vpx_codec_encode(&encoder, frame_avail ? &raw : NULL, frame_start,
                             next_frame_start - frame_start,
                             0, arg_deadline);
            vpx_usec_timer_mark(&timer);
            cx_time += vpx_usec_timer_elapsed(&timer);
            ctx_exit_on_error(&encoder, "Failed to encode frame");
            got_data = 0;

            while ((pkt = vpx_codec_get_cx_data(&encoder, &iter)))
            {
                got_data = 1;

                switch (pkt->kind)
                {
                case VPX_CODEC_CX_FRAME_PKT:
                    frames_out++;
                    fprintf(stderr, " %6luF",
                            (unsigned long)pkt->data.frame.sz);

                    if(write_webm)
                    {
                        /* Update the hash */
                        if(!ebml.debug)
                            hash = murmur(pkt->data.frame.buf,
                                          pkt->data.frame.sz, hash);

                        write_webm_block(&ebml, &cfg, pkt);
                    }
                    else
                    {
                        write_ivf_frame_header(outfile, pkt);
                        if(fwrite(pkt->data.frame.buf, 1,
                                  pkt->data.frame.sz, outfile));
                    }
                    nbytes += pkt->data.raw.sz;
                    break;
                case VPX_CODEC_STATS_PKT:
                    frames_out++;
                    fprintf(stderr, " %6luS",
                           (unsigned long)pkt->data.twopass_stats.sz);
                    stats_write(&stats,
                                pkt->data.twopass_stats.buf,
                                pkt->data.twopass_stats.sz);
                    nbytes += pkt->data.raw.sz;
                    break;
                case VPX_CODEC_PSNR_PKT:

                    if (show_psnr)
                    {
                        int i;

                        psnr_sse_total += pkt->data.psnr.sse[0];
                        psnr_samples_total += pkt->data.psnr.samples[0];
                        for (i = 0; i < 4; i++)
                        {
                            fprintf(stderr, "%.3lf ", pkt->data.psnr.psnr[i]);
                            psnr_totals[i] += pkt->data.psnr.psnr[i];
                        }
                        psnr_count++;
                    }

                    break;
                default:
                    break;
                }
            }

            fflush(stdout);
        }

        fprintf(stderr,
               "\rPass %d/%d frame %4d/%-4d %7ldB %7ldb/f %7"PRId64"b/s"
               " %7lu %s (%.2f fps)\033[K", pass + 1,
               arg_passes, frames_in, frames_out, nbytes, nbytes * 8 / frames_in,
               nbytes * 8 *(int64_t)arg_framerate.num / arg_framerate.den / frames_in,
               cx_time > 9999999 ? cx_time / 1000 : cx_time,
               cx_time > 9999999 ? "ms" : "us",
               (float)frames_in * 1000000.0 / (float)cx_time);

        if ( (show_psnr) && (psnr_count>0) )
        {
            int i;
            double ovpsnr = vp8_mse2psnr(psnr_samples_total, 255.0,
                                         psnr_sse_total);

            fprintf(stderr, "\nPSNR (Overall/Avg/Y/U/V)");

            fprintf(stderr, " %.3lf", ovpsnr);
            for (i = 0; i < 4; i++)
            {
                fprintf(stderr, " %.3lf", psnr_totals[i]/psnr_count);
            }
        }

        vpx_codec_destroy(&encoder);

        fclose(infile);

        if(write_webm)
        {
            write_webm_file_footer(&ebml, hash);
        }
        else
        {
            if (!fseek(outfile, 0, SEEK_SET))
                write_ivf_file_header(outfile, &cfg, codec->fourcc, frames_out);
        }

        fclose(outfile);
        stats_close(&stats, arg_passes-1);
        fprintf(stderr, "\n");

        if (one_pass_only)
            break;
    }

    vpx_img_free(&raw);
    free(argv);
    return EXIT_SUCCESS;
}
