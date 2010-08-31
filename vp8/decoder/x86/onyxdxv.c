/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     onyxdxv.c
*
*   Description  :     VP80 interface to DXV.
*
*****************************************************************************
*/
/****************************************************************************
*  Header Files
****************************************************************************/
#include <math.h>   // For Abs()
#include "pragmas.h"

#include "vpxdxv.h"
#include "vpxdxv_plugin.h"

#include "onyxd_int.h"
#include "onyx.h"
#include "codec_common_interface.h"
#include "vpx_scale/vpxscale.h"
#include "vpx_mem/vpx_mem.h"
#include "postproc.h"
#include "vpxblit.h"
#include "g_common.h"
#include "vpx_scale/yv12extend.h"

#include <limits.h>
#include <stdio.h>
#include "scale_mode.h"
#include "onyx_pb_interface.h"

/****************************************************************************
*  Macros
****************************************************************************/

#define VP8_FOURCC DXL_MKFOURCC( 'V', 'P', '8', '0')

extern void vp8_blit_text(const char *msg, unsigned char *address, const int pitch);


/****************************************************************************
*  Typedefs
****************************************************************************/

typedef struct  // YUV buffer configuration structure
{
    int   y_width;
    int   y_height;
    int   y_stride;

    int   uv_width;
    int   uv_height;
    int   uv_stride;

    char *y_buffer;
    char *u_buffer;
    char *v_buffer;

    char *uv_start;
    int   uv_dst_area;
    int   uv_used_area;

    unsigned char *y_ptr_scrn;
    unsigned char *u_ptr_scrn;
    unsigned char *v_ptr_scrn;


} DXV_YUV_BUFFER_CONFIG;


typedef void ((*vp8blit_func)(unsigned char *, int, YUV_BUFFER_CONFIG *));

/* define an x_image structure based on the core x_image struct */
typedef struct t_ximage_codec
{
    DXV_YUV_BUFFER_CONFIG frame_buffer;
    VP8D_COMP *my_pbi;
    VP8_COMMON *common;
    int owned;
    int decompressed_once;

    int sizeof_pixel;
    vp8blit_func blitter;

    unsigned int ppl_tag;
    unsigned int bd_tag;
    unsigned int *supported_output_format_list;

    int cpu_free;
    int postproc;
    int add_noise;
    int deinterlace;

    int post_proc2time;
    int post_proc4time;

    int hs;
    int hr;
    int vs;
    int vr;
    YV12_BUFFER_CONFIG this_buffer;
    YV12_BUFFER_CONFIG scaled_buffer;
    YV12_BUFFER_CONFIG *passed_in_buffer;

    int avgq;
    int ppcount;


} VP8_XIMAGE, *VP8_XIMAGE_HANDLE;


/****************************************************************************
*  Modul Statics
****************************************************************************/
static unsigned int g_vp8_preferred_output_format_list[] =
{
    VPXDXV_YUY2,
    VPXDXV_UYVY,
    VPXDXV_RGB8888,
    VPXDXV_RGB888,
    VPXDXV_RGB555,
    VPXDXV_RGB565,
    VPXDXV_YV12,
    VPXDXV_I420,

//    VPXDXV_YV12,
//    VPXDXV_YUY2,
//    VPXDXV_RGB565,
//    VPXDXV_UYVY,
    0
};

/****************************************************************************
*  Forward declarationss
****************************************************************************/
void onyx_set_parameter(XIMAGE_HANDLE src, int Command, unsigned int Parameter);

static int onyx_get_output_format(XIMAGE_HANDLE src, unsigned int *bd_tag);
static int onyx_set_output_format(XIMAGE_HANDLE src, unsigned int bd_tag);

static int vpx_get_size_of_pixel(unsigned int bd);

/****************************************************************************
*  Imports
****************************************************************************/

#define __Clamp255(x)   (unsigned char) ( (x) < 0 ? 0 : ( (x) <= 255 ? (x) : 255 ) )

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
void
convert_yv12_buffer_types(YV12_BUFFER_CONFIG *source, DXV_YUV_BUFFER_CONFIG *dest)
{
    dest->y_buffer = (char *)source->y_buffer;
    dest->u_buffer = (char *)source->u_buffer;
    dest->v_buffer = (char *)source->v_buffer;
    dest->y_width  = source->y_width;
    dest->y_height = source->y_height;
    dest->y_stride = source->y_stride;
    dest->uv_width  = source->uv_width;
    dest->uv_height = source->uv_height;
    dest->uv_stride = source->uv_stride;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


int onyx_blit
(
    XIMAGE_HANDLE src,
    VSCREEN_HANDLE v_screen,
    DXV_YUV_BUFFER_CONFIG *frame_buffer,
    int x,
    int y
)
{
    VP8_XIMAGE_HANDLE tab = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);
    VP8D_COMP *pbi;
    VP8_COMMON *common = tab->common;
    pbi = tab->my_pbi;

    if (v_screen) /* if there is a v_screen, blit to it */
    {
        unsigned char *ptr_scrn;
        int this_pitch, vs_height, vs_width;
        unsigned int start_tick, stop_tick;

        vpxdxv_get_vscreen_attributes(v_screen, (void **)&ptr_scrn,  &vs_width, &vs_height, &this_pitch);

        if (ptr_scrn)
        {
            int w, h;

            int p_size;
            int view_x, view_y, view_w;
            int hs, hr, vs, vr;
            int neww, newh;
            int cw, ch;
            int microseconds_available = (int)(1000000 / 30);

            microseconds_available = microseconds_available * tab->cpu_free / 100;

            if (pbi)
            {
                microseconds_available -= pbi->decode_microseconds;

                if (tab->cpu_free == 0)
                    microseconds_available = INT_MAX;

                if (tab->post_proc2time == 0)
                    tab->post_proc2time = pbi->decode_microseconds * 1 / 2;

                if (tab->post_proc4time == 0)
                    tab->post_proc4time = pbi->decode_microseconds;
            }


            if (tab->ppcount == 0)
            {
                tab->post_proc2time = 0;
                tab->post_proc4time = 0;
                tab->ppcount = 64;
            }
            else
            {
                tab->ppcount --;
            }

            vpxdxv_get_vscreen_view(v_screen, &view_x, &view_y, &view_w, NULL);

            Scale2Ratio(common->horiz_scale, &hr, &hs);
            Scale2Ratio(common->vert_scale, &vr, &vs);

            if (tab->postproc && tab->passed_in_buffer == 0)
            {
                int show_text = 0;

                unsigned char message[512];

                int pp = tab->postproc;
                int q = (tab->avgq + 4) / 8;
                int noise = 0;

                vp8_clear_system_state();

                if (pp >= 1000)
                {
                    pp -= 1000;
                    noise = pp / 100;
                    pp = pp - noise * 100;
                }

                if (pp >= 300)
                {
                    pp -= 300;
                    show_text = 3;
                }
                else if (pp >= 200)
                {
                    pp -= 200;
                    show_text = 2;
                }
                else if (pp >= 100)
                {
                    pp -= 100;
                    show_text = 1;
                }

                if (pbi && (pbi->mb.segmentation_enabled & SEGMENT_PF) && tab->deinterlace)
                {
                    de_interlace(common->frame_to_show->y_buffer, common->post_proc_buffer.y_buffer,
                                 common->post_proc_buffer.y_width, common->post_proc_buffer.y_height,
                                 common->post_proc_buffer.y_stride);

                    de_interlace(common->frame_to_show->u_buffer, common->post_proc_buffer.u_buffer,
                                 common->post_proc_buffer.uv_width, common->post_proc_buffer.uv_height,
                                 common->post_proc_buffer.uv_stride);
                    de_interlace(common->frame_to_show->v_buffer, common->post_proc_buffer.v_buffer,
                                 common->post_proc_buffer.uv_width, common->post_proc_buffer.uv_height,
                                 common->post_proc_buffer.uv_stride);
                }
                else
                {
                    if (pp >= 10 && pp <= 20)
                    {
                        q = q + (pp - 15) * 10;

                        if (q < 0)
                            q = 0;
                    }

                    start_tick = vp8_get_high_res_timer_tick();

                    if (pp > 3 && tab->post_proc4time < microseconds_available)
                    {
                        vp8_deblock_and_de_macro_block(common->frame_to_show, &common->post_proc_buffer, q, 1, 0);

                        stop_tick = vp8_get_high_res_timer_tick();

                        if (pbi)
                            tab->post_proc4time = vp8_get_time_in_micro_sec(start_tick, stop_tick);
                    }

                    else if (pp > 0 && tab->post_proc2time < microseconds_available)
                    {
                        vp8_deblock(common->frame_to_show, &common->post_proc_buffer, q , 1,  0);
                        stop_tick = vp8_get_high_res_timer_tick();

                        if (pbi)
                            tab->post_proc2time = vp8_get_time_in_micro_sec(start_tick, stop_tick);
                    }
                    else
                    {
                        vp8_yv12_copy_frame(common->frame_to_show, &common->post_proc_buffer);
                    }

                }

                vp8_clear_system_state();

                if (tab->add_noise == 1)
                {

                    vp8_plane_add_noise(common->post_proc_buffer.y_buffer,
                                        common->post_proc_buffer.y_width, common->post_proc_buffer.y_height,
                                        common->post_proc_buffer.y_stride, 63 - q, noise);
                }


                if (show_text == 1)
                {
#ifdef PACKET_TESTING
                    {
                        VP8_HEADER *oh2 = (VP8_HEADER *) pbi->Source;
                        sprintf(message, "%8d %d%d%d%d%d size:%d\n",
                        oh2->frame_number ,
                        oh2->update_gold  ,
                        oh2->update_last  ,
                        oh2->uses_gold    ,
                        oh2->uses_last    ,
                        oh2->type,
                        vpxdxv_get_ximage_csize(src));
                    }
#else
                    sprintf(message, "F:%1ldG:%1ldQ:%3ldF:%3ld,%3ldP:%d_s:%6ld,N:%d,",
                            (common->frame_type == KEY_FRAME),
                            common->refresh_golden_frame,
                            common->base_qindex,
                            common->filter_level,
                            q,
                            tab->postproc,
                            vpxdxv_get_ximage_csize(src), noise);
#endif

                    vp8_blit_text(message, common->post_proc_buffer.y_buffer, common->post_proc_buffer.y_stride);

                }
                else if (show_text == 2)
                {
                    int i, j;
                    unsigned char *y_ptr;
                    YV12_BUFFER_CONFIG *post = &common->post_proc_buffer;
                    int mb_rows = post->y_height >> 4;
                    int mb_cols = post->y_width  >> 4;
                    int mb_index = 0;
                    MODE_INFO *mi = common->mi;

                    y_ptr = post->y_buffer + 4 * post->y_stride + 4;

                    // vp8_filter each macro block
                    for (i = 0; i < mb_rows; i++)
                    {
                        for (j = 0; j < mb_cols; j++)
                        {
                            char zz[4];

                            if (pp == 4)
                                sprintf(zz, "%c", mi[mb_index].mbmi.mode + 'a');
                            else
                                sprintf(zz, "%c", mi[mb_index].mbmi.ref_frame + 'a');

                            vp8_blit_text(zz, y_ptr, post->y_stride);
                            mb_index ++;
                            y_ptr += 16;
                        }

                        mb_index ++; //border
                        y_ptr += post->y_stride  * 16 - post->y_width;

                    }
                }
                else if (show_text == 3)
                {
                    int i, j;
                    unsigned char *y_ptr;
                    YV12_BUFFER_CONFIG *post = &common->post_proc_buffer;
                    int mb_rows = post->y_height >> 4;
                    int mb_cols = post->y_width  >> 4;
                    int mb_index = 0;
                    MODE_INFO *mi = common->mi;

                    y_ptr = post->y_buffer + 4 * post->y_stride + 4;

                    // vp8_filter each macro block
                    for (i = 0; i < mb_rows; i++)
                    {
                        for (j = 0; j < mb_cols; j++)
                        {
                            char zz[4];

                            if (j == 0)
                                sprintf(zz, "%c", '0' + i % 10);
                            else
                                sprintf(zz, "%c", '0' + j % 10);

                            vp8_blit_text(zz, y_ptr, post->y_stride);
                            mb_index ++;
                            y_ptr += 16;
                        }

                        y_ptr += post->y_stride  * 16 - post->y_width;

                    }
                }

                vpx_memcpy(&tab->this_buffer, &common->post_proc_buffer, sizeof(YV12_BUFFER_CONFIG));
            }
            else
            {
                vpx_memcpy(&tab->this_buffer, common->frame_to_show, sizeof(YV12_BUFFER_CONFIG));
            }


            /* get a frame pointer to the scaled and postprocessed reconstructed buffer */
            if (tab->passed_in_buffer == 0)
            {
                if (common->horiz_scale != NORMAL || common->vert_scale != NORMAL)
                {
                    neww = hs * tab->this_buffer.y_width / hr;
                    newh = vs * tab->this_buffer.y_height / vr;

                    neww += neww & 1;

                    if (tab->hs != hs || tab->hr != hr || tab->vs != vs || tab->vr != vr)
                    {
                        vp8_yv12_alloc_frame_buffer(&tab->scaled_buffer, neww, newh , 8);
                    }

                    vp8_yv12_scale_or_center(&tab->this_buffer,
                                             &tab->scaled_buffer,
                                             neww, newh, SCALE_TO_FIT, hs, hr, vs, vr);

                    convert_yv12_buffer_types(&tab->scaled_buffer, frame_buffer);

                    cw = hs * common->Width / hr;
                    ch = vs * common->Height / vr;

                }
                else
                {
                    convert_yv12_buffer_types(&tab->this_buffer, frame_buffer);

                    cw = common->Width;
                    ch = common->Height;
                }
            }
            else
            {
                convert_yv12_buffer_types(tab->passed_in_buffer, frame_buffer);
                cw = common->Width;
                ch = common->Height;
                tab->passed_in_buffer = 0;
            }

            frame_buffer->y_width = cw;
            frame_buffer->y_height = ch;
            frame_buffer->uv_width = cw / 2;
            frame_buffer->uv_height = ch / 2;

            p_size = vpx_get_size_of_pixel(tab->bd_tag);

            /* remember to offset if requested */
            y += view_y;
            x += view_x ;

            /* for planar destinations */
            w = view_w;
            h = vs_height;

            if (w < frame_buffer->y_width)
            {
                frame_buffer->y_width = w;
                frame_buffer->uv_width = (w + 1) / 2;
            }

            if (h < frame_buffer->y_height)
            {
                frame_buffer->y_height = h;
                frame_buffer->uv_height = (h + 1) / 2;
            }

            if (frame_buffer->y_width < view_w)
                x += (view_w - frame_buffer->y_width) / 2;

            if (x & 1)
                x -= 1;

            if (frame_buffer->y_height < vs_height)
                y += (vs_height - frame_buffer->y_height) / 2;


            ptr_scrn += (x * p_size) + (y * this_pitch);

            frame_buffer->y_stride *= -1;
            frame_buffer->uv_stride *= -1;

            if (tab->bd_tag == VPXDXV_YV12 || tab->bd_tag == VPXDXV_I420)
            {
                if (this_pitch < 0)
                {
                    frame_buffer->uv_start = (char *)(ptr_scrn + abs(this_pitch) + abs(this_pitch) * h / 4 + this_pitch / 2);
                    frame_buffer->uv_dst_area = abs((this_pitch * h) / 4);
                    frame_buffer->uv_used_area = 0;
                }
                else
                {
                    frame_buffer->uv_start = (char *)(ptr_scrn + (this_pitch * h));
                    frame_buffer->uv_dst_area = (((this_pitch + 1) / 2) * ((h + 1) / 2));
                    frame_buffer->uv_used_area = (((this_pitch + 1) / 2) * frame_buffer->uv_height);
                }
            }

            if ((pbi->mb.segmentation_enabled & SEGMENT_PF) && (tab->bd_tag != VPXDXV_YV12 && tab->bd_tag != VPXDXV_I420))
            {
                int ypitch = frame_buffer->y_stride;
                int uvpitch = frame_buffer->uv_stride;

                frame_buffer->y_stride <<= 1;
                frame_buffer->y_height >>= 1;
                frame_buffer->uv_stride <<= 1;
                frame_buffer->uv_height >>= 1;

                ptr_scrn += this_pitch;
                frame_buffer->y_buffer -= ypitch;
                frame_buffer->u_buffer -= uvpitch;
                frame_buffer->v_buffer -= uvpitch;
                tab->blitter(ptr_scrn, 2 * this_pitch, (YUV_BUFFER_CONFIG *)(&tab->frame_buffer));

                ptr_scrn -= this_pitch;
                frame_buffer->y_buffer += ypitch;
                frame_buffer->u_buffer += uvpitch;
                frame_buffer->v_buffer += uvpitch;
                tab->blitter(ptr_scrn, 2 * this_pitch, (YUV_BUFFER_CONFIG *)(&tab->frame_buffer));

            }
            else
            {
                /* blit the screen */
                tab->blitter(ptr_scrn, this_pitch, (YUV_BUFFER_CONFIG *)(&tab->frame_buffer));
                vpx_log("Decoder: Frame shown \n");
            }

        }
        else
            vpx_log("Decoder: Frame not shown scrn pointer 0\n");
    }
    else
        vpx_log("Decoder: Frame not shown vscreen 0\n");

    return DXV_OK;
}
/****************************************************************************
 *
 *  ROUTINE       :     onyx_decompress
 *
 *  INPUTS        :     None
 *
 *  OUTPUTS       :     None
 *
 *  RETURNS       :     None.
 *
 *  FUNCTION      :
 *
 *  SPECIAL NOTES :
 *
 ****************************************************************************/
static
int onyx_decompress(XIMAGE_HANDLE src, VSCREEN_HANDLE v_screen)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);
    unsigned char *c_addr;
    unsigned int c_size;
    int w, h, x, y;
    int vp8_rv;

    c_addr = vpxdxv_get_ximage_cdata_addr(src);
    c_size = vpxdxv_get_ximage_csize(src);
    vpxdxv_get_ximage_xywh(src, &x, &y, &w, &h);

    // if we have a compressed frame decompress it ( otherwise we'll just redo
    // the scaling and postprocessing from the last frame )
    if (c_addr)
    {
        if (c_size != 0)
        {
            int flags;
            int ret_val;

            int f;

            // decode the frame
            ret_val = vp8d_decompress_frame((VP8D_PTR) this_algorithm_base->my_pbi,
                                            c_size,
                                            (char *) c_addr,
                                            &this_algorithm_base->this_buffer,
                                            &flags);


            f = this_algorithm_base->my_pbi->common.filter_level * 10 / 6;

            if (this_algorithm_base->my_pbi->common.frame_type == KEY_FRAME)
                this_algorithm_base->avgq = 8 * f;
            else
                this_algorithm_base->avgq = this_algorithm_base->avgq * 7 / 8 + f;



            if (ret_val != 0)
            {
                if (ret_val == -1)
                    return DXV_VERSION_CONFLICT;
                else
                    return DXV_BAD_DATA;
            }

        }
    }


    vp8_rv = onyx_blit(src, v_screen, &this_algorithm_base->frame_buffer, x, y);


    return vp8_rv;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static
int vp8_ximagedestroy(XIMAGE_HANDLE src)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    if (this_algorithm_base)
    {

        vp8_yv12_de_alloc_frame_buffer(&this_algorithm_base->scaled_buffer);

        /* safety check in case stopdecode was not called */
        if (this_algorithm_base->owned)
            vp8dx_remove_decompressor((VP8D_PTR)(this_algorithm_base->my_pbi));

        duck_free(this_algorithm_base);
    }

    return DXV_OK;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static int
onyx_get_post_proc(XIMAGE_HANDLE src, unsigned int *ppl)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    if (this_algorithm_base)
    {
        *ppl = this_algorithm_base->ppl_tag;

        return DXV_OK;
    }

    return DXV_NULL_BASE;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static int
onyx_set_post_proc(XIMAGE_HANDLE src, unsigned int ppl)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    if (this_algorithm_base)
    {
        this_algorithm_base->ppl_tag = ppl;

        return DXV_OK;
    }

    return DXV_NULL_BASE;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static
int vp8_ximagestop_decode(XIMAGE_HANDLE src)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    if (this_algorithm_base)
    {

        vp8_yv12_de_alloc_frame_buffer(&this_algorithm_base->scaled_buffer);

        if (this_algorithm_base->owned)
            vp8dx_remove_decompressor((VP8D_PTR)(this_algorithm_base->my_pbi));

        this_algorithm_base->owned = 0;
    }

    return DXV_OK;
}


/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static
int vp8_ximagestart_decode
(
    XIMAGE_HANDLE src
)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);
    XIMAGE_INFO_PTR xinfo = vpxdxv_get_ximage_info(src);
    VP8D_CONFIG ocf;

    if (xinfo)
    {
        ocf.Width = xinfo->width;
        ocf.Height = xinfo->height;
    }

    if (this_algorithm_base->common == 0)
    {
        this_algorithm_base->my_pbi = (VP8D_COMP *) vp8dx_create_decompressor(&ocf);
        this_algorithm_base->owned = 1;
        this_algorithm_base->common = &this_algorithm_base->my_pbi->common;
        this_algorithm_base->avgq = 0;

    }

    this_algorithm_base->passed_in_buffer = 0;
    this_algorithm_base->post_proc2time = 0;
    this_algorithm_base->post_proc4time = 0;
    this_algorithm_base->ppcount = 64;

    return DXV_OK;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static
DXV_HANDLE vp8_ximagecreate(XIMAGE_HANDLE src)
{
    VP8_XIMAGE_HANDLE this_algorithm_base;

    /* create a new algorithm base container */
    this_algorithm_base = (VP8_XIMAGE_HANDLE)duck_calloc(1, sizeof(VP8_XIMAGE), DMEM_GENERAL);

    if (this_algorithm_base == NULL)
        return NULL;

    vp8_scale_machine_specific_config();

    vpxdxv_register_ximage_start_decode(src, vp8_ximagestart_decode);

    vpxdxv_register_ximage_stop_decode(src, vp8_ximagestop_decode);

    vpxdxv_register_ximage_destroy(src, vp8_ximagedestroy);

    vpxdxv_register_ximage_dx(src, onyx_decompress);

    vpxdxv_register_ximage_set_parameter(src, onyx_set_parameter);

    vpxdxv_register_ximage_output_format_func(src,
            onyx_get_output_format,
            onyx_set_output_format);

    vpxdxv_register_ximage_post_proc_level_func(src,
            onyx_get_post_proc,
            onyx_set_post_proc);

    return (DXV_HANDLE)this_algorithm_base;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

static int store_output_list(unsigned int supported, int count,
                             unsigned int *outlist)
{
    int i = 0, j = 0,
        ret = DXV_OK;

    while (i < count)
    {
        while (supported && !(supported & 0x01))
        {
            supported >>= 1;
            ++j;
        }

        *(outlist + i) = g_vp8_preferred_output_format_list[j];
        ++i;
        ++j;
        supported >>= 1;
    }


    return ret;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static int onyx_get_output_list(XIMAGE_INFO_PTR xinfo, unsigned int *outlist,
                                unsigned int *size)
{
    int i,
        ret = DXV_INVALID_REQUEST;
    unsigned int supported = 0,
                 count = 0;
    (void)xinfo;

    if (size)
    {
        for (i = 0; i < sizeof(g_vp8_preferred_output_format_list) / sizeof(unsigned int) && i < 32; ++i)
        {
            if (vpx_get_blitter(g_vp8_preferred_output_format_list[i]) != (void *)0xffffffff)
            {
                supported |= (1 << i);
                ++count;
            }
        }

        if (outlist)
        {
            if (count && ((count + 1) == (*size / sizeof(int))))
                ret = store_output_list(supported, count, outlist);
            else
                *outlist = 0;
        }
        else
        {
            *size = (count + 1) * sizeof(int);
            ret = DXV_OK;
        }
    }

    return ret;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
int onyx_init(void)
{
    int vp8_rv;

    /* register VPX blitters based on cpu */
    vpx_set_blit();

    vp8_rv = vpxdxv_register_ximage(vp8_ximagecreate, onyx_get_output_list, VP8_FOURCC);
    return vp8_rv;

    return DXV_OK;
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
int onyx_exit(void)
{

    vpxdxv_un_register_ximage(VP8_FOURCC);

    return DXV_OK;
}
/****************************************************************************
 *
 *  ROUTINE       :  onyx_set_parameter
 *
 *  INPUTS        :  XIMAGE_HANDLE src   :
 *                   int Command             :
 *                   unsigned long Parameter :
 *
 *  OUTPUTS       :  None.
 *
 *  RETURNS       :  void
 *
 *  FUNCTION      :
 *
 *
 *  SPECIAL NOTES :  None.
 *
 ****************************************************************************/
void onyx_set_parameter(XIMAGE_HANDLE src, int Command, unsigned int Parameter)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    switch (Command)
    {
    case PBC_SET_CPUFREE:
        this_algorithm_base->cpu_free  = Parameter;
        break;
    case PBC_SET_POSTPROC:
        this_algorithm_base->postproc = Parameter;
        break;

    case PBC_SET_BLITBUFF:
        this_algorithm_base->passed_in_buffer = (YV12_BUFFER_CONFIG *) Parameter;
        break;

    case PBC_SET_REFERENCEFRAME:
    {
        VP8_XIMAGE_HANDLE tab = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);
        VP8D_COMP *pbi;
        pbi = tab->my_pbi;
        vp8_yv12_copy_frame((YV12_BUFFER_CONFIG *) Parameter, &pbi->common.last_frame);
    }
    break;

    case PBC_SET_COMMON:

        if (Parameter)
        {
            this_algorithm_base->common = (VP8_COMMON *)Parameter;
        }

        break;
    case PBC_SET_ADDNOISE:
        this_algorithm_base->add_noise = Parameter;
        break;
    case PBC_SET_DEINTERLACEMODE:
        this_algorithm_base->deinterlace = Parameter;
        break;

    }
}
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static int
onyx_get_output_format(XIMAGE_HANDLE src, unsigned int *format_tag)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);

    if (this_algorithm_base)
    {
        *format_tag = this_algorithm_base->bd_tag;
        return DXV_OK;
    }

    return DXV_NULL_BASE;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
static int
onyx_set_output_format(XIMAGE_HANDLE src, unsigned int bd_tag)
{
    VP8_XIMAGE_HANDLE this_algorithm_base = (VP8_XIMAGE_HANDLE)vpxdxv_get_algorithm_base_ptr(src);
    int i;
    unsigned int bd_tag_found;

    if (this_algorithm_base)
    {
        i = 0;
        bd_tag_found = 0;

        while (g_vp8_preferred_output_format_list[i] != 0)
        {
            if (g_vp8_preferred_output_format_list[i] == bd_tag)
            {
                bd_tag_found = 1;
                break;
            }

            i++;
        }

        if (bd_tag_found)
        {
            this_algorithm_base->blitter = (vp8blit_func)vpx_get_blitter(bd_tag);
            this_algorithm_base->bd_tag = bd_tag;
            return DXV_OK;
        }

        return DXV_INVALID_BLIT;
    }

    return DXV_NULL_BASE;
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
int
vpx_get_size_of_pixel(unsigned int bd)
{
    int vp8_rv;

    switch (bd)
    {
    case VPXDXV_YV12:
    case VPXDXV_I420:
        vp8_rv = 1;
        break;

#ifdef _ENABLE_SPLIT_PIXEL_
    case VPXDXV_SPLIT565:
#endif
    case VPXDXV_RGB555:
    case VPXDXV_RGB565:
    case VPXDXV_YUY2:
    case VPXDXV_UYVY:
    case VPXDXV_YVYU:
        vp8_rv = 2;
        break;

    case VPXDXV_RGB888:
        vp8_rv = 3;
        break;

    case VPXDXV_RGB8888:
        vp8_rv = 4;
        break;

    default:
        vp8_rv = -1;
        break;
    }

    return vp8_rv;
}
