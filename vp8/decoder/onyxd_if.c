/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8/common/onyxc_int.h"
#if CONFIG_POSTPROC
#include "vp8/common/postproc.h"
#endif
#include "vp8/common/onyxd.h"
#include "onyxd_int.h"
#include "vpx_mem/vpx_mem.h"
#include "vp8/common/alloccommon.h"
#include "vpx_scale/yv12extend.h"
#include "vp8/common/loopfilter.h"
#include "vp8/common/swapyv12buffer.h"
#include "vp8/common/g_common.h"
#include "vp8/common/threading.h"
#include "decoderthreading.h"
#include <stdio.h>

#include "vp8/common/quant_common.h"
#include "vpx_scale/vpxscale.h"
#include "vp8/common/systemdependent.h"
#include "vpx_ports/vpx_timer.h"
#include "detokenize.h"
#if ARCH_ARM
#include "vpx_ports/arm.h"
#endif

#include "vpx_config.h"
#if CONFIG_OPENCL
#include "vp8/common/opencl/blockd_cl.h"
#include "vp8/common/opencl/vp8_opencl.h"
#endif

extern void vp8_init_loop_filter(VP8_COMMON *cm);
extern void vp8cx_init_de_quantizer(VP8D_COMP *pbi);

#define PROFILE_OUTPUT 0
#if PROFILE_OUTPUT
struct vpx_usec_timer frame_timer;
struct vpx_usec_timer loop_filter_timer;
unsigned int total_mb = 0;
unsigned int total_loop_filter = 0;
#endif

void vp8dx_initialize()
{
    static int init_done = 0;

    if (!init_done)
    {
        vp8_initialize_common();
        vp8_scale_machine_specific_config();
        init_done = 1;
    }
}


VP8D_PTR vp8dx_create_decompressor(VP8D_CONFIG *oxcf)
{
    VP8D_COMP *pbi = vpx_memalign(32, sizeof(VP8D_COMP));

    if (!pbi)
        return NULL;

    vpx_memset(pbi, 0, sizeof(VP8D_COMP));

    if (setjmp(pbi->common.error.jmp))
    {
        pbi->common.error.setjmp = 0;
        vp8dx_remove_decompressor(pbi);
        return 0;
    }

    pbi->common.error.setjmp = 1;
    vp8dx_initialize();

    vp8_create_common(&pbi->common);
    vp8_dmachine_specific_config(pbi);

    pbi->common.current_video_frame = 0;
    pbi->ready_for_new_data = 1;

#if CONFIG_MULTITHREAD
    pbi->max_threads = oxcf->max_threads;
    vp8_decoder_create_threads(pbi);
#endif

    /* vp8cx_init_de_quantizer() is first called here. Add check in frame_init_dequantizer() to avoid
     *  unnecessary calling of vp8cx_init_de_quantizer() for every frame.
     */
    vp8cx_init_de_quantizer(pbi);

    {
        VP8_COMMON *cm = &pbi->common;

        vp8_init_loop_filter(cm);
        cm->last_frame_type = KEY_FRAME;
        cm->last_filter_type = cm->filter_type;
        cm->last_sharpness_level = cm->sharpness_level;
    }

    pbi->common.error.setjmp = 0;
    return (VP8D_PTR) pbi;
}


void vp8dx_remove_decompressor(VP8D_PTR ptr)
{
    VP8D_COMP *pbi = (VP8D_COMP *) ptr;

    if (!pbi)
        return;

#if CONFIG_MULTITHREAD
    if (pbi->b_multithreaded_rd)
        vp8mt_de_alloc_temp_buffers(pbi, pbi->common.mb_rows);
    vp8_decoder_remove_threads(pbi);
#endif
    vp8_remove_common(&pbi->common);

    vpx_free(pbi);
}


int vp8dx_get_reference(VP8D_PTR ptr, VP8_REFFRAME ref_frame_flag, YV12_BUFFER_CONFIG *sd)
{
    VP8D_COMP *pbi = (VP8D_COMP *) ptr;
    VP8_COMMON *cm = &pbi->common;
    int ref_fb_idx;

    if (ref_frame_flag == VP8_LAST_FLAG)
        ref_fb_idx = cm->lst_fb_idx;
    else if (ref_frame_flag == VP8_GOLD_FLAG)
        ref_fb_idx = cm->gld_fb_idx;
    else if (ref_frame_flag == VP8_ALT_FLAG)
        ref_fb_idx = cm->alt_fb_idx;
    else
        return -1;

    vp8_yv12_copy_frame_ptr(&cm->yv12_fb[ref_fb_idx], sd);

    return 0;
}


int vp8dx_set_reference(VP8D_PTR ptr, VP8_REFFRAME ref_frame_flag, YV12_BUFFER_CONFIG *sd)
{
    VP8D_COMP *pbi = (VP8D_COMP *) ptr;
    VP8_COMMON *cm = &pbi->common;
    int ref_fb_idx;

    if (ref_frame_flag == VP8_LAST_FLAG)
        ref_fb_idx = cm->lst_fb_idx;
    else if (ref_frame_flag == VP8_GOLD_FLAG)
        ref_fb_idx = cm->gld_fb_idx;
    else if (ref_frame_flag == VP8_ALT_FLAG)
        ref_fb_idx = cm->alt_fb_idx;
    else
        return -1;

    vp8_yv12_copy_frame_ptr(sd, &cm->yv12_fb[ref_fb_idx]);

    return 0;
}

/*For ARM NEON, d8-d15 are callee-saved registers, and need to be saved by us.*/
#if HAVE_ARMV7
extern void vp8_push_neon(INT64 *store);
extern void vp8_pop_neon(INT64 *store);
#endif

static int get_free_fb (VP8_COMMON *cm)
{
    int i;
    for (i = 0; i < NUM_YV12_BUFFERS; i++)
        if (cm->fb_idx_ref_cnt[i] == 0)
            break;

    cm->fb_idx_ref_cnt[i] = 1;
    return i;
}

static void ref_cnt_fb (int *buf, int *idx, int new_idx)
{
    if (buf[*idx] > 0)
        buf[*idx]--;

    *idx = new_idx;

    buf[new_idx]++;
}

/* If any buffer copy / swapping is signalled it should be done here. */
static int swap_frame_buffers (VP8_COMMON *cm)
{
    int err = 0;

    /* The alternate reference frame or golden frame can be updated
     *  using the new, last, or golden/alt ref frame.  If it
     *  is updated using the newly decoded frame it is a refresh.
     *  An update using the last or golden/alt ref frame is a copy.
     */
    if (cm->copy_buffer_to_arf)
    {
        int new_fb = 0;

        if (cm->copy_buffer_to_arf == 1)
            new_fb = cm->lst_fb_idx;
        else if (cm->copy_buffer_to_arf == 2)
            new_fb = cm->gld_fb_idx;
        else
            err = -1;

        ref_cnt_fb (cm->fb_idx_ref_cnt, &cm->alt_fb_idx, new_fb);
    }

    if (cm->copy_buffer_to_gf)
    {
        int new_fb = 0;

        if (cm->copy_buffer_to_gf == 1)
            new_fb = cm->lst_fb_idx;
        else if (cm->copy_buffer_to_gf == 2)
            new_fb = cm->alt_fb_idx;
        else
            err = -1;

        ref_cnt_fb (cm->fb_idx_ref_cnt, &cm->gld_fb_idx, new_fb);
    }

    if (cm->refresh_golden_frame)
        ref_cnt_fb (cm->fb_idx_ref_cnt, &cm->gld_fb_idx, cm->new_fb_idx);

    if (cm->refresh_alt_ref_frame)
        ref_cnt_fb (cm->fb_idx_ref_cnt, &cm->alt_fb_idx, cm->new_fb_idx);

    if (cm->refresh_last_frame)
    {
        ref_cnt_fb (cm->fb_idx_ref_cnt, &cm->lst_fb_idx, cm->new_fb_idx);

        cm->frame_to_show = &cm->yv12_fb[cm->lst_fb_idx];
    }
    else
        cm->frame_to_show = &cm->yv12_fb[cm->new_fb_idx];

    cm->fb_idx_ref_cnt[cm->new_fb_idx]--;

    return err;
}

int vp8dx_receive_compressed_data(VP8D_PTR ptr, unsigned long size, const unsigned char *source, INT64 time_stamp)
{
#if HAVE_ARMV7
    INT64 dx_store_reg[8];
#endif
    VP8D_COMP *pbi = (VP8D_COMP *) ptr;
    VP8_COMMON *cm = &pbi->common;
    int retcode = 0;

    /*if(pbi->ready_for_new_data == 0)
        return -1;*/

    if (ptr == 0)
    {
        return -1;
    }

    pbi->common.error.error_code = VPX_CODEC_OK;

    if (size == 0)
    {
       /* This is used to signal that we are missing frames.
        * We do not know if the missing frame(s) was supposed to update
        * any of the reference buffers, but we act conservative and
        * mark only the last buffer as corrupted.
        */
        cm->yv12_fb[cm->lst_fb_idx].corrupted = 1;

        /* Signal that we have no frame to show. */
        cm->show_frame = 0;

        /* Nothing more to do. */
        return 0;
    }


#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_push_neon(dx_store_reg);
    }
#endif

    cm->new_fb_idx = get_free_fb (cm);

    if (setjmp(pbi->common.error.jmp))
    {
#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
        if (cm->rtcd.flags & HAS_NEON)
#endif
        {
            vp8_pop_neon(dx_store_reg);
        }
#endif
        pbi->common.error.setjmp = 0;

       /* We do not know if the missing frame(s) was supposed to update
        * any of the reference buffers, but we act conservative and
        * mark only the last buffer as corrupted.
        */
        cm->yv12_fb[cm->lst_fb_idx].corrupted = 1;

        if (cm->fb_idx_ref_cnt[cm->new_fb_idx] > 0)
          cm->fb_idx_ref_cnt[cm->new_fb_idx]--;
        return -1;
    }

    pbi->common.error.setjmp = 1;

    /*cm->current_video_frame++;*/
    pbi->Source = source;
    pbi->source_sz = size;

#if CONFIG_OPENCL
    pbi->mb.cl_commands = NULL;
    if (cl_initialized == CL_SUCCESS){
        int err;
        //Create command queue for macroblock.
        pbi->mb.cl_commands = clCreateCommandQueue(cl_data.context, cl_data.device_id, 0, &err);
        if (!pbi->mb.cl_commands || err != CL_SUCCESS) {
            printf("Error: Failed to create a command queue!\n");
            cl_destroy(NULL, VP8_CL_TRIED_BUT_FAILED);
        }

        pbi->mb.cl_diff_mem = NULL;
        pbi->mb.cl_predictor_mem = NULL;
        pbi->mb.cl_qcoeff_mem = NULL;
        pbi->mb.cl_dqcoeff_mem = NULL;
        pbi->mb.cl_eobs_mem = NULL;

#define SET_ON_ALLOC 0
#if SET_ON_ALLOC
        
#if ENABLE_CL_SUBPIXEL || ENABLE_CL_IDCT_DEQUANT
            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_predictor_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                    sizeof(cl_uchar)*384, pbi->mb.predictor, goto BUF_DONE, -1);
#endif

#if ENABLE_CL_IDCT_DEQUANT
            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_diff_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                    sizeof(cl_short)*400, pbi->mb.diff, goto BUF_DONE, -1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_qcoeff_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                    sizeof(cl_short)*400, pbi->mb.qcoeff, goto BUF_DONE,-1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_dqcoeff_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                    sizeof(cl_short)*400, pbi->mb.dqcoeff, goto BUF_DONE,-1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_eobs_mem, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                    sizeof(cl_char)*25, pbi->mb.eobs, goto BUF_DONE,-1);
#endif
#else
#if ENABLE_CL_IDCT_DEQUANT || ENABLE_CL_SUBPIXEL
            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_predictor_mem, CL_MEM_READ_WRITE,
                    sizeof(cl_uchar)*384, NULL, goto BUF_DONE,-1);
#endif

#if ENABLE_CL_IDCT_DEQUANT
            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_diff_mem, CL_MEM_READ_WRITE,
                    sizeof(cl_short)*400, NULL, goto BUF_DONE,-1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_qcoeff_mem, CL_MEM_READ_WRITE,
                    sizeof(cl_short)*400, NULL, goto BUF_DONE,-1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_dqcoeff_mem, CL_MEM_READ_WRITE,
                    sizeof(cl_short)*400, NULL, goto BUF_DONE,-1);

            VP8_CL_CREATE_BUF(pbi->mb.cl_commands, pbi->mb.cl_eobs_mem, CL_MEM_READ_WRITE,
                    sizeof(cl_char) * 25, NULL, goto BUF_DONE,-1);
#endif
#endif
    }
#if ENABLE_CL_IDCT_DEQUANT || ENABLE_CL_SUBPIXEL
    BUF_DONE:
#endif
#endif

#if PROFILE_OUTPUT
    printf("Frame size = %d * %d\n", cm->Height, cm->Width);
    printf("Macroblocks = %d * %d\n", cm->mb_rows, cm->mb_cols);

    vpx_usec_timer_start(&frame_timer);
#endif
    retcode = vp8_decode_frame(pbi);

#if PROFILE_OUTPUT
    vpx_usec_timer_mark(&frame_timer);
    total_mb += vpx_usec_timer_elapsed(&frame_timer);
#endif

    if (retcode < 0)
    {
#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
        if (cm->rtcd.flags & HAS_NEON)
#endif
        {
            vp8_pop_neon(dx_store_reg);
        }
#endif
        pbi->common.error.error_code = VPX_CODEC_ERROR;
        pbi->common.error.setjmp = 0;
        if (cm->fb_idx_ref_cnt[cm->new_fb_idx] > 0)
          cm->fb_idx_ref_cnt[cm->new_fb_idx]--;
        return retcode;
    }

#if CONFIG_MULTITHREAD
    if (pbi->b_multithreaded_rd && cm->multi_token_partition != ONE_PARTITION)
    {
        if (swap_frame_buffers (cm))
        {
#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
            if (cm->rtcd.flags & HAS_NEON)
#endif
            {
                vp8_pop_neon(dx_store_reg);
            }
#endif
            pbi->common.error.error_code = VPX_CODEC_ERROR;
            pbi->common.error.setjmp = 0;
            return -1;
        }
    } else
#endif
    {
        if (swap_frame_buffers (cm))
        {
#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
            if (cm->rtcd.flags & HAS_NEON)
#endif
            {
                vp8_pop_neon(dx_store_reg);
            }
#endif
            pbi->common.error.error_code = VPX_CODEC_ERROR;
            pbi->common.error.setjmp = 0;
            return -1;
        }

        if(pbi->common.filter_level)
        {

#if PROFILE_OUTPUT
            struct vpx_usec_timer lpftimer;
            vpx_usec_timer_start(&lpftimer);
#endif
           
            /* Apply the loop filter if appropriate. */
            vp8_loop_filter_frame(cm, &pbi->mb, cm->filter_level);

#if PROFILE_OUTPUT
            vpx_usec_timer_mark(&lpftimer);
            pbi->time_loop_filtering += vpx_usec_timer_elapsed(&lpftimer);

            printf("Loop Filter\n");
            total_loop_filter += vpx_usec_timer_elapsed(&lpftimer);
#if 0
            if (pbi->common.filter_type == NORMAL_LOOPFILTER){
                printf("Normal LF Time (us): %d\n", vpx_usec_timer_elapsed(&lpftimer));
            } else {
                printf("Simple LF Time (us): %d\n", vpx_usec_timer_elapsed(&lpftimer));
            }
#endif
#endif

            cm->last_frame_type = cm->frame_type;
            cm->last_filter_type = cm->filter_type;
            cm->last_sharpness_level = cm->sharpness_level;
        }
#if PROFILE_OUTPUT
        else {
            printf("No Loop Filter\n");
        }
#endif
        vp8_yv12_extend_frame_borders_ptr(cm->frame_to_show);
    }

#if CONFIG_OPENCL
    if (cl_initialized == CL_SUCCESS){
        //Copy buffer_alloc to buffer_mem so YV12_BUFFER_CONFIG can be used as
        //a reference frame (e.g. YV12..buffer_mem contains same as buffer_alloc).
        vp8_cl_mb_prep(&pbi->mb, DST_BUF);

        if (pbi->mb.cl_commands != NULL)
            clReleaseCommandQueue(pbi->mb.cl_commands);
        pbi->mb.cl_commands = NULL;
    }
#endif

    vp8_clear_system_state();

    /*vp8_print_modes_and_motion_vectors( cm->mi, cm->mb_rows,cm->mb_cols, cm->current_video_frame);*/

    if (cm->show_frame)
        cm->current_video_frame++;

    pbi->ready_for_new_data = 0;
    pbi->last_time_stamp = time_stamp;

#if 0
    {
        int i;
        INT64 earliest_time = pbi->dr[0].time_stamp;
        INT64 latest_time = pbi->dr[0].time_stamp;
        INT64 time_diff = 0;
        int bytes = 0;

        pbi->dr[pbi->common.current_video_frame&0xf].size = pbi->bc.pos + pbi->bc2.pos + 4;;
        pbi->dr[pbi->common.current_video_frame&0xf].time_stamp = time_stamp;

        for (i = 0; i < 16; i++)
        {

            bytes += pbi->dr[i].size;

            if (pbi->dr[i].time_stamp < earliest_time)
                earliest_time = pbi->dr[i].time_stamp;

            if (pbi->dr[i].time_stamp > latest_time)
                latest_time = pbi->dr[i].time_stamp;
        }

        time_diff = latest_time - earliest_time;

        if (time_diff > 0)
        {
            pbi->common.bitrate = 80000.00 * bytes / time_diff  ;
            pbi->common.framerate = 160000000.00 / time_diff ;
        }

    }
#endif

#if HAVE_ARMV7
#if CONFIG_RUNTIME_CPU_DETECT
    if (cm->rtcd.flags & HAS_NEON)
#endif
    {
        vp8_pop_neon(dx_store_reg);
    }
#endif
    pbi->common.error.setjmp = 0;


#if PROFILE_OUTPUT
    //Dump the total MB/Loop Filter processing times.
    //This is cumulative between frames, so only use the last output value.
    printf("MB Time (us): %d, LF Time (us): %d\n", total_mb, total_loop_filter);
#endif


    return retcode;
}

int vp8dx_get_raw_frame(VP8D_PTR ptr, YV12_BUFFER_CONFIG *sd, INT64 *time_stamp, INT64 *time_end_stamp, vp8_ppflags_t *flags)
{
    int ret = -1;
    VP8D_COMP *pbi = (VP8D_COMP *) ptr;

    if (pbi->ready_for_new_data == 1)
        return ret;

    /* ie no raw frame to show!!! */
    if (pbi->common.show_frame == 0)
        return ret;

    pbi->ready_for_new_data = 1;
    *time_stamp = pbi->last_time_stamp;
    *time_end_stamp = 0;

    sd->clrtype = pbi->common.clr_type;
#if CONFIG_POSTPROC
    ret = vp8_post_proc_frame(&pbi->common, sd, flags);
#else

    if (pbi->common.frame_to_show)
    {
        *sd = *pbi->common.frame_to_show;
        sd->y_width = pbi->common.Width;
        sd->y_height = pbi->common.Height;
        sd->uv_height = pbi->common.Height / 2;
        ret = 0;
    }
    else
    {
        ret = -1;
    }

#endif /*!CONFIG_POSTPROC*/
    vp8_clear_system_state();
    return ret;
}
