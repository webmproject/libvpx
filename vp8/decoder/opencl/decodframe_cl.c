/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "../onyxd_int.h"
#include "vp8/common/header.h"
#include "vp8/common/reconintra.h"
#include "vp8/common/reconintra4x4.h"
#include "vp8/common/recon.h"
#include "vp8/common/reconinter.h"
//#include "../dequantize.h"
//#include "../detokenize.h"
//#include "vp8/common/alloccommon.h"
//#include "vp8/common/entropymode.h"
//#include "vp8/common/quant_common.h"
//#include "vpx_scale/vpxscale.h"
//#include "vpx_scale/yv12extend.h"
//#include "vp8/common/setupintrarecon.h"

//#include "../decodemv.h"
//#include "vp8/common/extend.h"
//#include "vpx_mem/vpx_mem.h"
//#include "vp8/common/idct.h"
//#include "../dequantize.h"
//#include "vp8/common/predictdc.h"
//#include "vp8/common/threading.h"
//#include "../decoderthreading.h"
//#include "../dboolhuff.h"
//#include "vp8/common/blockd.h"

#include <assert.h>
#include <stdio.h>

#include "vpx_config.h"
#if CONFIG_OPENCL
#include "vp8/common/opencl/vp8_opencl.h"
#include "vp8/common/opencl/blockd_cl.h"
#include "vp8/common/opencl/reconinter_cl.h"
#include "dequantize_cl.h"
#endif

#define PROFILE_OUTPUT 0

//Implemented in ../decodframe.c
extern void mb_init_dequantizer(VP8D_COMP *pbi, MACROBLOCKD *xd);

void mb_init_dequantizer_cl(MACROBLOCKD *xd){
    int i, err;
    //Set up per-block dequant CL memory. Eventually, might be able to set up
    //one large buffer containing the entire large dequant buffer.
    if (cl_initialized == CL_SUCCESS){
        for (i=0; i < 25; i++){

#if 1 //Initialize CL memory on allocation?
            VP8_CL_CREATE_BUF(xd->cl_commands, xd->block[i].cl_dequant_mem,
                ,
                16*sizeof(cl_short),
                xd->block[i].dequant,,
            );
#else
            VP8_CL_CREATE_BUF(xd->cl_commands, xd->block[i].cl_dequant_mem,
                ,
                16*sizeof(cl_short),
                NULL,,
            );
#endif
        }
    }
}

#if CONFIG_RUNTIME_CPU_DETECT
#define RTCD_VTABLE(x) (&(pbi)->common.rtcd.x)
#else
#define RTCD_VTABLE(x) NULL
#endif

/* skip_recon_mb() is Modified: Instead of writing the result to predictor buffer and then copying it
 *  to dst buffer, we can write the result directly to dst buffer. This eliminates unnecessary copy.
 */
static void skip_recon_mb_cl(VP8D_COMP *pbi, MACROBLOCKD *xd)
{
    if (xd->frame_type == KEY_FRAME  ||  xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME)
    {

        vp8_build_intra_predictors_mbuv_s(xd);
        RECON_INVOKE(&pbi->common.rtcd.recon,
                     build_intra_predictors_mby_s)(xd);

    }
    else
    {
#if ENABLE_CL_SUBPIXEL
        if (cl_initialized == CL_SUCCESS)
        {
            vp8_build_inter_predictors_mb_s_cl(xd);
        } else
#endif
        {
            vp8_build_inter_predictors_mb_s(xd);
        }
        VP8_CL_FINISH(xd->cl_commands);
#if !ONE_CQ_PER_MB
        VP8_CL_FINISH(xd->block[0].cl_commands);
        VP8_CL_FINISH(xd->block[16].cl_commands);
        VP8_CL_FINISH(xd->block[20].cl_commands);
#endif
    }
}

void vp8_decode_macroblock_cl(VP8D_COMP *pbi, MACROBLOCKD *xd, int eobtotal)
{
    int i;

    if (xd->mode_info_context->mbmi.mode != B_PRED && xd->mode_info_context->mbmi.mode != SPLITMV && eobtotal == 0)
    {
        xd->mode_info_context->mbmi.dc_diff = 0;
        skip_recon_mb_cl(pbi, xd);
        return;
    }

    if (xd->segmentation_enabled)
        mb_init_dequantizer(pbi, xd);

    /* do prediction */
    if (xd->frame_type == KEY_FRAME  ||  xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME)
    {
        vp8_build_intra_predictors_mbuv(xd);

        if (xd->mode_info_context->mbmi.mode != B_PRED)
        {
            RECON_INVOKE(&pbi->common.rtcd.recon,
                         build_intra_predictors_mby)(xd);
        } else {
            vp8_intra_prediction_down_copy(xd);
        }
    }
    else
    {
#if ENABLE_CL_SUBPIXEL
        vp8_build_inter_predictors_mb_cl(xd);
#else
        vp8_build_inter_predictors_mb(xd);
#endif

#if !ENABLE_CL_IDCT_DEQUANT
        //Wait for inter-predict if dequant/IDCT is being done on the CPU
        VP8_CL_FINISH(xd->cl_commands);
#endif
    }

    /* dequantization and idct */
    if (xd->mode_info_context->mbmi.mode != B_PRED && xd->mode_info_context->mbmi.mode != SPLITMV)
    {
        BLOCKD *b = &xd->block[24];
        short *qcoeff = b->qcoeff_base + b->qcoeff_offset;
        vp8_second_order_fn_t second_order;

#if ENABLE_CL_IDCT_DEQUANT
        if (cl_initialized == CL_SUCCESS){
            vp8_cl_block_prep(b, DEQUANT|QCOEFF);
            vp8_dequantize_b_cl(b);
            vp8_cl_block_finish(b, DQCOEFF);
            VP8_CL_FINISH(b->cl_commands); //Keep until qcoeff memset below is CL
        }
        else
#endif
        {
            DEQUANT_INVOKE(&pbi->dequant, block)(b);
        }


        /* do 2nd order transform on the dc block */
        if (xd->eobs[24] > 1){
            second_order = IDCT_INVOKE(RTCD_VTABLE(idct), iwalsh16);
            ((int *)qcoeff)[0] = 0;
            ((int *)qcoeff)[1] = 0;
            ((int *)qcoeff)[2] = 0;
            ((int *)qcoeff)[3] = 0;
            ((int *)qcoeff)[4] = 0;
            ((int *)qcoeff)[5] = 0;
            ((int *)qcoeff)[6] = 0;
            ((int *)qcoeff)[7] = 0;
        } else {
            second_order = IDCT_INVOKE(RTCD_VTABLE(idct), iwalsh1);
            ((int *)qcoeff)[0] = 0;
        }

#if ENABLE_CL_IDCT_DEQUANT
        if (cl_initialized == CL_SUCCESS){
            int y_off = xd->dst.y_buffer - xd->dst.buffer_alloc;
            vp8_cl_block_prep(b, DQCOEFF|DIFF);

            if (xd->eobs[24] > 1)
            {
                vp8_short_inv_walsh4x4_cl(b);
            } else {
                vp8_short_inv_walsh4x4_1_cl(b);
            }
            vp8_cl_block_finish(b, DIFF);

            vp8_dequant_dc_idct_add_y_block_cl(&xd->block[0], 
                    xd->dst.buffer_alloc, xd->dst.buffer_mem, y_off, xd->dst.y_stride, xd->eobs,
                    xd->block[24].diff_offset);
        }
        else
#endif
        {
            second_order(b->dqcoeff_base + b->dqcoeff_offset, &b->diff_base[b->diff_offset]);
            DEQUANT_INVOKE (&pbi->dequant, dc_idct_add_y_block)
                (xd->qcoeff, xd->block[0].dequant,
                 xd->predictor, xd->dst.y_buffer,
                 xd->dst.y_stride, xd->eobs, &xd->block[24].diff_base[xd->block[24].diff_offset]);
        }
    }
    else if ((xd->frame_type == KEY_FRAME  ||  xd->mode_info_context->mbmi.ref_frame == INTRA_FRAME) && xd->mode_info_context->mbmi.mode == B_PRED)
    {
#if ENABLE_CL_IDCT_DEQUANT
        if (cl_initialized == CL_SUCCESS)
            vp8_cl_mb_prep(xd, DST_BUF);
#endif
        for (i = 0; i < 16; i++)
        {
            BLOCKD *b = &xd->block[i];
            short *qcoeff = b->qcoeff_base + b->qcoeff_offset;
#if ENABLE_CL_IDCT_DEQUANT
            VP8_CL_FINISH(b->cl_commands);
#endif
            vp8_predict_intra4x4(b, b->bmi.mode, b->predictor_base + b->predictor_offset);

#if ENABLE_CL_IDCT_DEQUANT
            if (cl_initialized == CL_SUCCESS){
                size_t dst_size = (4*b->dst_stride + b->dst + 4);
                cl_mem dst_mem = xd->dst.buffer_mem;

                int dst_off = *(b->base_dst) - xd->dst.buffer_alloc;

                if (xd->eobs[i] > 1)
                {
                    vp8_cl_block_prep(b, QCOEFF|DEQUANT|PREDICTOR);
                    vp8_dequant_idct_add_cl(b, *(b->base_dst), dst_mem, dst_off+b->dst, dst_size, b->qcoeff_offset, b->predictor_offset, 16, b->dst_stride, DEQUANT_INVOKE(&pbi->dequant, idct_add));
                    vp8_cl_block_finish(b, QCOEFF);
                }
                else
                {
                    vp8_cl_block_prep(b, PREDICTOR|DIFF|QCOEFF|DEQUANT);
                    vp8_dc_only_idct_add_cl(b, CL_FALSE, 0, b->qcoeff_offset, b->predictor_offset,
                        *(b->base_dst), dst_mem, dst_off+b->dst, dst_size, 16, b->dst_stride);
                    VP8_CL_FINISH(b->cl_commands);
                    ((int *)(b->qcoeff_base + b->qcoeff_offset))[0] = 0; //Move into follow-up kernel?
                }
                vp8_cl_mb_finish(xd,DST_BUF);
            }
            else
#endif
            {
                if (xd->eobs[i] > 1)
                {
                    DEQUANT_INVOKE(&pbi->dequant, idct_add)
                        (qcoeff, b->dequant,  b->predictor_base + b->predictor_offset,
                        *(b->base_dst) + b->dst, 16, b->dst_stride);
                }
                else
                {
                    IDCT_INVOKE(RTCD_VTABLE(idct), idct1_scalar_add)
                        (qcoeff[0] * b->dequant[0], b->predictor_base + b->predictor_offset,
                        *(b->base_dst) + b->dst, 16, b->dst_stride);
                    ((int *)qcoeff)[0] = 0;
                }
            }
            
        }
    }
    else
    {
#if ENABLE_CL_IDCT_DEQUANT
        if (cl_initialized == CL_SUCCESS){
            vp8_cl_mb_prep(xd,DST_BUF);
            vp8_dequant_idct_add_y_block_cl(pbi, xd);
            vp8_cl_mb_finish(xd,DST_BUF);
        }
        else
#endif
        {
            DEQUANT_INVOKE (&pbi->dequant, idct_add_y_block)
                (xd->qcoeff, xd->block[0].dequant,
                 xd->predictor, xd->dst.y_buffer,
                 xd->dst.y_stride, xd->eobs);
        }
    }

#if ENABLE_CL_IDCT_DEQUANT
    if (cl_initialized == CL_SUCCESS){
        vp8_cl_mb_prep(xd,DST_BUF);
        vp8_dequant_idct_add_uv_block_cl(pbi, xd,  DEQUANT_INVOKE (&pbi->dequant, idct_add_uv_block));
        vp8_cl_mb_finish(xd,DST_BUF);
        VP8_CL_FINISH(xd->cl_commands);
    } else
#endif
    {
    DEQUANT_INVOKE (&pbi->dequant, idct_add_uv_block)
        (xd->qcoeff+16*16, xd->block[16].dequant,
         xd->predictor+16*16, xd->dst.u_buffer, xd->dst.v_buffer,
         xd->dst.uv_stride, xd->eobs+16);
    }
}

void vp8_decode_frame_cl_finish(VP8D_COMP *pbi){

    //If using OpenCL, free all of the GPU buffers we've allocated.
    if (cl_initialized == CL_SUCCESS){
#if ENABLE_CL_IDCT_DEQUANT
        int i;
#endif

        //Wait for stuff to finish, just in case
        clFinish(pbi->mb.cl_commands);

#if !ONE_CQ_PER_MB
        clFinish(pbi->mb.block[0].cl_commands);
        clFinish(pbi->mb.block[16].cl_commands);
        clFinish(pbi->mb.block[20].cl_commands);
        clReleaseCommandQueue(pbi->mb.block[0].cl_commands);
        clReleaseCommandQueue(pbi->mb.block[16].cl_commands);
        clReleaseCommandQueue(pbi->mb.block[20].cl_commands);
#endif

#if ENABLE_CL_IDCT_DEQUANT || ENABLE_CL_SUBPIXEL
        //Free Predictor CL buffer
        if (pbi->mb.cl_predictor_mem != NULL)
            clReleaseMemObject(pbi->mb.cl_predictor_mem);
#endif

#if ENABLE_CL_IDCT_DEQUANT
        //Free other CL Block/MBlock buffers
        if (pbi->mb.cl_diff_mem != NULL)
            clReleaseMemObject(pbi->mb.cl_diff_mem);
        if (pbi->mb.cl_qcoeff_mem != NULL)
            clReleaseMemObject(pbi->mb.cl_qcoeff_mem);
        if (pbi->mb.cl_dqcoeff_mem != NULL)
            clReleaseMemObject(pbi->mb.cl_dqcoeff_mem);
        if (pbi->mb.cl_eobs_mem != NULL)
            clReleaseMemObject(pbi->mb.cl_eobs_mem);

        for (i = 0; i < 25; i++){
            clReleaseMemObject(pbi->mb.block[i].cl_dequant_mem);
            pbi->mb.block[i].cl_dequant_mem = NULL;
        }
#endif
    }
}
