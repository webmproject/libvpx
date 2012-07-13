/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vpx/vpx_integer.h"
#include "recon.h"
#include "subpixel.h"
#include "blockd.h"
#include "reconinter.h"
#if CONFIG_RUNTIME_CPU_DETECT
#include "onyxc_int.h"
#endif

void vp8_copy_mem16x16_c(
  unsigned char *src,
  int src_stride,
  unsigned char *dst,
  int dst_stride) {

  int r;

  for (r = 0; r < 16; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
    dst[8] = src[8];
    dst[9] = src[9];
    dst[10] = src[10];
    dst[11] = src[11];
    dst[12] = src[12];
    dst[13] = src[13];
    dst[14] = src[14];
    dst[15] = src[15];

#else
    ((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((uint32_t *)src)[1];
    ((uint32_t *)dst)[2] = ((uint32_t *)src)[2];
    ((uint32_t *)dst)[3] = ((uint32_t *)src)[3];

#endif
    src += src_stride;
    dst += dst_stride;

  }

}

void vp8_avg_mem16x16_c(
  unsigned char *src,
  int src_stride,
  unsigned char *dst,
  int dst_stride) {
  int r;

  for (r = 0; r < 16; r++) {
    int n;

    for (n = 0; n < 16; n++) {
      dst[n] = (dst[n] + src[n] + 1) >> 1;
    }

    src += src_stride;
    dst += dst_stride;
  }
}

void vp8_copy_mem8x8_c(
  unsigned char *src,
  int src_stride,
  unsigned char *dst,
  int dst_stride) {
  int r;

  for (r = 0; r < 8; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
#else
    ((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((uint32_t *)src)[1];
#endif
    src += src_stride;
    dst += dst_stride;

  }

}

void vp8_avg_mem8x8_c(
  unsigned char *src,
  int src_stride,
  unsigned char *dst,
  int dst_stride) {
  int r;

  for (r = 0; r < 8; r++) {
    int n;

    for (n = 0; n < 8; n++) {
      dst[n] = (dst[n] + src[n] + 1) >> 1;
    }

    src += src_stride;
    dst += dst_stride;
  }
}

void vp8_copy_mem8x4_c(
  unsigned char *src,
  int src_stride,
  unsigned char *dst,
  int dst_stride) {
  int r;

  for (r = 0; r < 4; r++) {
#if !(CONFIG_FAST_UNALIGNED)
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
    dst[7] = src[7];
#else
    ((uint32_t *)dst)[0] = ((uint32_t *)src)[0];
    ((uint32_t *)dst)[1] = ((uint32_t *)src)[1];
#endif
    src += src_stride;
    dst += dst_stride;

  }

}



void vp8_build_inter_predictors_b(BLOCKD *d, int pitch, vp8_subpix_fn_t sppf) {
  int r;
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;

  if (mv.as_mv.row & 7 || mv.as_mv.col & 7) {
    ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);
#if CONFIG_SIXTEENTH_SUBPEL_UV
    sppf(ptr, d->pre_stride, (mv.as_mv.col & 7) << 1, (mv.as_mv.row & 7) << 1, pred_ptr, pitch);
#else
    sppf(ptr, d->pre_stride, mv.as_mv.col & 7, mv.as_mv.row & 7, pred_ptr, pitch);
#endif
  } else {
    ptr_base += d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);
    ptr = ptr_base;

    for (r = 0; r < 4; r++) {
#if !(CONFIG_FAST_UNALIGNED)
      pred_ptr[0]  = ptr[0];
      pred_ptr[1]  = ptr[1];
      pred_ptr[2]  = ptr[2];
      pred_ptr[3]  = ptr[3];
#else
      *(uint32_t *)pred_ptr = *(uint32_t *)ptr;
#endif
      pred_ptr     += pitch;
      ptr         += d->pre_stride;
    }
  }
}

/*
 * Similar to vp8_build_inter_predictors_b(), but instead of storing the
 * results in d->predictor, we average the contents of d->predictor (which
 * come from an earlier call to vp8_build_inter_predictors_b()) with the
 * predictor of the second reference frame / motion vector.
 */
void vp8_build_2nd_inter_predictors_b(BLOCKD *d, int pitch, vp8_subpix_fn_t sppf) {
  int r;
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_second_pre);
  mv.as_int = d->bmi.as_mv.second.as_int;

  if (mv.as_mv.row & 7 || mv.as_mv.col & 7) {
    ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);
#if CONFIG_SIXTEENTH_SUBPEL_UV
    sppf(ptr, d->pre_stride, (mv.as_mv.col & 7) << 1, (mv.as_mv.row & 7) << 1, pred_ptr, pitch);
#else
    sppf(ptr, d->pre_stride, mv.as_mv.col & 7, mv.as_mv.row & 7, pred_ptr, pitch);
#endif
  } else {
    ptr_base += d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);
    ptr = ptr_base;

    for (r = 0; r < 4; r++) {
      pred_ptr[0]  = (pred_ptr[0] + ptr[0] + 1) >> 1;
      pred_ptr[1]  = (pred_ptr[1] + ptr[1] + 1) >> 1;
      pred_ptr[2]  = (pred_ptr[2] + ptr[2] + 1) >> 1;
      pred_ptr[3]  = (pred_ptr[3] + ptr[3] + 1) >> 1;
      pred_ptr    += pitch;
      ptr         += d->pre_stride;
    }
  }
}

static void build_inter_predictors4b(MACROBLOCKD *x, BLOCKD *d, int pitch) {
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);

  if (mv.as_mv.row & 7 || mv.as_mv.col & 7) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
    x->subpixel_predict8x8(ptr, d->pre_stride, (mv.as_mv.col & 7) << 1, (mv.as_mv.row & 7) << 1, pred_ptr, pitch);
#else
    x->subpixel_predict8x8(ptr, d->pre_stride, mv.as_mv.col & 7, mv.as_mv.row & 7, pred_ptr, pitch);
#endif
  } else {
    RECON_INVOKE(&x->rtcd->recon, copy8x8)(ptr, d->pre_stride, pred_ptr, pitch);
  }
}

/*
 * Similar to build_inter_predictors_4b(), but instead of storing the
 * results in d->predictor, we average the contents of d->predictor (which
 * come from an earlier call to build_inter_predictors_4b()) with the
 * predictor of the second reference frame / motion vector.
 */
static void build_2nd_inter_predictors4b(MACROBLOCKD *x, BLOCKD *d, int pitch) {
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_second_pre);
  mv.as_int = d->bmi.as_mv.second.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);

  if (mv.as_mv.row & 7 || mv.as_mv.col & 7) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
    x->subpixel_predict_avg8x8(ptr, d->pre_stride, (mv.as_mv.col & 7) << 1, (mv.as_mv.row & 7) << 1, pred_ptr, pitch);
#else
    x->subpixel_predict_avg8x8(ptr, d->pre_stride, mv.as_mv.col & 7, mv.as_mv.row & 7, pred_ptr, pitch);
#endif
  } else {
    RECON_INVOKE(&x->rtcd->recon, avg8x8)(ptr, d->pre_stride, pred_ptr, pitch);
  }
}

static void build_inter_predictors2b(MACROBLOCKD *x, BLOCKD *d, int pitch) {
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = d->predictor;
  int_mv mv;

  ptr_base = *(d->base_pre);
  mv.as_int = d->bmi.as_mv.first.as_int;
  ptr = ptr_base + d->pre + (mv.as_mv.row >> 3) * d->pre_stride + (mv.as_mv.col >> 3);

  if (mv.as_mv.row & 7 || mv.as_mv.col & 7) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
    x->subpixel_predict8x4(ptr, d->pre_stride, (mv.as_mv.col & 7) << 1, (mv.as_mv.row & 7) << 1, pred_ptr, pitch);
#else
    x->subpixel_predict8x4(ptr, d->pre_stride, mv.as_mv.col & 7, mv.as_mv.row & 7, pred_ptr, pitch);
#endif
  } else {
    RECON_INVOKE(&x->rtcd->recon, copy8x4)(ptr, d->pre_stride, pred_ptr, pitch);
  }
}


/*encoder only*/
#if CONFIG_PRED_FILTER

// Select the thresholded or non-thresholded filter
#define USE_THRESH_FILTER 0

#define PRED_FILT_LEN 5

static const int filt_shift = 4;
static const int pred_filter[PRED_FILT_LEN] = {1, 2, 10, 2, 1};
// Alternative filter {1, 1, 4, 1, 1}

#if !USE_THRESH_FILTER
void filter_mb(unsigned char *src, int src_stride,
               unsigned char *dst, int dst_stride,
               int width, int height) {
  int i, j, k;
  unsigned int Temp[32 * 32];
  unsigned int  *pTmp = Temp;
  unsigned char *pSrc = src - (1 + src_stride) * (PRED_FILT_LEN / 2);

  // Horizontal
  for (i = 0; i < height + PRED_FILT_LEN - 1; i++) {
    for (j = 0; j < width; j++) {
      int sum = 0;
      for (k = 0; k < PRED_FILT_LEN; k++)
        sum += pSrc[j + k] * pred_filter[k];
      pTmp[j] = sum;
    }

    pSrc += src_stride;
    pTmp += width;
  }

  // Vertical
  pTmp = Temp;
  for (i = 0; i < width; i++) {
    unsigned char *pDst = dst + i;
    for (j = 0; j < height; j++) {
      int sum = 0;
      for (k = 0; k < PRED_FILT_LEN; k++)
        sum += pTmp[(j + k) * width] * pred_filter[k];
      // Round
      sum = (sum + ((1 << (filt_shift << 1)) >> 1)) >> (filt_shift << 1);
      pDst[j * dst_stride] = (sum < 0 ? 0 : sum > 255 ? 255 : sum);
    }
    ++pTmp;
  }
}
#else
// Based on vp8_post_proc_down_and_across_c (postproc.c)
void filter_mb(unsigned char *src, int src_stride,
               unsigned char *dst, int dst_stride,
               int width, int height) {
  unsigned char *pSrc, *pDst;
  int row;
  int col;
  int i;
  int v;
  unsigned char d[8];

  /* TODO flimit should be linked to the quantizer value */
  int flimit = 7;

  for (row = 0; row < height; row++) {
    /* post_proc_down for one row */
    pSrc = src;
    pDst = dst;

    for (col = 0; col < width; col++) {
      int kernel = (1 << (filt_shift - 1));
      int v = pSrc[col];

      for (i = -2; i <= 2; i++) {
        if (abs(v - pSrc[col + i * src_stride]) > flimit)
          goto down_skip_convolve;

        kernel += pred_filter[2 + i] * pSrc[col + i * src_stride];
      }

      v = (kernel >> filt_shift);
    down_skip_convolve:
      pDst[col] = v;
    }

    /* now post_proc_across */
    pSrc = dst;
    pDst = dst;

    for (i = 0; i < 8; i++)
      d[i] = pSrc[i];

    for (col = 0; col < width; col++) {
      int kernel = (1 << (filt_shift - 1));
      v = pSrc[col];

      d[col & 7] = v;

      for (i = -2; i <= 2; i++) {
        if (abs(v - pSrc[col + i]) > flimit)
          goto across_skip_convolve;

        kernel += pred_filter[2 + i] * pSrc[col + i];
      }

      d[col & 7] = (kernel >> filt_shift);
    across_skip_convolve:

      if (col >= 2)
        pDst[col - 2] = d[(col - 2) & 7];
    }

    /* handle the last two pixels */
    pDst[col - 2] = d[(col - 2) & 7];
    pDst[col - 1] = d[(col - 1) & 7];

    /* next row */
    src += src_stride;
    dst += dst_stride;
  }
}
#endif  // !USE_THRESH_FILTER

#endif  // CONFIG_PRED_FILTER

void vp8_build_inter16x16_predictors_mbuv(MACROBLOCKD *xd) {
  unsigned char *uptr, *vptr;
  unsigned char *upred_ptr = &xd->predictor[256];
  unsigned char *vpred_ptr = &xd->predictor[320];

  int omv_row = xd->mode_info_context->mbmi.mv.as_mv.row;
  int omv_col = xd->mode_info_context->mbmi.mv.as_mv.col;
  int mv_row  = omv_row;
  int mv_col  = omv_col;
  int offset;
  int pre_stride = xd->block[16].pre_stride;

  /* calc uv motion vectors */
  if (mv_row < 0)
    mv_row -= 1;
  else
    mv_row += 1;

  if (mv_col < 0)
    mv_col -= 1;
  else
    mv_col += 1;

  mv_row /= 2;
  mv_col /= 2;

  mv_row &= xd->fullpixel_mask;
  mv_col &= xd->fullpixel_mask;

  offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
  uptr = xd->pre.u_buffer + offset;
  vptr = xd->pre.v_buffer + offset;

#if CONFIG_PRED_FILTER
  if (xd->mode_info_context->mbmi.pred_filter_enabled) {
    int i;
#if CONFIG_ENHANCED_INTERP
    int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
    int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
    int len = 7 + (Interp_Extend << 1);
    unsigned char Temp[32 * 32]; // Input data required by sub-pel filter
    unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);
    unsigned char *pSrc = uptr;
    unsigned char *pDst = upred_ptr;

    // U & V
    for (i = 0; i < 2; i++) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      if ((omv_row | omv_col) & 15) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel interpolation
        xd->subpixel_predict8x8(pTemp, len, omv_col & 15,
                                omv_row & 15, pDst, 8);
      }
#else   /* CONFIG_SIXTEENTH_SUBPEL_UV */
      if ((mv_row | mv_col) & 7) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel interpolation
        xd->subpixel_predict8x8(pTemp, len, mv_col & 7,
                                mv_row & 7, pDst, 8);
      }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
      else {
        // Apply prediction filter as we copy from source to destination
        filter_mb(pSrc, pre_stride, pDst, 8, 8, 8);
      }

      // V
      pSrc = vptr;
      pDst = vpred_ptr;
    }
  } else
#endif
#if CONFIG_SIXTEENTH_SUBPEL_UV
    if ((omv_row | omv_col) & 15) {
      xd->subpixel_predict8x8(uptr, pre_stride, omv_col & 15, omv_row & 15, upred_ptr, 8);
      xd->subpixel_predict8x8(vptr, pre_stride, omv_col & 15, omv_row & 15, vpred_ptr, 8);
    }
#else   /* CONFIG_SIXTEENTH_SUBPEL_UV */
    if ((mv_row | mv_col) & 7) {
      xd->subpixel_predict8x8(uptr, pre_stride, mv_col & 7, mv_row & 7, upred_ptr, 8);
      xd->subpixel_predict8x8(vptr, pre_stride, mv_col & 7, mv_row & 7, vpred_ptr, 8);
    }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
    else {
      RECON_INVOKE(&xd->rtcd->recon, copy8x8)(uptr, pre_stride, upred_ptr, 8);
      RECON_INVOKE(&xd->rtcd->recon, copy8x8)(vptr, pre_stride, vpred_ptr, 8);
    }
}

/*encoder only*/
void vp8_build_inter4x4_predictors_mbuv(MACROBLOCKD *x) {
  int i, j;

  /* build uv mvs */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int yoffset = i * 8 + j * 2;
      int uoffset = 16 + i * 2 + j;
      int voffset = 20 + i * 2 + j;
      int temp;

      temp = x->block[yoffset  ].bmi.as_mv.first.as_mv.row
             + x->block[yoffset + 1].bmi.as_mv.first.as_mv.row
             + x->block[yoffset + 4].bmi.as_mv.first.as_mv.row
             + x->block[yoffset + 5].bmi.as_mv.first.as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      x->block[uoffset].bmi.as_mv.first.as_mv.row = (temp / 8) & x->fullpixel_mask;

      temp = x->block[yoffset  ].bmi.as_mv.first.as_mv.col
             + x->block[yoffset + 1].bmi.as_mv.first.as_mv.col
             + x->block[yoffset + 4].bmi.as_mv.first.as_mv.col
             + x->block[yoffset + 5].bmi.as_mv.first.as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      x->block[uoffset].bmi.as_mv.first.as_mv.col = (temp / 8) & x->fullpixel_mask;

      x->block[voffset].bmi.as_mv.first.as_mv.row =
        x->block[uoffset].bmi.as_mv.first.as_mv.row;
      x->block[voffset].bmi.as_mv.first.as_mv.col =
        x->block[uoffset].bmi.as_mv.first.as_mv.col;

      if (x->mode_info_context->mbmi.second_ref_frame) {
        temp = x->block[yoffset  ].bmi.as_mv.second.as_mv.row
               + x->block[yoffset + 1].bmi.as_mv.second.as_mv.row
               + x->block[yoffset + 4].bmi.as_mv.second.as_mv.row
               + x->block[yoffset + 5].bmi.as_mv.second.as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        x->block[uoffset].bmi.as_mv.second.as_mv.row = (temp / 8) & x->fullpixel_mask;

        temp = x->block[yoffset  ].bmi.as_mv.second.as_mv.col
               + x->block[yoffset + 1].bmi.as_mv.second.as_mv.col
               + x->block[yoffset + 4].bmi.as_mv.second.as_mv.col
               + x->block[yoffset + 5].bmi.as_mv.second.as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        x->block[uoffset].bmi.as_mv.second.as_mv.col = (temp / 8) & x->fullpixel_mask;

        x->block[voffset].bmi.as_mv.second.as_mv.row =
          x->block[uoffset].bmi.as_mv.second.as_mv.row;
        x->block[voffset].bmi.as_mv.second.as_mv.col =
          x->block[uoffset].bmi.as_mv.second.as_mv.col;
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *d0 = &x->block[i];
    BLOCKD *d1 = &x->block[i + 1];

    if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
      build_inter_predictors2b(x, d0, 8);
    else {
      vp8_build_inter_predictors_b(d0, 8, x->subpixel_predict);
      vp8_build_inter_predictors_b(d1, 8, x->subpixel_predict);
    }

    if (x->mode_info_context->mbmi.second_ref_frame) {
      vp8_build_2nd_inter_predictors_b(d0, 8, x->subpixel_predict_avg);
      vp8_build_2nd_inter_predictors_b(d1, 8, x->subpixel_predict_avg);
    }
  }
}


/*encoder only*/
void vp8_build_inter16x16_predictors_mby(MACROBLOCKD *xd) {
  unsigned char *ptr_base;
  unsigned char *ptr;
  unsigned char *pred_ptr = xd->predictor;
  int mv_row = xd->mode_info_context->mbmi.mv.as_mv.row;
  int mv_col = xd->mode_info_context->mbmi.mv.as_mv.col;
  int pre_stride = xd->block[0].pre_stride;

  ptr_base = xd->pre.y_buffer;
  ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

#if CONFIG_PRED_FILTER
  if (xd->mode_info_context->mbmi.pred_filter_enabled) {
    // Produce predictor from the filtered source
    if ((mv_row | mv_col) & 7) {
      // Sub-pel filter needs extended input
#if CONFIG_ENHANCED_INTERP
      int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
      int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
      int len = 15 + (Interp_Extend << 1);
      unsigned char Temp[32 * 32]; // Data required by sub-pel filter
      unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);

      // Copy extended MB into Temp array, applying the spatial filter
      filter_mb(ptr - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                Temp, len, len, len);

      // Sub-pel interpolation
#if CONFIG_SIXTEENTH_SUBPEL_UV
      xd->subpixel_predict16x16(pTemp, len, (mv_col & 7) << 1,
                                (mv_row & 7) << 1, pred_ptr, 16);
#else
      xd->subpixel_predict16x16(pTemp, len, mv_col & 7,
                                mv_row & 7, pred_ptr, 16);
#endif
    } else {
      // Apply spatial filter to create the prediction directly
      filter_mb(ptr, pre_stride, pred_ptr, 16, 16, 16);
    }
  } else
#endif
    if ((mv_row | mv_col) & 7) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      xd->subpixel_predict16x16(ptr, pre_stride, (mv_col & 7) << 1,
                                (mv_row & 7) << 1, pred_ptr, 16);
#else
      xd->subpixel_predict16x16(ptr, pre_stride, mv_col & 7,
                                mv_row & 7, pred_ptr, 16);
#endif
    } else {
      RECON_INVOKE(&xd->rtcd->recon, copy16x16)(ptr, pre_stride, pred_ptr, 16);
    }
}

static void clamp_mv_to_umv_border(MV *mv, const MACROBLOCKD *xd) {
  /* If the MV points so far into the UMV border that no visible pixels
   * are used for reconstruction, the subpel part of the MV can be
   * discarded and the MV limited to 16 pixels with equivalent results.
   *
   * This limit kicks in at 19 pixels for the top and left edges, for
   * the 16 pixels plus 3 taps right of the central pixel when subpel
   * filtering. The bottom and right edges use 16 pixels plus 2 pixels
   * left of the central pixel when filtering.
   */
  if (mv->col < (xd->mb_to_left_edge - ((16 + INTERP_EXTEND) << 3)))
    mv->col = xd->mb_to_left_edge - (16 << 3);
  else if (mv->col > xd->mb_to_right_edge + ((15 + INTERP_EXTEND) << 3))
    mv->col = xd->mb_to_right_edge + (16 << 3);

  if (mv->row < (xd->mb_to_top_edge - ((16 + INTERP_EXTEND) << 3)))
    mv->row = xd->mb_to_top_edge - (16 << 3);
  else if (mv->row > xd->mb_to_bottom_edge + ((15 + INTERP_EXTEND) << 3))
    mv->row = xd->mb_to_bottom_edge + (16 << 3);
}

/* A version of the above function for chroma block MVs.*/
static void clamp_uvmv_to_umv_border(MV *mv, const MACROBLOCKD *xd) {
  mv->col = (2 * mv->col < (xd->mb_to_left_edge - ((16 + INTERP_EXTEND) << 3))) ?
            (xd->mb_to_left_edge - (16 << 3)) >> 1 : mv->col;
  mv->col = (2 * mv->col > xd->mb_to_right_edge + ((15 + INTERP_EXTEND) << 3)) ?
            (xd->mb_to_right_edge + (16 << 3)) >> 1 : mv->col;

  mv->row = (2 * mv->row < (xd->mb_to_top_edge - ((16 + INTERP_EXTEND) << 3))) ?
            (xd->mb_to_top_edge - (16 << 3)) >> 1 : mv->row;
  mv->row = (2 * mv->row > xd->mb_to_bottom_edge + ((15 + INTERP_EXTEND) << 3)) ?
            (xd->mb_to_bottom_edge + (16 << 3)) >> 1 : mv->row;
}



void vp8_build_inter16x16_predictors_mb(MACROBLOCKD *x,
                                        unsigned char *dst_y,
                                        unsigned char *dst_u,
                                        unsigned char *dst_v,
                                        int dst_ystride,
                                        int dst_uvstride) {
  int offset;
  unsigned char *ptr;
  unsigned char *uptr, *vptr;

  int_mv _o16x16mv;
  int_mv _16x16mv;

  unsigned char *ptr_base = x->pre.y_buffer;
  int pre_stride = x->block[0].pre_stride;

  _16x16mv.as_int = x->mode_info_context->mbmi.mv.as_int;

  if (x->mode_info_context->mbmi.need_to_clamp_mvs) {
    clamp_mv_to_umv_border(&_16x16mv.as_mv, x);
  }

  ptr = ptr_base + (_16x16mv.as_mv.row >> 3) * pre_stride +
        (_16x16mv.as_mv.col >> 3);

#if CONFIG_PRED_FILTER
  if (x->mode_info_context->mbmi.pred_filter_enabled) {
    if (_16x16mv.as_int & 0x00070007) {
      // Sub-pel filter needs extended input
#if CONFIG_ENHANCED_INTERP
      int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
      int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
      int len = 15 + (Interp_Extend << 1);
      unsigned char Temp[32 * 32]; // Data required by the sub-pel filter
      unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);

      // Copy extended MB into Temp array, applying the spatial filter
      filter_mb(ptr - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                Temp, len, len, len);

      // Sub-pel filter
#if CONFIG_SIXTEENTH_SUBPEL_UV
      x->subpixel_predict16x16(pTemp, len,
                               (_16x16mv.as_mv.col & 7) << 1,
                               (_16x16mv.as_mv.row & 7) << 1,
                               dst_y, dst_ystride);
#else
      x->subpixel_predict16x16(pTemp, len,
                               _16x16mv.as_mv.col & 7,
                               _16x16mv.as_mv.row & 7,
                               dst_y, dst_ystride);
#endif
    } else {
      // Apply spatial filter to create the prediction directly
      filter_mb(ptr, pre_stride, dst_y, dst_ystride, 16, 16);
    }
  } else
#endif
    if (_16x16mv.as_int & 0x00070007) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      x->subpixel_predict16x16(ptr, pre_stride, (_16x16mv.as_mv.col & 7) << 1,
                               (_16x16mv.as_mv.row & 7) << 1,
                               dst_y, dst_ystride);
#else
      x->subpixel_predict16x16(ptr, pre_stride, _16x16mv.as_mv.col & 7,
                               _16x16mv.as_mv.row & 7, dst_y, dst_ystride);
#endif
    } else {
      RECON_INVOKE(&x->rtcd->recon, copy16x16)(ptr, pre_stride, dst_y,
                                               dst_ystride);
    }

  _o16x16mv = _16x16mv;
  /* calc uv motion vectors */
  if (_16x16mv.as_mv.row < 0)
    _16x16mv.as_mv.row -= 1;
  else
    _16x16mv.as_mv.row += 1;

  if (_16x16mv.as_mv.col < 0)
    _16x16mv.as_mv.col -= 1;
  else
    _16x16mv.as_mv.col += 1;

  _16x16mv.as_mv.row /= 2;
  _16x16mv.as_mv.col /= 2;

  _16x16mv.as_mv.row &= x->fullpixel_mask;
  _16x16mv.as_mv.col &= x->fullpixel_mask;

  pre_stride >>= 1;
  offset = (_16x16mv.as_mv.row >> 3) * pre_stride + (_16x16mv.as_mv.col >> 3);
  uptr = x->pre.u_buffer + offset;
  vptr = x->pre.v_buffer + offset;

#if CONFIG_PRED_FILTER
  if (x->mode_info_context->mbmi.pred_filter_enabled) {
    int i;
    unsigned char *pSrc = uptr;
    unsigned char *pDst = dst_u;
#if CONFIG_ENHANCED_INTERP
    int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
    int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
    int len = 7 + (Interp_Extend << 1);
    unsigned char Temp[32 * 32]; // Data required by the sub-pel filter
    unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);

    // U & V
    for (i = 0; i < 2; i++) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      if (_o16x16mv.as_int & 0x000f000f) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel filter
        x->subpixel_predict8x8(pTemp, len,
                               _o16x16mv.as_mv.col & 15,
                               _o16x16mv.as_mv.row & 15,
                               pDst, dst_uvstride);
      }
#else  /* CONFIG_SIXTEENTH_SUBPEL_UV */
      if (_16x16mv.as_int & 0x00070007) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel filter
        x->subpixel_predict8x8(pTemp, len,
                               _16x16mv.as_mv.col & 7,
                               _16x16mv.as_mv.row & 7,
                               pDst, dst_uvstride);
      }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
      else {
        filter_mb(pSrc, pre_stride, pDst, dst_uvstride, 8, 8);
      }

      // V
      pSrc = vptr;
      pDst = dst_v;
    }
  } else
#endif
#if CONFIG_SIXTEENTH_SUBPEL_UV
    if (_o16x16mv.as_int & 0x000f000f) {
      x->subpixel_predict8x8(uptr, pre_stride, _o16x16mv.as_mv.col & 15,  _o16x16mv.as_mv.row & 15, dst_u, dst_uvstride);
      x->subpixel_predict8x8(vptr, pre_stride, _o16x16mv.as_mv.col & 15,  _o16x16mv.as_mv.row & 15, dst_v, dst_uvstride);
    }
#else  /* CONFIG_SIXTEENTH_SUBPEL_UV */
    if (_16x16mv.as_int & 0x00070007) {
      x->subpixel_predict8x8(uptr, pre_stride, _16x16mv.as_mv.col & 7,  _16x16mv.as_mv.row & 7, dst_u, dst_uvstride);
      x->subpixel_predict8x8(vptr, pre_stride, _16x16mv.as_mv.col & 7,  _16x16mv.as_mv.row & 7, dst_v, dst_uvstride);
    }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
    else {
      RECON_INVOKE(&x->rtcd->recon, copy8x8)(uptr, pre_stride, dst_u, dst_uvstride);
      RECON_INVOKE(&x->rtcd->recon, copy8x8)(vptr, pre_stride, dst_v, dst_uvstride);
    }

}

/*
 * This function should be called after an initial call to
 * vp8_build_inter16x16_predictors_mb() or _mby()/_mbuv().
 * It will run a second sixtap filter on a (different) ref
 * frame and average the result with the output of the
 * first sixtap filter. The second reference frame is stored
 * in x->second_pre (the reference frame index is in
 * x->mode_info_context->mbmi.second_ref_frame). The second
 * motion vector is x->mode_info_context->mbmi.second_mv.
 *
 * This allows blending prediction from two reference frames
 * which sometimes leads to better prediction than from a
 * single reference framer.
 */
void vp8_build_2nd_inter16x16_predictors_mb(MACROBLOCKD *x,
                                            unsigned char *dst_y,
                                            unsigned char *dst_u,
                                            unsigned char *dst_v,
                                            int dst_ystride,
                                            int dst_uvstride) {
  int offset;
  unsigned char *ptr;
  unsigned char *uptr, *vptr;

  int_mv _16x16mv;
  int mv_row;
  int mv_col;

  int omv_row, omv_col;

  unsigned char *ptr_base = x->second_pre.y_buffer;
  int pre_stride = x->block[0].pre_stride;

  _16x16mv.as_int = x->mode_info_context->mbmi.second_mv.as_int;

  if (x->mode_info_context->mbmi.need_to_clamp_secondmv) {
    clamp_mv_to_umv_border(&_16x16mv.as_mv, x);
  }

  mv_row = _16x16mv.as_mv.row;
  mv_col = _16x16mv.as_mv.col;

  ptr = ptr_base + (mv_row >> 3) * pre_stride + (mv_col >> 3);

#if CONFIG_PRED_FILTER
  if (x->mode_info_context->mbmi.pred_filter_enabled) {
    if ((mv_row | mv_col) & 7) {
      // Sub-pel filter needs extended input
#if CONFIG_ENHANCED_INTERP
      int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
      int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
      int len = 15 + (Interp_Extend << 1);
      unsigned char Temp[32 * 32]; // Data required by sub-pel filter
      unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);

      // Copy extended MB into Temp array, applying the spatial filter
      filter_mb(ptr - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                Temp, len, len, len);

      // Sub-pel filter
#if CONFIG_SIXTEENTH_SUBPEL_UV
      x->subpixel_predict_avg16x16(pTemp, len, (mv_col & 7) << 1,
                                   (mv_row & 7) << 1, dst_y, dst_ystride);
#else
      x->subpixel_predict_avg16x16(pTemp, len, mv_col & 7,
                                   mv_row & 7, dst_y, dst_ystride);
#endif
    } else {
      // TODO Needs to AVERAGE with the dst_y
      // For now, do not apply the prediction filter in these cases!
      RECON_INVOKE(&x->rtcd->recon, avg16x16)(ptr, pre_stride, dst_y,
                                              dst_ystride);
    }
  } else
#endif  // CONFIG_PRED_FILTER
  {
    if ((mv_row | mv_col) & 7) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      x->subpixel_predict_avg16x16(ptr, pre_stride, (mv_col & 7) << 1,
                                   (mv_row & 7) << 1, dst_y, dst_ystride);
#else
      x->subpixel_predict_avg16x16(ptr, pre_stride, mv_col & 7,
                                   mv_row & 7, dst_y, dst_ystride);
#endif
    } else {
      RECON_INVOKE(&x->rtcd->recon, avg16x16)(ptr, pre_stride, dst_y,
                                              dst_ystride);
    }
  }

  /* calc uv motion vectors */
  omv_row = mv_row;
  omv_col = mv_col;
  mv_row = (mv_row + (mv_row > 0)) >> 1;
  mv_col = (mv_col + (mv_col > 0)) >> 1;

  mv_row &= x->fullpixel_mask;
  mv_col &= x->fullpixel_mask;

  pre_stride >>= 1;
  offset = (mv_row >> 3) * pre_stride + (mv_col >> 3);
  uptr = x->second_pre.u_buffer + offset;
  vptr = x->second_pre.v_buffer + offset;

#if CONFIG_PRED_FILTER
  if (x->mode_info_context->mbmi.pred_filter_enabled) {
    int i;
#if CONFIG_ENHANCED_INTERP
    int Interp_Extend = 4;  // 8-tap filter needs 3+4 pels extension
#else
    int Interp_Extend = 3;  // 6-tap filter needs 2+3 pels extension
#endif
    int len = 7 + (Interp_Extend << 1);
    unsigned char Temp[32 * 32]; // Data required by sub-pel filter
    unsigned char *pTemp = Temp + (Interp_Extend - 1) * (len + 1);
    unsigned char *pSrc = uptr;
    unsigned char *pDst = dst_u;

    // U & V
    for (i = 0; i < 2; i++) {
#if CONFIG_SIXTEENTH_SUBPEL_UV
      if ((omv_row | omv_col) & 15) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel filter
        x->subpixel_predict_avg8x8(pTemp, len, omv_col & 15,
                                   omv_row & 15, pDst, dst_uvstride);
      }
#else  /* CONFIG_SIXTEENTH_SUBPEL_UV */
      if ((mv_row | mv_col) & 7) {
        // Copy extended MB into Temp array, applying the spatial filter
        filter_mb(pSrc - (Interp_Extend - 1) * (pre_stride + 1), pre_stride,
                  Temp, len, len, len);

        // Sub-pel filter
        x->subpixel_predict_avg8x8(pTemp, len, mv_col & 7, mv_row & 7,
                                   pDst, dst_uvstride);
      }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
      else {
        // TODO Needs to AVERAGE with the dst_[u|v]
        // For now, do not apply the prediction filter here!
        RECON_INVOKE(&x->rtcd->recon, avg8x8)(pSrc, pre_stride, pDst,
                                              dst_uvstride);
      }

      // V
      pSrc = vptr;
      pDst = dst_v;
    }
  } else
#endif  // CONFIG_PRED_FILTER
#if CONFIG_SIXTEENTH_SUBPEL_UV
    if ((omv_row | omv_col) & 15) {
      x->subpixel_predict_avg8x8(uptr, pre_stride, omv_col & 15, omv_row & 15, dst_u, dst_uvstride);
      x->subpixel_predict_avg8x8(vptr, pre_stride, omv_col & 15, omv_row & 15, dst_v, dst_uvstride);
    }
#else  /* CONFIG_SIXTEENTH_SUBPEL_UV */
    if ((mv_row | mv_col) & 7) {
      x->subpixel_predict_avg8x8(uptr, pre_stride, mv_col & 7, mv_row & 7, dst_u, dst_uvstride);
      x->subpixel_predict_avg8x8(vptr, pre_stride, mv_col & 7, mv_row & 7, dst_v, dst_uvstride);
    }
#endif  /* CONFIG_SIXTEENTH_SUBPEL_UV */
    else {
      RECON_INVOKE(&x->rtcd->recon, avg8x8)(uptr, pre_stride, dst_u, dst_uvstride);
      RECON_INVOKE(&x->rtcd->recon, avg8x8)(vptr, pre_stride, dst_v, dst_uvstride);
    }
}

static void build_inter4x4_predictors_mb(MACROBLOCKD *x) {
  int i;

  if (x->mode_info_context->mbmi.partitioning < 3) {
    x->block[ 0].bmi = x->mode_info_context->bmi[ 0];
    x->block[ 2].bmi = x->mode_info_context->bmi[ 2];
    x->block[ 8].bmi = x->mode_info_context->bmi[ 8];
    x->block[10].bmi = x->mode_info_context->bmi[10];

    if (x->mode_info_context->mbmi.need_to_clamp_mvs) {
      clamp_mv_to_umv_border(&x->block[ 0].bmi.as_mv.first.as_mv, x);
      clamp_mv_to_umv_border(&x->block[ 2].bmi.as_mv.first.as_mv, x);
      clamp_mv_to_umv_border(&x->block[ 8].bmi.as_mv.first.as_mv, x);
      clamp_mv_to_umv_border(&x->block[10].bmi.as_mv.first.as_mv, x);
      if (x->mode_info_context->mbmi.second_ref_frame) {
        clamp_mv_to_umv_border(&x->block[ 0].bmi.as_mv.second.as_mv, x);
        clamp_mv_to_umv_border(&x->block[ 2].bmi.as_mv.second.as_mv, x);
        clamp_mv_to_umv_border(&x->block[ 8].bmi.as_mv.second.as_mv, x);
        clamp_mv_to_umv_border(&x->block[10].bmi.as_mv.second.as_mv, x);
      }
    }


    build_inter_predictors4b(x, &x->block[ 0], 16);
    build_inter_predictors4b(x, &x->block[ 2], 16);
    build_inter_predictors4b(x, &x->block[ 8], 16);
    build_inter_predictors4b(x, &x->block[10], 16);

    if (x->mode_info_context->mbmi.second_ref_frame) {
      build_2nd_inter_predictors4b(x, &x->block[ 0], 16);
      build_2nd_inter_predictors4b(x, &x->block[ 2], 16);
      build_2nd_inter_predictors4b(x, &x->block[ 8], 16);
      build_2nd_inter_predictors4b(x, &x->block[10], 16);
    }
  } else {
    for (i = 0; i < 16; i += 2) {
      BLOCKD *d0 = &x->block[i];
      BLOCKD *d1 = &x->block[i + 1];

      x->block[i + 0].bmi = x->mode_info_context->bmi[i + 0];
      x->block[i + 1].bmi = x->mode_info_context->bmi[i + 1];

      if (x->mode_info_context->mbmi.need_to_clamp_mvs) {
        clamp_mv_to_umv_border(&x->block[i + 0].bmi.as_mv.first.as_mv, x);
        clamp_mv_to_umv_border(&x->block[i + 1].bmi.as_mv.first.as_mv, x);
        if (x->mode_info_context->mbmi.second_ref_frame) {
          clamp_mv_to_umv_border(&x->block[i + 0].bmi.as_mv.second.as_mv, x);
          clamp_mv_to_umv_border(&x->block[i + 1].bmi.as_mv.second.as_mv, x);
        }
      }

      if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
        build_inter_predictors2b(x, d0, 16);
      else {
        vp8_build_inter_predictors_b(d0, 16, x->subpixel_predict);
        vp8_build_inter_predictors_b(d1, 16, x->subpixel_predict);
      }

      if (x->mode_info_context->mbmi.second_ref_frame) {
        vp8_build_2nd_inter_predictors_b(d0, 16, x->subpixel_predict_avg);
        vp8_build_2nd_inter_predictors_b(d1, 16, x->subpixel_predict_avg);
      }
    }
  }

  for (i = 16; i < 24; i += 2) {
    BLOCKD *d0 = &x->block[i];
    BLOCKD *d1 = &x->block[i + 1];

    if (d0->bmi.as_mv.first.as_int == d1->bmi.as_mv.first.as_int)
      build_inter_predictors2b(x, d0, 8);
    else {
      vp8_build_inter_predictors_b(d0, 8, x->subpixel_predict);
      vp8_build_inter_predictors_b(d1, 8, x->subpixel_predict);
    }

    if (x->mode_info_context->mbmi.second_ref_frame) {
      vp8_build_2nd_inter_predictors_b(d0, 8, x->subpixel_predict_avg);
      vp8_build_2nd_inter_predictors_b(d1, 8, x->subpixel_predict_avg);
    }
  }
}

static
void build_4x4uvmvs(MACROBLOCKD *x) {
  int i, j;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int yoffset = i * 8 + j * 2;
      int uoffset = 16 + i * 2 + j;
      int voffset = 20 + i * 2 + j;

      int temp;

      temp = x->mode_info_context->bmi[yoffset + 0].as_mv.first.as_mv.row
             + x->mode_info_context->bmi[yoffset + 1].as_mv.first.as_mv.row
             + x->mode_info_context->bmi[yoffset + 4].as_mv.first.as_mv.row
             + x->mode_info_context->bmi[yoffset + 5].as_mv.first.as_mv.row;

      if (temp < 0) temp -= 4;
      else temp += 4;

      x->block[uoffset].bmi.as_mv.first.as_mv.row = (temp / 8) & x->fullpixel_mask;

      temp = x->mode_info_context->bmi[yoffset + 0].as_mv.first.as_mv.col
             + x->mode_info_context->bmi[yoffset + 1].as_mv.first.as_mv.col
             + x->mode_info_context->bmi[yoffset + 4].as_mv.first.as_mv.col
             + x->mode_info_context->bmi[yoffset + 5].as_mv.first.as_mv.col;

      if (temp < 0) temp -= 4;
      else temp += 4;

      x->block[uoffset].bmi.as_mv.first.as_mv.col = (temp / 8) & x->fullpixel_mask;

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&x->block[uoffset].bmi.as_mv.first.as_mv, x);

      // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
      clamp_uvmv_to_umv_border(&x->block[uoffset].bmi.as_mv.first.as_mv, x);

      x->block[voffset].bmi.as_mv.first.as_mv.row =
        x->block[uoffset].bmi.as_mv.first.as_mv.row;
      x->block[voffset].bmi.as_mv.first.as_mv.col =
        x->block[uoffset].bmi.as_mv.first.as_mv.col;

      if (x->mode_info_context->mbmi.second_ref_frame) {
        temp = x->mode_info_context->bmi[yoffset + 0].as_mv.second.as_mv.row
               + x->mode_info_context->bmi[yoffset + 1].as_mv.second.as_mv.row
               + x->mode_info_context->bmi[yoffset + 4].as_mv.second.as_mv.row
               + x->mode_info_context->bmi[yoffset + 5].as_mv.second.as_mv.row;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        x->block[uoffset].bmi.as_mv.second.as_mv.row = (temp / 8) & x->fullpixel_mask;

        temp = x->mode_info_context->bmi[yoffset + 0].as_mv.second.as_mv.col
               + x->mode_info_context->bmi[yoffset + 1].as_mv.second.as_mv.col
               + x->mode_info_context->bmi[yoffset + 4].as_mv.second.as_mv.col
               + x->mode_info_context->bmi[yoffset + 5].as_mv.second.as_mv.col;

        if (temp < 0) {
          temp -= 4;
        } else {
          temp += 4;
        }

        x->block[uoffset].bmi.as_mv.second.as_mv.col = (temp / 8) & x->fullpixel_mask;

        // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(&x->block[uoffset].bmi.as_mv.second.as_mv, x);

        // if (x->mode_info_context->mbmi.need_to_clamp_mvs)
        clamp_uvmv_to_umv_border(&x->block[uoffset].bmi.as_mv.second.as_mv, x);

        x->block[voffset].bmi.as_mv.second.as_mv.row =
          x->block[uoffset].bmi.as_mv.second.as_mv.row;
        x->block[voffset].bmi.as_mv.second.as_mv.col =
          x->block[uoffset].bmi.as_mv.second.as_mv.col;
      }
    }
  }
}

void vp8_build_inter_predictors_mb(MACROBLOCKD *x) {
  if (x->mode_info_context->mbmi.mode != SPLITMV) {
    vp8_build_inter16x16_predictors_mb(x, x->predictor, &x->predictor[256],
                                       &x->predictor[320], 16, 8);

    if (x->mode_info_context->mbmi.second_ref_frame) {
      /* 256 = offset of U plane in Y+U+V buffer;
       * 320 = offset of V plane in Y+U+V buffer.
       * (256=16x16, 320=16x16+8x8). */
      vp8_build_2nd_inter16x16_predictors_mb(x, x->predictor,
                                             &x->predictor[256],
                                             &x->predictor[320], 16, 8);
    }
  } else {
    build_4x4uvmvs(x);
    build_inter4x4_predictors_mb(x);
  }
}
