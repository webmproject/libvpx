/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "invtrans.h"



static void recon_dcblock(MACROBLOCKD *x) {
  BLOCKD *b = &x->block[24];
  int i;

  for (i = 0; i < 16; i++) {
    x->block[i].dqcoeff[0] = b->diff[i];
  }

}
static void recon_dcblock_8x8(MACROBLOCKD *x) {
  BLOCKD *b = &x->block[24]; // for coeff 0, 2, 8, 10
  x->block[0].dqcoeff[0] = b->diff[0];
  x->block[4].dqcoeff[0] = b->diff[1];
  x->block[8].dqcoeff[0] = b->diff[4];
  x->block[12].dqcoeff[0] = b->diff[8];

}

#if CONFIG_HYBRIDTRANSFORM
void vp8_inverse_htransform_b(const vp8_idct_rtcd_vtable_t *rtcd, BLOCKD *b, int pitch) {
  vp8_iht4x4llm_c(b->dqcoeff, b->diff, pitch, b->bmi.as_mode.tx_type);
}
#endif

void vp8_inverse_transform_b(const vp8_idct_rtcd_vtable_t *rtcd, BLOCKD *b, int pitch) {
  if (b->eob <= 1)
    IDCT_INVOKE(rtcd, idct1)(b->dqcoeff, b->diff, pitch);
  else
    IDCT_INVOKE(rtcd, idct16)(b->dqcoeff, b->diff, pitch);
}


void vp8_inverse_transform_mby(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  /* do 2nd order transform on the dc block */
  IDCT_INVOKE(rtcd, iwalsh16)(x->block[24].dqcoeff, x->block[24].diff);

  recon_dcblock(x);

  for (i = 0; i < 16; i++) {
    vp8_inverse_transform_b(rtcd, &x->block[i], 32);
  }

}
void vp8_inverse_transform_mbuv(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  for (i = 16; i < 24; i++) {
    vp8_inverse_transform_b(rtcd, &x->block[i], 16);
  }

}


void vp8_inverse_transform_mb(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  if (x->mode_info_context->mbmi.mode != B_PRED &&
      x->mode_info_context->mbmi.mode != I8X8_PRED &&
      x->mode_info_context->mbmi.mode != SPLITMV) {
    /* do 2nd order transform on the dc block */

    IDCT_INVOKE(rtcd, iwalsh16)(&x->block[24].dqcoeff[0], x->block[24].diff);
    recon_dcblock(x);
  }

  for (i = 0; i < 16; i++) {
    vp8_inverse_transform_b(rtcd, &x->block[i], 32);
  }


  for (i = 16; i < 24; i++) {
    vp8_inverse_transform_b(rtcd, &x->block[i], 16);
  }

}


void vp8_inverse_transform_b_8x8(const vp8_idct_rtcd_vtable_t *rtcd, short *input_dqcoeff, short *output_coeff, int pitch) { // pay attention to use when 8x8
  // int b,i;
  // if (b->eob > 1)
  IDCT_INVOKE(rtcd, idct8)(input_dqcoeff, output_coeff, pitch);
  // else
  // IDCT_INVOKE(rtcd, idct8_1)(b->dqcoeff, b->diff, pitch);//pitch

}


void vp8_inverse_transform_mby_8x8(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  // do 2nd order transform on the dc block
  IDCT_INVOKE(rtcd, ihaar2)(x->block[24].dqcoeff, x->block[24].diff, 8);

  recon_dcblock_8x8(x); // need to change for 8x8
  for (i = 0; i < 9; i += 8) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 32);
  }
  for (i = 2; i < 11; i += 8) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i + 2].dqcoeff[0], &x->block[i].diff[0], 32);
  }

}
void vp8_inverse_transform_mbuv_8x8(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  for (i = 16; i < 24; i += 4) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 16);
  }

}


void vp8_inverse_transform_mb_8x8(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  if (x->mode_info_context->mbmi.mode != B_PRED &&
      x->mode_info_context->mbmi.mode != SPLITMV) {
    // do 2nd order transform on the dc block

    IDCT_INVOKE(rtcd, ihaar2)(&x->block[24].dqcoeff[0], x->block[24].diff, 8);// dqcoeff[0]
    recon_dcblock_8x8(x); // need to change for 8x8

  }

  for (i = 0; i < 9; i += 8) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 32);
  }
  for (i = 2; i < 11; i += 8) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i + 2].dqcoeff[0], &x->block[i].diff[0], 32);
  }


  for (i = 16; i < 24; i += 4) {
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 16);
  }

}

#if CONFIG_TX16X16
void vp8_inverse_transform_b_16x16(const vp8_idct_rtcd_vtable_t *rtcd,
                                   short *input_dqcoeff,
                                   short *output_coeff, int pitch) {
  IDCT_INVOKE(rtcd, idct16x16)(input_dqcoeff, output_coeff, pitch);
}

void vp8_inverse_transform_mby_16x16(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
    vp8_inverse_transform_b_16x16(rtcd, &x->block[0].dqcoeff[0], &x->block[0].diff[0], 32);
}

// U,V blocks are 8x8 per macroblock, so just run 8x8
void vp8_inverse_transform_mbuv_16x16(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;
  for (i = 16; i < 24; i += 4)
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 16);
}

void vp8_inverse_transform_mb_16x16(const vp8_idct_rtcd_vtable_t *rtcd, MACROBLOCKD *x) {
  int i;

  // Luma
  vp8_inverse_transform_b_16x16(rtcd, &x->block[0].dqcoeff[0], &x->block[0].diff[0], 32);

  // U, V
  // Chroma blocks are downscaled, so run an 8x8 on them.
  for (i = 16; i < 24; i+= 4)
    vp8_inverse_transform_b_8x8(rtcd, &x->block[i].dqcoeff[0], &x->block[i].diff[0], 16);
}
#endif
