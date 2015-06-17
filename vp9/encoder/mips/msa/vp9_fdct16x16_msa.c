/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./vp9_rtcd.h"
#include "vp9/encoder/mips/msa/vp9_fdct_msa.h"

static void fdct8x16_1d_column(const int16_t *input, int16_t *tmp_ptr,
                               int32_t src_stride) {
  v8i16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;
  v8i16 in8, in9, in10, in11, in12, in13, in14, in15;
  v8i16 stp21, stp22, stp23, stp24, stp25, stp26, stp30;
  v8i16 stp31, stp32, stp33, stp34, stp35, stp36, stp37;
  v8i16 vec0, vec1, vec2, vec3, vec4, vec5, cnst0, cnst1, cnst4, cnst5;
  v8i16 coeff = { cospi_16_64, -cospi_16_64, cospi_8_64, cospi_24_64,
                 -cospi_8_64, -cospi_24_64, cospi_12_64, cospi_20_64 };
  v8i16 coeff1 = { cospi_2_64, cospi_30_64, cospi_14_64, cospi_18_64,
                   cospi_10_64, cospi_22_64, cospi_6_64, cospi_26_64 };
  v8i16 coeff2 = { -cospi_2_64, -cospi_10_64, -cospi_18_64, -cospi_26_64,
                   0, 0, 0, 0 };

  LD_SH16(input, src_stride,
          in0, in1, in2, in3, in4, in5, in6, in7,
          in8, in9, in10, in11, in12, in13, in14, in15);
  SLLI_4V(in0, in1, in2, in3, 2);
  SLLI_4V(in4, in5, in6, in7, 2);
  SLLI_4V(in8, in9, in10, in11, 2);
  SLLI_4V(in12, in13, in14, in15, 2);
  ADD4(in0, in15, in1, in14, in2, in13, in3, in12, tmp0, tmp1, tmp2, tmp3);
  ADD4(in4, in11, in5, in10, in6, in9, in7, in8, tmp4, tmp5, tmp6, tmp7);
  VP9_FDCT8x16_EVEN(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,
                    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7);
  ST_SH8(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp_ptr, 32);
  SUB4(in0, in15, in1, in14, in2, in13, in3, in12, in15, in14, in13, in12);
  SUB4(in4, in11, in5, in10, in6, in9, in7, in8, in11, in10, in9, in8);

  tmp_ptr += 16;

  /* stp 1 */
  ILVL_H2_SH(in10, in13, in11, in12, vec2, vec4);
  ILVR_H2_SH(in10, in13, in11, in12, vec3, vec5);

  cnst4 = __msa_splati_h(coeff, 0);
  stp25 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2, vec3, cnst4);

  cnst5 = __msa_splati_h(coeff, 1);
  cnst5 = __msa_ilvev_h(cnst5, cnst4);
  stp22 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2, vec3, cnst5);
  stp24 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4, vec5, cnst4);
  stp23 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4, vec5, cnst5);

  /* stp2 */
  BUTTERFLY_4(in8, in9, stp22, stp23, stp30, stp31, stp32, stp33);
  BUTTERFLY_4(in15, in14, stp25, stp24, stp37, stp36, stp35, stp34);
  ILVL_H2_SH(stp36, stp31, stp35, stp32, vec2, vec4);
  ILVR_H2_SH(stp36, stp31, stp35, stp32, vec3, vec5);
  SPLATI_H2_SH(coeff, 2, 3, cnst0, cnst1);
  cnst0 = __msa_ilvev_h(cnst0, cnst1);
  stp26 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2, vec3, cnst0);

  cnst0 = __msa_splati_h(coeff, 4);
  cnst1 = __msa_ilvev_h(cnst1, cnst0);
  stp21 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec2, vec3, cnst1);

  BUTTERFLY_4(stp30, stp37, stp26, stp21, in8, in15, in14, in9);
  ILVRL_H2_SH(in15, in8, vec1, vec0);
  SPLATI_H2_SH(coeff1, 0, 1, cnst0, cnst1);
  cnst0 = __msa_ilvev_h(cnst0, cnst1);

  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr);

  cnst0 = __msa_splati_h(coeff2, 0);
  cnst0 = __msa_ilvev_h(cnst1, cnst0);
  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr + 224);

  ILVRL_H2_SH(in14, in9, vec1, vec0);
  SPLATI_H2_SH(coeff1, 2, 3, cnst0, cnst1);
  cnst1 = __msa_ilvev_h(cnst1, cnst0);

  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst1);
  ST_SH(in8, tmp_ptr + 128);

  cnst1 = __msa_splati_h(coeff2, 2);
  cnst0 = __msa_ilvev_h(cnst0, cnst1);
  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr + 96);

  SPLATI_H2_SH(coeff, 2, 5, cnst0, cnst1);
  cnst1 = __msa_ilvev_h(cnst1, cnst0);

  stp25 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4, vec5, cnst1);

  cnst1 = __msa_splati_h(coeff, 3);
  cnst1 = __msa_ilvev_h(cnst0, cnst1);
  stp22 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec4, vec5, cnst1);

  /* stp4 */
  ADD2(stp34, stp25, stp33, stp22, in13, in10);

  ILVRL_H2_SH(in13, in10, vec1, vec0);
  SPLATI_H2_SH(coeff1, 4, 5, cnst0, cnst1);
  cnst0 = __msa_ilvev_h(cnst0, cnst1);
  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr + 64);

  cnst0 = __msa_splati_h(coeff2, 1);
  cnst0 = __msa_ilvev_h(cnst1, cnst0);
  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr + 160);

  SUB2(stp34, stp25, stp33, stp22, in12, in11);
  ILVRL_H2_SH(in12, in11, vec1, vec0);
  SPLATI_H2_SH(coeff1, 6, 7, cnst0, cnst1);
  cnst1 = __msa_ilvev_h(cnst1, cnst0);

  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst1);
  ST_SH(in8, tmp_ptr + 192);

  cnst1 = __msa_splati_h(coeff2, 3);
  cnst0 = __msa_ilvev_h(cnst0, cnst1);
  in8 = VP9_DOT_SHIFT_RIGHT_PCK_H(vec0, vec1, cnst0);
  ST_SH(in8, tmp_ptr + 32);
}

static void fdct16x8_1d_row(int16_t *input, int16_t *output) {
  v8i16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7;
  v8i16 in8, in9, in10, in11, in12, in13, in14, in15;

  LD_SH8(input, 16, in0, in1, in2, in3, in4, in5, in6, in7);
  LD_SH8((input + 8), 16, in8, in9, in10, in11, in12, in13, in14, in15);
  TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(in8, in9, in10, in11, in12, in13, in14, in15,
                     in8, in9, in10, in11, in12, in13, in14, in15);
  ADD4(in0, 1, in1, 1, in2, 1, in3, 1, in0, in1, in2, in3);
  ADD4(in4, 1, in5, 1, in6, 1, in7, 1, in4, in5, in6, in7);
  ADD4(in8, 1, in9, 1, in10, 1, in11, 1, in8, in9, in10, in11);
  ADD4(in12, 1, in13, 1, in14, 1, in15, 1, in12, in13, in14, in15);
  SRA_4V(in0, in1, in2, in3, 2);
  SRA_4V(in4, in5, in6, in7, 2);
  SRA_4V(in8, in9, in10, in11, 2);
  SRA_4V(in12, in13, in14, in15, 2);
  BUTTERFLY_16(in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11,
               in12, in13, in14, in15, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5,
               tmp6, tmp7, in8, in9, in10, in11, in12, in13, in14, in15);
  ST_SH8(in8, in9, in10, in11, in12, in13, in14, in15, input, 16);
  VP9_FDCT8x16_EVEN(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,
                    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7);
  LD_SH8(input, 16, in8, in9, in10, in11, in12, in13, in14, in15);
  VP9_FDCT8x16_ODD(in8, in9, in10, in11, in12, in13, in14, in15,
                   in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3,
                     tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3);
  ST_SH8(tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3, output, 16);
  TRANSPOSE8x8_SH_SH(tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7,
                     tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7);
  ST_SH8(tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7, output + 8, 16);
}

void vp9_fdct16x16_msa(const int16_t *input, int16_t *output,
                       int32_t src_stride) {
  int32_t i;
  DECLARE_ALIGNED(32, int16_t, tmp_buf[16 * 16]);

  /* column transform */
  for (i = 0; i < 2; ++i) {
    fdct8x16_1d_column((input + 8 * i), (&tmp_buf[0] + 8 * i), src_stride);
  }

  /* row transform */
  for (i = 0; i < 2; ++i) {
    fdct16x8_1d_row((&tmp_buf[0] + (128 * i)), (output + (128 * i)));
  }
}

void vp9_fdct16x16_1_msa(const int16_t *input, int16_t *out, int32_t stride) {
  out[1] = 0;

  out[0] = VP9_LD_HADD(input, stride);
  out[0] += VP9_LD_HADD(input + 8, stride);
  out[0] += VP9_LD_HADD(input + 16 * 8, stride);
  out[0] += VP9_LD_HADD(input + 16 * 8 + 8, stride);
  out[0] >>= 1;
}

static void fadst16_cols_step1_msa(const int16_t *input, int32_t stride,
                                   const int32_t *const0, int16_t *int_buf) {
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 tp0, tp1, tp2, tp3, g0, g1, g2, g3, g8, g9, g10, g11, h0, h1, h2, h3;
  v4i32 k0, k1, k2, k3;

  /* load input data */
  r0 = LD_SH(input);
  r15 = LD_SH(input + 15 * stride);
  r7 = LD_SH(input + 7 * stride);
  r8 = LD_SH(input + 8 * stride);
  SLLI_4V(r0, r15, r7, r8, 2);

  /* stage 1 */
  LD_SW2(const0, 4, k0, k1);
  LD_SW2(const0 + 8, 4, k2, k3);
  VP9_MADD_BF(r15, r0, r7, r8, k0, k1, k2, k3, g0, g1, g2, g3);

  r3 = LD_SH(input + 3 * stride);
  r4 = LD_SH(input + 4 * stride);
  r11 = LD_SH(input + 11 * stride);
  r12 = LD_SH(input + 12 * stride);
  SLLI_4V(r3, r4, r11, r12, 2);

  LD_SW2(const0 + 4 * 4, 4, k0, k1);
  LD_SW2(const0 + 4 * 6, 4, k2, k3);
  VP9_MADD_BF(r11, r4, r3, r12, k0, k1, k2, k3, g8, g9, g10, g11);

  /* stage 2 */
  BUTTERFLY_4(g0, g2, g10, g8, tp0, tp2, tp3, tp1);
  ST_SH2(tp0, tp2, int_buf, 8);
  ST_SH2(tp1, tp3, int_buf + 4 * 8, 8);

  LD_SW2(const0 + 4 * 8, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 10);
  VP9_MADD_BF(g1, g3, g9, g11, k0, k1, k2, k0, h0, h1, h2, h3);

  ST_SH2(h0, h1, int_buf + 8 * 8, 8);
  ST_SH2(h3, h2, int_buf + 12 * 8, 8);

  r9 = LD_SH(input + 9 * stride);
  r6 = LD_SH(input + 6 * stride);
  r1 = LD_SH(input + stride);
  r14 = LD_SH(input + 14 * stride);
  SLLI_4V(r9, r6, r1, r14, 2);

  LD_SW2(const0 + 4 * 11, 4, k0, k1);
  LD_SW2(const0 + 4 * 13, 4, k2, k3);
  VP9_MADD_BF(r9, r6, r1, r14, k0, k1, k2, k3, g0, g1, g2, g3);

  ST_SH2(g1, g3, int_buf + 3 * 8, 4 * 8);

  r13 = LD_SH(input + 13 * stride);
  r2 = LD_SH(input + 2 * stride);
  r5 = LD_SH(input + 5 * stride);
  r10 = LD_SH(input + 10 * stride);
  SLLI_4V(r13, r2, r5, r10, 2);

  LD_SW2(const0 + 4 * 15, 4, k0, k1);
  LD_SW2(const0 + 4 * 17, 4, k2, k3);
  VP9_MADD_BF(r13, r2, r5, r10, k0, k1, k2, k3, h0, h1, h2, h3);

  ST_SH2(h1, h3, int_buf + 11 * 8, 4 * 8);

  BUTTERFLY_4(h0, h2, g2, g0, tp0, tp1, tp2, tp3);
  ST_SH4(tp0, tp1, tp2, tp3, int_buf + 2 * 8, 4 * 8);
}

static void fadst16_cols_step2_msa(int16_t *int_buf, const int32_t *const0,
                                   int16_t *out) {
  int16_t *out_ptr = out + 128;
  v8i16 tp0, tp1, tp2, tp3, g5, g7, g13, g15;
  v8i16 h0, h1, h2, h3, h4, h5, h6, h7, h10, h11;
  v8i16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8i16 out8, out9, out10, out11, out12, out13, out14, out15;
  v4i32 k0, k1, k2, k3;

  LD_SH2(int_buf + 3 * 8, 4 * 8, g13, g15);
  LD_SH2(int_buf + 11 * 8, 4 * 8, g5, g7);
  LD_SW2(const0 + 4 * 19, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 21);
  VP9_MADD_BF(g7, g5, g15, g13, k0, k1, k2, k0, h4, h5, h6, h7);

  tp0 = LD_SH(int_buf + 4 * 8);
  tp1 = LD_SH(int_buf + 5 * 8);
  tp3 = LD_SH(int_buf + 10 * 8);
  tp2 = LD_SH(int_buf + 14 * 8);
  LD_SW2(const0 + 4 * 22, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 24);
  VP9_MADD_BF(tp0, tp1, tp2, tp3, k0, k1, k2, k0, out4, out6, out5, out7);
  out4 = -out4;
  ST_SH(out4, (out + 3 * 16));
  ST_SH(out5, (out_ptr + 4 * 16));

  h1 = LD_SH(int_buf + 9 * 8);
  h3 = LD_SH(int_buf + 12 * 8);
  VP9_MADD_BF(h1, h3, h5, h7, k0, k1, k2, k0, out12, out14, out13, out15);
  out13 = -out13;
  ST_SH(out12, (out + 2 * 16));
  ST_SH(out13, (out_ptr + 5 * 16));

  tp0 = LD_SH(int_buf);
  tp1 = LD_SH(int_buf + 8);
  tp2 = LD_SH(int_buf + 2 * 8);
  tp3 = LD_SH(int_buf + 6 * 8);

  BUTTERFLY_4(tp0, tp1, tp3, tp2, out0, out1, h11, h10);
  out1 = -out1;
  ST_SH(out0, (out));
  ST_SH(out1, (out_ptr + 7 * 16));

  h0 = LD_SH(int_buf + 8 * 8);
  h2 = LD_SH(int_buf + 13 * 8);

  BUTTERFLY_4(h0, h2, h6, h4, out8, out9, out11, out10);
  out8 = -out8;
  ST_SH(out8, (out + 16));
  ST_SH(out9, (out_ptr + 6 * 16));

  /* stage 4 */
  LD_SW2(const0 + 4 * 25, 4, k0, k1);
  LD_SW2(const0 + 4 * 27, 4, k2, k3);
  VP9_MADD_SHORT(h10, h11, k1, k2, out2, out3);
  ST_SH(out2, (out + 7 * 16));
  ST_SH(out3, (out_ptr));

  VP9_MADD_SHORT(out6, out7, k0, k3, out6, out7);
  ST_SH(out6, (out + 4 * 16));
  ST_SH(out7, (out_ptr + 3 * 16));

  VP9_MADD_SHORT(out10, out11, k0, k3, out10, out11);
  ST_SH(out10, (out + 6 * 16));
  ST_SH(out11, (out_ptr + 16));

  VP9_MADD_SHORT(out14, out15, k1, k2, out14, out15);
  ST_SH(out14, (out + 5 * 16));
  ST_SH(out15, (out_ptr + 2 * 16));
}

static void fadst16_transpose_postproc_msa(int16_t *input, int16_t *out) {
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;

  /* load input data */
  LD_SH8(input, 16, l0, l1, l2, l3, l4, l5, l6, l7);
  TRANSPOSE8x8_SH_SH(l0, l1, l2, l3, l4, l5, l6, l7,
                     r0, r1, r2, r3, r4, r5, r6, r7);
  VP9_FDCT_POSTPROC_2V_NEG_H(r0, r1);
  VP9_FDCT_POSTPROC_2V_NEG_H(r2, r3);
  VP9_FDCT_POSTPROC_2V_NEG_H(r4, r5);
  VP9_FDCT_POSTPROC_2V_NEG_H(r6, r7);
  ST_SH8(r0, r1, r2, r3, r4, r5, r6, r7, out, 8);
  out += 64;

  LD_SH8(input + 8, 16, l8, l9, l10, l11, l12, l13, l14, l15);
  TRANSPOSE8x8_SH_SH(l8, l9, l10, l11, l12, l13, l14, l15,
                     r8, r9, r10, r11, r12, r13, r14, r15);
  VP9_FDCT_POSTPROC_2V_NEG_H(r8, r9);
  VP9_FDCT_POSTPROC_2V_NEG_H(r10, r11);
  VP9_FDCT_POSTPROC_2V_NEG_H(r12, r13);
  VP9_FDCT_POSTPROC_2V_NEG_H(r14, r15);
  ST_SH8(r8, r9, r10, r11, r12, r13, r14, r15, out, 8);
  out += 64;

  /* load input data */
  input += 128;
  LD_SH8(input, 16, l0, l1, l2, l3, l4, l5, l6, l7);
  TRANSPOSE8x8_SH_SH(l0, l1, l2, l3, l4, l5, l6, l7,
                     r0, r1, r2, r3, r4, r5, r6, r7);
  VP9_FDCT_POSTPROC_2V_NEG_H(r0, r1);
  VP9_FDCT_POSTPROC_2V_NEG_H(r2, r3);
  VP9_FDCT_POSTPROC_2V_NEG_H(r4, r5);
  VP9_FDCT_POSTPROC_2V_NEG_H(r6, r7);
  ST_SH8(r0, r1, r2, r3, r4, r5, r6, r7, out, 8);
  out += 64;

  LD_SH8(input + 8, 16, l8, l9, l10, l11, l12, l13, l14, l15);
  TRANSPOSE8x8_SH_SH(l8, l9, l10, l11, l12, l13, l14, l15,
                     r8, r9, r10, r11, r12, r13, r14, r15);
  VP9_FDCT_POSTPROC_2V_NEG_H(r8, r9);
  VP9_FDCT_POSTPROC_2V_NEG_H(r10, r11);
  VP9_FDCT_POSTPROC_2V_NEG_H(r12, r13);
  VP9_FDCT_POSTPROC_2V_NEG_H(r14, r15);
  ST_SH8(r8, r9, r10, r11, r12, r13, r14, r15, out, 8);
}

static void fadst16_rows_step1_msa(int16_t *input, const int32_t *const0,
                                   int16_t *int_buf) {
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 tp0, tp1, tp2, tp3, g0, g1, g2, g3, g8, g9, g10, g11, h0, h1, h2, h3;
  v4i32 k0, k1, k2, k3;

  /* load input data */
  r0 = LD_SH(input);
  r7 = LD_SH(input + 7 * 8);
  r8 = LD_SH(input + 8 * 8);
  r15 = LD_SH(input + 15 * 8);

  /* stage 1 */
  LD_SW2(const0, 4, k0, k1);
  LD_SW2(const0 + 4 * 2, 4, k2, k3);
  VP9_MADD_BF(r15, r0, r7, r8, k0, k1, k2, k3, g0, g1, g2, g3);

  r3 = LD_SH(input + 3 * 8);
  r4 = LD_SH(input + 4 * 8);
  r11 = LD_SH(input + 11 * 8);
  r12 = LD_SH(input + 12 * 8);

  LD_SW2(const0 + 4 * 4, 4, k0, k1);
  LD_SW2(const0 + 4 * 6, 4, k2, k3);
  VP9_MADD_BF(r11, r4, r3, r12, k0, k1, k2, k3, g8, g9, g10, g11);

  /* stage 2 */
  BUTTERFLY_4(g0, g2, g10, g8, tp0, tp2, tp3, tp1);
  ST_SH2(tp0, tp1, int_buf, 4 * 8);
  ST_SH2(tp2, tp3, int_buf + 8, 4 * 8);

  LD_SW2(const0 + 4 * 8, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 10);
  VP9_MADD_BF(g1, g3, g9, g11, k0, k1, k2, k0, h0, h1, h2, h3);
  ST_SH2(h0, h3, int_buf + 8 * 8, 4 * 8);
  ST_SH2(h1, h2, int_buf + 9 * 8, 4 * 8);

  r1 = LD_SH(input + 8);
  r6 = LD_SH(input + 6 * 8);
  r9 = LD_SH(input + 9 * 8);
  r14 = LD_SH(input + 14 * 8);

  LD_SW2(const0 + 4 * 11, 4, k0, k1);
  LD_SW2(const0 + 4 * 13, 4, k2, k3);
  VP9_MADD_BF(r9, r6, r1, r14, k0, k1, k2, k3, g0, g1, g2, g3);
  ST_SH2(g1, g3, int_buf + 3 * 8, 4 * 8);

  r2 = LD_SH(input + 2 * 8);
  r5 = LD_SH(input + 5 * 8);
  r10 = LD_SH(input + 10 * 8);
  r13 = LD_SH(input + 13 * 8);

  LD_SW2(const0 + 4 * 15, 4, k0, k1);
  LD_SW2(const0 + 4 * 17, 4, k2, k3);
  VP9_MADD_BF(r13, r2, r5, r10, k0, k1, k2, k3, h0, h1, h2, h3);
  ST_SH2(h1, h3, int_buf + 11 * 8, 4 * 8);
  BUTTERFLY_4(h0, h2, g2, g0, tp0, tp1, tp2, tp3);
  ST_SH4(tp0, tp1, tp2, tp3, int_buf + 2 * 8, 4 * 8);
}

static void fadst16_rows_step2_msa(int16_t *int_buf, const int32_t *const0,
                                   int16_t *out) {
  int16_t *out_ptr = out + 8;
  v8i16 tp0, tp1, tp2, tp3, g5, g7, g13, g15;
  v8i16 h0, h1, h2, h3, h4, h5, h6, h7, h10, h11;
  v8i16 out0, out1, out2, out3, out4, out5, out6, out7;
  v8i16 out8, out9, out10, out11, out12, out13, out14, out15;
  v4i32 k0, k1, k2, k3;

  g13 = LD_SH(int_buf + 3 * 8);
  g15 = LD_SH(int_buf + 7 * 8);
  g5 = LD_SH(int_buf + 11 * 8);
  g7 = LD_SH(int_buf + 15 * 8);

  LD_SW2(const0 + 4 * 19, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 21);
  VP9_MADD_BF(g7, g5, g15, g13, k0, k1, k2, k0, h4, h5, h6, h7);

  tp0 = LD_SH(int_buf + 4 * 8);
  tp1 = LD_SH(int_buf + 5 * 8);
  tp3 = LD_SH(int_buf + 10 * 8);
  tp2 = LD_SH(int_buf + 14 * 8);

  LD_SW2(const0 + 4 * 22, 4, k0, k1);
  k2 = LD_SW(const0 + 4 * 24);
  VP9_MADD_BF(tp0, tp1, tp2, tp3, k0, k1, k2, k0, out4, out6, out5, out7);
  out4 = -out4;
  ST_SH(out4, (out + 3 * 16));
  ST_SH(out5, (out_ptr + 4 * 16));

  h1 = LD_SH(int_buf + 9 * 8);
  h3 = LD_SH(int_buf + 12 * 8);
  VP9_MADD_BF(h1, h3, h5, h7, k0, k1, k2, k0, out12, out14, out13, out15);
  out13 = -out13;
  ST_SH(out12, (out + 2 * 16));
  ST_SH(out13, (out_ptr + 5 * 16));

  tp0 = LD_SH(int_buf);
  tp1 = LD_SH(int_buf + 8);
  tp2 = LD_SH(int_buf + 2 * 8);
  tp3 = LD_SH(int_buf + 6 * 8);

  BUTTERFLY_4(tp0, tp1, tp3, tp2, out0, out1, h11, h10);
  out1 = -out1;
  ST_SH(out0, (out));
  ST_SH(out1, (out_ptr + 7 * 16));

  h0 = LD_SH(int_buf + 8 * 8);
  h2 = LD_SH(int_buf + 13 * 8);
  BUTTERFLY_4(h0, h2, h6, h4, out8, out9, out11, out10);
  out8 = -out8;
  ST_SH(out8, (out + 16));
  ST_SH(out9, (out_ptr + 6 * 16));

  /* stage 4 */
  LD_SW2(const0 + 4 * 25, 4, k0, k1);
  LD_SW2(const0 + 4 * 27, 4, k2, k3);
  VP9_MADD_SHORT(h10, h11, k1, k2, out2, out3);
  ST_SH(out2, (out + 7 * 16));
  ST_SH(out3, (out_ptr));

  VP9_MADD_SHORT(out6, out7, k0, k3, out6, out7);
  ST_SH(out6, (out + 4 * 16));
  ST_SH(out7, (out_ptr + 3 * 16));

  VP9_MADD_SHORT(out10, out11, k0, k3, out10, out11);
  ST_SH(out10, (out + 6 * 16));
  ST_SH(out11, (out_ptr + 16));

  VP9_MADD_SHORT(out14, out15, k1, k2, out14, out15);
  ST_SH(out14, (out + 5 * 16));
  ST_SH(out15, (out_ptr + 2 * 16));
}

static void fadst16_transpose_msa(int16_t *input, int16_t *out) {
  v8i16 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
  v8i16 l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;

  /* load input data */
  LD_SH16(input, 8, l0, l8, l1, l9, l2, l10, l3, l11,
          l4, l12, l5, l13, l6, l14, l7, l15);
  TRANSPOSE8x8_SH_SH(l0, l1, l2, l3, l4, l5, l6, l7,
                     r0, r1, r2, r3, r4, r5, r6, r7);
  TRANSPOSE8x8_SH_SH(l8, l9, l10, l11, l12, l13, l14, l15,
                     r8, r9, r10, r11, r12, r13, r14, r15);
  ST_SH8(r0, r8, r1, r9, r2, r10, r3, r11, out, 8);
  ST_SH8(r4, r12, r5, r13, r6, r14, r7, r15, (out + 64), 8);
  out += 16 * 8;

  /* load input data */
  input += 128;
  LD_SH16(input, 8, l0, l8, l1, l9, l2, l10, l3, l11,
          l4, l12, l5, l13, l6, l14, l7, l15);
  TRANSPOSE8x8_SH_SH(l0, l1, l2, l3, l4, l5, l6, l7,
                     r0, r1, r2, r3, r4, r5, r6, r7);
  TRANSPOSE8x8_SH_SH(l8, l9, l10, l11, l12, l13, l14, l15,
                     r8, r9, r10, r11, r12, r13, r14, r15);
  ST_SH8(r0, r8, r1, r9, r2, r10, r3, r11, out, 8);
  ST_SH8(r4, r12, r5, r13, r6, r14, r7, r15, (out + 64), 8);
}

static void postproc_fdct16x8_1d_row(int16_t *intermediate, int16_t *output) {
  int16_t *temp = intermediate;
  int16_t *out = output;
  v8i16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  v8i16 in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11;
  v8i16 in12, in13, in14, in15;

  LD_SH8(temp, 16, in0, in1, in2, in3, in4, in5, in6, in7);
  temp = intermediate + 8;
  LD_SH8(temp, 16, in8, in9, in10, in11, in12, in13, in14, in15);
  TRANSPOSE8x8_SH_SH(in0, in1, in2, in3, in4, in5, in6, in7,
                     in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(in8, in9, in10, in11, in12, in13, in14, in15,
                     in8, in9, in10, in11, in12, in13, in14, in15);
  VP9_FDCT_POSTPROC_2V_NEG_H(in0, in1);
  VP9_FDCT_POSTPROC_2V_NEG_H(in2, in3);
  VP9_FDCT_POSTPROC_2V_NEG_H(in4, in5);
  VP9_FDCT_POSTPROC_2V_NEG_H(in6, in7);
  VP9_FDCT_POSTPROC_2V_NEG_H(in8, in9);
  VP9_FDCT_POSTPROC_2V_NEG_H(in10, in11);
  VP9_FDCT_POSTPROC_2V_NEG_H(in12, in13);
  VP9_FDCT_POSTPROC_2V_NEG_H(in14, in15);
  BUTTERFLY_16(in0, in1, in2, in3, in4, in5, in6, in7,
               in8, in9, in10, in11, in12, in13, in14, in15,
               tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,
               in8, in9, in10, in11, in12, in13, in14, in15);
  temp = intermediate;
  ST_SH8(in8, in9, in10, in11, in12, in13, in14, in15, temp, 16);
  VP9_FDCT8x16_EVEN(tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,
                    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7);
  temp = intermediate;
  LD_SH8(temp, 16, in8, in9, in10, in11, in12, in13, in14, in15);
  VP9_FDCT8x16_ODD(in8, in9, in10, in11, in12, in13, in14, in15,
                   in0, in1, in2, in3, in4, in5, in6, in7);
  TRANSPOSE8x8_SH_SH(tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3,
                     tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3);
  ST_SH8(tmp0, in0, tmp1, in1, tmp2, in2, tmp3, in3, out, 16);
  TRANSPOSE8x8_SH_SH(tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7,
                     tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7);
  out = output + 8;
  ST_SH8(tmp4, in4, tmp5, in5, tmp6, in6, tmp7, in7, out, 16);
}

void vp9_fht16x16_msa(const int16_t *input, int16_t *output,
                      int32_t stride, int32_t tx_type) {
  DECLARE_ALIGNED(32, int16_t, tmp[256]);
  DECLARE_ALIGNED(32, int16_t, trans_buf[256]);
  DECLARE_ALIGNED(32, int16_t, tmp_buf[128]);
  int32_t i;
  int16_t *ptmpbuf = &tmp_buf[0];
  int16_t *trans = &trans_buf[0];
  const int32_t const_arr[29 * 4] = {
    52707308, 52707308, 52707308, 52707308,
    -1072430300, -1072430300, -1072430300, -1072430300,
    795618043, 795618043, 795618043, 795618043,
    -721080468, -721080468, -721080468, -721080468,
    459094491, 459094491, 459094491, 459094491,
    -970646691, -970646691, -970646691, -970646691,
    1010963856, 1010963856, 1010963856, 1010963856,
    -361743294, -361743294, -361743294, -361743294,
    209469125, 209469125, 209469125, 209469125,
    -1053094788, -1053094788, -1053094788, -1053094788,
    1053160324, 1053160324, 1053160324, 1053160324,
    639644520, 639644520, 639644520, 639644520,
    -862444000, -862444000, -862444000, -862444000,
    1062144356, 1062144356, 1062144356, 1062144356,
    -157532337, -157532337, -157532337, -157532337,
    260914709, 260914709, 260914709, 260914709,
    -1041559667, -1041559667, -1041559667, -1041559667,
    920985831, 920985831, 920985831, 920985831,
    -551995675, -551995675, -551995675, -551995675,
    596522295, 596522295, 596522295, 596522295,
    892853362, 892853362, 892853362, 892853362,
    -892787826, -892787826, -892787826, -892787826,
    410925857, 410925857, 410925857, 410925857,
    -992012162, -992012162, -992012162, -992012162,
    992077698, 992077698, 992077698, 992077698,
    759246145, 759246145, 759246145, 759246145,
    -759180609, -759180609, -759180609, -759180609,
    -759222975, -759222975, -759222975, -759222975,
    759288511, 759288511, 759288511, 759288511 };

  switch (tx_type) {
    case DCT_DCT:
      /* column transform */
      for (i = 0; i < 2; ++i) {
        fdct8x16_1d_column(input + 8 * i, tmp + 8 * i, stride);
      }

      /* row transform */
      for (i = 0; i < 2; ++i) {
        fdct16x8_1d_row(tmp + (128 * i), output + (128 * i));
      }
      break;
    case ADST_DCT:
      /* column transform */
      for (i = 0; i < 2; ++i) {
        fadst16_cols_step1_msa(input + (i << 3), stride, const_arr, ptmpbuf);
        fadst16_cols_step2_msa(ptmpbuf, const_arr, tmp + (i << 3));
      }

      /* row transform */
      for (i = 0; i < 2; ++i) {
        postproc_fdct16x8_1d_row(tmp + (128 * i), output + (128 * i));
      }
      break;
    case DCT_ADST:
      /* column transform */
      for (i = 0; i < 2; ++i) {
        fdct8x16_1d_column(input + 8 * i, tmp + 8 * i, stride);
      }

      fadst16_transpose_postproc_msa(tmp, trans);

      /* row transform */
      for (i = 0; i < 2; ++i) {
        fadst16_rows_step1_msa(trans + (i << 7), const_arr, ptmpbuf);
        fadst16_rows_step2_msa(ptmpbuf, const_arr, tmp + (i << 7));
      }

      fadst16_transpose_msa(tmp, output);
      break;
    case ADST_ADST:
      /* column transform */
      for (i = 0; i < 2; ++i) {
        fadst16_cols_step1_msa(input + (i << 3), stride, const_arr, ptmpbuf);
        fadst16_cols_step2_msa(ptmpbuf, const_arr, tmp + (i << 3));
      }

      fadst16_transpose_postproc_msa(tmp, trans);

      /* row transform */
      for (i = 0; i < 2; ++i) {
        fadst16_rows_step1_msa(trans + (i << 7), const_arr, ptmpbuf);
        fadst16_rows_step2_msa(ptmpbuf, const_arr, tmp + (i << 7));
      }

      fadst16_transpose_msa(tmp, output);
      break;
    default:
      assert(0);
      break;
  }
}
