/*
 * Copyright (c) 2001-2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/* clang-format off */

#if !defined(_pvq_encoder_H)
# define _pvq_encoder_H (1)
# include "aom_dsp/entenc.h"
# include "av1/common/blockd.h"
# include "av1/common/pvq.h"
# include "av1/encoder/encint.h"

#define PVQ_CHROMA_RD 1

void od_encode_band_pvq_splits(od_ec_enc *ec, od_pvq_codeword_ctx *adapt,
 const int *y, int n, int k, int level);

void od_laplace_encode_special(od_ec_enc *enc, int x, unsigned decay, int max);
void od_laplace_encode(od_ec_enc *enc, int x, int ex_q8, int k);
void od_laplace_encode_vector(od_ec_enc *enc, const od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means);

#if OD_SIGNAL_Q_SCALING
void od_encode_quantizer_scaling(daala_enc_ctx *enc, int q_scaling, int bx,
 int by, int skip);
#endif

void pvq_encode_partition(od_ec_enc *ec,
                                 int qg,
                                 int theta,
                                 int max_theta,
                                 const od_coeff *in,
                                 int n,
                                 int k,
                                 generic_encoder model[3],
                                 od_adapt_ctx *adapt,
                                 int *exg,
                                 int *ext,
                                 int nodesync,
                                 int cdf_ctx,
                                 int is_keyframe,
                                 int code_skip,
                                 int skip_rest,
                                 int encode_flip,
                                 int flip);

int od_pvq_encode(daala_enc_ctx *enc, od_coeff *ref, const od_coeff *in,
 od_coeff *out, int q_dc, int q_ac, int pli, int bs, const od_val16 *beta, int robust,
 int is_keyframe, int q_scaling, int bx, int by, const int16_t *qm,
 const int16_t *qm_inv, int speed, PVQ_INFO *pvq_info);

#endif
