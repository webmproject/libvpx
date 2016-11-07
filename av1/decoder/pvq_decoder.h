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

#if !defined(_pvq_decoder_H)
# define _pvq_decoder_H (1)
# include "aom_dsp/entdec.h"
# include "av1/common/pvq.h"
# include "av1/decoder/decint.h"

void od_decode_band_pvq_splits(od_ec_dec *ec, od_pvq_codeword_ctx *adapt,
 od_coeff *y, int n, int k, int level);

#if OD_ACCOUNTING
# define laplace_decode_special(dec, decay, max, str) od_laplace_decode_special_(dec, decay, max, str)
# define laplace_decode(dec, ex_q8, k, str) od_laplace_decode_(dec, ex_q8, k, str)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) od_laplace_decode_vector_(dec, y, n, k, curr, means, str)
#else
# define laplace_decode_special(dec, decay, max, str) od_laplace_decode_special_(dec, decay, max)
# define laplace_decode(dec, ex_q8, k, str) od_laplace_decode_(dec, ex_q8, k)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) od_laplace_decode_vector_(dec, y, n, k, curr, means)
#endif

int od_laplace_decode_special_(od_ec_dec *dec, unsigned decay, int max OD_ACC_STR);
int od_laplace_decode_(od_ec_dec *dec, unsigned ex_q8, int k OD_ACC_STR);
void od_laplace_decode_vector_(od_ec_dec *dec, od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means
                                  OD_ACC_STR);


void od_pvq_decode(daala_dec_ctx *dec, od_coeff *ref, od_coeff *out, int q0,
 int pli, int bs, const od_val16 *beta, int robust, int is_keyframe,
 unsigned int *flags, int block_skip, const int16_t *qm,
 const int16_t *qm_inv);

#endif
