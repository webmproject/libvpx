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

#if !defined(_generic_code_H)
# define _generic_code_H

# include "aom_dsp/entdec.h"
# include "aom_dsp/entenc.h"

# define GENERIC_TABLES 12

#if OD_ACCOUNTING
# define generic_decode(dec, model, max, ex_q16, integration, str) generic_decode_(dec, model, max, ex_q16, integration, str)
# define od_decode_cdf_adapt_q15(ec, cdf, n, count, rate, str) od_decode_cdf_adapt_q15_(ec, cdf, n, count, rate, str)
# define od_decode_cdf_adapt(ec, cdf, n, increment, str) od_decode_cdf_adapt_(ec, cdf, n, increment, str)
#else
# define generic_decode(dec, model, max, ex_q16, integration, str) generic_decode_(dec, model, max, ex_q16, integration)
# define od_decode_cdf_adapt_q15(ec, cdf, n, count, rate, str) od_decode_cdf_adapt_q15_(ec, cdf, n, count, rate)
# define od_decode_cdf_adapt(ec, cdf, n, increment, str) od_decode_cdf_adapt_(ec, cdf, n, increment)
#endif

typedef struct {
  /** cdf for multiple expectations of x */
  uint16_t cdf[GENERIC_TABLES][16];
  /** Frequency increment for learning the cdfs */
  int increment;
} generic_encoder;

#define OD_IIR_DIADIC(y, x, shift) ((y) += ((x) - (y)) >> (shift))

void generic_model_init(generic_encoder *model);

#define OD_CDFS_INIT(cdf, val) od_cdf_init(&cdf[0][0],\
 sizeof(cdf)/sizeof(cdf[0]), sizeof(cdf[0])/sizeof(cdf[0][0]), val, val)

#define OD_CDFS_INIT_FIRST(cdf, val, first) od_cdf_init(&cdf[0][0],\
 sizeof(cdf)/sizeof(cdf[0]), sizeof(cdf[0])/sizeof(cdf[0][0]), val, first)

#define OD_SINGLE_CDF_INIT(cdf, val) od_cdf_init(cdf,\
 1, sizeof(cdf)/sizeof(cdf[0]), val, val)

#define OD_SINGLE_CDF_INIT_FIRST(cdf, val, first) od_cdf_init(cdf,\
 1, sizeof(cdf)/sizeof(cdf[0]), val, first)

void od_cdf_init(uint16_t *cdf, int ncdfs, int nsyms, int val, int first);

void od_cdf_adapt_q15(int val, uint16_t *cdf, int n, int *count, int rate);

void od_encode_cdf_adapt_q15(od_ec_enc *ec, int val, uint16_t *cdf, int n,
 int *count, int rate);

void od_encode_cdf_adapt(od_ec_enc *ec, int val, uint16_t *cdf, int n,
 int increment);

int od_decode_cdf_adapt_(od_ec_dec *ec, uint16_t *cdf, int n,
 int increment OD_ACC_STR);

void generic_encode(od_ec_enc *enc, generic_encoder *model, int x, int max,
 int *ex_q16, int integration);
double generic_encode_cost(generic_encoder *model, int x, int max,
 int *ex_q16);

double od_encode_cdf_cost(int val, uint16_t *cdf, int n);

int od_decode_cdf_adapt_q15_(od_ec_dec *ec, uint16_t *cdf, int n,
 int *count, int rate OD_ACC_STR);

int generic_decode_(od_ec_dec *dec, generic_encoder *model, int max,
 int *ex_q16, int integration OD_ACC_STR);

int log_ex(int ex_q16);

void generic_model_update(generic_encoder *model, int *ex_q16, int x, int xs,
 int id, int integration);

#endif
