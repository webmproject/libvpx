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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#include "aom_dsp/entdec.h"
#include "av1/common/pvq.h"
#include "pvq_decoder.h"

#if OD_ACCOUNTING
# define od_decode_pvq_split(ec, adapt, sum, ctx, str) od_decode_pvq_split_(ec, adapt, sum, ctx, str)
#else
# define od_decode_pvq_split(ec, adapt, sum, ctx, str) od_decode_pvq_split_(ec, adapt, sum, ctx)
#endif

static int od_decode_pvq_split_(od_ec_dec *ec, od_pvq_codeword_ctx *adapt,
 int sum, int ctx OD_ACC_STR) {
  int shift;
  int count;
  int msbs;
  int fctx;
  count = 0;
  if (sum == 0) return 0;
  shift = OD_MAXI(0, OD_ILOG(sum) - 3);
  fctx = 7*ctx + (sum >> shift) - 1;
  msbs = od_decode_cdf_adapt(ec, adapt->pvq_split_cdf[fctx],
   (sum >> shift) + 1, adapt->pvq_split_increment, acc_str);
  if (shift) count = od_ec_dec_bits(ec, shift, acc_str);
  count += msbs << shift;
  if (count > sum) {
    count = sum;
    ec->error = 1;
  }
  return count;
}

void od_decode_band_pvq_splits(od_ec_dec *ec, od_pvq_codeword_ctx *adapt,
 od_coeff *y, int n, int k, int level) {
  int mid;
  int count_right;
  if (n == 1) {
    y[0] = k;
  }
  else if (k == 0) {
    OD_CLEAR(y, n);
  }
  else if (k == 1 && n <= 16) {
    int cdf_id;
    int pos;
    cdf_id = od_pvq_k1_ctx(n, level == 0);
    OD_CLEAR(y, n);
    pos = od_decode_cdf_adapt(ec, adapt->pvq_k1_cdf[cdf_id], n,
     adapt->pvq_k1_increment, "pvq:k1");
    y[pos] = 1;
  }
  else {
    mid = n >> 1;
    count_right = od_decode_pvq_split(ec, adapt, k, od_pvq_size_ctx(n),
     "pvq:split");
    od_decode_band_pvq_splits(ec, adapt, y, mid, k - count_right, level + 1);
    od_decode_band_pvq_splits(ec, adapt, y + mid, n - mid, count_right,
     level + 1);
  }
}

/** Decodes the tail of a Laplace-distributed variable, i.e. it doesn't
 * do anything special for the zero case.
 *
 * @param [dec] range decoder
 * @param [decay] decay factor of the distribution, i.e. pdf ~= decay^x
 * @param [max] maximum possible value of x (used to truncate the pdf)
 *
 * @retval decoded variable x
 */
int od_laplace_decode_special_(od_ec_dec *dec, unsigned decay, int max OD_ACC_STR) {
  int pos;
  int shift;
  int xs;
  int ms;
  int sym;
  const uint16_t *cdf;
  shift = 0;
  if (max == 0) return 0;
  /* We don't want a large decay value because that would require too many
     symbols. However, it's OK if the max is below 15. */
  while (((max >> shift) >= 15 || max == -1) && decay > 235) {
    decay = (decay*decay + 128) >> 8;
    shift++;
  }
  decay = OD_MINI(decay, 254);
  decay = OD_MAXI(decay, 2);
  ms = max >> shift;
  cdf = EXP_CDF_TABLE[(decay + 1) >> 1];
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "decay = %d\n", decay));
  xs = 0;
  do {
    sym = OD_MINI(xs, 15);
    {
      int i;
      OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d %d %d", xs, shift, sym, max));
      for (i = 0; i < 16; i++) {
        OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", cdf[i]));
      }
      OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "\n"));
    }
    if (ms > 0 && ms < 15) {
      /* Simple way of truncating the pdf when we have a bound. */
      sym = od_ec_decode_cdf_unscaled(dec, cdf, ms + 1);
    }
    else sym = od_ec_decode_cdf_q15(dec, cdf, 16);
    xs += sym;
    ms -= 15;
  }
  while (sym >= 15 && ms != 0);
  if (shift) pos = (xs << shift) + od_ec_dec_bits(dec, shift, acc_str);
  else pos = xs;
  OD_ASSERT(pos >> shift <= max >> shift || max == -1);
  if (max != -1 && pos > max) {
    pos = max;
    dec->error = 1;
  }
  OD_ASSERT(pos <= max || max == -1);
  return pos;
}

/** Decodes a Laplace-distributed variable for use in PVQ.
 *
 * @param [in,out] dec  range decoder
 * @param [in]     ExQ8 expectation of the absolute value of x
 * @param [in]     K    maximum value of |x|
 *
 * @retval decoded variable (including sign)
 */
int od_laplace_decode_(od_ec_dec *dec, unsigned ex_q8, int k OD_ACC_STR) {
  int j;
  int shift;
  uint16_t cdf[16];
  int sym;
  int lsb;
  int decay;
  int offset;
  lsb = 0;
  /* Shift down x if expectation is too high. */
  shift = OD_ILOG(ex_q8) - 11;
  if (shift < 0) shift = 0;
  /* Apply the shift with rounding to Ex, K and xs. */
  ex_q8 = (ex_q8 + (1 << shift >> 1)) >> shift;
  k = (k + (1 << shift >> 1)) >> shift;
  decay = OD_MINI(254, OD_DIVU(256*ex_q8, (ex_q8 + 256)));
  offset = LAPLACE_OFFSET[(decay + 1) >> 1];
  for (j = 0; j < 16; j++) {
    cdf[j] = EXP_CDF_TABLE[(decay + 1) >> 1][j] - offset;
  }
  /* Simple way of truncating the pdf when we have a bound */
  if (k == 0) sym = 0;
  else sym = od_ec_decode_cdf_unscaled(dec, cdf, OD_MINI(k + 1, 16));
  if (shift) {
    int special;
    /* Because of the rounding, there's only half the number of possibilities
       for xs=0 */
    special = (sym == 0);
    if (shift - special > 0) lsb = od_ec_dec_bits(dec, shift - special, acc_str);
    lsb -= (!special << (shift - 1));
  }
  /* Handle the exponentially-decaying tail of the distribution */
  if (sym == 15) sym += laplace_decode_special(dec, decay, k - 15, acc_str);
  return (sym << shift) + lsb;
}

#if OD_ACCOUNTING
# define laplace_decode_vector_delta(dec, y, n, k, curr, means, str) laplace_decode_vector_delta_(dec, y, n, k, curr, means, str)
#else
# define laplace_decode_vector_delta(dec, y, n, k, curr, means, str) laplace_decode_vector_delta_(dec, y, n, k, curr, means)
#endif

static void laplace_decode_vector_delta_(od_ec_dec *dec, od_coeff *y, int n, int k,
                                        int32_t *curr, const int32_t *means
                                        OD_ACC_STR) {
  int i;
  int prev;
  int sum_ex;
  int sum_c;
  int coef;
  int pos;
  int k0;
  int sign;
  int first;
  int k_left;
  prev = 0;
  sum_ex = 0;
  sum_c = 0;
  coef = 256*means[OD_ADAPT_COUNT_Q8]/
   (1 + means[OD_ADAPT_COUNT_EX_Q8]);
  pos = 0;
  sign = 0;
  first = 1;
  k_left = k;
  for (i = 0; i < n; i++) y[i] = 0;
  k0 = k_left;
  coef = OD_MAXI(coef, 1);
  for (i = 0; i < k0; i++) {
    int count;
    if (first) {
      int decay;
      int ex = coef*(n - prev)/k_left;
      if (ex > 65280) decay = 255;
      else {
        decay = OD_MINI(255,
         (int)((256*ex/(ex + 256) + (ex>>5)*ex/((n + 1)*(n - 1)*(n - 1)))));
      }
      /*Update mean position.*/
      count = laplace_decode_special(dec, decay, n - 1, acc_str);
      first = 0;
    }
    else count = laplace_decode(dec, coef*(n - prev)/k_left, n - prev - 1, acc_str);
    sum_ex += 256*(n - prev);
    sum_c += count*k_left;
    pos += count;
    OD_ASSERT(pos < n);
    if (y[pos] == 0)
      sign = od_ec_dec_bits(dec, 1, acc_str);
    y[pos] += sign ? -1 : 1;
    prev = pos;
    k_left--;
    if (k_left == 0) break;
  }
  if (k > 0) {
    curr[OD_ADAPT_COUNT_Q8] = 256*sum_c;
    curr[OD_ADAPT_COUNT_EX_Q8] = sum_ex;
  }
  else {
    curr[OD_ADAPT_COUNT_Q8] = -1;
    curr[OD_ADAPT_COUNT_EX_Q8] = 0;
  }
  curr[OD_ADAPT_K_Q8] = 0;
  curr[OD_ADAPT_SUM_EX_Q8] = 0;
}

/** Decodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] dec range decoder
 * @param [in]     y     decoded vector
 * @param [in]     N     dimension of the vector
 * @param [in]     K     sum of the absolute value of components of y
 * @param [out]    curr  Adaptation context output, may alias means.
 * @param [in]     means Adaptation context input.
 */
void od_laplace_decode_vector_(od_ec_dec *dec, od_coeff *y, int n, int k,
                           int32_t *curr, const int32_t *means OD_ACC_STR) {
  int i;
  int sum_ex;
  int kn;
  int exp_q8;
  int mean_k_q8;
  int mean_sum_ex_q8;
  int ran_delta;
  ran_delta = 0;
  if (k <= 1) {
    laplace_decode_vector_delta(dec, y, n, k, curr, means, acc_str);
    return;
  }
  if (k == 0) {
    curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
    curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
    curr[OD_ADAPT_K_Q8] = 0;
    curr[OD_ADAPT_SUM_EX_Q8] = 0;
    for (i = 0; i < n; i++) y[i] = 0;
    return;
  }
  sum_ex = 0;
  kn = k;
  /* Estimates the factor relating pulses_left and positions_left to E(|x|).*/
  mean_k_q8 = means[OD_ADAPT_K_Q8];
  mean_sum_ex_q8 = means[OD_ADAPT_SUM_EX_Q8];
  if (mean_k_q8 < 1 << 23) exp_q8 = 256*mean_k_q8/(1 + mean_sum_ex_q8);
  else exp_q8 = mean_k_q8/(1 + (mean_sum_ex_q8 >> 8));
  for (i = 0; i < n; i++) {
    int ex;
    int x;
    if (kn == 0) break;
    if (kn <= 1 && i != n - 1) {
      laplace_decode_vector_delta(dec, y + i, n - i, kn, curr, means, acc_str);
      ran_delta = 1;
      i = n;
      break;
    }
    /* Expected value of x (round-to-nearest) is
       expQ8*pulses_left/positions_left. */
    ex = (2*exp_q8*kn + (n - i))/(2*(n - i));
    if (ex > kn*256) ex = kn*256;
    sum_ex += (2*256*kn + (n - i))/(2*(n - i));
    /* No need to encode the magnitude for the last bin. */
    if (i != n - 1) x = laplace_decode(dec, ex, kn, acc_str);
    else x = kn;
    if (x != 0) {
      if (od_ec_dec_bits(dec, 1, acc_str)) x = -x;
    }
    y[i] = x;
    kn -= abs(x);
  }
  /* Adapting the estimates for expQ8. */
  if (!ran_delta) {
    curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
    curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
  }
  curr[OD_ADAPT_K_Q8] = k - kn;
  curr[OD_ADAPT_SUM_EX_Q8] = sum_ex;
  for (; i < n; i++) y[i] = 0;
}
