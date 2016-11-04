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

#include "generic_code.h"

void od_cdf_init(uint16_t *cdf, int ncdfs, int nsyms, int val, int first) {
  int i;
  int j;
  for (i = 0; i < ncdfs; i++) {
    for (j = 0; j < nsyms; j++) {
      cdf[i*nsyms + j] = val*j + first;
    }
  }
}

/** Adapts a Q15 cdf after encoding/decoding a symbol. */
void od_cdf_adapt_q15(int val, uint16_t *cdf, int n, int *count, int rate) {
  int i;
  *count = OD_MINI(*count + 1, 1 << rate);
  OD_ASSERT(cdf[n - 1] == 32768);
  if (*count >= 1 << rate) {
    /* Steady-state adaptation based on a simple IIR with dyadic rate. */
    for (i = 0; i < n; i++) {
      int tmp;
      /* When (i < val), we want the adjustment ((cdf[i] - tmp) >> rate) to be
         positive so long as (cdf[i] > i + 1), and 0 when (cdf[i] == i + 1),
         to ensure we don't drive any probabilities to 0. Replacing cdf[i] with
         (i + 2) and solving ((i + 2 - tmp) >> rate == 1) for tmp produces
         tmp == i + 2 - (1 << rate). Using this value of tmp with
         cdf[i] == i + 1 instead gives an adjustment of 0 as desired.

         When (i >= val), we want ((cdf[i] - tmp) >> rate) to be negative so
         long as cdf[i] < 32768 - (n - 1 - i), and 0 when
         cdf[i] == 32768 - (n - 1 - i), again to ensure we don't drive any
         probabilities to 0. Since right-shifting any negative value is still
         negative, we can solve (32768 - (n - 1 - i) - tmp == 0) for tmp,
         producing tmp = 32769 - n + i. Using this value of tmp with smaller
         values of cdf[i] instead gives negative adjustments, as desired.

         Combining the two cases gives the expression below. These could be
         stored in a lookup table indexed by n and rate to avoid the
         arithmetic. */
      tmp = 2 - (1<<rate) + i + (32767 + (1<<rate) - n)*(i >= val);
      cdf[i] -= (cdf[i] - tmp) >> rate;
    }
  }
  else {
    int alpha;
    /* Initial adaptation for the first symbols. The adaptation rate is
       computed to be equivalent to what od_{en,de}code_cdf_adapt() does
       when the initial cdf is set to increment/4. */
    alpha = 4*32768/(n + 4**count);
    for (i = 0; i < n; i++) {
      int tmp;
      tmp = (32768 - n)*(i >= val) + i + 1;
      cdf[i] -= ((cdf[i] - tmp)*alpha) >> 15;
    }
  }
  OD_ASSERT(cdf[n - 1] == 32768);
}

/** Initializes the cdfs and freq counts for a model.
 *
 * @param [out] model model being initialized
 */
void generic_model_init(generic_encoder *model) {
  int i;
  int j;
  model->increment = 64;
  for (i = 0; i < GENERIC_TABLES; i++) {
    for (j = 0; j < 16; j++) {
      /* Do flat initialization equivalent to a single symbol in each bin. */
      model->cdf[i][j] = (j + 1) * model->increment;
    }
  }
}

/** Takes the base-2 log of E(x) in Q1.
 *
 * @param [in] ExQ16 expectation of x in Q16
 *
 * @retval 2*log2(ExQ16/2^16)
 */
int log_ex(int ex_q16) {
  int lg;
  int lg_q1;
  int odd;
  lg = OD_ILOG(ex_q16);
  if (lg < 15) {
    odd = ex_q16*ex_q16 > 2 << 2*lg;
  }
  else {
    int tmp;
    tmp = ex_q16 >> (lg - 8);
    odd = tmp*tmp > (1 << 15);
  }
  lg_q1 = OD_MAXI(0, 2*lg - 33 + odd);
  return lg_q1;
}

/** Updates the probability model based on the encoded/decoded value
 *
 * @param [in,out] model generic prob model
 * @param [in,out] ExQ16 expectation of x
 * @param [in]     x     variable encoded/decoded (used for ExQ16)
 * @param [in]     xs    variable x after shift (used for the model)
 * @param [in]     id    id of the icdf to adapt
 * @param [in]     integration integration period of ExQ16 (leaky average over
 * 1<<integration samples)
 */
void generic_model_update(generic_encoder *model, int *ex_q16, int x, int xs,
 int id, int integration) {
  int i;
  int xenc;
  uint16_t *cdf;
  cdf = model->cdf[id];
  /* Renormalize if we cannot add increment */
  if (cdf[15] + model->increment > 32767) {
    for (i = 0; i < 16; i++) {
      /* Second term ensures that the pdf is non-null */
      cdf[i] = (cdf[i] >> 1) + i + 1;
    }
  }
  /* Update freq count */
  xenc = OD_MINI(15, xs);
  /* This can be easily vectorized */
  for (i = xenc; i < 16; i++) cdf[i] += model->increment;
  /* We could have saturated ExQ16 directly, but this is safe and simpler */
  x = OD_MINI(x, 32767);
  OD_IIR_DIADIC(*ex_q16, x << 16, integration);
}
