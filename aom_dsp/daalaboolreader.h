/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_DSP_DAALABOOLREADER_H_
#define AOM_DSP_DAALABOOLREADER_H_

#include "aom/aom_integer.h"
#include "aom_dsp/entdec.h"
#include "aom_dsp/prob.h"
#if CONFIG_ACCOUNTING
#include "av1/common/accounting.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct daala_reader {
  const uint8_t *buffer;
  const uint8_t *buffer_end;
  od_ec_dec ec;
#if CONFIG_ACCOUNTING
  Accounting *accounting;
#endif
};

typedef struct daala_reader daala_reader;

int aom_daala_reader_init(daala_reader *r, const uint8_t *buffer, int size);
const uint8_t *aom_daala_reader_find_end(daala_reader *r);
uint32_t aom_daala_reader_tell(const daala_reader *r);
uint32_t aom_daala_reader_tell_frac(const daala_reader *r);

static INLINE int aom_daala_read(daala_reader *r, int prob) {
  if (prob == 128) {
    return od_ec_dec_bits(&r->ec, 1, "aom_bits");
  } else {
    int p = ((prob << 15) + (256 - prob)) >> 8;
    return od_ec_decode_bool_q15(&r->ec, p);
  }
}

static INLINE int aom_daala_read_bit(daala_reader *r) {
  return aom_daala_read(r, 128);
}

static INLINE int aom_daala_reader_has_error(daala_reader *r) {
  return r->ec.error;
}

static INLINE int daala_read_tree_bits(daala_reader *r,
                                       const aom_tree_index *tree,
                                       const aom_prob *probs) {
  aom_tree_index i = 0;
  do {
    aom_cdf_prob cdf[16];
    aom_tree_index index[16];
    int path[16];
    int dist[16];
    int nsymbs;
    int symb;
    nsymbs = tree_to_cdf(tree, probs, i, cdf, index, path, dist);
    symb = od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs);
    OD_ASSERT(symb >= 0 && symb < nsymbs);
    i = index[symb];
  } while (i > 0);
  return -i;
}

static INLINE int daala_read_symbol(daala_reader *r, const aom_cdf_prob *cdf,
                                    int nsymbs) {
  return od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
