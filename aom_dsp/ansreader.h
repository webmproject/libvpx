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

#ifndef AOM_DSP_ANSREADER_H_
#define AOM_DSP_ANSREADER_H_
// A uABS and rANS decoder implementation of Asymmetric Numeral Systems
// http://arxiv.org/abs/1311.2540v2

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/prob.h"
#include "aom_dsp/ans.h"
#include "aom_ports/mem_ops.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

struct AnsDecoder {
  const uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

static INLINE int uabs_read(struct AnsDecoder *ans, AnsP8 p0) {
  AnsP8 p = ANS_P8_PRECISION - p0;
  int s;
  unsigned xp, sp;
  unsigned state = ans->state;
  while (state < L_BASE && ans->buf_offset > 0) {
    state = state * IO_BASE + ans->buf[--ans->buf_offset];
  }
  sp = state * p;
  xp = sp / ANS_P8_PRECISION;
  s = (sp & 0xFF) >= p0;
  if (s)
    ans->state = xp;
  else
    ans->state = state - xp;
  return s;
}

static INLINE int uabs_read_bit(struct AnsDecoder *ans) {
  int s;
  unsigned state = ans->state;
  while (state < L_BASE && ans->buf_offset > 0) {
    state = state * IO_BASE + ans->buf[--ans->buf_offset];
  }
  s = (int)(state & 1);
  ans->state = state >> 1;
  return s;
}

struct rans_dec_sym {
  uint8_t val;
  AnsP10 prob;
  AnsP10 cum_prob;  // not-inclusive
};

static INLINE void fetch_sym(struct rans_dec_sym *out, const rans_lut cdf,
                             AnsP10 rem) {
  int i = 0;
  // TODO(skal): if critical, could be a binary search.
  // Or, better, an O(1) alias-table.
  while (rem >= cdf[i]) {
    ++i;
  }
  out->val = i - 1;
  out->prob = (AnsP10)(cdf[i] - cdf[i - 1]);
  out->cum_prob = (AnsP10)cdf[i - 1];
}

static INLINE int rans_read(struct AnsDecoder *ans, const rans_lut tab) {
  unsigned rem;
  unsigned quo;
  struct rans_dec_sym sym;
  while (ans->state < L_BASE && ans->buf_offset > 0) {
    ans->state = ans->state * IO_BASE + ans->buf[--ans->buf_offset];
  }
  quo = ans->state / RANS_PRECISION;
  rem = ans->state % RANS_PRECISION;
  fetch_sym(&sym, tab, rem);
  ans->state = quo * sym.prob + rem - sym.cum_prob;
  return sym.val;
}

static INLINE int ans_read_init(struct AnsDecoder *const ans,
                                const uint8_t *const buf, int offset) {
  unsigned x;
  if (offset < 1) return 1;
  ans->buf = buf;
  x = buf[offset - 1] >> 6;
  if (x == 0) {
    ans->buf_offset = offset - 1;
    ans->state = buf[offset - 1] & 0x3F;
  } else if (x == 1) {
    if (offset < 2) return 1;
    ans->buf_offset = offset - 2;
    ans->state = mem_get_le16(buf + offset - 2) & 0x3FFF;
  } else if (x == 2) {
    if (offset < 3) return 1;
    ans->buf_offset = offset - 3;
    ans->state = mem_get_le24(buf + offset - 3) & 0x3FFFFF;
  } else {
    // x == 3 implies this byte is a superframe marker
    return 1;
  }
  ans->state += L_BASE;
  if (ans->state >= L_BASE * IO_BASE) return 1;
  return 0;
}

static INLINE int ans_read_end(struct AnsDecoder *const ans) {
  return ans->state == L_BASE;
}

static INLINE int ans_reader_has_error(const struct AnsDecoder *const ans) {
  return ans->state < L_BASE && ans->buf_offset == 0;
}
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // AOM_DSP_ANSREADER_H_
