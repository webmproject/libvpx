/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_ANS_H_
#define VP10_COMMON_ANS_H_
// An implementation of Asymmetric Numeral Systems
// http://arxiv.org/abs/1311.2540v2

#include <assert.h>
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/prob.h"
#include "vpx_ports/mem_ops.h"

#define ANS_DIVIDE_BY_MULTIPLY 1
#if ANS_DIVIDE_BY_MULTIPLY
#include "vp10/common/divide.h"
#define ANS_DIVREM(quotient, remainder, dividend, divisor) \
  do { \
    quotient = fastdiv(dividend, divisor); \
    remainder = dividend - quotient * divisor; \
  } while (0)
#define ANS_DIV(dividend, divisor) \
  fastdiv(dividend, divisor)
#else
#define ANS_DIVREM(quotient, remainder, dividend, divisor) \
  do { \
    quotient = dividend / divisor; \
    remainder = dividend % divisor; \
  } while (0)
#define ANS_DIV(dividend, divisor) \
    ((dividend) / (divisor))
#endif

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

struct AnsCoder {
  uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

struct AnsDecoder {
  const uint8_t *buf;
  int buf_offset;
  uint32_t state;
};

typedef uint8_t AnsP8;
#define ans_p8_precision 256u
#define ans_p8_shift 8
typedef uint16_t AnsP10;
#define ans_p10_precision 1024u

#define rans_precision ans_p10_precision

#define l_base (ans_p10_precision * 4)  // l_base % precision must be 0
#define io_base 256
// Range I = { l_base, l_base + 1, ..., l_base * io_base - 1 }

static INLINE void ans_write_init(struct AnsCoder *const ans,
                                  uint8_t *const buf) {
  ans->buf = buf;
  ans->buf_offset = 0;
  ans->state = l_base;
}

static INLINE int ans_write_end(struct AnsCoder *const ans) {
  uint32_t state;
  assert(ans->state >= l_base);
  assert(ans->state < l_base * io_base);
  state = ans->state - l_base;
  if (state < (1 << 6)) {
    ans->buf[ans->buf_offset] = (0x00 << 6) + state;
    return ans->buf_offset + 1;
  } else if (state < (1 << 14)) {
    mem_put_le16(ans->buf + ans->buf_offset, (0x01 << 14) + state);
    return ans->buf_offset + 2;
  } else if (state < (1 << 22)) {
    mem_put_le24(ans->buf + ans->buf_offset, (0x02 << 22) + state);
    return ans->buf_offset + 3;
  } else {
    assert(0 && "State is too large to be serialized");
    return ans->buf_offset;
  }
}

// rABS with descending spread
// p or p0 takes the place of l_s from the paper
// ans_p8_precision is m
static INLINE void rabs_desc_write(struct AnsCoder *ans, int val, AnsP8 p0) {
  const AnsP8 p = ans_p8_precision - p0;
  const unsigned l_s = val ? p : p0;
  unsigned quot, rem;
  if (ans->state >= l_base / ans_p8_precision * io_base * l_s) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  ANS_DIVREM(quot, rem, ans->state, l_s);
  ans->state = quot * ans_p8_precision + rem + (val ? 0 : p);
}

#define ANS_IMPL1 0
#define UNPREDICTABLE(x) x
static INLINE int rabs_desc_read(struct AnsDecoder *ans, AnsP8 p0) {
  int val;
#if ANS_IMPL1
  unsigned l_s;
#else
  unsigned quot, rem, x, xn;
#endif
  const AnsP8 p = ans_p8_precision - p0;
  if (ans->state < l_base) {
    ans->state = ans->state * io_base + ans->buf[--ans->buf_offset];
  }
#if ANS_IMPL1
  val = ans->state % ans_p8_precision < p;
  l_s = val ? p : p0;
  ans->state = (ans->state / ans_p8_precision) * l_s +
               ans->state % ans_p8_precision - (!val * p);
#else
  x = ans->state;
  quot = x / ans_p8_precision;
  rem = x % ans_p8_precision;
  xn = quot * p;
  val = rem < p;
  if (UNPREDICTABLE(val)) {
    ans->state = xn + rem;
  } else {
    // ans->state = quot * p0 + rem - p;
    ans->state = x - xn - p;
  }
#endif
  return val;
}

// rABS with ascending spread
// p or p0 takes the place of l_s from the paper
// ans_p8_precision is m
static INLINE void rabs_asc_write(struct AnsCoder *ans, int val, AnsP8 p0) {
  const AnsP8 p = ans_p8_precision - p0;
  const unsigned l_s = val ? p : p0;
  unsigned quot, rem;
  if (ans->state >= l_base / ans_p8_precision * io_base * l_s) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  ANS_DIVREM(quot, rem, ans->state, l_s);
  ans->state = quot * ans_p8_precision + rem + (val ? p0 : 0);
}

static INLINE int rabs_asc_read(struct AnsDecoder *ans, AnsP8 p0) {
  int val;
#if ANS_IMPL1
  unsigned l_s;
#else
  unsigned quot, rem, x, xn;
#endif
  const AnsP8 p = ans_p8_precision - p0;
  if (ans->state < l_base) {
    ans->state = ans->state * io_base + ans->buf[--ans->buf_offset];
  }
#if ANS_IMPL1
  val = ans->state % ans_p8_precision < p;
  l_s = val ? p : p0;
  ans->state = (ans->state / ans_p8_precision) * l_s +
               ans->state % ans_p8_precision - (!val * p);
#else
  x = ans->state;
  quot = x / ans_p8_precision;
  rem = x % ans_p8_precision;
  xn = quot * p;
  val = rem >= p0;
  if (UNPREDICTABLE(val)) {
    ans->state = xn + rem - p0;
  } else {
    // ans->state = quot * p0 + rem - p0;
    ans->state = x - xn;
  }
#endif
  return val;
}

#define rabs_read rabs_desc_read
#define rabs_write rabs_desc_write

// uABS with normalization
static INLINE void uabs_write(struct AnsCoder *ans, int val, AnsP8 p0) {
  AnsP8 p = ans_p8_precision - p0;
  const unsigned l_s = val ? p : p0;
  while (ans->state >= l_base / ans_p8_precision * io_base * l_s) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  if (!val)
    ans->state = ANS_DIV(ans->state * ans_p8_precision, p0);
  else
    ans->state = ANS_DIV((ans->state + 1) * ans_p8_precision + p - 1, p) - 1;
}

static INLINE int uabs_read(struct AnsDecoder *ans, AnsP8 p0) {
  AnsP8 p = ans_p8_precision - p0;
  int s;
  // unsigned int xp1;
  unsigned xp, sp;
  unsigned state = ans->state;
  while (state < l_base && ans->buf_offset > 0) {
    state = state * io_base + ans->buf[--ans->buf_offset];
  }
  sp = state * p;
  // xp1 = (sp + p) / ans_p8_precision;
  xp = sp / ans_p8_precision;
  // s = xp1 - xp;
  s = (sp & 0xFF) >= p0;
  if (UNPREDICTABLE(s))
    ans->state = xp;
  else
    ans->state = state - xp;
  return s;
}

static INLINE int uabs_read_bit(struct AnsDecoder *ans) {
  int s;
  unsigned state = ans->state;
  while (state < l_base && ans->buf_offset > 0) {
    state = state * io_base + ans->buf[--ans->buf_offset];
  }
  s = (int)(state & 1);
  ans->state = state >> 1;
  return s;
}

static INLINE int uabs_read_literal(struct AnsDecoder *ans, int bits) {
  int literal = 0, bit;
  assert(bits < 31);

  // TODO(aconverse): Investigate ways to read/write literals faster,
  // e.g. 8-bit chunks.
  for (bit = bits - 1; bit >= 0; bit--)
    literal |= uabs_read_bit(ans) << bit;

  return literal;
}

// TODO(aconverse): Replace trees with tokensets.
static INLINE int uabs_read_tree(struct AnsDecoder *ans,
                                 const vpx_tree_index *tree,
                                 const AnsP8 *probs) {
  vpx_tree_index i = 0;

  while ((i = tree[i + uabs_read(ans, probs[i >> 1])]) > 0)
    continue;

  return -i;
}

struct rans_sym {
  AnsP10 prob;
  AnsP10 cum_prob;  // not-inclusive
};

struct rans_dec_sym {
  uint8_t val;
  AnsP10 prob;
  AnsP10 cum_prob;  // not-inclusive
};

// This is now just a boring cdf. It starts with an explicit zero.
// TODO(aconverse): Remove starting zero.
typedef uint16_t rans_dec_lut[16];

static INLINE void rans_build_cdf_from_pdf(const AnsP10 token_probs[],
                                           rans_dec_lut cdf_tab) {
  int i;
  cdf_tab[0] = 0;
  for (i = 1; cdf_tab[i - 1] < rans_precision; ++i) {
    cdf_tab[i] = cdf_tab[i - 1] + token_probs[i - 1];
  }
  assert(cdf_tab[i - 1] == rans_precision);
}

static INLINE int ans_find_largest(const AnsP10 *const pdf_tab,
                                   int num_syms) {
  int largest_idx = -1;
  int largest_p = -1;
  int i;
  for (i = 0; i < num_syms; ++i) {
    int p = pdf_tab[i];
    if (p > largest_p) {
      largest_p = p;
      largest_idx = i;
    }
  }
  return largest_idx;
}

static INLINE void rans_merge_prob8_pdf(AnsP10 *const out_pdf,
                                        const AnsP8 node_prob,
                                        const AnsP10 *const src_pdf,
                                        int in_syms) {
  int i;
  int adjustment = rans_precision;
  const int round_fact = ans_p8_precision >> 1;
  const AnsP8 p1 = ans_p8_precision - node_prob;
  const int out_syms = in_syms + 1;
  assert(src_pdf != out_pdf);

  out_pdf[0] = node_prob << (10 - 8);
  adjustment -= out_pdf[0];
  for (i = 0; i < in_syms; ++i) {
    int p = (p1 * src_pdf[i] + round_fact) >> ans_p8_shift;
    p = VPXMIN(p, (int)rans_precision - in_syms);
    p = VPXMAX(p, 1);
    out_pdf[i + 1] = p;
    adjustment -= p;
  }

  // Adjust probabilities so they sum to the total probability
  if (adjustment > 0) {
    i = ans_find_largest(out_pdf, out_syms);
    out_pdf[i] += adjustment;
  } else {
    while (adjustment < 0) {
      i = ans_find_largest(out_pdf, out_syms);
      --out_pdf[i];
      assert(out_pdf[i] > 0);
      adjustment++;
    }
  }
}

// rANS with normalization
// sym->prob takes the place of l_s from the paper
// ans_p10_precision is m
static INLINE void rans_write(struct AnsCoder *ans,
                              const struct rans_sym *const sym) {
  const AnsP10 p = sym->prob;
  while (ans->state >= l_base / rans_precision * io_base * p) {
    ans->buf[ans->buf_offset++] = ans->state % io_base;
    ans->state /= io_base;
  }
  ans->state =
      (ans->state / p) * rans_precision + ans->state % p + sym->cum_prob;
}

static INLINE void fetch_sym(struct rans_dec_sym *out, const rans_dec_lut cdf,
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

static INLINE int rans_read(struct AnsDecoder *ans,
                            const rans_dec_lut tab) {
  unsigned rem;
  unsigned quo;
  struct rans_dec_sym sym;
  while (ans->state < l_base && ans->buf_offset > 0) {
    ans->state = ans->state * io_base + ans->buf[--ans->buf_offset];
  }
  quo = ans->state / rans_precision;
  rem = ans->state % rans_precision;
  fetch_sym(&sym, tab, rem);
  ans->state = quo * sym.prob + rem - sym.cum_prob;
  return sym.val;
}

static INLINE int ans_read_init(struct AnsDecoder *const ans,
                                const uint8_t *const buf,
                                int offset) {
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
  ans->state += l_base;
  if (ans->state >= l_base * io_base)
    return 1;
  return 0;
}

static INLINE int ans_read_end(struct AnsDecoder *const ans) {
  return ans->state == l_base;
}

static INLINE int ans_reader_has_error(const struct AnsDecoder *const ans) {
  return ans->state < l_base && ans->buf_offset == 0;
}
#undef ANS_DIVREM
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // VP10_COMMON_ANS_H_
