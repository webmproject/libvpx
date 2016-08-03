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

#ifndef AOM_DSP_BUF_ANS_H_
#define AOM_DSP_BUF_ANS_H_
// Buffered forward ANS writer.
// Symbols are written to the writer in forward (decode) order and serialized
// backwards due to ANS's stack like behavior.

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/ans.h"
#include "aom_dsp/answriter.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

#define ANS_METHOD_UABS 0
#define ANS_METHOD_RANS 1

struct buffered_ans_symbol {
  unsigned int method : 1;  // one of ANS_METHOD_UABS or ANS_METHOD_RANS
  // TODO(aconverse): Should be possible to write this in terms of start for ABS
  unsigned int val_start : RANS_PROB_BITS;  // Boolean value for ABS
                                            // start in symbol cycle for Rans
  unsigned int prob : RANS_PROB_BITS;       // Probability of this symbol
};

struct BufAnsCoder {
  struct aom_internal_error_info *error;
  struct buffered_ans_symbol *buf;
  int size;
  int offset;
};

void aom_buf_ans_alloc(struct BufAnsCoder *c,
                       struct aom_internal_error_info *error, int size_hint);

void aom_buf_ans_free(struct BufAnsCoder *c);

void aom_buf_ans_grow(struct BufAnsCoder *c);

static INLINE void buf_ans_write_reset(struct BufAnsCoder *const c) {
  c->offset = 0;
}

static INLINE void buf_uabs_write(struct BufAnsCoder *const c, uint8_t val,
                                  AnsP8 prob) {
  assert(c->offset <= c->size);
  if (c->offset == c->size) {
    aom_buf_ans_grow(c);
  }
  c->buf[c->offset].method = ANS_METHOD_UABS;
  c->buf[c->offset].val_start = val;
  c->buf[c->offset].prob = prob;
  ++c->offset;
}

static INLINE void buf_rans_write(struct BufAnsCoder *const c,
                                  const struct rans_sym *const sym) {
  assert(c->offset <= c->size);
  if (c->offset == c->size) {
    aom_buf_ans_grow(c);
  }
  c->buf[c->offset].method = ANS_METHOD_RANS;
  c->buf[c->offset].val_start = sym->cum_prob;
  c->buf[c->offset].prob = sym->prob;
  ++c->offset;
}

static INLINE void buf_ans_flush(const struct BufAnsCoder *const c,
                                 struct AnsCoder *ans) {
  int offset;
  for (offset = c->offset - 1; offset >= 0; --offset) {
    if (c->buf[offset].method == ANS_METHOD_RANS) {
      struct rans_sym sym;
      sym.prob = c->buf[offset].prob;
      sym.cum_prob = c->buf[offset].val_start;
      rans_write(ans, &sym);
    } else {
      uabs_write(ans, (uint8_t)c->buf[offset].val_start,
                 (AnsP8)c->buf[offset].prob);
    }
  }
}

static INLINE void buf_uabs_write_bit(struct BufAnsCoder *c, int bit) {
  buf_uabs_write(c, bit, 128);
}

static INLINE void buf_uabs_write_literal(struct BufAnsCoder *c, int literal,
                                          int bits) {
  int bit;

  assert(bits < 31);
  for (bit = bits - 1; bit >= 0; bit--)
    buf_uabs_write_bit(c, 1 & (literal >> bit));
}
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // AOM_DSP_BUF_ANS_H_
