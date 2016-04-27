/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_BUF_ANS_H_
#define VP10_ENCODER_BUF_ANS_H_
// Buffered forward ANS writer.
// Symbols are written to the writer in forward (decode) order and serialzed
// backwards due to ANS's stack like behavior.

#include <assert.h>
#include "./vpx_config.h"
#include "vpx/vpx_integer.h"
#include "vp10/common/ans.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

#define ANS_METHOD_UABS 0
#define ANS_METHOD_RANS 1

struct buffered_ans_symbol {
  uint8_t method;    // one of ANS_METHOD_UABS or ANS_METHOD_RANS
  // TODO(aconverse): Should be possible to write this interms of start for ABS
  AnsP10 val_start;  // Boolean value for ABS, start in symbol cycle for Rans
  AnsP10 prob;  // Probability of this symbol
};

struct BufAnsCoder {
  struct VP10Common *cm;
  struct buffered_ans_symbol *buf;
  int size;
  int offset;
};

void vp10_buf_ans_alloc(struct BufAnsCoder *c, struct VP10Common *cm,
                        int size_hint);

void vp10_buf_ans_free(struct BufAnsCoder *c);

void vp10_buf_ans_grow(struct BufAnsCoder *c);

static INLINE void buf_ans_write_reset(struct BufAnsCoder *const c) {
  c->offset = 0;
}

static INLINE void buf_uabs_write(struct BufAnsCoder *const c,
                                  uint8_t val, AnsP8 prob) {
  assert(c->offset <= c->size);
  if (c->offset == c->size) {
    vp10_buf_ans_grow(c);
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
    vp10_buf_ans_grow(c);
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

static INLINE void buf_uabs_write_literal(struct BufAnsCoder *c,
                                          int literal, int bits) {
  int bit;

  assert(bits < 31);
  for (bit = bits - 1; bit >= 0; bit--)
    buf_uabs_write_bit(c, 1 & (literal >> bit));
}
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // VP10_ENCODER_BUF_ANS_H_
