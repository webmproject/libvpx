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

#ifndef AOM_DSP_BITREADER_H_
#define AOM_DSP_BITREADER_H_

#include <assert.h>

#include "./aom_config.h"
#include "aom/aomdx.h"
#include "aom/aom_integer.h"
#if CONFIG_ANS
#include "aom_dsp/ans.h"
#else
#include "aom_dsp/dkboolreader.h"
#endif
#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_ANS
typedef struct AnsDecoder aom_reader;
#else
typedef struct aom_dk_reader aom_reader;
#endif

static INLINE int aom_reader_init(aom_reader *r, const uint8_t *buffer,
                                  size_t size, aom_decrypt_cb decrypt_cb,
                                  void *decrypt_state) {
#if CONFIG_ANS
  (void)decrypt_cb;
  (void)decrypt_state;
  assert(size <= INT_MAX);
  return ans_read_init(r, buffer, size);
#else
  return aom_dk_reader_init(r, buffer, size, decrypt_cb, decrypt_state);
#endif
}

static INLINE const uint8_t *aom_reader_find_end(aom_reader *r) {
#if CONFIG_ANS
  (void)r;
  assert(0 && "Use the raw buffer size with ANS");
  return NULL;
#else
  return aom_dk_reader_find_end(r);
#endif
}

static INLINE int aom_reader_has_error(aom_reader *r) {
#if CONFIG_ANS
  return ans_reader_has_error(r);
#else
  return aom_dk_reader_has_error(r);
#endif
}

static INLINE int aom_read(aom_reader *r, int prob) {
#if CONFIG_ANS
  return uabs_read(r, prob);
#else
  return aom_dk_read(r, prob);
#endif
}

static INLINE int aom_read_bit(aom_reader *r) {
#if CONFIG_ANS
  return uabs_read_bit(r);  // Non trivial optimization at half probability
#else
  return aom_read(r, 128);  // aom_prob_half
#endif
}

static INLINE int aom_read_literal(aom_reader *r, int bits) {
  int literal = 0, bit;

  for (bit = bits - 1; bit >= 0; bit--) literal |= aom_read_bit(r) << bit;

  return literal;
}

static INLINE int aom_read_tree(aom_reader *r, const aom_tree_index *tree,
                                const aom_prob *probs) {
  aom_tree_index i = 0;

  while ((i = tree[i + aom_read(r, probs[i >> 1])]) > 0) continue;

  return -i;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITREADER_H_
