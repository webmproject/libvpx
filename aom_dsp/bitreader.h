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

#include "./aom_config.h"
#include "aom/aomdx.h"
#include "aom/aom_integer.h"
#include "aom_dsp/dkboolreader.h"
#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct aom_dk_reader aom_reader;

static INLINE int aom_reader_init(aom_reader *r, const uint8_t *buffer,
                                  size_t size, aom_decrypt_cb decrypt_cb,
                                  void *decrypt_state) {
  return aom_dk_reader_init(r, buffer, size, decrypt_cb, decrypt_state);
}

static INLINE const uint8_t *aom_reader_find_end(aom_reader *r) {
  return aom_dk_reader_find_end(r);
}

static INLINE int aom_reader_has_error(aom_reader *r) {
  return aom_dk_reader_has_error(r);
}

static INLINE int aom_read(aom_reader *r, int prob) {
  return aom_dk_read(r, prob);
}

static INLINE int aom_read_bit(aom_reader *r) { return aom_dk_read_bit(r); }

static INLINE int aom_read_literal(aom_reader *r, int bits) {
  return aom_dk_read_literal(r, bits);
}

static INLINE int aom_read_tree(aom_reader *r, const aom_tree_index *tree,
                                const aom_prob *probs) {
  return aom_dk_read_tree(r, tree, probs);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITREADER_H_
