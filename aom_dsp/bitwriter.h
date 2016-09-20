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

#ifndef AOM_DSP_BITWRITER_H_
#define AOM_DSP_BITWRITER_H_

#include <assert.h>
#include "./aom_config.h"
#if CONFIG_ANS
#include "aom_dsp/buf_ans.h"
#elif CONFIG_DAALA_EC
#include "aom_dsp/daalaboolwriter.h"
#else
#include "aom_dsp/dkboolwriter.h"
#endif
#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_ANS
typedef struct BufAnsCoder aom_writer;
#else
typedef struct aom_dk_writer aom_writer;
#endif

static INLINE void aom_start_encode(aom_writer *bc, uint8_t *buffer) {
#if CONFIG_ANS
  (void)bc;
  (void)buffer;
  assert(0 && "buf_ans requires a more complicated startup procedure");
#else
  aom_dk_start_encode(bc, buffer);
#endif
}

static INLINE void aom_stop_encode(aom_writer *bc) {
#if CONFIG_ANS
  (void)bc;
  assert(0 && "buf_ans requires a more complicated shutdown procedure");
#else
  aom_dk_stop_encode(bc);
#endif
}

static INLINE void aom_write(aom_writer *br, int bit, int probability) {
#if CONFIG_ANS
  buf_uabs_write(br, bit, probability);
#else
  aom_dk_write(br, bit, probability);
#endif
}

static INLINE void aom_write_bit(aom_writer *w, int bit) {
  aom_write(w, bit, 128);  // aom_prob_half
}

static INLINE void aom_write_literal(aom_writer *w, int data, int bits) {
  int bit;

  for (bit = bits - 1; bit >= 0; bit--) aom_write_bit(w, 1 & (data >> bit));
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITWRITER_H_
