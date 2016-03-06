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

#ifndef AOM_DSP_DAALABOOLWRITER_H_
#define AOM_DSP_DAALABOOLWRITER_H_

#include "aom_dsp/entenc.h"

#ifdef __cplusplus
extern "C" {
#endif

struct daala_writer {
  unsigned int pos;
  uint8_t *buffer;
  od_ec_enc ec;
};

typedef struct daala_writer daala_writer;

void aom_daala_start_encode(daala_writer *w, uint8_t *buffer);
void aom_daala_stop_encode(daala_writer *w);

static INLINE void aom_daala_write(daala_writer *w, int bit, int prob) {
  if (prob == 128) {
    od_ec_enc_bits(&w->ec, bit, 1);
  } else {
    int p = ((prob << 15) + (256 - prob)) >> 8;
    od_ec_encode_bool_q15(&w->ec, bit, p);
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
