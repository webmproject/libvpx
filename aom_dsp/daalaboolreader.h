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

#include "aom_dsp/entdec.h"
#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif

struct daala_reader {
  const uint8_t *buffer;
  const uint8_t *buffer_end;
  od_ec_dec ec;
};

typedef struct daala_reader daala_reader;

int aom_daala_reader_init(daala_reader *r, const uint8_t *buffer, int size);
const uint8_t *aom_daala_reader_find_end(daala_reader *r);

static INLINE int aom_daala_read(daala_reader *r, int prob) {
  if (prob == 128) {
    return od_ec_dec_bits(&r->ec, 1, "aom_bits");
  } else {
    int p = ((prob << 15) + (256 - prob)) >> 8;
    return od_ec_decode_bool_q15(&r->ec, p, "aom");
  }
}

static INLINE int aom_daala_read_bit(daala_reader *r) {
  return aom_daala_read(r, 128);
}

static INLINE int aom_daala_reader_has_error(daala_reader *r) {
  return r->ec.error;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
