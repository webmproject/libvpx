/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./dkboolwriter.h"

static INLINE void aom_dk_write_bit(aom_dk_writer *w, int bit) {
  aom_dk_write(w, bit, 128);  // aom_prob_half
}

void aom_dk_start_encode(aom_dk_writer *br, uint8_t *source) {
  br->lowvalue = 0;
  br->range = 255;
  br->count = -24;
  br->buffer = source;
  br->pos = 0;
  aom_dk_write_bit(br, 0);
}

void aom_dk_stop_encode(aom_dk_writer *br) {
  int i;

#if CONFIG_BITSTREAM_DEBUG
  bitstream_queue_set_skip_write(1);
#endif  // CONFIG_BITSTREAM_DEBUG

  for (i = 0; i < 32; i++) aom_dk_write_bit(br, 0);

#if CONFIG_BITSTREAM_DEBUG
  bitstream_queue_set_skip_write(0);
#endif  // CONFIG_BITSTREAM_DEBUG

  // Ensure there's no ambigous collision with any index marker bytes
  if ((br->buffer[br->pos - 1] & 0xe0) == 0xc0) br->buffer[br->pos++] = 0;
}
