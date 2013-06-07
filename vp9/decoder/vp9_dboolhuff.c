/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_ports/mem.h"
#include "vpx_mem/vpx_mem.h"

#include "vp9/decoder/vp9_dboolhuff.h"

int vp9_reader_init(vp9_reader *r, const uint8_t *buffer, size_t size) {
  int marker_bit;

  r->buffer_end = buffer + size;
  r->buffer = buffer;
  r->value = 0;
  r->count = -8;
  r->range = 255;

  if (size && !buffer)
    return 1;

  vp9_reader_fill(r);
  marker_bit = vp9_read_bit(r);
  return marker_bit != 0;
}

void vp9_reader_fill(vp9_reader *r) {
  const uint8_t *const buffer_end = r->buffer_end;
  const uint8_t *buffer = r->buffer;
  VP9_BD_VALUE value = r->value;
  int count = r->count;
  int shift = VP9_BD_VALUE_SIZE - 8 - (count + 8);
  int loop_end = 0;
  const int bits_left = (int)((buffer_end - buffer)*CHAR_BIT);
  const int x = shift + CHAR_BIT - bits_left;

  if (x >= 0) {
    count += VP9_LOTS_OF_BITS;
    loop_end = x;
  }

  if (x < 0 || bits_left) {
    while (shift >= loop_end) {
      count += CHAR_BIT;
      value |= (VP9_BD_VALUE)*buffer++ << shift;
      shift -= CHAR_BIT;
    }
  }

  r->buffer = buffer;
  r->value = value;
  r->count = count;
}

const uint8_t *vp9_reader_find_end(vp9_reader *r) {
  // Find the end of the coded buffer
  while (r->count > CHAR_BIT && r->count < VP9_BD_VALUE_SIZE) {
    r->count -= CHAR_BIT;
    r->buffer--;
  }
  return r->buffer;
}

