/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_BIT_WRITE_BUFFER_H_
#define VP9_BIT_WRITE_BUFFER_H_

#include "vpx/vpx_integer.h"

struct vp9_write_bit_buffer {
  uint8_t *const bit_buffer;
  size_t bit_offset;
};

static void vp9_wb_write_bit(struct vp9_write_bit_buffer *wb, int bit) {
  const int off = wb->bit_offset;
  const int p = off / CHAR_BIT;
  const int q = /*CHAR_BIT - 1 -*/ off % CHAR_BIT;
  wb->bit_buffer[p] &= ~(1 << q);
  wb->bit_buffer[p] |= bit << q;
  wb->bit_offset = off + 1;
}

static void vp9_wb_write_literal(struct vp9_write_bit_buffer *wb,
                              int data, int bits) {
  int bit;
  for (bit = bits - 1; bit >= 0; bit--)
    vp9_wb_write_bit(wb, (data >> bit) & 1);
}


#endif  // VP9_BIT_WRITE_BUFFER_H_
