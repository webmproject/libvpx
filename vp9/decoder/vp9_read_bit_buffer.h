/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_READ_BIT_BUFFER_
#define VP9_READ_BIT_BUFFER_

struct vp9_read_bit_buffer {
  const uint8_t *const bit_buffer;
  size_t bit_offset;
};

static int vp9_rb_read_bit(struct vp9_read_bit_buffer *rb) {
  const int off = rb->bit_offset;
  const int p = off / CHAR_BIT;
  const int q = /*CHAR_BIT - 1 -*/ off % CHAR_BIT;
  const int bit = (rb->bit_buffer[p] & (1 << q)) >> q;
  rb->bit_offset = off + 1;
  return bit;
}

static int vp9_rb_read_literal(struct vp9_read_bit_buffer *rb, int bits) {
  int value = 0, bit;
  for (bit = bits - 1; bit >= 0; bit--)
    value |= vp9_rb_read_bit(rb) << bit;
  return value;
}

#endif  // VP9_READ_BIT_BUFFER_
