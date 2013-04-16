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

int vp9_start_decode(BOOL_DECODER *br, const uint8_t *buffer, size_t size) {
  br->buffer_end = buffer + size;
  br->buffer = buffer;
  br->value = 0;
  br->count = -8;
  br->range = 255;

  if (size && !buffer)
    return 1;

  vp9_reader_fill(br);
  return 0;
}

void vp9_reader_fill(BOOL_DECODER *br) {
  const uint8_t *const buffer_end = br->buffer_end;
  const uint8_t *buffer = br->buffer;
  VP9_BD_VALUE value = br->value;
  int count = br->count;
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

  br->buffer = buffer;
  br->value = value;
  br->count = count;
}


static int get_unsigned_bits(unsigned int num_values) {
  int cat = 0;
  if (num_values <= 1)
    return 0;
  num_values--;
  while (num_values > 0) {
    cat++;
    num_values >>= 1;
  }
  return cat;
}

int vp9_inv_recenter_nonneg(int v, int m) {
  if (v > (m << 1))
    return v;
  else if ((v & 1) == 0)
    return (v >> 1) + m;
  else
    return m - ((v + 1) >> 1);
}

int vp9_decode_uniform(BOOL_DECODER *br, int n) {
  int v;
  const int l = get_unsigned_bits(n);
  const int m = (1 << l) - n;
  if (!l)
    return 0;

  v = vp9_read_literal(br, l - 1);
  return v < m ?  v : (v << 1) - m + vp9_read_bit(br);
}

int vp9_decode_term_subexp(BOOL_DECODER *br, int k, int num_syms) {
  int i = 0, mk = 0, word;
  while (1) {
    const int b = i ? k + i - 1 : k;
    const int a = 1 << b;
    if (num_syms <= mk + 3 * a) {
      word = vp9_decode_uniform(br, num_syms - mk) + mk;
      break;
    } else {
      if (vp9_read_bit(br)) {
        i++;
        mk += a;
      } else {
        word = vp9_read_literal(br, b) + mk;
        break;
      }
    }
  }
  return word;
}

int vp9_decode_unsigned_max(BOOL_DECODER *br, int max) {
  int data = 0, bit = 0, lmax = max;

  while (lmax) {
    data |= vp9_read_bit(br) << bit++;
    lmax >>= 1;
  }
  return data > max ? max : data;
}
