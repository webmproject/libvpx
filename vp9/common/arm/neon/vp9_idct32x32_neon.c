/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_common.h"

// defined in vp9/common/arm/neon/vp9_short_idct32x32_add_neon.asm
extern void idct32_transpose_and_transform(int16_t *transpose_buffer,
                                           int16_t *output, int16_t *input);
extern void idct32_combine_add(uint8_t *dest, int16_t *out, int dest_stride);


// defined in vp9/common/arm/neon/vp9_short_idct16x16_add_neon.asm
extern void save_neon_registers();
extern void restore_neon_registers();

void vp9_short_idct32x32_add_neon(int16_t *input, uint8_t *dest,
                                  int dest_stride) {
  // TODO(cd): move the creation of these buffers within the ASM file
  // internal buffer used to transpose 8 lines into before transforming them
  int16_t transpose_buffer[32 * 8];
  // results of the first pass (transpose and transform rows)
  int16_t pass1[32 * 32];
  // results of the second pass (transpose and transform columns)
  int16_t pass2[32 * 32];

  // save register we need to preserve
  save_neon_registers();
  // process rows
  idct32_transpose_and_transform(transpose_buffer, pass1, input);
  // process columns
  // TODO(cd): do these two steps/passes within the ASM file
  idct32_transpose_and_transform(transpose_buffer, pass2, pass1);
  // combine and add to dest
  // TODO(cd): integrate this within the last storage step of the second pass
  idct32_combine_add(dest, pass2, dest_stride);
  // restore register we need to preserve
  restore_neon_registers();
}

// TODO(cd): Eliminate this file altogether when everything is in ASM file
