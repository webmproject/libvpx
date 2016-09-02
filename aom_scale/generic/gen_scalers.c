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

#include "./aom_scale_rtcd.h"
#include "aom_scale/aom_scale.h"
#include "aom_mem/aom_mem.h"
/****************************************************************************
*  Imports
****************************************************************************/

/****************************************************************************
 *
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width   : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width     : Stride of destination
 *                                                (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 4 to 5.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void aom_horizontal_line_5_4_scale_c(const unsigned char *source,
                                     unsigned int source_width,
                                     unsigned char *dest,
                                     unsigned int dest_width) {
  unsigned i;
  unsigned int a, b, c, d, e;
  unsigned char *des = dest;
  const unsigned char *src = source;

  (void)dest_width;

  for (i = 0; i < source_width; i += 5) {
    a = src[0];
    b = src[1];
    c = src[2];
    d = src[3];
    e = src[4];

    des[0] = (unsigned char)a;
    des[1] = (unsigned char)((b * 192 + c * 64 + 128) >> 8);
    des[2] = (unsigned char)((c * 128 + d * 128 + 128) >> 8);
    des[3] = (unsigned char)((d * 64 + e * 192 + 128) >> 8);

    src += 5;
    des += 4;
  }
}

void aom_vertical_band_5_4_scale_c(unsigned char *source,
                                   unsigned int src_pitch, unsigned char *dest,
                                   unsigned int dest_pitch,
                                   unsigned int dest_width) {
  unsigned int i;
  unsigned int a, b, c, d, e;
  unsigned char *des = dest;
  unsigned char *src = source;

  for (i = 0; i < dest_width; i++) {
    a = src[0 * src_pitch];
    b = src[1 * src_pitch];
    c = src[2 * src_pitch];
    d = src[3 * src_pitch];
    e = src[4 * src_pitch];

    des[0 * dest_pitch] = (unsigned char)a;
    des[1 * dest_pitch] = (unsigned char)((b * 192 + c * 64 + 128) >> 8);
    des[2 * dest_pitch] = (unsigned char)((c * 128 + d * 128 + 128) >> 8);
    des[3 * dest_pitch] = (unsigned char)((d * 64 + e * 192 + 128) >> 8);

    src++;
    des++;
  }
}

/*7***************************************************************************
 *
 *  ROUTINE       : aom_horizontal_line_3_5_scale_c
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width   : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width     : Stride of destination
 *                                                (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 3 to 5.
 *
 *  SPECIAL NOTES : None.
 *
 *
 ****************************************************************************/
void aom_horizontal_line_5_3_scale_c(const unsigned char *source,
                                     unsigned int source_width,
                                     unsigned char *dest,
                                     unsigned int dest_width) {
  unsigned int i;
  unsigned int a, b, c, d, e;
  unsigned char *des = dest;
  const unsigned char *src = source;

  (void)dest_width;

  for (i = 0; i < source_width; i += 5) {
    a = src[0];
    b = src[1];
    c = src[2];
    d = src[3];
    e = src[4];

    des[0] = (unsigned char)a;
    des[1] = (unsigned char)((b * 85 + c * 171 + 128) >> 8);
    des[2] = (unsigned char)((d * 171 + e * 85 + 128) >> 8);

    src += 5;
    des += 3;
  }
}

void aom_vertical_band_5_3_scale_c(unsigned char *source,
                                   unsigned int src_pitch, unsigned char *dest,
                                   unsigned int dest_pitch,
                                   unsigned int dest_width) {
  unsigned int i;
  unsigned int a, b, c, d, e;
  unsigned char *des = dest;
  unsigned char *src = source;

  for (i = 0; i < dest_width; i++) {
    a = src[0 * src_pitch];
    b = src[1 * src_pitch];
    c = src[2 * src_pitch];
    d = src[3 * src_pitch];
    e = src[4 * src_pitch];

    des[0 * dest_pitch] = (unsigned char)a;
    des[1 * dest_pitch] = (unsigned char)((b * 85 + c * 171 + 128) >> 8);
    des[2 * dest_pitch] = (unsigned char)((d * 171 + e * 85 + 128) >> 8);

    src++;
    des++;
  }
}

/****************************************************************************
 *
 *  ROUTINE       : aom_horizontal_line_1_2_scale_c
 *
 *  INPUTS        : const unsigned char *source : Pointer to source data.
 *                  unsigned int source_width   : Stride of source.
 *                  unsigned char *dest         : Pointer to destination data.
 *                  unsigned int dest_width     : Stride of destination
 *                                                (NOT USED).
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Copies horizontal line of pixels from source to
 *                  destination scaling up by 1 to 2.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void aom_horizontal_line_2_1_scale_c(const unsigned char *source,
                                     unsigned int source_width,
                                     unsigned char *dest,
                                     unsigned int dest_width) {
  unsigned int i;
  unsigned int a;
  unsigned char *des = dest;
  const unsigned char *src = source;

  (void)dest_width;

  for (i = 0; i < source_width; i += 2) {
    a = src[0];
    des[0] = (unsigned char)(a);
    src += 2;
    des += 1;
  }
}

void aom_vertical_band_2_1_scale_c(unsigned char *source,
                                   unsigned int src_pitch, unsigned char *dest,
                                   unsigned int dest_pitch,
                                   unsigned int dest_width) {
  (void)dest_pitch;
  (void)src_pitch;
  memcpy(dest, source, dest_width);
}

void aom_vertical_band_2_1_scale_i_c(unsigned char *source,
                                     unsigned int src_pitch,
                                     unsigned char *dest,
                                     unsigned int dest_pitch,
                                     unsigned int dest_width) {
  int i;
  int temp;
  int width = dest_width;

  (void)dest_pitch;

  for (i = 0; i < width; i++) {
    temp = 8;
    temp += source[i - (int)src_pitch] * 3;
    temp += source[i] * 10;
    temp += source[i + src_pitch] * 3;
    temp >>= 4;
    dest[i] = (unsigned char)(temp);
  }
}
