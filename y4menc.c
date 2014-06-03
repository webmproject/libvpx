/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./y4menc.h"

int y4m_write_file_header(char *buf, size_t len, int width, int height,
                          const struct VpxRational *framerate,
                          vpx_img_fmt_t fmt, vpx_bit_depth_t bit_depth) {
  char *color;
  switch (bit_depth) {
    default:
      color = fmt == VPX_IMG_FMT_444A ? "C444alpha\n" :
              fmt == VPX_IMG_FMT_I444 ? "C444\n" :
              fmt == VPX_IMG_FMT_I422 ? "C422\n" :
              "C420jpeg\n";
      break;
    case VPX_BITS_10:
      color = fmt == VPX_IMG_FMT_I44416 ? "C444p10 XYSCSS=444P10\n" :
              fmt == VPX_IMG_FMT_I42216 ? "C422p10 XYSCSS=422P10\n" :
              "C420p10 XYSCSS=420P10\n";
      break;
    case VPX_BITS_12:
      color = fmt == VPX_IMG_FMT_I44416 ? "C444p12 XYSCSS=444P12\n" :
              fmt == VPX_IMG_FMT_I42216 ? "C422p12 XYSCSS=422P12\n" :
              "C420p12 XYSCSS=420P12\n";
      break;
  }

  return snprintf(buf, len, "YUV4MPEG2 W%u H%u F%u:%u I%c %s", width, height,
                  framerate->numerator, framerate->denominator, 'p', color);
}

int y4m_write_frame_header(char *buf, size_t len) {
  return snprintf(buf, len, "FRAME\n");
}
