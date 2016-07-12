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

#include "vp10/common/filter.h"

DECLARE_ALIGNED(256, static const InterpKernel,
                bilinear_filters[SUBPEL_SHIFTS]) = {
  { 0, 0, 0, 128,   0, 0, 0, 0 },
  { 0, 0, 0, 120,   8, 0, 0, 0 },
  { 0, 0, 0, 112,  16, 0, 0, 0 },
  { 0, 0, 0, 104,  24, 0, 0, 0 },
  { 0, 0, 0,  96,  32, 0, 0, 0 },
  { 0, 0, 0,  88,  40, 0, 0, 0 },
  { 0, 0, 0,  80,  48, 0, 0, 0 },
  { 0, 0, 0,  72,  56, 0, 0, 0 },
  { 0, 0, 0,  64,  64, 0, 0, 0 },
  { 0, 0, 0,  56,  72, 0, 0, 0 },
  { 0, 0, 0,  48,  80, 0, 0, 0 },
  { 0, 0, 0,  40,  88, 0, 0, 0 },
  { 0, 0, 0,  32,  96, 0, 0, 0 },
  { 0, 0, 0,  24, 104, 0, 0, 0 },
  { 0, 0, 0,  16, 112, 0, 0, 0 },
  { 0, 0, 0,   8, 120, 0, 0, 0 }
};

#if USE_TEMPORALFILTER_12TAP
DECLARE_ALIGNED(16, static const int16_t,
                sub_pel_filters_temporalfilter_12[SUBPEL_SHIFTS][12]) = {
  // intfilt 0.8
  {0,   0,   0,   0,   0, 128,   0,   0,   0,   0,   0, 0},
  {0,   1,  -1,   3,  -7, 127,   8,  -4,   2,  -1,   0, 0},
  {0,   1,  -3,   5, -12, 124,  18,  -8,   4,  -2,   1, 0},
  {-1,   2,  -4,   8, -17, 120,  28, -11,   6,  -3,   1, -1},
  {-1,   2,  -4,  10, -21, 114,  38, -15,   8,  -4,   2, -1},
  {-1,   3,  -5,  11, -23, 107,  49, -18,   9,  -5,   2, -1},
  {-1,   3,  -6,  12, -25,  99,  60, -21,  11,  -6,   3, -1},
  {-1,   3,  -6,  12, -25,  90,  70, -23,  12,  -6,   3, -1},
  {-1,   3,  -6,  12, -24,  80,  80, -24,  12,  -6,   3, -1},
  {-1,   3,  -6,  12, -23,  70,  90, -25,  12,  -6,   3, -1},
  {-1,   3,  -6,  11, -21,  60,  99, -25,  12,  -6,   3, -1},
  {-1,   2,  -5,   9, -18,  49, 107, -23,  11,  -5,   3, -1},
  {-1,   2,  -4,   8, -15,  38, 114, -21,  10,  -4,   2, -1},
  {-1,   1,  -3,   6, -11,  28, 120, -17,   8,  -4,   2, -1},
  {0,   1,  -2,   4,  -8,  18, 124, -12,   5,  -3,   1, 0},
  {0,   0,  -1,   2,  -4,   8, 127,  -7,   3,  -1,   1, 0},
};
#endif  // USE_TEMPORALFILTER_12TAP

#if CONFIG_EXT_INTERP
DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8[SUBPEL_SHIFTS]) = {
  // intfilt 0.575
  {0,   0,   0, 128,   0,   0,   0, 0},
  {0,   1,  -5, 126,   8,  -3,   1, 0},
  {-1,   3, -10, 123,  18,  -6,   2, -1},
  {-1,   4, -14, 118,  27,  -9,   3, 0},
  {-1,   5, -16, 112,  37, -12,   4, -1},
  {-1,   5, -18, 105,  48, -14,   4, -1},
  {-1,   6, -19,  97,  58, -17,   5, -1},
  {-1,   6, -20,  88,  68, -18,   6, -1},
  {-1,   6, -19,  78,  78, -19,   6, -1},
  {-1,   6, -18,  68,  88, -20,   6, -1},
  {-1,   5, -17,  58,  97, -19,   6, -1},
  {-1,   4, -14,  48, 105, -18,   5, -1},
  {-1,   4, -12,  37, 112, -16,   5, -1},
  {0,   3,  -9,  27, 118, -14,   4, -1},
  {-1,   2,  -6,  18, 123, -10,   3, -1},
  {0,   1,  -3,   8, 126,  -5,   1, 0},
};

#if CONFIG_EXT_INTRA
DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8sharp[SUBPEL_SHIFTS]) = {
  // intfilt 0.8
  {0,   0,   0, 128,   0,   0,   0, 0},
  {-1,   2,  -6, 127,   9,  -4,   2, -1},
  {-2,   5, -12, 124,  18,  -7,   4, -2},
  {-2,   7, -16, 119,  28, -11,   5, -2},
  {-3,   8, -19, 114,  38, -14,   7, -3},
  {-3,   9, -22, 107,  49, -17,   8, -3},
  {-4,  10, -23,  99,  60, -20,  10, -4},
  {-4,  11, -23,  90,  70, -22,  10, -4},
  {-4,  11, -23,  80,  80, -23,  11, -4},
  {-4,  10, -22,  70,  90, -23,  11, -4},
  {-4,  10, -20,  60,  99, -23,  10, -4},
  {-3,   8, -17,  49, 107, -22,   9, -3},
  {-3,   7, -14,  38, 114, -19,   8, -3},
  {-2,   5, -11,  28, 119, -16,   7, -2},
  {-2,   4,  -7,  18, 124, -12,   5, -2},
  {-1,   2,  -4,   9, 127,  -6,   2, -1},
};
#endif  // CONFIG_EXT_INTRA

DECLARE_ALIGNED(256, static const int16_t,
                sub_pel_filters_10sharp[SUBPEL_SHIFTS][10]) = {
  // intfilt 0.77
  {0,   0,   0,   0, 128,   0,   0,   0,   0, 0},
  {0,  -1,   3,  -6, 127,   8,  -4,   2,  -1, 0},
  {1,  -2,   5, -12, 124,  18,  -7,   3,  -2, 0},
  {1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1},
  {1,  -4,   8, -20, 114,  38, -14,   7,  -3, 1},
  {1,  -4,   9, -22, 107,  49, -17,   8,  -4, 1},
  {2,  -5,  10, -24,  99,  59, -20,   9,  -4, 2},
  {2,  -5,  10, -24,  90,  70, -22,  10,  -5, 2},
  {2,  -5,  10, -23,  80,  80, -23,  10,  -5, 2},
  {2,  -5,  10, -22,  70,  90, -24,  10,  -5, 2},
  {2,  -4,   9, -20,  59,  99, -24,  10,  -5, 2},
  {1,  -4,   8, -17,  49, 107, -22,   9,  -4, 1},
  {1,  -3,   7, -14,  38, 114, -20,   8,  -4, 1},
  {1,  -2,   5, -11,  28, 119, -17,   7,  -3, 1},
  {0,  -2,   3,  -7,  18, 124, -12,   5,  -2, 1},
  {0,  -1,   2,  -4,   8, 127,  -6,   3,  -1, 0},
};

DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8smooth2[SUBPEL_SHIFTS]) = {
// freqmultiplier = 0.35
  {0,  0,  0, 128,  0,  0,  0,  0},
  {-1,  8, 31, 47, 34, 10,  0, -1},
  {-1,  7, 29, 46, 36, 12,  0, -1},
  {-1,  6, 28, 46, 37, 13,  0, -1},
  {-1,  5, 26, 46, 38, 14,  1, -1},
  {-1,  4, 25, 45, 39, 16,  1, -1},
  {-1,  4, 23, 44, 41, 17,  1, -1},
  {-1,  3, 21, 44, 42, 18,  2, -1},
  {-1,  2, 20, 43, 43, 20,  2, -1},
  {-1,  2, 18, 42, 44, 21,  3, -1},
  {-1,  1, 17, 41, 44, 23,  4, -1},
  {-1,  1, 16, 39, 45, 25,  4, -1},
  {-1,  1, 14, 38, 46, 26,  5, -1},
  {-1,  0, 13, 37, 46, 28,  6, -1},
  {-1,  0, 12, 36, 46, 29,  7, -1},
  {-1,  0, 10, 34, 47, 31,  8, -1},
};

DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8smooth[SUBPEL_SHIFTS]) = {
// freqmultiplier = 0.75
  {0,  0,  0, 128,  0,  0,  0,  0},
  {2, -10,  19,  95,  31, -11,   2, 0},
  {2,  -9,  14,  94,  37, -12,   2, 0},
  {2,  -8,   9,  92,  43, -12,   1, 1},
  {2,  -7,   5,  90,  49, -12,   1, 0},
  {2,  -5,   1,  86,  55, -12,   0, 1},
  {1,  -4,  -2,  82,  61, -11,   0, 1},
  {1, -3, -5, 77, 67, -9, -1, 1},
  {1, -2, -7, 72, 72, -7, -2, 1},
  {1, -1, -9, 67, 77, -5, -3, 1},
  {1,   0, -11,  61,  82,  -2,  -4, 1},
  {1,   0, -12,  55,  86,   1,  -5, 2},
  {0,   1, -12,  49,  90,   5,  -7, 2},
  {1,   1, -12,  43,  92,   9,  -8, 2},
  {0,   2, -12,  37,  94,  14,  -9, 2},
  {0,   2, -11,  31,  95,  19, -10, 2},
};

DECLARE_ALIGNED(16, static const int16_t,
                sub_pel_filters_12sharp[SUBPEL_SHIFTS][12]) = {
  // intfilt 0.85
  {0,   0,   0,   0,   0, 128,   0,   0,   0,   0,   0, 0},
  {0,   1,  -2,   3,  -7, 127,   8,  -4,   2,  -1,   1, 0},
  {-1,   2,  -3,   6, -13, 124,  18,  -8,   4,  -2,   2, -1},
  {-1,   3,  -4,   8, -18, 120,  28, -12,   7,  -4,   2, -1},
  {-1,   3,  -6,  10, -21, 115,  38, -15,   8,  -5,   3, -1},
  {-2,   4,  -6,  12, -24, 108,  49, -18,  10,  -6,   3, -2},
  {-2,   4,  -7,  13, -25, 100,  60, -21,  11,  -7,   4, -2},
  {-2,   4,  -7,  13, -26,  91,  71, -24,  13,  -7,   4, -2},
  {-2,   4,  -7,  13, -25,  81,  81, -25,  13,  -7,   4, -2},
  {-2,   4,  -7,  13, -24,  71,  91, -26,  13,  -7,   4, -2},
  {-2,   4,  -7,  11, -21,  60, 100, -25,  13,  -7,   4, -2},
  {-2,   3,  -6,  10, -18,  49, 108, -24,  12,  -6,   4, -2},
  {-1,   3,  -5,   8, -15,  38, 115, -21,  10,  -6,   3, -1},
  {-1,   2,  -4,   7, -12,  28, 120, -18,   8,  -4,   3, -1},
  {-1,   2,  -2,   4,  -8,  18, 124, -13,   6,  -3,   2, -1},
  {0,   1,  -1,   2,  -4,   8, 127,  -7,   3,  -2,   1, 0},
};
#else  // CONFIG_EXT_INTERP

DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8[SUBPEL_SHIFTS]) = {
  // Lagrangian interpolation filter
  { 0,   0,   0, 128,   0,   0,   0,  0},
  { 0,   1,  -5, 126,   8,  -3,   1,  0},
  { -1,   3, -10, 122,  18,  -6,   2,  0},
  { -1,   4, -13, 118,  27,  -9,   3, -1},
  { -1,   4, -16, 112,  37, -11,   4, -1},
  { -1,   5, -18, 105,  48, -14,   4, -1},
  { -1,   5, -19,  97,  58, -16,   5, -1},
  { -1,   6, -19,  88,  68, -18,   5, -1},
  { -1,   6, -19,  78,  78, -19,   6, -1},
  { -1,   5, -18,  68,  88, -19,   6, -1},
  { -1,   5, -16,  58,  97, -19,   5, -1},
  { -1,   4, -14,  48, 105, -18,   5, -1},
  { -1,   4, -11,  37, 112, -16,   4, -1},
  { -1,   3,  -9,  27, 118, -13,   4, -1},
  { 0,   2,  -6,  18, 122, -10,   3, -1},
  { 0,   1,  -3,   8, 126,  -5,   1,  0}
};

DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8sharp[SUBPEL_SHIFTS]) = {
  // DCT based filter
  {0,   0,   0, 128,   0,   0,   0, 0},
  {-1,   3,  -7, 127,   8,  -3,   1, 0},
  {-2,   5, -13, 125,  17,  -6,   3, -1},
  {-3,   7, -17, 121,  27, -10,   5, -2},
  {-4,   9, -20, 115,  37, -13,   6, -2},
  {-4,  10, -23, 108,  48, -16,   8, -3},
  {-4,  10, -24, 100,  59, -19,   9, -3},
  {-4,  11, -24,  90,  70, -21,  10, -4},
  {-4,  11, -23,  80,  80, -23,  11, -4},
  {-4,  10, -21,  70,  90, -24,  11, -4},
  {-3,   9, -19,  59, 100, -24,  10, -4},
  {-3,   8, -16,  48, 108, -23,  10, -4},
  {-2,   6, -13,  37, 115, -20,   9, -4},
  {-2,   5, -10,  27, 121, -17,   7, -3},
  {-1,   3,  -6,  17, 125, -13,   5, -2},
  {0,   1,  -3,   8, 127,  -7,   3, -1}
};

DECLARE_ALIGNED(256, static const InterpKernel,
                sub_pel_filters_8smooth[SUBPEL_SHIFTS]) = {
// freqmultiplier = 0.5
  { 0,  0,  0, 128,  0,  0,  0,  0},
  {-3, -1, 32,  64, 38,  1, -3,  0},
  {-2, -2, 29,  63, 41,  2, -3,  0},
  {-2, -2, 26,  63, 43,  4, -4,  0},
  {-2, -3, 24,  62, 46,  5, -4,  0},
  {-2, -3, 21,  60, 49,  7, -4,  0},
  {-1, -4, 18,  59, 51,  9, -4,  0},
  {-1, -4, 16,  57, 53, 12, -4, -1},
  {-1, -4, 14,  55, 55, 14, -4, -1},
  {-1, -4, 12,  53, 57, 16, -4, -1},
  { 0, -4,  9,  51, 59, 18, -4, -1},
  { 0, -4,  7,  49, 60, 21, -3, -2},
  { 0, -4,  5,  46, 62, 24, -3, -2},
  { 0, -4,  4,  43, 63, 26, -2, -2},
  { 0, -3,  2,  41, 63, 29, -2, -2},
  { 0, -3,  1,  38, 64, 32, -1, -3}
};
#endif  // CONFIG_EXT_INTERP

#if CONFIG_EXT_INTRA
const InterpKernel *vp10_intra_filter_kernels[INTRA_FILTERS] = {
  bilinear_filters,         // INTRA_FILTER_LINEAR
  sub_pel_filters_8,        // INTRA_FILTER_8TAP
  sub_pel_filters_8sharp,   // INTRA_FILTER_8TAP_SHARP
  sub_pel_filters_8smooth,  // INTRA_FILTER_8TAP_SMOOTH
};
#endif  // CONFIG_EXT_INTRA

#if CONFIG_EXT_INTERP
static const InterpFilterParams
vp10_interp_filter_params_list[SWITCHABLE_FILTERS + 1] = {
  {(const int16_t*)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_10sharp, 10, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_8smooth2, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_12sharp, 12, SUBPEL_SHIFTS},
  {(const int16_t*)bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS}
};
#else
static const InterpFilterParams
vp10_interp_filter_params_list[SWITCHABLE_FILTERS + 1] = {
  {(const int16_t*)sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)sub_pel_filters_8sharp, SUBPEL_TAPS, SUBPEL_SHIFTS},
  {(const int16_t*)bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS}
};
#endif  // CONFIG_EXT_INTERP

#if USE_TEMPORALFILTER_12TAP
static const InterpFilterParams vp10_interp_temporalfilter_12tap = {
  (const int16_t*)sub_pel_filters_temporalfilter_12, 12, SUBPEL_SHIFTS
};
#endif  // USE_TEMPORALFILTER_12TAP

InterpFilterParams vp10_get_interp_filter_params(
    const INTERP_FILTER interp_filter) {
#if USE_TEMPORALFILTER_12TAP
  if (interp_filter == TEMPORALFILTER_12TAP)
    return vp10_interp_temporalfilter_12tap;
#endif  // USE_TEMPORALFILTER_12TAP
  return vp10_interp_filter_params_list[interp_filter];
}

const int16_t *vp10_get_interp_filter_kernel(
    const INTERP_FILTER interp_filter) {
#if USE_TEMPORALFILTER_12TAP
  if (interp_filter == TEMPORALFILTER_12TAP)
    return vp10_interp_temporalfilter_12tap.filter_ptr;
#endif  // USE_TEMPORALFILTER_12TAP
  return (const int16_t*)
      vp10_interp_filter_params_list[interp_filter].filter_ptr;
}

SubpelFilterCoeffs vp10_get_subpel_filter_signal_dir(
    const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}

SubpelFilterCoeffs vp10_get_subpel_filter_ver_signal_dir(
    const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_ver_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_ver_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSSE3
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_ver_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}

#if CONFIG_VP9_HIGHBITDEPTH
HbdSubpelFilterCoeffs vp10_hbd_get_subpel_filter_ver_signal_dir(
    const InterpFilterParams p, int index) {
#if CONFIG_EXT_INTERP && HAVE_SSE4_1
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_12sharp) {
    return &sub_pel_filters_12sharp_highbd_ver_signal_dir[index][0];
  }
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_10sharp) {
    return &sub_pel_filters_10sharp_highbd_ver_signal_dir[index][0];
  }
#endif
#if USE_TEMPORALFILTER_12TAP && HAVE_SSE4_1
  if (p.filter_ptr == (const int16_t *)sub_pel_filters_temporalfilter_12) {
    return &sub_pel_filters_temporalfilter_12_highbd_ver_signal_dir[index][0];
  }
#endif
  (void)p;
  (void)index;
  return NULL;
}
#endif
