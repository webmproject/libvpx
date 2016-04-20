/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_MV_H_
#define VP9_COMMON_VP9_MV_H_

#include "vpx/vpx_integer.h"

#include "vp9/common/vp9_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mv {
  int16_t row;
  int16_t col;
} MV;

typedef union int_mv {
  uint32_t as_int;
  MV as_mv;
} int_mv; /* facilitates faster equality tests and copies */

typedef struct mv32 {
  int32_t row;
  int32_t col;
} MV32;

static INLINE int is_zero_mv(const MV *mv) {
  return *((const uint32_t *)mv) == 0;
}

static INLINE int is_equal_mv(const MV *a, const MV *b) {
  return  *((const uint32_t *)a) == *((const uint32_t *)b);
}

static INLINE void clamp_mv(MV *mv, int min_col, int max_col,
                            int min_row, int max_row) {
  mv->col = clamp(mv->col, min_col, max_col);
  mv->row = clamp(mv->row, min_row, max_row);
}

#if CONFIG_GLOBAL_MOTION
#define MAX_GLOBAL_MOTION_MODELS  1

#define ZOOM_PRECISION_BITS       11
#define ROTATION_PRECISION_BITS   11

#define ABS_ZOOM_BITS             11
#define ABS_ROTATION_BITS         11
#define ABS_TRANSLATION_BITS      11

typedef enum {
  GLOBAL_ZERO = 0,
  GLOBAL_TRANSLATION = 1,
  GLOBAL_ROTZOOM = 2,
  GLOBAL_MOTION_TYPES
} GLOBAL_MOTION_TYPE;

// Currently this is specialized for rotzoom model only
typedef struct {
  GLOBAL_MOTION_TYPE gmtype;
  int rotation;   // positive or negative rotation angle in degrees
  int zoom;       // this is actually the zoom multiplier minus 1
  int_mv mv;
} Global_Motion_Params;

static INLINE GLOBAL_MOTION_TYPE get_gmtype(const Global_Motion_Params *gm) {
  if (gm->rotation == 0 && gm->zoom == 0) {
    return (gm->mv.as_int == 0 ? GLOBAL_ZERO : GLOBAL_TRANSLATION);
  } else {
    return GLOBAL_ROTZOOM;
  }
}
#endif  // CONFIG_GLOBAL_MOTION

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_COMMON_VP9_MV_H_
