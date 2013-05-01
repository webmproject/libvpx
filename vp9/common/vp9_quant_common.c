/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp9/common/vp9_common.h"
#include "vp9/common/vp9_quant_common.h"
#include "vp9/common/vp9_seg_common.h"

static int16_t dc_qlookup[QINDEX_RANGE];
static int16_t ac_qlookup[QINDEX_RANGE];

#define ACDC_MIN 4

// TODO(dkovalev) move to common and reuse
static double poly3(double a, double b, double c, double d, double x) {
  return a*x*x*x + b*x*x + c*x + d;
}

void vp9_init_quant_tables() {
  int i, val = 4;

  for (i = 0; i < QINDEX_RANGE; i++) {
    const int ac_val = val;

    val = (int)(val * 1.02);
    if (val == ac_val)
      ++val;

    ac_qlookup[i] = (int16_t)ac_val;
    dc_qlookup[i] = (int16_t)MAX(ACDC_MIN, poly3(0.000000305, -0.00065, 0.9,
                                                 0.5, ac_val));
  }
}

int16_t vp9_dc_quant(int qindex, int delta) {
  return dc_qlookup[clamp(qindex + delta, 0, MAXQ)];
}

int16_t vp9_ac_quant(int qindex, int delta) {
  return ac_qlookup[clamp(qindex + delta, 0, MAXQ)];
}


int vp9_get_qindex(MACROBLOCKD *xd, int segment_id, int base_qindex) {
  if (vp9_segfeature_active(xd, segment_id, SEG_LVL_ALT_Q)) {
    const int data = vp9_get_segdata(xd, segment_id, SEG_LVL_ALT_Q);
    return xd->mb_segment_abs_delta == SEGMENT_ABSDATA ?
               data :  // Abs value
               clamp(base_qindex + data, 0, MAXQ);  // Delta value
  } else {
    return base_qindex;
  }
}

