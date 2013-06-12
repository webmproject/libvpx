/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_blockd.h"

void vp9_setup_block_dptrs(MACROBLOCKD *mb,
                           int subsampling_x, int subsampling_y) {
  int i;

  for (i = 0; i < MAX_MB_PLANE; i++) {
    mb->plane[i].plane_type = i ? PLANE_TYPE_UV : PLANE_TYPE_Y_WITH_DC;
    mb->plane[i].subsampling_x = i ? subsampling_x : 0;
    mb->plane[i].subsampling_y = i ? subsampling_y : 0;
  }
#if CONFIG_ALPHA
  // TODO(jkoleszar): Using the Y w/h for now
  mb->plane[3].subsampling_x = 0;
  mb->plane[3].subsampling_y = 0;
#endif
}
