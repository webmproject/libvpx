/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vp9/encoder/vp9_variance.h"
#include "vp9/encoder/vp9_onyx_int.h"


void vp9_cmachine_specific_config(VP9_COMP *cpi) {
#if CONFIG_RUNTIME_CPU_DETECT
  cpi->rtcd.common                    = &cpi->common.rtcd;
#endif
}
