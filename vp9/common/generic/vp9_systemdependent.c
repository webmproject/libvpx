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
#include "vp9_rtcd.h"
#include "vp9/common/vp9_subpixel.h"
#include "vp9/common/vp9_loopfilter.h"
#include "vp9/common/vp9_onyxc_int.h"

extern void vp9_arch_x86_common_init(VP9_COMMON *ctx);
extern void vp9_arch_arm_common_init(VP9_COMMON *ctx);

void vp9_machine_specific_config(VP9_COMMON *ctx) {
#if CONFIG_RUNTIME_CPU_DETECT
  VP9_COMMON_RTCD *rtcd = &ctx->rtcd;

#endif

#if ARCH_X86 || ARCH_X86_64
  vp9_arch_x86_common_init(ctx);
#endif

#if ARCH_ARM
  vp9_arch_arm_common_init(ctx);
#endif

  vp9_rtcd();
}
