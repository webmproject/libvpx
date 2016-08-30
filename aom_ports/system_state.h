/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AOM_PORTS_SYSTEM_STATE_H_
#define AOM_PORTS_SYSTEM_STATE_H_

#include "./aom_config.h"

#if ARCH_X86 || ARCH_X86_64
void aom_reset_mmx_state(void);
#define aom_clear_system_state() aom_reset_mmx_state()
#else
#define aom_clear_system_state()
#endif  // ARCH_X86 || ARCH_X86_64
#endif  // AOM_PORTS_SYSTEM_STATE_H_
