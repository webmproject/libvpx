/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "./aom_config.h"
#define RTCD_C
#include "./aom_dsp_rtcd.h"
#include "aom_ports/aom_once.h"

void aom_dsp_rtcd() { once(setup_rtcd_internal); }
