/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "./aom_config.h"
#define RTCD_C
#include "./av1_rtcd.h"
#include "aom_ports/aom_once.h"

void av1_rtcd() {
  // TODO(JBB): Remove this once, by insuring that both the encoder and
  // decoder setup functions are protected by once();
  once(setup_rtcd_internal);
}
