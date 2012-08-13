/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "vpx_config.h"
#define RTCD_C
#include "vpx_rtcd.h"

/* No-op version that performs no synchronization. vpx_rtcd() is idempotent,
 * so as long as your platform provides atomic loads/stores of pointers
 * no synchronization is strictly necessary.
 */

static void once(void (*func)(void)) {
  static int done;

  if(!done) {
    func();
    done = 1;
  }
}

void vpx_rtcd() {
  once(setup_rtcd_internal);
}
