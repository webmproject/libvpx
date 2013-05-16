/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp9/common/vp9_entropy.h"

const int vp9_default_mode_contexts[INTER_MODE_CONTEXTS][4] = {
  {2,       173,   34,   229},  // 0 = both zero mv
  {7,       145,   85,   225},  // 1 = one zero mv + one a predicted mv
  {7,       166,   63,   231},  // 2 = two predicted mvs
  {7,       94,    66,   219},  // 3 = one predicted/zero and one new mv
  {8,       64,    46,   213},  // 4 = two new mvs
  {17,      81,    31,   231},  // 5 = one intra neighbour + x
  {25,      29,    30,   246},  // 6 = two intra neighbours
};
