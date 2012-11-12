/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "entropy.h"

const int vp9_default_mode_contexts[6][4] = {
  {117,     1,     1,    141},
  {234,     1,     1,    213},
  {128,     90,    22,   145},
  {30,      104,   61,   159},
  {13,      169,   18,   206},
  {15,      76,    24,   166}
};
const int vp9_default_mode_contexts_a[6][4] = {
  {117,     1,     1,    141},
  {234,     1,     1,    213},
  {128,     90,    22,   145},
  {30,      104,   61,   159},
  {13,      169,   18,   206},
  {15,      76,    24,   166}
};
