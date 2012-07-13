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
const int default_vp8_mode_contexts[6][4] = {
  {
    /* 0 */
    7,     1,     1,   183
  },
  {
    /* 1 */
    14,    18,    14,   147
  },
  {
    /* 2 */
    135,    64,    57,    68
  },
  {
    /* 3 */
    60,    56,   128,   65
  },
  {
    /* 4 */
    159,   134,   128,   34
  },
  {
    /* 5 */
    234,   188,   128,   28
  },
};
const int default_vp8_mode_contexts_a[6][4] = {
  {
    /* 0 */
    4,     1,    1,   143
  },
  {
    /* 1 */
    7,     9,    7,   107
  },
  {
    /* 2 */
    95,    34,   57,    68
  },
  {
    /* 3 */
    95,    56,   128,   65
  },
  {
    /* 4 */
    159,   67,   128,   34
  },
  {
    /* 5 */
    234,   94,   128,   28
  },
};
