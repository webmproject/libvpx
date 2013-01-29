/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "./vpx_config.h"
#include "vpx_ports/mem.h"
///////////////////////////////////////////////////////////////////////////
// the mmx function that does the bilinear filtering and var calculation //
// int one pass                                                          //
///////////////////////////////////////////////////////////////////////////
DECLARE_ALIGNED(16, const short, vp9_bilinear_filters_mmx[16][8]) = {
  { 128, 128, 128, 128,  0,  0,  0,  0 },
  { 120, 120, 120, 120,  8,  8,  8,  8 },
  { 112, 112, 112, 112, 16, 16, 16, 16 },
  { 104, 104, 104, 104, 24, 24, 24, 24 },
  {  96, 96, 96, 96, 32, 32, 32, 32 },
  {  88, 88, 88, 88, 40, 40, 40, 40 },
  {  80, 80, 80, 80, 48, 48, 48, 48 },
  {  72, 72, 72, 72, 56, 56, 56, 56 },
  {  64, 64, 64, 64, 64, 64, 64, 64 },
  {  56, 56, 56, 56, 72, 72, 72, 72 },
  {  48, 48, 48, 48, 80, 80, 80, 80 },
  {  40, 40, 40, 40, 88, 88, 88, 88 },
  {  32, 32, 32, 32, 96, 96, 96, 96 },
  {  24, 24, 24, 24, 104, 104, 104, 104 },
  {  16, 16, 16, 16, 112, 112, 112, 112 },
  {   8,  8,  8,  8, 120, 120, 120, 120 }
};
