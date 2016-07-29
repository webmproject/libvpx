/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_TXFM_TEST_H_
#define VP10_TXFM_TEST_H_

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "vp10/common/enums.h"
#include "vp10/common/vp10_txfm.h"
#include "./vp10_rtcd.h"

namespace libvpx_test {
typedef enum {
  TYPE_DCT = 0,
  TYPE_ADST,
  TYPE_IDCT,
  TYPE_IADST,
  TYPE_LAST
} TYPE_TXFM;

int get_txfm1d_size(TX_SIZE tx_size);

void get_txfm1d_type(TX_TYPE txfm2d_type, TYPE_TXFM* type0,
                     TYPE_TXFM* type1);

void reference_dct_1d(const double* in, double* out, int size);

void reference_adst_1d(const double* in, double* out, int size);

void reference_hybrid_1d(double* in, double* out, int size, int type);

void reference_hybrid_2d(double* in, double* out, int size,
                         int type0, int type1);
template <typename Type1, typename Type2>
static double compute_avg_abs_error(const Type1* a, const Type2* b,
                                    const int size) {
  double error = 0;
  for (int i = 0; i < size; i++) {
    error += fabs(static_cast<double>(a[i]) - static_cast<double>(b[i]));
  }
  error = error / size;
  return error;
}

template<typename Type>
void fliplr(Type *dest, int stride, int length);

template<typename Type>
void flipud(Type *dest, int stride, int length);

template<typename Type>
void fliplrud(Type *dest, int stride, int length);

typedef void (*TxfmFunc)(const int32_t* in, int32_t* out, const int8_t* cos_bit,
                         const int8_t* range_bit);

typedef void (*Fwd_Txfm2d_Func)(const int16_t*, int32_t*, int, int, int);
typedef void (*Inv_Txfm2d_Func)(const int32_t*, uint16_t*, int, int, int);

static const int bd = 10;
static const int input_base = (1 << bd);

#if CONFIG_VPX_HIGHBITDEPTH
static const Fwd_Txfm2d_Func fwd_txfm_func_ls[TX_SIZES] = {
    vp10_fwd_txfm2d_4x4_c, vp10_fwd_txfm2d_8x8_c, vp10_fwd_txfm2d_16x16_c,
    vp10_fwd_txfm2d_32x32_c};

static const Inv_Txfm2d_Func inv_txfm_func_ls[TX_SIZES] = {
    vp10_inv_txfm2d_add_4x4_c, vp10_inv_txfm2d_add_8x8_c,
    vp10_inv_txfm2d_add_16x16_c, vp10_inv_txfm2d_add_32x32_c};
#endif  // CONFIG_VPX_HIGHBITDEPTH

}  // namespace libvpx_test
#endif  // VP10_TXFM_TEST_H_
