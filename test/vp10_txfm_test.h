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
#include <math.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "vp10/common/vp10_txfm.h"

typedef enum {
  TYPE_DCT = 0,
  TYPE_ADST,
  TYPE_IDCT,
  TYPE_IADST,
  TYPE_LAST
} TYPE_TXFM;

static double invSqrt2 = 1 / pow(2, 0.5);

static void reference_dct_1d(const double* in, double* out, int size) {
  for (int k = 0; k < size; ++k) {
    out[k] = 0;
    for (int n = 0; n < size; ++n) {
      out[k] += in[n] * cos(M_PI * (2 * n + 1) * k / (2 * size));
    }
    if (k == 0) out[k] = out[k] * invSqrt2;
  }
}

static void reference_adst_1d(const double* in, double* out, int size) {
  for (int k = 0; k < size; ++k) {
    out[k] = 0;
    for (int n = 0; n < size; ++n) {
      out[k] += in[n] * sin(M_PI * (2 * n + 1) * (2 * k + 1) / (4 * size));
    }
  }
}

static void reference_hybrid_1d(double* in, double* out, int size, int type) {
  if (type == TYPE_DCT)
    reference_dct_1d(in, out, size);
  else
    reference_adst_1d(in, out, size);
}

static void reference_hybrid_2d(double* in, double* out, int size, int type0,
                                int type1) {
  double* tempOut = new double[size * size];

  for (int r = 0; r < size; r++) {
    // out ->tempOut
    for (int c = 0; c < size; c++) {
      tempOut[r * size + c] = in[c * size + r];
    }
  }

  // dct each row: in -> out
  for (int r = 0; r < size; r++) {
    reference_hybrid_1d(tempOut + r * size, out + r * size, size, type0);
  }

  for (int r = 0; r < size; r++) {
    // out ->tempOut
    for (int c = 0; c < size; c++) {
      tempOut[r * size + c] = out[c * size + r];
    }
  }

  for (int r = 0; r < size; r++) {
    reference_hybrid_1d(tempOut + r * size, out + r * size, size, type1);
  }
  delete[] tempOut;
}

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

typedef void (*TxfmFunc)(const int32_t* in, int32_t* out, const int8_t* cos_bit,
                         const int8_t* range_bit);

typedef void (*Fwd_Txfm2d_Func)(const int16_t*, int32_t*, const int,
                                const TXFM_2D_CFG*, const int);
typedef void (*Inv_Txfm2d_Func)(const int32_t*, uint16_t*, const int,
                                const TXFM_2D_CFG*, const int);

static const int bd = 10;
static const int base = (1 << bd);

#endif  // VP10_TXFM_TEST_H_
