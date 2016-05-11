/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdio.h>
#include "test/vp10_txfm_test.h"

namespace libvpx_test {

int get_txfm1d_size(TX_SIZE tx_size) {
  return 1 << (tx_size + 2);
}

void get_txfm1d_type(TX_TYPE txfm2d_type, TYPE_TXFM* type0,
                     TYPE_TXFM* type1) {
  switch (txfm2d_type) {
    case DCT_DCT:
      *type0 = TYPE_DCT;
      *type1 = TYPE_DCT;
      break;
    case ADST_DCT:
      *type0 = TYPE_ADST;
      *type1 = TYPE_DCT;
      break;
    case DCT_ADST:
      *type0 = TYPE_DCT;
      *type1 = TYPE_ADST;
      break;
    case ADST_ADST:
      *type0 = TYPE_ADST;
      *type1 = TYPE_ADST;
      break;
    default:
      *type0 = TYPE_DCT;
      *type1 = TYPE_DCT;
      assert(0);
      break;
  }
}

double invSqrt2 = 1 / pow(2, 0.5);

void reference_dct_1d(const double* in, double* out, int size) {
  for (int k = 0; k < size; ++k) {
    out[k] = 0;
    for (int n = 0; n < size; ++n) {
      out[k] += in[n] * cos(M_PI * (2 * n + 1) * k / (2 * size));
    }
    if (k == 0) out[k] = out[k] * invSqrt2;
  }
}

void reference_adst_1d(const double* in, double* out, int size) {
  for (int k = 0; k < size; ++k) {
    out[k] = 0;
    for (int n = 0; n < size; ++n) {
      out[k] += in[n] * sin(M_PI * (2 * n + 1) * (2 * k + 1) / (4 * size));
    }
  }
}

void reference_hybrid_1d(double* in, double* out, int size, int type) {
  if (type == TYPE_DCT)
    reference_dct_1d(in, out, size);
  else
    reference_adst_1d(in, out, size);
}

void reference_hybrid_2d(double* in, double* out, int size,
                         int type0, int type1) {
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
}  // namespace libvpx_test
