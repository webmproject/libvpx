/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./vp9_rtcd.h"
#include "./vpx_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/transform_test_base.h"
#include "test/util.h"
#include "vpx/vpx_codec.h"
#include "vpx/vpx_integer.h"
#include "vpx_ports/mem.h"

using libvpx_test::ACMRandom;

namespace {
typedef void (*FdctFunc)(const int16_t *in, tran_low_t *out, int stride);
typedef void (*IdctFunc)(const tran_low_t *in, uint8_t *out, int stride);
typedef void (*IhtFunc)(const tran_low_t *in, uint8_t *out, int stride,
                        int tx_type);
using libvpx_test::FhtFunc;

typedef std::tr1::tuple<FdctFunc, IdctFunc, int, vpx_bit_depth_t, int>
Dct4x4Param;
typedef std::tr1::tuple<FhtFunc, IhtFunc, int, vpx_bit_depth_t, int>
Ht4x4Param;

void fdct4x4_ref(const int16_t *in, tran_low_t *out, int stride,
                 int /*tx_type*/) {
  vpx_fdct4x4_c(in, out, stride);
}

void fht4x4_ref(const int16_t *in, tran_low_t *out, int stride, int tx_type) {
  vp9_fht4x4_c(in, out, stride, tx_type);
}

void fwht4x4_ref(const int16_t *in, tran_low_t *out, int stride,
                 int /*tx_type*/) {
  vp9_fwht4x4_c(in, out, stride);
}

#if CONFIG_VP9_HIGHBITDEPTH
void idct4x4_10(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_idct4x4_16_add_c(in, out, stride, 10);
}

void idct4x4_12(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_idct4x4_16_add_c(in, out, stride, 12);
}

void iht4x4_10(const tran_low_t *in, uint8_t *out, int stride, int tx_type) {
  vp9_highbd_iht4x4_16_add_c(in, out, stride, tx_type, 10);
}

void iht4x4_12(const tran_low_t *in, uint8_t *out, int stride, int tx_type) {
  vp9_highbd_iht4x4_16_add_c(in, out, stride, tx_type, 12);
}

void iwht4x4_10(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_iwht4x4_16_add_c(in, out, stride, 10);
}

void iwht4x4_12(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_iwht4x4_16_add_c(in, out, stride, 12);
}

#if HAVE_SSE2
void idct4x4_10_sse2(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_idct4x4_16_add_sse2(in, out, stride, 10);
}

void idct4x4_12_sse2(const tran_low_t *in, uint8_t *out, int stride) {
  vpx_highbd_idct4x4_16_add_sse2(in, out, stride, 12);
}
#endif  // HAVE_SSE2
#endif  // CONFIG_VP9_HIGHBITDEPTH


class Trans4x4DCT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<Dct4x4Param> {
 public:
  virtual ~Trans4x4DCT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 4;
    fwd_txfm_ref = fdct4x4_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride);
  }
  void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride);
  }

  FdctFunc fwd_txfm_;
  IdctFunc inv_txfm_;
};

TEST_P(Trans4x4DCT, AccuracyCheck) {
  RunAccuracyCheck(1);
}

TEST_P(Trans4x4DCT, CoeffCheck) {
  RunCoeffCheck();
}

TEST_P(Trans4x4DCT, MemCheck) {
  RunMemCheck();
}

TEST_P(Trans4x4DCT, InvAccuracyCheck) {
  RunInvAccuracyCheck(1);
}

class Trans4x4HT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<Ht4x4Param> {
 public:
  virtual ~Trans4x4HT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 4;
    fwd_txfm_ref = fht4x4_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride, tx_type_);
  }

  void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride, tx_type_);
  }

  FhtFunc fwd_txfm_;
  IhtFunc inv_txfm_;
};

TEST_P(Trans4x4HT, AccuracyCheck) {
  RunAccuracyCheck(1);
}

TEST_P(Trans4x4HT, CoeffCheck) {
  RunCoeffCheck();
}

TEST_P(Trans4x4HT, MemCheck) {
  RunMemCheck();
}

TEST_P(Trans4x4HT, InvAccuracyCheck) {
  RunInvAccuracyCheck(1);
}

class Trans4x4WHT
    : public libvpx_test::TransformTestBase,
      public ::testing::TestWithParam<Dct4x4Param> {
 public:
  virtual ~Trans4x4WHT() {}

  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    inv_txfm_ = GET_PARAM(1);
    tx_type_  = GET_PARAM(2);
    pitch_    = 4;
    fwd_txfm_ref = fwht4x4_ref;
    bit_depth_ = GET_PARAM(3);
    mask_ = (1 << bit_depth_) - 1;
    num_coeffs_ = GET_PARAM(4);
  }
  virtual void TearDown() { libvpx_test::ClearSystemState(); }

 protected:
  void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) {
    fwd_txfm_(in, out, stride);
  }
  void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) {
    inv_txfm_(out, dst, stride);
  }

  FdctFunc fwd_txfm_;
  IdctFunc inv_txfm_;
};

TEST_P(Trans4x4WHT, AccuracyCheck) {
  RunAccuracyCheck(0);
}

TEST_P(Trans4x4WHT, CoeffCheck) {
  RunCoeffCheck();
}

TEST_P(Trans4x4WHT, MemCheck) {
  RunMemCheck();
}

TEST_P(Trans4x4WHT, InvAccuracyCheck) {
  RunInvAccuracyCheck(0);
}
using std::tr1::make_tuple;

#if CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4DCT,
    ::testing::Values(
        make_tuple(&vpx_highbd_fdct4x4_c, &idct4x4_10, 0, VPX_BITS_10, 16),
        make_tuple(&vpx_highbd_fdct4x4_c, &idct4x4_12, 0, VPX_BITS_12, 16),
        make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c, 0, VPX_BITS_8, 16)));
#else
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4DCT,
    ::testing::Values(
         make_tuple(&vpx_fdct4x4_c, &vpx_idct4x4_16_add_c, 0, VPX_BITS_8, 16)));
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_10, 0, VPX_BITS_10, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_10, 1, VPX_BITS_10, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_10, 2, VPX_BITS_10, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_10, 3, VPX_BITS_10, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_12, 0, VPX_BITS_12, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_12, 1, VPX_BITS_12, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_12, 2, VPX_BITS_12, 16),
        make_tuple(&vp9_highbd_fht4x4_c, &iht4x4_12, 3, VPX_BITS_12, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 0, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 1, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 2, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 3, VPX_BITS_8, 16)));
#else
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 0, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 1, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 2, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_c, 3, VPX_BITS_8, 16)));
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if CONFIG_VP9_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4WHT,
    ::testing::Values(
        make_tuple(&vp9_highbd_fwht4x4_c, &iwht4x4_10, 0, VPX_BITS_10, 16),
        make_tuple(&vp9_highbd_fwht4x4_c, &iwht4x4_12, 0, VPX_BITS_12, 16),
        make_tuple(&vp9_fwht4x4_c, &vpx_iwht4x4_16_add_c, 0, VPX_BITS_8, 16)));
#else
INSTANTIATE_TEST_CASE_P(
    C, Trans4x4WHT,
    ::testing::Values(
        make_tuple(&vp9_fwht4x4_c, &vpx_iwht4x4_16_add_c, 0, VPX_BITS_8, 16)));
#endif  // CONFIG_VP9_HIGHBITDEPTH

#if HAVE_NEON_ASM && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    NEON, Trans4x4DCT,
    ::testing::Values(
        make_tuple(&vpx_fdct4x4_c,
                   &vpx_idct4x4_16_add_neon, 0, VPX_BITS_8, 16)));
#endif  // HAVE_NEON_ASM && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

#if HAVE_NEON && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    NEON, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_neon, 0, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_neon, 1, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_neon, 2, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_c, &vp9_iht4x4_16_add_neon, 3, VPX_BITS_8, 16)));
#endif  // HAVE_NEON && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

#if CONFIG_USE_X86INC && HAVE_SSE2 && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSE2, Trans4x4WHT,
    ::testing::Values(
        make_tuple(&vp9_fwht4x4_sse2, &vpx_iwht4x4_16_add_c, 0, VPX_BITS_8, 16),
        make_tuple(&vp9_fwht4x4_c, &vpx_iwht4x4_16_add_sse2, 0, VPX_BITS_8, 16)));
#endif

#if HAVE_SSE2 && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSE2, Trans4x4DCT,
    ::testing::Values(
        make_tuple(&vpx_fdct4x4_sse2,
                   &vpx_idct4x4_16_add_sse2, 0, VPX_BITS_8, 16)));
INSTANTIATE_TEST_CASE_P(
    SSE2, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_sse2, 0,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_sse2, 1,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_sse2, 2,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_sse2, 3,
                   VPX_BITS_8, 16)));
#endif  // HAVE_SSE2 && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

#if HAVE_SSE2 && CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    SSE2, Trans4x4DCT,
    ::testing::Values(
        make_tuple(&vpx_highbd_fdct4x4_c,    &idct4x4_10_sse2, 0,
                   VPX_BITS_10, 16),
        make_tuple(&vpx_highbd_fdct4x4_sse2, &idct4x4_10_sse2, 0,
                   VPX_BITS_10, 16),
        make_tuple(&vpx_highbd_fdct4x4_c,    &idct4x4_12_sse2, 0,
                   VPX_BITS_12, 16),
        make_tuple(&vpx_highbd_fdct4x4_sse2, &idct4x4_12_sse2, 0,
                   VPX_BITS_12, 16),
        make_tuple(&vpx_fdct4x4_sse2,      &vpx_idct4x4_16_add_c, 0,
                   VPX_BITS_8, 16)));

INSTANTIATE_TEST_CASE_P(
    SSE2, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_c, 0, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_c, 1, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_c, 2, VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_sse2, &vp9_iht4x4_16_add_c, 3, VPX_BITS_8, 16)));
#endif  // HAVE_SSE2 && CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE

#if HAVE_MSA && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
INSTANTIATE_TEST_CASE_P(
    MSA, Trans4x4DCT,
    ::testing::Values(
        make_tuple(&vpx_fdct4x4_msa, &vpx_idct4x4_16_add_msa, 0,
                   VPX_BITS_8, 16)));
INSTANTIATE_TEST_CASE_P(
    MSA, Trans4x4HT,
    ::testing::Values(
        make_tuple(&vp9_fht4x4_msa, &vp9_iht4x4_16_add_msa, 0,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_msa, &vp9_iht4x4_16_add_msa, 1,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_msa, &vp9_iht4x4_16_add_msa, 2,
                   VPX_BITS_8, 16),
        make_tuple(&vp9_fht4x4_msa, &vp9_iht4x4_16_add_msa, 3,
                   VPX_BITS_8, 16)));
#endif  // HAVE_MSA && !CONFIG_VP9_HIGHBITDEPTH && !CONFIG_EMULATE_HARDWARE
}  // namespace
