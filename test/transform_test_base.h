/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef TEST_TRANSFORM_TEST_BASE_H_
#define TEST_TRANSFORM_TEST_BASE_H_

#include "./vpx_config.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx/vpx_codec.h"


namespace libvpx_test {

//  Note:
//   Same constant are defined in vp9/common/vp9_entropy.h and
//   vp10/common/entropy.h.  Goal is to make this base class
//   to use for future codec transform testing.  But including
//   either of them would lead to compiling error when we do
//   unit test for another codec. Suggest to move the definition
//   to a vpx header file.
const int kDctMaxValue = 16384;

typedef void (*FhtFunc)(const int16_t *in, tran_low_t *out, int stride,
                        int tx_type);

class TransformTestBase {
 public:
  virtual ~TransformTestBase() {}

 protected:
  virtual void RunFwdTxfm(const int16_t *in, tran_low_t *out, int stride) = 0;

  virtual void RunInvTxfm(const tran_low_t *out, uint8_t *dst, int stride) = 0;

  void RunAccuracyCheck(int limit) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    uint32_t max_error = 0;
    int64_t total_error = 0;
    const int count_test_block = 10000;

    int16_t *test_input_block = reinterpret_cast<int16_t *>
        (vpx_memalign(16, sizeof(int16_t) * num_coeffs_));
    tran_low_t *test_temp_block = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));
    uint8_t *dst = reinterpret_cast<uint8_t *>
        (vpx_memalign(16, sizeof(uint8_t) * num_coeffs_));
    uint8_t *src = reinterpret_cast<uint8_t *>
        (vpx_memalign(16, sizeof(uint8_t) * num_coeffs_));
#if CONFIG_VPX_HIGHBITDEPTH
    uint16_t *dst16 = reinterpret_cast<uint16_t *>
        (vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
    uint16_t *src16 = reinterpret_cast<uint16_t *>
        (vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
#endif

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-255, 255].
      for (int j = 0; j < num_coeffs_; ++j) {
        if (bit_depth_ == VPX_BITS_8) {
          src[j] = rnd.Rand8();
          dst[j] = rnd.Rand8();
          test_input_block[j] = src[j] - dst[j];
#if CONFIG_VPX_HIGHBITDEPTH
        } else {
          src16[j] = rnd.Rand16() & mask_;
          dst16[j] = rnd.Rand16() & mask_;
          test_input_block[j] = src16[j] - dst16[j];
#endif
        }
      }

      ASM_REGISTER_STATE_CHECK(RunFwdTxfm(test_input_block,
                                          test_temp_block, pitch_));
      if (bit_depth_ == VPX_BITS_8) {
        ASM_REGISTER_STATE_CHECK(RunInvTxfm(test_temp_block, dst, pitch_));
#if CONFIG_VPX_HIGHBITDEPTH
      } else {
        ASM_REGISTER_STATE_CHECK(RunInvTxfm(test_temp_block,
                                            CONVERT_TO_BYTEPTR(dst16), pitch_));
#endif
      }

      for (int j = 0; j < num_coeffs_; ++j) {
#if CONFIG_VPX_HIGHBITDEPTH
        const uint32_t diff =
            bit_depth_ == VPX_BITS_8 ? dst[j] - src[j] : dst16[j] - src16[j];
#else
        ASSERT_EQ(VPX_BITS_8, bit_depth_);
        const uint32_t diff = dst[j] - src[j];
#endif
        const uint32_t error = diff * diff;
        if (max_error < error)
          max_error = error;
        total_error += error;
      }
    }

    EXPECT_GE(static_cast<uint32_t>(limit), max_error)
        << "Error: 4x4 FHT/IHT has an individual round trip error > "
        << limit;

    EXPECT_GE(count_test_block * limit, total_error)
        << "Error: 4x4 FHT/IHT has average round trip error > " << limit
        << " per block";

    vpx_free(test_input_block);
    vpx_free(test_temp_block);
    vpx_free(dst);
    vpx_free(src);
#if CONFIG_VPX_HIGHBITDEPTH
    vpx_free(dst16);
    vpx_free(src16);
#endif
  }

  void RunCoeffCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 5000;

    int16_t *input_block = reinterpret_cast<int16_t *>
        (vpx_memalign(16, sizeof(int16_t) * num_coeffs_));
    tran_low_t *output_ref_block = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));
    tran_low_t *output_block = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-mask_, mask_].
      for (int j = 0; j < num_coeffs_; ++j)
        input_block[j] = (rnd.Rand16() & mask_) - (rnd.Rand16() & mask_);

      fwd_txfm_ref(input_block, output_ref_block, pitch_, tx_type_);
      ASM_REGISTER_STATE_CHECK(RunFwdTxfm(input_block, output_block, pitch_));

      // The minimum quant value is 4.
      for (int j = 0; j < num_coeffs_; ++j) {
        EXPECT_EQ(output_block[j], output_ref_block[j])
            << "Error: not bit-exact result at index: " << j
            << " at test block: " << i;
      }
    }
    vpx_free(input_block);
    vpx_free(output_ref_block);
    vpx_free(output_block);
  }

  void RunMemCheck() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 5000;

    int16_t *input_extreme_block = reinterpret_cast<int16_t *>
        (vpx_memalign(16, sizeof(int16_t) * num_coeffs_));
    tran_low_t *output_ref_block = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));
    tran_low_t *output_block = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-mask_, mask_].
      for (int j = 0; j < num_coeffs_; ++j) {
        input_extreme_block[j] = rnd.Rand8() % 2 ? mask_ : -mask_;
      }
      if (i == 0) {
        for (int j = 0; j < num_coeffs_; ++j)
          input_extreme_block[j] = mask_;
      } else if (i == 1) {
        for (int j = 0; j < num_coeffs_; ++j)
          input_extreme_block[j] = -mask_;
      }

      fwd_txfm_ref(input_extreme_block, output_ref_block, pitch_, tx_type_);
      ASM_REGISTER_STATE_CHECK(RunFwdTxfm(input_extreme_block,
                                          output_block, pitch_));

      int row_length = FindRowLength();
      // The minimum quant value is 4.
      for (int j = 0; j < num_coeffs_; ++j) {
        EXPECT_EQ(output_block[j], output_ref_block[j]);
        EXPECT_GE(row_length * kDctMaxValue << (bit_depth_ - 8),
                  abs(output_block[j]))
            << "Error: NxN FDCT has coefficient larger than N*DCT_MAX_VALUE";
      }
    }
    vpx_free(input_extreme_block);
    vpx_free(output_ref_block);
    vpx_free(output_block);
  }

  void RunInvAccuracyCheck(int limit) {
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 1000;

    int16_t *in = reinterpret_cast<int16_t *>
        (vpx_memalign(16, sizeof(int16_t) * num_coeffs_));
    tran_low_t *coeff = reinterpret_cast<tran_low_t *>
        (vpx_memalign(16, sizeof(tran_low_t) * num_coeffs_));
    uint8_t *dst = reinterpret_cast<uint8_t *>
        (vpx_memalign(16, sizeof(uint8_t) * num_coeffs_));
    uint8_t *src = reinterpret_cast<uint8_t *>
        (vpx_memalign(16, sizeof(uint8_t) * num_coeffs_));

#if CONFIG_VPX_HIGHBITDEPTH
    uint16_t *dst16 = reinterpret_cast<uint16_t *>
        (vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
    uint16_t *src16 = reinterpret_cast<uint16_t *>
        (vpx_memalign(16, sizeof(uint16_t) * num_coeffs_));
#endif

    for (int i = 0; i < count_test_block; ++i) {
      // Initialize a test block with input range [-mask_, mask_].
      for (int j = 0; j < num_coeffs_; ++j) {
        if (bit_depth_ == VPX_BITS_8) {
          src[j] = rnd.Rand8();
          dst[j] = rnd.Rand8();
          in[j] = src[j] - dst[j];
#if CONFIG_VPX_HIGHBITDEPTH
        } else {
          src16[j] = rnd.Rand16() & mask_;
          dst16[j] = rnd.Rand16() & mask_;
          in[j] = src16[j] - dst16[j];
#endif
        }
      }

      fwd_txfm_ref(in, coeff, pitch_, tx_type_);

      if (bit_depth_ == VPX_BITS_8) {
        ASM_REGISTER_STATE_CHECK(RunInvTxfm(coeff, dst, pitch_));
#if CONFIG_VPX_HIGHBITDEPTH
      } else {
        ASM_REGISTER_STATE_CHECK(RunInvTxfm(coeff, CONVERT_TO_BYTEPTR(dst16),
                                            pitch_));
#endif
      }

      for (int j = 0; j < num_coeffs_; ++j) {
#if CONFIG_VPX_HIGHBITDEPTH
        const uint32_t diff =
            bit_depth_ == VPX_BITS_8 ? dst[j] - src[j] : dst16[j] - src16[j];
#else
        const uint32_t diff = dst[j] - src[j];
#endif
        const uint32_t error = diff * diff;
        EXPECT_GE(static_cast<uint32_t>(limit), error)
            << "Error: 4x4 IDCT has error " << error
            << " at index " << j;
      }
    }
    vpx_free(in);
    vpx_free(coeff);
    vpx_free(dst);
    vpx_free(src);
#if CONFIG_VPX_HIGHBITDEPTH
    vpx_free(src16);
    vpx_free(dst16);
#endif
  }

  int pitch_;
  int tx_type_;
  FhtFunc fwd_txfm_ref;
  vpx_bit_depth_t bit_depth_;
  int mask_;
  int num_coeffs_;

 private:
  //  Assume transform size is 4x4, 8x8, 16x16,...
  int FindRowLength() const {
    int row = 4;
    if (16 == num_coeffs_) {
      row = 4;
    } else if (64 == num_coeffs_) {
      row = 8;
    } else if (256 == num_coeffs_) {
      row = 16;
    } else if (1024 == num_coeffs_) {
      row = 32;
    }
    return row;
  }
};

}  // namespace libvpx_test

#endif  // TEST_TRANSFORM_TEST_BASE_H_
