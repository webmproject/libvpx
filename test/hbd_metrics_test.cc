/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <stdlib.h>
#include <new>

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/util.h"
#include "./vpx_config.h"
#include "vpx_dsp/psnr.h"
#include "vpx_dsp/ssim.h"
#include "vpx_ports/mem.h"
#include "vpx_ports/msvc.h"
#include "vpx_scale/yv12config.h"


using libvpx_test::ACMRandom;

namespace {

typedef double (*LBDMetricFunc)(const YV12_BUFFER_CONFIG *source,
                                const YV12_BUFFER_CONFIG *dest);
typedef double (*HBDMetricFunc)(const YV12_BUFFER_CONFIG *source,
                                const YV12_BUFFER_CONFIG *dest,
                                uint32_t in_bd, uint32_t bd);

double compute_hbd_psnr(const YV12_BUFFER_CONFIG *source,
                        const YV12_BUFFER_CONFIG *dest,
                        uint32_t in_bd, uint32_t bd) {
  PSNR_STATS psnr;
  vpx_calc_highbd_psnr(source, dest, &psnr, bd, in_bd);
  return psnr.psnr[0];
}

double compute_psnr(const YV12_BUFFER_CONFIG *source,
  const YV12_BUFFER_CONFIG *dest) {
  PSNR_STATS psnr;
  vpx_calc_psnr(source, dest, &psnr);
  return psnr.psnr[0];
}

double compute_hbd_psnrhvs(const YV12_BUFFER_CONFIG *source,
                           const YV12_BUFFER_CONFIG *dest,
                           uint32_t in_bd, uint32_t bd) {
  double tempy, tempu, tempv;
  return vpx_psnrhvs(source, dest,
                     &tempy, &tempu, &tempv, bd, in_bd);
}

double compute_psnrhvs(const YV12_BUFFER_CONFIG *source,
  const YV12_BUFFER_CONFIG *dest) {
  double tempy, tempu, tempv;
  return vpx_psnrhvs(source, dest,
                     &tempy, &tempu, &tempv, 8, 8);
}

double compute_hbd_fastssim(const YV12_BUFFER_CONFIG *source,
                            const YV12_BUFFER_CONFIG *dest,
                            uint32_t in_bd, uint32_t bd) {
  double tempy, tempu, tempv;
  return vpx_calc_fastssim(source, dest,
                               &tempy, &tempu, &tempv, bd, in_bd);
}

double compute_fastssim(const YV12_BUFFER_CONFIG *source,
                        const YV12_BUFFER_CONFIG *dest) {
  double tempy, tempu, tempv;
  return vpx_calc_fastssim(source, dest,
                           &tempy, &tempu, &tempv, 8, 8);
}

double compute_hbd_vpxssim(const YV12_BUFFER_CONFIG *source,
                           const YV12_BUFFER_CONFIG *dest,
                           uint32_t in_bd, uint32_t bd) {
  double ssim, weight;
  ssim = vpx_highbd_calc_ssim(source, dest, &weight, bd, in_bd);
  return 100 * pow(ssim / weight, 8.0);
}

double compute_vpxssim(const YV12_BUFFER_CONFIG *source,
  const YV12_BUFFER_CONFIG *dest) {
  double ssim, weight;
  ssim = vpx_calc_ssim(source, dest, &weight);
  return 100 * pow(ssim / weight, 8.0);
}


class HBDMetricsTestBase {
 public:
  virtual ~HBDMetricsTestBase() {}

 protected:
  void RunAccuracyCheck() {
    const int width = 1920;
    const int height = 1080;
    int i = 0;
    const uint8_t kPixFiller = 128;
    YV12_BUFFER_CONFIG lbd_src, lbd_dst;
    YV12_BUFFER_CONFIG hbd_src, hbd_dst;
    ACMRandom rnd(ACMRandom::DeterministicSeed());
    double lbd_db, hbd_db;

    memset(&lbd_src, 0, sizeof(lbd_src));
    memset(&lbd_dst, 0, sizeof(lbd_dst));
    memset(&hbd_src, 0, sizeof(hbd_src));
    memset(&hbd_dst, 0, sizeof(hbd_dst));

    vpx_alloc_frame_buffer(&lbd_src, width, height, 1, 1, 0, 32, 16);
    vpx_alloc_frame_buffer(&lbd_dst, width, height, 1, 1, 0, 32, 16);
    vpx_alloc_frame_buffer(&hbd_src, width, height, 1, 1, 1, 32, 16);
    vpx_alloc_frame_buffer(&hbd_dst, width, height, 1, 1, 1, 32, 16);

    memset(lbd_src.buffer_alloc, kPixFiller, lbd_src.buffer_alloc_sz);
    while (i < lbd_src.buffer_alloc_sz) {
      uint16_t spel, dpel;
      spel = lbd_src.buffer_alloc[i];
      // Create some distortion for dst buffer.
      dpel = rnd.Rand8();
      lbd_dst.buffer_alloc[i] = (uint8_t)dpel;
      ((uint16_t*)(hbd_src.buffer_alloc))[i] = spel << (bit_depth_ - 8);
      ((uint16_t*)(hbd_dst.buffer_alloc))[i] = dpel << (bit_depth_ - 8);
      i++;
    }

    lbd_db = lbd_metric_(&lbd_src, &lbd_dst);
    hbd_db = hbd_metric_(&hbd_src, &hbd_dst, input_bit_depth_, bit_depth_);
    EXPECT_LE(fabs(lbd_db - hbd_db), threshold_);

    i = 0;
    while (i < lbd_src.buffer_alloc_sz) {
      uint16_t dpel;
      // Create some small distortion for dst buffer.
      dpel = 120 + (rnd.Rand8() >> 4);
      lbd_dst.buffer_alloc[i] = (uint8_t)dpel;
      ((uint16_t*)(hbd_dst.buffer_alloc))[i] = dpel << (bit_depth_ - 8);
      i++;
    }

    lbd_db = lbd_metric_(&lbd_src, &lbd_dst);
    hbd_db = hbd_metric_(&hbd_src, &hbd_dst, input_bit_depth_, bit_depth_);
    EXPECT_LE(fabs(lbd_db - hbd_db), threshold_);

    i = 0;
    while (i < lbd_src.buffer_alloc_sz) {
      uint16_t dpel;
      // Create some small distortion for dst buffer.
      dpel = 126 + (rnd.Rand8() >> 6);
      lbd_dst.buffer_alloc[i] = (uint8_t)dpel;
      ((uint16_t*)(hbd_dst.buffer_alloc))[i] = dpel << (bit_depth_ - 8);
      i++;
    }

    lbd_db = lbd_metric_(&lbd_src, &lbd_dst);
    hbd_db = hbd_metric_(&hbd_src, &hbd_dst, input_bit_depth_, bit_depth_);
    EXPECT_LE(fabs(lbd_db - hbd_db), threshold_);

    vpx_free_frame_buffer(&lbd_src);
    vpx_free_frame_buffer(&lbd_dst);
    vpx_free_frame_buffer(&hbd_src);
    vpx_free_frame_buffer(&hbd_dst);
  }

  int input_bit_depth_;
  int bit_depth_;
  double threshold_;
  LBDMetricFunc lbd_metric_;
  HBDMetricFunc hbd_metric_;
};

typedef std::tr1::tuple<LBDMetricFunc,
                        HBDMetricFunc, int, int, double> MetricTestTParam;
class HBDMetricsTest
    : public HBDMetricsTestBase,
      public ::testing::TestWithParam<MetricTestTParam> {
 public:
  virtual void SetUp() {
    lbd_metric_ = GET_PARAM(0);
    hbd_metric_ = GET_PARAM(1);
    input_bit_depth_ = GET_PARAM(2);
    bit_depth_ = GET_PARAM(3);
    threshold_ = GET_PARAM(4);
  }
  virtual void TearDown() {}
};

TEST_P(HBDMetricsTest, RunAccuracyCheck) {
  RunAccuracyCheck();
}

// Allow small variation due to floating point operations.
static const double kSsim_thresh = 0.001;
// Allow some additional errors accumulated in floating point operations.
static const double kFSsim_thresh = 0.03;
// Allow some extra variation due to rounding error accumulated in dct.
static const double kPhvs_thresh = 0.3;

INSTANTIATE_TEST_CASE_P(
    VPXSSIM, HBDMetricsTest,
    ::testing::Values(
        MetricTestTParam(&compute_vpxssim, &compute_hbd_vpxssim, 8, 10,
                         kSsim_thresh),
        MetricTestTParam(&compute_vpxssim, &compute_hbd_vpxssim, 10, 10,
                         kPhvs_thresh),
        MetricTestTParam(&compute_vpxssim, &compute_hbd_vpxssim, 8, 12,
                         kSsim_thresh),
        MetricTestTParam(&compute_vpxssim, &compute_hbd_vpxssim, 12, 12,
                         kPhvs_thresh)));
INSTANTIATE_TEST_CASE_P(
    FASTSSIM, HBDMetricsTest,
    ::testing::Values(
        MetricTestTParam(&compute_fastssim, &compute_hbd_fastssim, 8, 10,
                         kFSsim_thresh),
        MetricTestTParam(&compute_fastssim, &compute_hbd_fastssim, 10, 10,
                         kFSsim_thresh),
        MetricTestTParam(&compute_fastssim, &compute_hbd_fastssim, 8, 12,
                         kFSsim_thresh),
        MetricTestTParam(&compute_fastssim, &compute_hbd_fastssim, 12, 12,
                         kFSsim_thresh)));
INSTANTIATE_TEST_CASE_P(
    PSNRHVS, HBDMetricsTest,
    ::testing::Values(
        MetricTestTParam(&compute_psnrhvs, &compute_hbd_psnrhvs, 8, 10,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnrhvs, &compute_hbd_psnrhvs, 10, 10,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnrhvs, &compute_hbd_psnrhvs, 8, 12,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnrhvs, &compute_hbd_psnrhvs, 12, 12,
                         kPhvs_thresh)));
INSTANTIATE_TEST_CASE_P(
    PSNR, HBDMetricsTest,
    ::testing::Values(
        MetricTestTParam(&compute_psnr, &compute_hbd_psnr, 8, 10,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnr, &compute_hbd_psnr, 10, 10,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnr, &compute_hbd_psnr, 8, 12,
                         kPhvs_thresh),
        MetricTestTParam(&compute_psnr, &compute_hbd_psnr, 12, 12,
                         kPhvs_thresh)));
}  // namespace

