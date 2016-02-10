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
#include "vpx_dsp/ssim.h"
#include "vpx_ports/mem.h"
#include "vpx_ports/msvc.h"
#include "vpx_scale/yv12config.h"


using libvpx_test::ACMRandom;

namespace {

typedef double (*LBDMetricFunc)(const YV12_BUFFER_CONFIG *source,
                                const YV12_BUFFER_CONFIG *dest,
                                double *weight);
typedef double (*HBDMetricFunc)(const YV12_BUFFER_CONFIG *source,
                                const YV12_BUFFER_CONFIG *dest,
                                double *weight, unsigned int bd);

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
    double lbd_score, hbd_score, lbd_db, hbd_db, lbd_w, hbd_w;

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
      lbd_dst.buffer_alloc[i] = rnd.Rand8();
      dpel = lbd_dst.buffer_alloc[i];
      ((uint16_t*)(hbd_src.buffer_alloc))[i] = spel << (bit_depth_ - 8);
      ((uint16_t*)(hbd_dst.buffer_alloc))[i] = dpel << (bit_depth_ - 8);
      i++;
    }

    lbd_score = lbd_metric_(&lbd_src, &lbd_dst, &lbd_w);
    hbd_score = hbd_metric_(&hbd_src, &hbd_dst, &hbd_w, bit_depth_);

    lbd_db = 100 * pow(lbd_score / lbd_w, 8.0);
    hbd_db = 100 * pow(hbd_score / hbd_w, 8.0);

    vpx_free_frame_buffer(&lbd_src);
    vpx_free_frame_buffer(&lbd_dst);
    vpx_free_frame_buffer(&hbd_src);
    vpx_free_frame_buffer(&hbd_dst);

    EXPECT_LE(fabs(lbd_db - hbd_db), threshold_);
  }

  int bit_depth_;
  double threshold_;
  LBDMetricFunc lbd_metric_;
  HBDMetricFunc hbd_metric_;
};

typedef std::tr1::tuple<LBDMetricFunc,
                        HBDMetricFunc, int, double> MetricTestTParam;
class HBDMetricsTest
    : public HBDMetricsTestBase,
      public ::testing::TestWithParam<MetricTestTParam> {
 public:
  virtual void SetUp() {
    lbd_metric_ = GET_PARAM(0);
    hbd_metric_ = GET_PARAM(1);
    bit_depth_ = GET_PARAM(2);
    threshold_ = GET_PARAM(3);
  }
  virtual void TearDown() {}
};

TEST_P(HBDMetricsTest, RunAccuracyCheck) {
  RunAccuracyCheck();
}

// Allow small variation due to floating point operations.
static const double kSsim_thresh = 0.001;

INSTANTIATE_TEST_CASE_P(
    C, HBDMetricsTest,
    ::testing::Values(
        MetricTestTParam(&vpx_calc_ssim, &vpx_highbd_calc_ssim, 10,
                         kSsim_thresh),
        MetricTestTParam(&vpx_calc_ssim, &vpx_highbd_calc_ssim, 12,
                         kSsim_thresh)));
}  // namespace

