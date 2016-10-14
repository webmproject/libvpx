/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <utility>
#include <vector>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "aom_dsp/ansreader.h"
#include "aom_dsp/answriter.h"

namespace {
typedef std::vector<std::pair<uint8_t, bool> > PvVec;

const int kPrintStats = 0;

PvVec abs_encode_build_vals(int iters) {
  PvVec ret;
  libaom_test::ACMRandom gen(0x30317076);
  double entropy = 0;
  for (int i = 0; i < iters; ++i) {
    uint8_t p;
    do {
      p = gen.Rand8();
    } while (p == 0);  // zero is not a valid coding probability
    bool b = gen.Rand8() < p;
    ret.push_back(std::make_pair(static_cast<uint8_t>(p), b));
    if (kPrintStats) {
      double d = p / 256.;
      entropy += -d * log2(d) - (1 - d) * log2(1 - d);
    }
  }
  if (kPrintStats) printf("entropy %f\n", entropy);
  return ret;
}

bool check_uabs(const PvVec &pv_vec, uint8_t *buf) {
  AnsCoder a;
  ans_write_init(&a, buf);

  std::clock_t start = std::clock();
  for (PvVec::const_reverse_iterator it = pv_vec.rbegin(); it != pv_vec.rend();
       ++it) {
    uabs_write(&a, it->second, 256 - it->first);
  }
  std::clock_t enc_time = std::clock() - start;
  int offset = ans_write_end(&a);
  bool okay = true;
  AnsDecoder d;
  if (ans_read_init(&d, buf, offset)) return false;
  start = std::clock();
  for (PvVec::const_iterator it = pv_vec.begin(); it != pv_vec.end(); ++it) {
    okay &= uabs_read(&d, 256 - it->first) == it->second;
  }
  std::clock_t dec_time = std::clock() - start;
  if (!okay) return false;
  if (kPrintStats)
    printf("uABS size %d enc_time %f dec_time %f\n", offset,
           static_cast<float>(enc_time) / CLOCKS_PER_SEC,
           static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
}

const aom_cdf_prob spareto65[] = { 8320, 6018, 4402, 3254, 4259,
                                   3919, 2057, 492,  45,   2 };

const int kRansSymbols =
    static_cast<int>(sizeof(spareto65) / sizeof(spareto65[0]));

std::vector<int> ans_encode_build_vals(rans_sym *const tab, int iters) {
  aom_cdf_prob sum = 0;
  for (int i = 0; i < kRansSymbols; ++i) {
    tab[i].cum_prob = sum;
    tab[i].prob = spareto65[i];
    sum += spareto65[i];
  }
  std::vector<int> p_to_sym;
  for (int i = 0; i < kRansSymbols; ++i) {
    p_to_sym.insert(p_to_sym.end(), tab[i].prob, i);
  }
  assert(p_to_sym.size() == RANS_PRECISION);
  std::vector<int> ret;
  libaom_test::ACMRandom gen(18543637);
  for (int i = 0; i < iters; ++i) {
    int sym =
        p_to_sym[((gen.Rand8() << 8) + gen.Rand8()) & (RANS_PRECISION - 1)];
    ret.push_back(sym);
  }
  return ret;
}

void rans_build_dec_tab(const struct rans_sym sym_tab[],
                        aom_cdf_prob *dec_tab) {
  unsigned int sum = 0;
  for (int i = 0; sum < RANS_PRECISION; ++i) {
    dec_tab[i] = sum += sym_tab[i].prob;
  }
}

bool check_rans(const std::vector<int> &sym_vec, const rans_sym *const tab,
                uint8_t *buf) {
  AnsCoder a;
  ans_write_init(&a, buf);
  aom_cdf_prob dec_tab[kRansSymbols];
  rans_build_dec_tab(tab, dec_tab);

  std::clock_t start = std::clock();
  for (std::vector<int>::const_reverse_iterator it = sym_vec.rbegin();
       it != sym_vec.rend(); ++it) {
    rans_write(&a, &tab[*it]);
  }
  std::clock_t enc_time = std::clock() - start;
  int offset = ans_write_end(&a);
  bool okay = true;
  AnsDecoder d;
  if (ans_read_init(&d, buf, offset)) return false;
  start = std::clock();
  for (std::vector<int>::const_iterator it = sym_vec.begin();
       it != sym_vec.end(); ++it) {
    okay &= rans_read(&d, dec_tab) == *it;
  }
  std::clock_t dec_time = std::clock() - start;
  if (!okay) return false;
  if (kPrintStats)
    printf("rANS size %d enc_time %f dec_time %f\n", offset,
           static_cast<float>(enc_time) / CLOCKS_PER_SEC,
           static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
}

class AbsTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() { pv_vec_ = abs_encode_build_vals(kNumBools); }
  virtual void SetUp() { buf_ = new uint8_t[kNumBools / 8]; }
  virtual void TearDown() { delete[] buf_; }
  static const int kNumBools = 100000000;
  static PvVec pv_vec_;
  uint8_t *buf_;
};
PvVec AbsTest::pv_vec_;

class AnsTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    sym_vec_ = ans_encode_build_vals(rans_sym_tab_, kNumSyms);
  }
  virtual void SetUp() { buf_ = new uint8_t[kNumSyms / 2]; }
  virtual void TearDown() { delete[] buf_; }
  static const int kNumSyms = 25000000;
  static std::vector<int> sym_vec_;
  static rans_sym rans_sym_tab_[kRansSymbols];
  uint8_t *buf_;
};
std::vector<int> AnsTest::sym_vec_;
rans_sym AnsTest::rans_sym_tab_[kRansSymbols];

TEST_F(AbsTest, Uabs) { EXPECT_TRUE(check_uabs(pv_vec_, buf_)); }
TEST_F(AnsTest, Rans) {
  EXPECT_TRUE(check_rans(sym_vec_, rans_sym_tab_, buf_));
}
}  // namespace
