/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#define AV1_FORCE_AOMBOOL_TREEWRITER

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <utility>
#include <vector>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "test/acm_random.h"
#include "av1/common/ans.h"
#include "av1/encoder/treewriter.h"
#include "aom_dsp/bitreader.h"
#include "aom_dsp/bitwriter.h"

namespace {
typedef std::vector<std::pair<uint8_t, bool> > PvVec;

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
    double d = p / 256.;
    entropy += -d * log2(d) - (1 - d) * log2(1 - d);
  }
  printf("entropy %f\n", entropy);
  return ret;
}

bool check_rabs(const PvVec &pv_vec, uint8_t *buf) {
  AnsCoder a;
  ans_write_init(&a, buf);

  std::clock_t start = std::clock();
  for (PvVec::const_reverse_iterator it = pv_vec.rbegin(); it != pv_vec.rend();
       ++it) {
    rabs_write(&a, it->second, 256 - it->first);
  }
  std::clock_t enc_time = std::clock() - start;
  int offset = ans_write_end(&a);
  bool okay = true;
  AnsDecoder d;
  if (ans_read_init(&d, buf, offset)) return false;
  start = std::clock();
  for (PvVec::const_iterator it = pv_vec.begin(); it != pv_vec.end(); ++it) {
    okay &= rabs_read(&d, 256 - it->first) == it->second;
  }
  std::clock_t dec_time = std::clock() - start;
  if (!okay) return false;
  printf("rABS size %d enc_time %f dec_time %f\n", offset,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
}

bool check_rabs_asc(const PvVec &pv_vec, uint8_t *buf) {
  AnsCoder a;
  ans_write_init(&a, buf);

  std::clock_t start = std::clock();
  for (PvVec::const_reverse_iterator it = pv_vec.rbegin(); it != pv_vec.rend();
       ++it) {
    rabs_asc_write(&a, it->second, 256 - it->first);
  }
  std::clock_t enc_time = std::clock() - start;
  int offset = ans_write_end(&a);
  bool okay = true;
  AnsDecoder d;
  if (ans_read_init(&d, buf, offset)) return false;
  start = std::clock();
  for (PvVec::const_iterator it = pv_vec.begin(); it != pv_vec.end(); ++it) {
    okay &= rabs_asc_read(&d, 256 - it->first) == it->second;
  }
  std::clock_t dec_time = std::clock() - start;
  if (!okay) return false;
  printf("rABS (asc) size %d enc_time %f dec_time %f\n", offset,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
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
  printf("uABS size %d enc_time %f dec_time %f\n", offset,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
}

bool check_aombool(const PvVec &pv_vec, uint8_t *buf) {
  aom_writer w;
  aom_reader r;
  aom_start_encode(&w, buf);

  std::clock_t start = std::clock();
  for (PvVec::const_iterator it = pv_vec.begin(); it != pv_vec.end(); ++it) {
    aom_write(&w, it->second, 256 - it->first);
  }
  std::clock_t enc_time = std::clock() - start;
  aom_stop_encode(&w);
  bool okay = true;
  aom_reader_init(&r, buf, w.pos, NULL, NULL);
  start = std::clock();
  for (PvVec::const_iterator it = pv_vec.begin(); it != pv_vec.end(); ++it) {
    okay &= aom_read(&r, 256 - it->first) == it->second;
  }
  std::clock_t dec_time = std::clock() - start;
  printf("AOM size %d enc_time %f dec_time %f\n", w.pos,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return okay;
}

// TODO(aconverse): replace this with a more representative distribution from
// the codec.
const rans_sym rans_sym_tab[] = {
  { 16 * 4, 0 * 4 },
  { 100 * 4, 16 * 4 },
  { 70 * 4, 116 * 4 },
  { 70 * 4, 186 * 4 },
};
const int kDistinctSyms = sizeof(rans_sym_tab) / sizeof(rans_sym_tab[0]);

std::vector<int> ans_encode_build_vals(const rans_sym *tab, int iters) {
  std::vector<int> p_to_sym;
  int i = 0;
  while (p_to_sym.size() < rans_precision) {
    p_to_sym.insert(p_to_sym.end(), tab[i].prob, i);
    ++i;
  }
  assert(p_to_sym.size() == rans_precision);
  std::vector<int> ret;
  libaom_test::ACMRandom gen(18543637);
  for (int i = 0; i < iters; ++i) {
    int sym = p_to_sym[gen.Rand8() * 4];
    ret.push_back(sym);
  }
  return ret;
}

void rans_build_dec_tab(const struct rans_sym sym_tab[], rans_dec_lut dec_tab) {
  dec_tab[0] = 0;
  for (int i = 1; dec_tab[i - 1] < rans_precision; ++i) {
    dec_tab[i] = dec_tab[i - 1] + sym_tab[i - 1].prob;
  }
}

bool check_rans(const std::vector<int> &sym_vec, const rans_sym *const tab,
                uint8_t *buf) {
  AnsCoder a;
  ans_write_init(&a, buf);
  rans_dec_lut dec_tab;
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
  printf("rANS size %d enc_time %f dec_time %f\n", offset,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return ans_read_end(&d);
}

void build_tree(aom_tree_index *tree, int num_syms) {
  aom_tree_index i;
  int sym = 0;
  for (i = 0; i < num_syms - 1; ++i) {
    tree[2 * i] = sym--;
    tree[2 * i + 1] = 2 * (i + 1);
  }
  tree[2 * i - 1] = sym;
}

/* The treep array contains the probabilities of nodes of a tree structured
 * like:
 *          *
 *         / \
 *    -sym0   *
 *           / \
 *       -sym1  *
 *             / \
 *        -sym2  -sym3
 */
void tab2tree(const rans_sym *tab, int tab_size, aom_prob *treep) {
  const unsigned basep = rans_precision;
  unsigned pleft = basep;
  for (int i = 0; i < tab_size - 1; ++i) {
    unsigned prob = (tab[i].prob * basep + basep * 2) / (pleft * 4);
    assert(prob > 0 && prob < 256);
    treep[i] = prob;
    pleft -= tab[i].prob;
  }
}

struct sym_bools {
  unsigned bits;
  int len;
};

static void make_tree_bits_tab(sym_bools *tab, int num_syms) {
  unsigned bits = 0;
  int len = 0;
  int i;
  for (i = 0; i < num_syms - 1; ++i) {
    bits *= 2;
    ++len;
    tab[i].bits = bits;
    tab[i].len = len;
    ++bits;
  }
  tab[i].bits = bits;
  tab[i].len = len;
}

void build_tpb(aom_prob probs[/*num_syms*/],
               aom_tree_index tree[/*2*num_syms*/],
               sym_bools bit_len[/*num_syms*/],
               const rans_sym sym_tab[/*num_syms*/], int num_syms) {
  tab2tree(sym_tab, num_syms, probs);
  build_tree(tree, num_syms);
  make_tree_bits_tab(bit_len, num_syms);
}

bool check_aomtree(const std::vector<int> &sym_vec, const rans_sym *sym_tab,
                   uint8_t *buf) {
  aom_writer w;
  aom_reader r;
  aom_start_encode(&w, buf);

  aom_prob probs[kDistinctSyms];
  aom_tree_index tree[2 * kDistinctSyms];
  sym_bools bit_len[kDistinctSyms];
  build_tpb(probs, tree, bit_len, sym_tab, kDistinctSyms);

  std::clock_t start = std::clock();
  for (std::vector<int>::const_iterator it = sym_vec.begin();
       it != sym_vec.end(); ++it) {
    av1_write_tree(&w, tree, probs, bit_len[*it].bits, bit_len[*it].len, 0);
  }
  std::clock_t enc_time = std::clock() - start;
  aom_stop_encode(&w);
  aom_reader_init(&r, buf, w.pos, NULL, NULL);
  start = std::clock();
  for (std::vector<int>::const_iterator it = sym_vec.begin();
       it != sym_vec.end(); ++it) {
    if (aom_read_tree(&r, tree, probs) != *it) return false;
  }
  std::clock_t dec_time = std::clock() - start;
  printf("AOMtree size %u enc_time %f dec_time %f\n", w.pos,
         static_cast<float>(enc_time) / CLOCKS_PER_SEC,
         static_cast<float>(dec_time) / CLOCKS_PER_SEC);
  return true;
}

class Av1AbsTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() { pv_vec_ = abs_encode_build_vals(kNumBools); }
  virtual void SetUp() { buf_ = new uint8_t[kNumBools / 8]; }
  virtual void TearDown() { delete[] buf_; }
  static const int kNumBools = 100000000;
  static PvVec pv_vec_;
  uint8_t *buf_;
};
PvVec Av1AbsTest::pv_vec_;

class Av1AnsTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    sym_vec_ = ans_encode_build_vals(rans_sym_tab, kNumSyms);
  }
  virtual void SetUp() { buf_ = new uint8_t[kNumSyms / 2]; }
  virtual void TearDown() { delete[] buf_; }
  static const int kNumSyms = 25000000;
  static std::vector<int> sym_vec_;
  uint8_t *buf_;
};
std::vector<int> Av1AnsTest::sym_vec_;

TEST_F(Av1AbsTest, Avxbool) { EXPECT_TRUE(check_aombool(pv_vec_, buf_)); }
TEST_F(Av1AbsTest, Rabs) { EXPECT_TRUE(check_rabs(pv_vec_, buf_)); }
TEST_F(Av1AbsTest, RabsAsc) { EXPECT_TRUE(check_rabs_asc(pv_vec_, buf_)); }
TEST_F(Av1AbsTest, Uabs) { EXPECT_TRUE(check_uabs(pv_vec_, buf_)); }

TEST_F(Av1AnsTest, Rans) {
  EXPECT_TRUE(check_rans(sym_vec_, rans_sym_tab, buf_));
}
TEST_F(Av1AnsTest, Avxtree) {
  EXPECT_TRUE(check_aomtree(sym_vec_, rans_sym_tab, buf_));
}
}  // namespace
