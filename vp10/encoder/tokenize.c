/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "vpx_mem/vpx_mem.h"

#include "vp10/common/entropy.h"
#include "vp10/common/pred_common.h"
#include "vp10/common/scan.h"
#include "vp10/common/seg_common.h"

#include "vp10/encoder/cost.h"
#include "vp10/encoder/encoder.h"
#include "vp10/encoder/tokenize.h"

static const TOKENVALUE dct_cat_lt_10_value_tokens[] = {
  {9, 63}, {9, 61}, {9, 59}, {9, 57}, {9, 55}, {9, 53}, {9, 51}, {9, 49},
  {9, 47}, {9, 45}, {9, 43}, {9, 41}, {9, 39}, {9, 37}, {9, 35}, {9, 33},
  {9, 31}, {9, 29}, {9, 27}, {9, 25}, {9, 23}, {9, 21}, {9, 19}, {9, 17},
  {9, 15}, {9, 13}, {9, 11}, {9, 9}, {9, 7}, {9, 5}, {9, 3}, {9, 1},
  {8, 31}, {8, 29}, {8, 27}, {8, 25}, {8, 23}, {8, 21},
  {8, 19}, {8, 17}, {8, 15}, {8, 13}, {8, 11}, {8, 9},
  {8, 7}, {8, 5}, {8, 3}, {8, 1},
  {7, 15}, {7, 13}, {7, 11}, {7, 9}, {7, 7}, {7, 5}, {7, 3}, {7, 1},
  {6, 7}, {6, 5}, {6, 3}, {6, 1}, {5, 3}, {5, 1},
  {4, 1}, {3, 1}, {2, 1}, {1, 1}, {0, 0},
  {1, 0},  {2, 0}, {3, 0}, {4, 0},
  {5, 0}, {5, 2}, {6, 0}, {6, 2}, {6, 4}, {6, 6},
  {7, 0}, {7, 2}, {7, 4}, {7, 6}, {7, 8}, {7, 10}, {7, 12}, {7, 14},
  {8, 0}, {8, 2}, {8, 4}, {8, 6}, {8, 8}, {8, 10}, {8, 12},
  {8, 14}, {8, 16}, {8, 18}, {8, 20}, {8, 22}, {8, 24},
  {8, 26}, {8, 28}, {8, 30}, {9, 0}, {9, 2},
  {9, 4}, {9, 6}, {9, 8}, {9, 10}, {9, 12}, {9, 14}, {9, 16},
  {9, 18}, {9, 20}, {9, 22}, {9, 24}, {9, 26}, {9, 28},
  {9, 30}, {9, 32}, {9, 34}, {9, 36}, {9, 38}, {9, 40},
  {9, 42}, {9, 44}, {9, 46}, {9, 48}, {9, 50}, {9, 52},
  {9, 54}, {9, 56}, {9, 58}, {9, 60}, {9, 62}
};
const TOKENVALUE *vp10_dct_cat_lt_10_value_tokens = dct_cat_lt_10_value_tokens +
    (sizeof(dct_cat_lt_10_value_tokens) / sizeof(*dct_cat_lt_10_value_tokens))
    / 2;
// The corresponding costs of the extrabits for the tokens in the above table
// are stored in the table below. The values are obtained from looking up the
// entry for the specified extrabits in the table corresponding to the token
// (as defined in cost element vp10_extra_bits)
// e.g. {9, 63} maps to cat5_cost[63 >> 1], {1, 1} maps to sign_cost[1 >> 1]
static const int dct_cat_lt_10_value_cost[] = {
  3773, 3750, 3704, 3681, 3623, 3600, 3554, 3531,
  3432, 3409, 3363, 3340, 3282, 3259, 3213, 3190,
  3136, 3113, 3067, 3044, 2986, 2963, 2917, 2894,
  2795, 2772, 2726, 2703, 2645, 2622, 2576, 2553,
  3197, 3116, 3058, 2977, 2881, 2800,
  2742, 2661, 2615, 2534, 2476, 2395,
  2299, 2218, 2160, 2079,
  2566, 2427, 2334, 2195, 2023, 1884, 1791, 1652,
  1893, 1696, 1453, 1256, 1229, 864,
  512, 512, 512, 512, 0,
  512, 512, 512, 512,
  864, 1229, 1256, 1453, 1696, 1893,
  1652, 1791, 1884, 2023, 2195, 2334, 2427, 2566,
  2079, 2160, 2218, 2299, 2395, 2476, 2534, 2615,
  2661, 2742, 2800, 2881, 2977, 3058, 3116, 3197,
  2553, 2576, 2622, 2645, 2703, 2726, 2772, 2795,
  2894, 2917, 2963, 2986, 3044, 3067, 3113, 3136,
  3190, 3213, 3259, 3282, 3340, 3363, 3409, 3432,
  3531, 3554, 3600, 3623, 3681, 3704, 3750, 3773,
};
const int *vp10_dct_cat_lt_10_value_cost = dct_cat_lt_10_value_cost +
    (sizeof(dct_cat_lt_10_value_cost) / sizeof(*dct_cat_lt_10_value_cost))
    / 2;

// Array indices are identical to previously-existing CONTEXT_NODE indices
const vpx_tree_index vp10_coef_tree[TREE_SIZE(ENTROPY_TOKENS)] = {
  -EOB_TOKEN, 2,                       // 0  = EOB
  -ZERO_TOKEN, 4,                      // 1  = ZERO
  -ONE_TOKEN, 6,                       // 2  = ONE
  8, 12,                               // 3  = LOW_VAL
  -TWO_TOKEN, 10,                      // 4  = TWO
  -THREE_TOKEN, -FOUR_TOKEN,           // 5  = THREE
  14, 16,                              // 6  = HIGH_LOW
  -CATEGORY1_TOKEN, -CATEGORY2_TOKEN,  // 7  = CAT_ONE
  18, 20,                              // 8  = CAT_THREEFOUR
  -CATEGORY3_TOKEN, -CATEGORY4_TOKEN,  // 9  = CAT_THREE
  -CATEGORY5_TOKEN, -CATEGORY6_TOKEN   // 10 = CAT_FIVE
};

static const vpx_tree_index cat1[2] = {0, 0};
static const vpx_tree_index cat2[4] = {2, 2, 0, 0};
static const vpx_tree_index cat3[6] = {2, 2, 4, 4, 0, 0};
static const vpx_tree_index cat4[8] = {2, 2, 4, 4, 6, 6, 0, 0};
static const vpx_tree_index cat5[10] = {2, 2, 4, 4, 6, 6, 8, 8, 0, 0};
static const vpx_tree_index cat6[28] = {2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12,
    14, 14, 16, 16, 18, 18, 20, 20, 22, 22, 24, 24, 26, 26, 0, 0};

static const int16_t zero_cost[] = {0};
static const int16_t sign_cost[1] = {512};
static const int16_t cat1_cost[1 << 1] = {864, 1229};
static const int16_t cat2_cost[1 << 2] = {1256, 1453, 1696, 1893};
static const int16_t cat3_cost[1 << 3] = {1652, 1791, 1884, 2023,
                                          2195, 2334, 2427, 2566};
static const int16_t cat4_cost[1 << 4] = {2079, 2160, 2218, 2299, 2395, 2476,
                                          2534, 2615, 2661, 2742, 2800, 2881,
                                          2977, 3058, 3116, 3197};
static const int16_t cat5_cost[1 << 5] = {
    2553, 2576, 2622, 2645, 2703, 2726, 2772, 2795, 2894, 2917, 2963,
    2986, 3044, 3067, 3113, 3136, 3190, 3213, 3259, 3282, 3340, 3363,
    3409, 3432, 3531, 3554, 3600, 3623, 3681, 3704, 3750, 3773};
const int16_t vp10_cat6_low_cost[256] = {
    3378, 3390, 3401, 3413, 3435, 3447, 3458, 3470, 3517, 3529, 3540, 3552,
    3574, 3586, 3597, 3609, 3671, 3683, 3694, 3706, 3728, 3740, 3751, 3763,
    3810, 3822, 3833, 3845, 3867, 3879, 3890, 3902, 3973, 3985, 3996, 4008,
    4030, 4042, 4053, 4065, 4112, 4124, 4135, 4147, 4169, 4181, 4192, 4204,
    4266, 4278, 4289, 4301, 4323, 4335, 4346, 4358, 4405, 4417, 4428, 4440,
    4462, 4474, 4485, 4497, 4253, 4265, 4276, 4288, 4310, 4322, 4333, 4345,
    4392, 4404, 4415, 4427, 4449, 4461, 4472, 4484, 4546, 4558, 4569, 4581,
    4603, 4615, 4626, 4638, 4685, 4697, 4708, 4720, 4742, 4754, 4765, 4777,
    4848, 4860, 4871, 4883, 4905, 4917, 4928, 4940, 4987, 4999, 5010, 5022,
    5044, 5056, 5067, 5079, 5141, 5153, 5164, 5176, 5198, 5210, 5221, 5233,
    5280, 5292, 5303, 5315, 5337, 5349, 5360, 5372, 4988, 5000, 5011, 5023,
    5045, 5057, 5068, 5080, 5127, 5139, 5150, 5162, 5184, 5196, 5207, 5219,
    5281, 5293, 5304, 5316, 5338, 5350, 5361, 5373, 5420, 5432, 5443, 5455,
    5477, 5489, 5500, 5512, 5583, 5595, 5606, 5618, 5640, 5652, 5663, 5675,
    5722, 5734, 5745, 5757, 5779, 5791, 5802, 5814, 5876, 5888, 5899, 5911,
    5933, 5945, 5956, 5968, 6015, 6027, 6038, 6050, 6072, 6084, 6095, 6107,
    5863, 5875, 5886, 5898, 5920, 5932, 5943, 5955, 6002, 6014, 6025, 6037,
    6059, 6071, 6082, 6094, 6156, 6168, 6179, 6191, 6213, 6225, 6236, 6248,
    6295, 6307, 6318, 6330, 6352, 6364, 6375, 6387, 6458, 6470, 6481, 6493,
    6515, 6527, 6538, 6550, 6597, 6609, 6620, 6632, 6654, 6666, 6677, 6689,
    6751, 6763, 6774, 6786, 6808, 6820, 6831, 6843, 6890, 6902, 6913, 6925,
    6947, 6959, 6970, 6982};
const int vp10_cat6_high_cost[64] = {
    88,    2251,  2727,  4890,  3148,  5311,  5787,  7950,  3666,  5829,  6305,
    8468,  6726,  8889,  9365,  11528, 3666,  5829,  6305,  8468,  6726,  8889,
    9365,  11528, 7244,  9407,  9883,  12046, 10304, 12467, 12943, 15106, 3666,
    5829,  6305,  8468,  6726,  8889,  9365,  11528, 7244,  9407,  9883,  12046,
    10304, 12467, 12943, 15106, 7244,  9407,  9883,  12046, 10304, 12467, 12943,
    15106, 10822, 12985, 13461, 15624, 13882, 16045, 16521, 18684};

#if CONFIG_VP9_HIGHBITDEPTH
const int vp10_cat6_high10_high_cost[256] = {
    94,    2257,  2733,  4896,  3154,  5317,  5793,  7956,  3672,  5835,  6311,
    8474,  6732,  8895,  9371,  11534, 3672,  5835,  6311,  8474,  6732,  8895,
    9371,  11534, 7250,  9413,  9889,  12052, 10310, 12473, 12949, 15112, 3672,
    5835,  6311,  8474,  6732,  8895,  9371,  11534, 7250,  9413,  9889,  12052,
    10310, 12473, 12949, 15112, 7250,  9413,  9889,  12052, 10310, 12473, 12949,
    15112, 10828, 12991, 13467, 15630, 13888, 16051, 16527, 18690, 4187,  6350,
    6826,  8989,  7247,  9410,  9886,  12049, 7765,  9928,  10404, 12567, 10825,
    12988, 13464, 15627, 7765,  9928,  10404, 12567, 10825, 12988, 13464, 15627,
    11343, 13506, 13982, 16145, 14403, 16566, 17042, 19205, 7765,  9928,  10404,
    12567, 10825, 12988, 13464, 15627, 11343, 13506, 13982, 16145, 14403, 16566,
    17042, 19205, 11343, 13506, 13982, 16145, 14403, 16566, 17042, 19205, 14921,
    17084, 17560, 19723, 17981, 20144, 20620, 22783, 4187,  6350,  6826,  8989,
    7247,  9410,  9886,  12049, 7765,  9928,  10404, 12567, 10825, 12988, 13464,
    15627, 7765,  9928,  10404, 12567, 10825, 12988, 13464, 15627, 11343, 13506,
    13982, 16145, 14403, 16566, 17042, 19205, 7765,  9928,  10404, 12567, 10825,
    12988, 13464, 15627, 11343, 13506, 13982, 16145, 14403, 16566, 17042, 19205,
    11343, 13506, 13982, 16145, 14403, 16566, 17042, 19205, 14921, 17084, 17560,
    19723, 17981, 20144, 20620, 22783, 8280,  10443, 10919, 13082, 11340, 13503,
    13979, 16142, 11858, 14021, 14497, 16660, 14918, 17081, 17557, 19720, 11858,
    14021, 14497, 16660, 14918, 17081, 17557, 19720, 15436, 17599, 18075, 20238,
    18496, 20659, 21135, 23298, 11858, 14021, 14497, 16660, 14918, 17081, 17557,
    19720, 15436, 17599, 18075, 20238, 18496, 20659, 21135, 23298, 15436, 17599,
    18075, 20238, 18496, 20659, 21135, 23298, 19014, 21177, 21653, 23816, 22074,
    24237, 24713, 26876};
const int vp10_cat6_high12_high_cost[1024] = {
    100,   2263,  2739,  4902,  3160,  5323,  5799,  7962,  3678,  5841,  6317,
    8480,  6738,  8901,  9377,  11540, 3678,  5841,  6317,  8480,  6738,  8901,
    9377,  11540, 7256,  9419,  9895,  12058, 10316, 12479, 12955, 15118, 3678,
    5841,  6317,  8480,  6738,  8901,  9377,  11540, 7256,  9419,  9895,  12058,
    10316, 12479, 12955, 15118, 7256,  9419,  9895,  12058, 10316, 12479, 12955,
    15118, 10834, 12997, 13473, 15636, 13894, 16057, 16533, 18696, 4193,  6356,
    6832,  8995,  7253,  9416,  9892,  12055, 7771,  9934,  10410, 12573, 10831,
    12994, 13470, 15633, 7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633,
    11349, 13512, 13988, 16151, 14409, 16572, 17048, 19211, 7771,  9934,  10410,
    12573, 10831, 12994, 13470, 15633, 11349, 13512, 13988, 16151, 14409, 16572,
    17048, 19211, 11349, 13512, 13988, 16151, 14409, 16572, 17048, 19211, 14927,
    17090, 17566, 19729, 17987, 20150, 20626, 22789, 4193,  6356,  6832,  8995,
    7253,  9416,  9892,  12055, 7771,  9934,  10410, 12573, 10831, 12994, 13470,
    15633, 7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633, 11349, 13512,
    13988, 16151, 14409, 16572, 17048, 19211, 7771,  9934,  10410, 12573, 10831,
    12994, 13470, 15633, 11349, 13512, 13988, 16151, 14409, 16572, 17048, 19211,
    11349, 13512, 13988, 16151, 14409, 16572, 17048, 19211, 14927, 17090, 17566,
    19729, 17987, 20150, 20626, 22789, 8286,  10449, 10925, 13088, 11346, 13509,
    13985, 16148, 11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726, 11864,
    14027, 14503, 16666, 14924, 17087, 17563, 19726, 15442, 17605, 18081, 20244,
    18502, 20665, 21141, 23304, 11864, 14027, 14503, 16666, 14924, 17087, 17563,
    19726, 15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304, 15442, 17605,
    18081, 20244, 18502, 20665, 21141, 23304, 19020, 21183, 21659, 23822, 22080,
    24243, 24719, 26882, 4193,  6356,  6832,  8995,  7253,  9416,  9892,  12055,
    7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633, 7771,  9934,  10410,
    12573, 10831, 12994, 13470, 15633, 11349, 13512, 13988, 16151, 14409, 16572,
    17048, 19211, 7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633, 11349,
    13512, 13988, 16151, 14409, 16572, 17048, 19211, 11349, 13512, 13988, 16151,
    14409, 16572, 17048, 19211, 14927, 17090, 17566, 19729, 17987, 20150, 20626,
    22789, 8286,  10449, 10925, 13088, 11346, 13509, 13985, 16148, 11864, 14027,
    14503, 16666, 14924, 17087, 17563, 19726, 11864, 14027, 14503, 16666, 14924,
    17087, 17563, 19726, 15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304,
    11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726, 15442, 17605, 18081,
    20244, 18502, 20665, 21141, 23304, 15442, 17605, 18081, 20244, 18502, 20665,
    21141, 23304, 19020, 21183, 21659, 23822, 22080, 24243, 24719, 26882, 8286,
    10449, 10925, 13088, 11346, 13509, 13985, 16148, 11864, 14027, 14503, 16666,
    14924, 17087, 17563, 19726, 11864, 14027, 14503, 16666, 14924, 17087, 17563,
    19726, 15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304, 11864, 14027,
    14503, 16666, 14924, 17087, 17563, 19726, 15442, 17605, 18081, 20244, 18502,
    20665, 21141, 23304, 15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304,
    19020, 21183, 21659, 23822, 22080, 24243, 24719, 26882, 12379, 14542, 15018,
    17181, 15439, 17602, 18078, 20241, 15957, 18120, 18596, 20759, 19017, 21180,
    21656, 23819, 15957, 18120, 18596, 20759, 19017, 21180, 21656, 23819, 19535,
    21698, 22174, 24337, 22595, 24758, 25234, 27397, 15957, 18120, 18596, 20759,
    19017, 21180, 21656, 23819, 19535, 21698, 22174, 24337, 22595, 24758, 25234,
    27397, 19535, 21698, 22174, 24337, 22595, 24758, 25234, 27397, 23113, 25276,
    25752, 27915, 26173, 28336, 28812, 30975, 4193,  6356,  6832,  8995,  7253,
    9416,  9892,  12055, 7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633,
    7771,  9934,  10410, 12573, 10831, 12994, 13470, 15633, 11349, 13512, 13988,
    16151, 14409, 16572, 17048, 19211, 7771,  9934,  10410, 12573, 10831, 12994,
    13470, 15633, 11349, 13512, 13988, 16151, 14409, 16572, 17048, 19211, 11349,
    13512, 13988, 16151, 14409, 16572, 17048, 19211, 14927, 17090, 17566, 19729,
    17987, 20150, 20626, 22789, 8286,  10449, 10925, 13088, 11346, 13509, 13985,
    16148, 11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726, 11864, 14027,
    14503, 16666, 14924, 17087, 17563, 19726, 15442, 17605, 18081, 20244, 18502,
    20665, 21141, 23304, 11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726,
    15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304, 15442, 17605, 18081,
    20244, 18502, 20665, 21141, 23304, 19020, 21183, 21659, 23822, 22080, 24243,
    24719, 26882, 8286,  10449, 10925, 13088, 11346, 13509, 13985, 16148, 11864,
    14027, 14503, 16666, 14924, 17087, 17563, 19726, 11864, 14027, 14503, 16666,
    14924, 17087, 17563, 19726, 15442, 17605, 18081, 20244, 18502, 20665, 21141,
    23304, 11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726, 15442, 17605,
    18081, 20244, 18502, 20665, 21141, 23304, 15442, 17605, 18081, 20244, 18502,
    20665, 21141, 23304, 19020, 21183, 21659, 23822, 22080, 24243, 24719, 26882,
    12379, 14542, 15018, 17181, 15439, 17602, 18078, 20241, 15957, 18120, 18596,
    20759, 19017, 21180, 21656, 23819, 15957, 18120, 18596, 20759, 19017, 21180,
    21656, 23819, 19535, 21698, 22174, 24337, 22595, 24758, 25234, 27397, 15957,
    18120, 18596, 20759, 19017, 21180, 21656, 23819, 19535, 21698, 22174, 24337,
    22595, 24758, 25234, 27397, 19535, 21698, 22174, 24337, 22595, 24758, 25234,
    27397, 23113, 25276, 25752, 27915, 26173, 28336, 28812, 30975, 8286,  10449,
    10925, 13088, 11346, 13509, 13985, 16148, 11864, 14027, 14503, 16666, 14924,
    17087, 17563, 19726, 11864, 14027, 14503, 16666, 14924, 17087, 17563, 19726,
    15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304, 11864, 14027, 14503,
    16666, 14924, 17087, 17563, 19726, 15442, 17605, 18081, 20244, 18502, 20665,
    21141, 23304, 15442, 17605, 18081, 20244, 18502, 20665, 21141, 23304, 19020,
    21183, 21659, 23822, 22080, 24243, 24719, 26882, 12379, 14542, 15018, 17181,
    15439, 17602, 18078, 20241, 15957, 18120, 18596, 20759, 19017, 21180, 21656,
    23819, 15957, 18120, 18596, 20759, 19017, 21180, 21656, 23819, 19535, 21698,
    22174, 24337, 22595, 24758, 25234, 27397, 15957, 18120, 18596, 20759, 19017,
    21180, 21656, 23819, 19535, 21698, 22174, 24337, 22595, 24758, 25234, 27397,
    19535, 21698, 22174, 24337, 22595, 24758, 25234, 27397, 23113, 25276, 25752,
    27915, 26173, 28336, 28812, 30975, 12379, 14542, 15018, 17181, 15439, 17602,
    18078, 20241, 15957, 18120, 18596, 20759, 19017, 21180, 21656, 23819, 15957,
    18120, 18596, 20759, 19017, 21180, 21656, 23819, 19535, 21698, 22174, 24337,
    22595, 24758, 25234, 27397, 15957, 18120, 18596, 20759, 19017, 21180, 21656,
    23819, 19535, 21698, 22174, 24337, 22595, 24758, 25234, 27397, 19535, 21698,
    22174, 24337, 22595, 24758, 25234, 27397, 23113, 25276, 25752, 27915, 26173,
    28336, 28812, 30975, 16472, 18635, 19111, 21274, 19532, 21695, 22171, 24334,
    20050, 22213, 22689, 24852, 23110, 25273, 25749, 27912, 20050, 22213, 22689,
    24852, 23110, 25273, 25749, 27912, 23628, 25791, 26267, 28430, 26688, 28851,
    29327, 31490, 20050, 22213, 22689, 24852, 23110, 25273, 25749, 27912, 23628,
    25791, 26267, 28430, 26688, 28851, 29327, 31490, 23628, 25791, 26267, 28430,
    26688, 28851, 29327, 31490, 27206, 29369, 29845, 32008, 30266, 32429, 32905,
    35068};
#endif

#if CONFIG_VP9_HIGHBITDEPTH
static const vpx_tree_index cat1_high10[2] = {0, 0};
static const vpx_tree_index cat2_high10[4] = {2, 2, 0, 0};
static const vpx_tree_index cat3_high10[6] = {2, 2, 4, 4, 0, 0};
static const vpx_tree_index cat4_high10[8] = {2, 2, 4, 4, 6, 6, 0, 0};
static const vpx_tree_index cat5_high10[10] = {2, 2, 4, 4, 6, 6, 8, 8, 0, 0};
static const vpx_tree_index cat6_high10[32] = {2, 2, 4, 4, 6, 6, 8, 8, 10, 10,
  12, 12, 14, 14, 16, 16, 18, 18, 20, 20, 22, 22, 24, 24, 26, 26, 28, 28,
  30, 30, 0, 0};
static const vpx_tree_index cat1_high12[2] = {0, 0};
static const vpx_tree_index cat2_high12[4] = {2, 2, 0, 0};
static const vpx_tree_index cat3_high12[6] = {2, 2, 4, 4, 0, 0};
static const vpx_tree_index cat4_high12[8] = {2, 2, 4, 4, 6, 6, 0, 0};
static const vpx_tree_index cat5_high12[10] = {2, 2, 4, 4, 6, 6, 8, 8, 0, 0};
static const vpx_tree_index cat6_high12[36] = {2, 2, 4, 4, 6, 6, 8, 8, 10, 10,
  12, 12, 14, 14, 16, 16, 18, 18, 20, 20, 22, 22, 24, 24, 26, 26, 28, 28,
  30, 30, 32, 32, 34, 34, 0, 0};
#endif

const vp10_extra_bit vp10_extra_bits[ENTROPY_TOKENS] = {
  {0, 0, 0, 0, zero_cost},                              // ZERO_TOKEN
  {0, 0, 0, 1, sign_cost},                              // ONE_TOKEN
  {0, 0, 0, 2, sign_cost},                              // TWO_TOKEN
  {0, 0, 0, 3, sign_cost},                              // THREE_TOKEN
  {0, 0, 0, 4, sign_cost},                              // FOUR_TOKEN
  {cat1, vp10_cat1_prob, 1,  CAT1_MIN_VAL, cat1_cost},  // CATEGORY1_TOKEN
  {cat2, vp10_cat2_prob, 2,  CAT2_MIN_VAL, cat2_cost},  // CATEGORY2_TOKEN
  {cat3, vp10_cat3_prob, 3,  CAT3_MIN_VAL, cat3_cost},  // CATEGORY3_TOKEN
  {cat4, vp10_cat4_prob, 4,  CAT4_MIN_VAL, cat4_cost},  // CATEGORY4_TOKEN
  {cat5, vp10_cat5_prob, 5,  CAT5_MIN_VAL, cat5_cost},  // CATEGORY5_TOKEN
  {cat6, vp10_cat6_prob, 14, CAT6_MIN_VAL, 0},          // CATEGORY6_TOKEN
  {0, 0, 0, 0, zero_cost}                               // EOB_TOKEN
};

#if CONFIG_VP9_HIGHBITDEPTH
const vp10_extra_bit vp10_extra_bits_high10[ENTROPY_TOKENS] = {
  {0, 0, 0, 0, zero_cost},                                            // ZERO
  {0, 0, 0, 1, sign_cost},                                            // ONE
  {0, 0, 0, 2, sign_cost},                                            // TWO
  {0, 0, 0, 3, sign_cost},                                            // THREE
  {0, 0, 0, 4, sign_cost},                                            // FOUR
  {cat1_high10, vp10_cat1_prob_high10, 1,  CAT1_MIN_VAL, cat1_cost},  // CAT1
  {cat2_high10, vp10_cat2_prob_high10, 2,  CAT2_MIN_VAL, cat2_cost},  // CAT2
  {cat3_high10, vp10_cat3_prob_high10, 3,  CAT3_MIN_VAL, cat3_cost},  // CAT3
  {cat4_high10, vp10_cat4_prob_high10, 4,  CAT4_MIN_VAL, cat4_cost},  // CAT4
  {cat5_high10, vp10_cat5_prob_high10, 5,  CAT5_MIN_VAL, cat5_cost},  // CAT5
  {cat6_high10, vp10_cat6_prob_high10, 16, CAT6_MIN_VAL, 0},          // CAT6
  {0, 0, 0, 0, zero_cost}                                             // EOB
};
const vp10_extra_bit vp10_extra_bits_high12[ENTROPY_TOKENS] = {
  {0, 0, 0, 0, zero_cost},                                            // ZERO
  {0, 0, 0, 1, sign_cost},                                            // ONE
  {0, 0, 0, 2, sign_cost},                                            // TWO
  {0, 0, 0, 3, sign_cost},                                            // THREE
  {0, 0, 0, 4, sign_cost},                                            // FOUR
  {cat1_high12, vp10_cat1_prob_high12, 1,  CAT1_MIN_VAL, cat1_cost},  // CAT1
  {cat2_high12, vp10_cat2_prob_high12, 2,  CAT2_MIN_VAL, cat2_cost},  // CAT2
  {cat3_high12, vp10_cat3_prob_high12, 3,  CAT3_MIN_VAL, cat3_cost},  // CAT3
  {cat4_high12, vp10_cat4_prob_high12, 4,  CAT4_MIN_VAL, cat4_cost},  // CAT4
  {cat5_high12, vp10_cat5_prob_high12, 5,  CAT5_MIN_VAL, cat5_cost},  // CAT5
  {cat6_high12, vp10_cat6_prob_high12, 18, CAT6_MIN_VAL, 0},          // CAT6
  {0, 0, 0, 0, zero_cost}                                             // EOB
};
#endif

#if !CONFIG_ANS
const struct vp10_token vp10_coef_encodings[ENTROPY_TOKENS] = {
  {2, 2}, {6, 3}, {28, 5}, {58, 6}, {59, 6}, {60, 6}, {61, 6}, {124, 7},
  {125, 7}, {126, 7}, {127, 7}, {0, 1}
};
#endif  // !CONFIG_ANS

struct tokenize_b_args {
  VP10_COMP *cpi;
  ThreadData *td;
  TOKENEXTRA **tp;
};

static void set_entropy_context_b(int plane, int block,
                                  int blk_row, int blk_col,
                                  BLOCK_SIZE plane_bsize,
                                  TX_SIZE tx_size, void *arg) {
  struct tokenize_b_args* const args = arg;
  ThreadData *const td = args->td;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *p = &x->plane[plane];
  struct macroblockd_plane *pd = &xd->plane[plane];
  vp10_set_contexts(xd, pd, plane_bsize, tx_size, p->eobs[block] > 0,
                    blk_col, blk_row);
}

static INLINE void add_token(TOKENEXTRA **t, const vpx_prob *context_tree,
#if CONFIG_ANS
                             const rans_dec_lut *token_cdf,
#endif  // CONFIG_ANS
                             int32_t extra, uint8_t token,
                             uint8_t skip_eob_node, unsigned int *counts) {
  (*t)->token = token;
  (*t)->extra = extra;
  (*t)->context_tree = context_tree;
#if CONFIG_ANS
  (*t)->token_cdf = token_cdf;
#endif  // CONFIG_ANS
  (*t)->skip_eob_node = skip_eob_node;
  (*t)++;
  ++counts[token];
}

static INLINE void add_token_no_extra(TOKENEXTRA **t,
                                      const vpx_prob *context_tree,
                                      uint8_t token,
                                      uint8_t skip_eob_node,
                                      unsigned int *counts) {
  (*t)->token = token;
  (*t)->context_tree = context_tree;
  (*t)->skip_eob_node = skip_eob_node;
  (*t)++;
  ++counts[token];
}

static INLINE int get_tx_eob(const struct segmentation *seg, int segment_id,
                             TX_SIZE tx_size) {
  const int eob_max = num_4x4_blocks_txsize_lookup[tx_size] << 4;
  return segfeature_active(seg, segment_id, SEG_LVL_SKIP) ? 0 : eob_max;
}

void vp10_tokenize_palette_sb(struct ThreadData *const td,
                              BLOCK_SIZE bsize, int plane,
                              TOKENEXTRA **t) {
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  uint8_t *color_map = xd->plane[plane != 0].color_index_map;
  PALETTE_MODE_INFO *pmi = &mbmi->palette_mode_info;
  int n = pmi->palette_size[plane != 0];
  int i, j, k;
  int color_new_idx = -1, color_ctx, color_order[PALETTE_MAX_SIZE];
  const int rows = (4 * num_4x4_blocks_high_lookup[bsize]) >>
      (xd->plane[plane != 0].subsampling_y);
  const int cols = (4 * num_4x4_blocks_wide_lookup[bsize]) >>
      (xd->plane[plane != 0].subsampling_x);
  const vpx_prob (* const probs)[PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1] =
      plane == 0 ? vp10_default_palette_y_color_prob :
          vp10_default_palette_uv_color_prob;

  for (i = 0; i < rows; ++i) {
    for (j = (i == 0 ? 1 : 0); j < cols; ++j) {
      color_ctx = vp10_get_palette_color_context(color_map, cols, i, j, n,
                                                 color_order);
      for (k = 0; k < n; ++k)
        if (color_map[i * cols + j] == color_order[k]) {
          color_new_idx = k;
          break;
        }
      assert(color_new_idx >= 0 && color_new_idx < n);
      (*t)->token = color_new_idx;
      (*t)->context_tree = probs[n - 2][color_ctx];
      (*t)->skip_eob_node = 0;
      ++(*t);
    }
  }
}

static void tokenize_b(int plane, int block, int blk_row, int blk_col,
                       BLOCK_SIZE plane_bsize,
                       TX_SIZE tx_size, void *arg) {
  struct tokenize_b_args* const args = arg;
  VP10_COMP *cpi = args->cpi;
  ThreadData *const td = args->td;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  TOKENEXTRA **tp = args->tp;
  uint8_t token_cache[MAX_TX_SQUARE];
  struct macroblock_plane *p = &x->plane[plane];
  struct macroblockd_plane *pd = &xd->plane[plane];
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  int pt; /* near block/prev token context index */
  int c;
  TOKENEXTRA *t = *tp;        /* store tokens starting here */
  int eob = p->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const tran_low_t *qcoeff = BLOCK_OFFSET(p->qcoeff, block);
#if CONFIG_SUPERTX
  const int segment_id = VPXMIN(mbmi->segment_id, mbmi->segment_id_supertx);
#else
  const int segment_id = mbmi->segment_id;
#endif  // CONFIG_SUEPRTX
  const int16_t *scan, *nb;
  const TX_TYPE tx_type = get_tx_type(type, xd, block, tx_size);
  const scan_order *const so = get_scan(tx_size, tx_type, is_inter_block(mbmi));
  const int ref = is_inter_block(mbmi);
  unsigned int (*const counts)[COEFF_CONTEXTS][ENTROPY_TOKENS] =
      td->rd_counts.coef_counts[txsize_sqr_map[tx_size]][type][ref];
#if CONFIG_ENTROPY
  vpx_prob (*coef_probs)[COEFF_CONTEXTS][UNCONSTRAINED_NODES] =
      cpi->subframe_stats.coef_probs_buf[cpi->common.coef_probs_update_idx]
                                        [txsize_sqr_map[tx_size]][type][ref];
#else
  vpx_prob (*const coef_probs)[COEFF_CONTEXTS][UNCONSTRAINED_NODES] =
      cpi->common.fc->coef_probs[txsize_sqr_map[tx_size]][type][ref];
#endif  // CONFIG_ENTROPY
#if CONFIG_ANS
  rans_dec_lut (*const coef_cdfs)[COEFF_CONTEXTS] =
      cpi->common.fc->coef_cdfs[txsize_sqr_map[tx_size]][type][ref];
#endif  // CONFIG_ANS
  unsigned int (*const eob_branch)[COEFF_CONTEXTS] =
      td->counts->eob_branch[txsize_sqr_map[tx_size]][type][ref];
  const uint8_t *const band = get_band_translate(tx_size);
  const int seg_eob = get_tx_eob(&cpi->common.seg, segment_id, tx_size);
  int skip_eob = 0;
  int16_t token;
  EXTRABIT extra;
  pt = get_entropy_context(tx_size, pd->above_context + blk_col,
                           pd->left_context + blk_row);
  scan = so->scan;
  nb = so->neighbors;
  c = 0;

  while (c < eob) {
    const int v = qcoeff[scan[c]];
    eob_branch[band[c]][pt] += !skip_eob;

    vp10_get_token_extra(v, &token, &extra);

    add_token(&t, coef_probs[band[c]][pt],
#if CONFIG_ANS
              (const rans_dec_lut*)&coef_cdfs[band[c]][pt],
#endif  // CONFIG_ANS
              extra, (uint8_t)token, (uint8_t)skip_eob, counts[band[c]][pt]);

    token_cache[scan[c]] = vp10_pt_energy_class[token];
    ++c;
    pt = get_coef_context(nb, token_cache, c);
    skip_eob = (token == ZERO_TOKEN);
  }
  if (c < seg_eob) {
    add_token_no_extra(&t, coef_probs[band[c]][pt], EOB_TOKEN, 0,
                       counts[band[c]][pt]);
    ++eob_branch[band[c]][pt];
  }

  *tp = t;

  vp10_set_contexts(xd, pd, plane_bsize, tx_size, c > 0, blk_col, blk_row);
}

struct is_skippable_args {
  uint16_t *eobs;
  int *skippable;
};
static void is_skippable(int plane, int block, int blk_row, int blk_col,
                         BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                         void *argv) {
  struct is_skippable_args *args = argv;
  (void)plane;
  (void)plane_bsize;
  (void)tx_size;
  (void)blk_row;
  (void)blk_col;
  args->skippable[0] &= (!args->eobs[block]);
}

// TODO(yaowu): rewrite and optimize this function to remove the usage of
//              vp10_foreach_transform_block() and simplify is_skippable().
int vp10_is_skippable_in_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  int result = 1;
  struct is_skippable_args args = {x->plane[plane].eobs, &result};
  vp10_foreach_transformed_block_in_plane(&x->e_mbd, bsize, plane, is_skippable,
                                          &args);
  return result;
}

static void has_high_freq_coeff(int plane, int block, int blk_row, int blk_col,
                                BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                                void *argv) {
  struct is_skippable_args *args = argv;
  int eobs = (tx_size == TX_4X4) ? 3 : 10;
  (void) plane;
  (void) plane_bsize;
  (void) blk_row;
  (void) blk_col;

  *(args->skippable) |= (args->eobs[block] > eobs);
}

int vp10_has_high_freq_in_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  int result = 0;
  struct is_skippable_args args = {x->plane[plane].eobs, &result};
  vp10_foreach_transformed_block_in_plane(&x->e_mbd, bsize, plane,
                                          has_high_freq_coeff, &args);
  return result;
}

#if CONFIG_VAR_TX
void tokenize_tx(ThreadData *td, TOKENEXTRA **t,
                 int dry_run, TX_SIZE tx_size, BLOCK_SIZE plane_bsize,
                 int blk_row, int blk_col, int block, int plane,
                 void *arg) {
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  const BLOCK_SIZE bsize = txsize_to_bsize[tx_size];
  const int tx_row = blk_row >> (1 - pd->subsampling_y);
  const int tx_col = blk_col >> (1 - pd->subsampling_x);
  const TX_SIZE plane_tx_size = plane ?
      get_uv_tx_size_impl(mbmi->inter_tx_size[tx_row][tx_col], bsize, 0, 0) :
      mbmi->inter_tx_size[tx_row][tx_col];

  int max_blocks_high = num_4x4_blocks_high_lookup[plane_bsize];
  int max_blocks_wide = num_4x4_blocks_wide_lookup[plane_bsize];

  assert(tx_size < TX_SIZES);

  if (xd->mb_to_bottom_edge < 0)
    max_blocks_high += xd->mb_to_bottom_edge >> (5 + pd->subsampling_y);
  if (xd->mb_to_right_edge < 0)
    max_blocks_wide += xd->mb_to_right_edge >> (5 + pd->subsampling_x);

  if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
    return;

  if (tx_size == plane_tx_size) {
    const struct macroblockd_plane *const pd = &xd->plane[plane];
    BLOCK_SIZE plane_bsize = get_plane_block_size(mbmi->sb_type, pd);
    if (!dry_run)
      tokenize_b(plane, block, blk_row, blk_col, plane_bsize, tx_size, arg);
    else
      set_entropy_context_b(plane, block, blk_row, blk_col,
                            plane_bsize, tx_size, arg);
  } else {
    int bsl = b_width_log2_lookup[bsize];
    int i;

    assert(bsl > 0);
    --bsl;

    for (i = 0; i < 4; ++i) {
      const int offsetr = blk_row + ((i >> 1) << bsl);
      const int offsetc = blk_col + ((i & 0x01) << bsl);
      int step = num_4x4_blocks_txsize_lookup[tx_size - 1];

      if (offsetr >= max_blocks_high || offsetc >= max_blocks_wide)
        continue;

      tokenize_tx(td, t, dry_run, tx_size - 1, plane_bsize,
                  offsetr, offsetc, block + i * step, plane, arg);
    }
  }
}

void vp10_tokenize_sb_inter(VP10_COMP *cpi, ThreadData *td, TOKENEXTRA **t,
                            int dry_run, int mi_row, int mi_col,
                            BLOCK_SIZE bsize) {
  VP10_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  TOKENEXTRA *t_backup = *t;
  const int ctx = vp10_get_skip_context(xd);
  const int skip_inc = !segfeature_active(&cm->seg, mbmi->segment_id,
                                          SEG_LVL_SKIP);
  struct tokenize_b_args arg = {cpi, td, t};
  int plane;
  if (mi_row >= cm->mi_rows || mi_col >= cm->mi_cols)
    return;

  if (mbmi->skip) {
    if (!dry_run)
      td->counts->skip[ctx][1] += skip_inc;
    reset_skip_context(xd, bsize);
    if (dry_run)
      *t = t_backup;
    return;
  }

  if (!dry_run)
    td->counts->skip[ctx][0] += skip_inc;
  else
    *t = t_backup;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    const struct macroblockd_plane *const pd = &xd->plane[plane];
    const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
    const int mi_width = num_4x4_blocks_wide_lookup[plane_bsize];
    const int mi_height = num_4x4_blocks_high_lookup[plane_bsize];
    const TX_SIZE max_tx_size = max_txsize_lookup[plane_bsize];
    const BLOCK_SIZE txb_size = txsize_to_bsize[max_tx_size];
    int bh = num_4x4_blocks_wide_lookup[txb_size];
    int idx, idy;
    int block = 0;
    int step = num_4x4_blocks_txsize_lookup[max_tx_size];
    for (idy = 0; idy < mi_height; idy += bh) {
      for (idx = 0; idx < mi_width; idx += bh) {
        tokenize_tx(td, t, dry_run, max_tx_size, plane_bsize, idy, idx,
                    block, plane, &arg);
        block += step;
      }
    }

    if (!dry_run) {
      (*t)->token = EOSB_TOKEN;
      (*t)++;
    }
  }
}
#endif  // CONFIG_VAR_TX

void vp10_tokenize_sb(VP10_COMP *cpi, ThreadData *td, TOKENEXTRA **t,
                     int dry_run, BLOCK_SIZE bsize) {
  VP10_COMMON *const cm = &cpi->common;
  MACROBLOCK *const x = &td->mb;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  const int ctx = vp10_get_skip_context(xd);
  const int skip_inc = !segfeature_active(&cm->seg, mbmi->segment_id,
                                          SEG_LVL_SKIP);
  struct tokenize_b_args arg = {cpi, td, t};
  if (mbmi->skip) {
    if (!dry_run)
      td->counts->skip[ctx][1] += skip_inc;
    reset_skip_context(xd, bsize);
    return;
  }

  if (!dry_run) {
    int plane;

    td->counts->skip[ctx][0] += skip_inc;

    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      vp10_foreach_transformed_block_in_plane(xd, bsize, plane, tokenize_b,
                                              &arg);
      (*t)->token = EOSB_TOKEN;
      (*t)++;
    }
  } else {
    vp10_foreach_transformed_block(xd, bsize, set_entropy_context_b, &arg);
  }
}

#if CONFIG_SUPERTX
void vp10_tokenize_sb_supertx(VP10_COMP *cpi, ThreadData *td, TOKENEXTRA **t,
                              int dry_run, BLOCK_SIZE bsize) {
  VP10_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &td->mb.e_mbd;
  MB_MODE_INFO *const mbmi = &xd->mi[0]->mbmi;
  TOKENEXTRA *t_backup = *t;
  const int ctx = vp10_get_skip_context(xd);
  const int skip_inc = !segfeature_active(&cm->seg, mbmi->segment_id_supertx,
                                          SEG_LVL_SKIP);
  struct tokenize_b_args arg = {cpi, td, t};
  if (mbmi->skip) {
    if (!dry_run)
      td->counts->skip[ctx][1] += skip_inc;
    reset_skip_context(xd, bsize);
    if (dry_run)
      *t = t_backup;
    return;
  }

  if (!dry_run) {
    int plane;
    td->counts->skip[ctx][0] += skip_inc;

    for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
      vp10_foreach_transformed_block_in_plane(xd, bsize, plane, tokenize_b,
                                              &arg);
      (*t)->token = EOSB_TOKEN;
      (*t)++;
    }
  } else {
    vp10_foreach_transformed_block(xd, bsize, set_entropy_context_b, &arg);
    *t = t_backup;
  }
}
#endif  // CONFIG_SUPERTX
