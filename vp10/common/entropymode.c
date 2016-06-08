/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_mem/vpx_mem.h"

#include "vp10/common/reconinter.h"
#include "vp10/common/onyxc_int.h"
#include "vp10/common/seg_common.h"

const vpx_prob vp10_kf_y_mode_prob[INTRA_MODES][INTRA_MODES][INTRA_MODES - 1] = {
  {  // above = dc
    { 137,  30,  42, 148, 151, 207,  70,  52,  91 },  // left = dc
    {  92,  45, 102, 136, 116, 180,  74,  90, 100 },  // left = v
    {  73,  32,  19, 187, 222, 215,  46,  34, 100 },  // left = h
    {  91,  30,  32, 116, 121, 186,  93,  86,  94 },  // left = d45
    {  72,  35,  36, 149,  68, 206,  68,  63, 105 },  // left = d135
    {  73,  31,  28, 138,  57, 124,  55, 122, 151 },  // left = d117
    {  67,  23,  21, 140, 126, 197,  40,  37, 171 },  // left = d153
    {  86,  27,  28, 128, 154, 212,  45,  43,  53 },  // left = d207
    {  74,  32,  27, 107,  86, 160,  63, 134, 102 },  // left = d63
    {  59,  67,  44, 140, 161, 202,  78,  67, 119 }   // left = tm
  }, {  // above = v
    {  63,  36, 126, 146, 123, 158,  60,  90,  96 },  // left = dc
    {  43,  46, 168, 134, 107, 128,  69, 142,  92 },  // left = v
    {  44,  29,  68, 159, 201, 177,  50,  57,  77 },  // left = h
    {  58,  38,  76, 114,  97, 172,  78, 133,  92 },  // left = d45
    {  46,  41,  76, 140,  63, 184,  69, 112,  57 },  // left = d135
    {  38,  32,  85, 140,  46, 112,  54, 151, 133 },  // left = d117
    {  39,  27,  61, 131, 110, 175,  44,  75, 136 },  // left = d153
    {  52,  30,  74, 113, 130, 175,  51,  64,  58 },  // left = d207
    {  47,  35,  80, 100,  74, 143,  64, 163,  74 },  // left = d63
    {  36,  61, 116, 114, 128, 162,  80, 125,  82 }   // left = tm
  }, {  // above = h
    {  82,  26,  26, 171, 208, 204,  44,  32, 105 },  // left = dc
    {  55,  44,  68, 166, 179, 192,  57,  57, 108 },  // left = v
    {  42,  26,  11, 199, 241, 228,  23,  15,  85 },  // left = h
    {  68,  42,  19, 131, 160, 199,  55,  52,  83 },  // left = d45
    {  58,  50,  25, 139, 115, 232,  39,  52, 118 },  // left = d135
    {  50,  35,  33, 153, 104, 162,  64,  59, 131 },  // left = d117
    {  44,  24,  16, 150, 177, 202,  33,  19, 156 },  // left = d153
    {  55,  27,  12, 153, 203, 218,  26,  27,  49 },  // left = d207
    {  53,  49,  21, 110, 116, 168,  59,  80,  76 },  // left = d63
    {  38,  72,  19, 168, 203, 212,  50,  50, 107 }   // left = tm
  }, {  // above = d45
    { 103,  26,  36, 129, 132, 201,  83,  80,  93 },  // left = dc
    {  59,  38,  83, 112, 103, 162,  98, 136,  90 },  // left = v
    {  62,  30,  23, 158, 200, 207,  59,  57,  50 },  // left = h
    {  67,  30,  29,  84,  86, 191, 102,  91,  59 },  // left = d45
    {  60,  32,  33, 112,  71, 220,  64,  89, 104 },  // left = d135
    {  53,  26,  34, 130,  56, 149,  84, 120, 103 },  // left = d117
    {  53,  21,  23, 133, 109, 210,  56,  77, 172 },  // left = d153
    {  77,  19,  29, 112, 142, 228,  55,  66,  36 },  // left = d207
    {  61,  29,  29,  93,  97, 165,  83, 175, 162 },  // left = d63
    {  47,  47,  43, 114, 137, 181, 100,  99,  95 }   // left = tm
  }, {  // above = d135
    {  69,  23,  29, 128,  83, 199,  46,  44, 101 },  // left = dc
    {  53,  40,  55, 139,  69, 183,  61,  80, 110 },  // left = v
    {  40,  29,  19, 161, 180, 207,  43,  24,  91 },  // left = h
    {  60,  34,  19, 105,  61, 198,  53,  64,  89 },  // left = d45
    {  52,  31,  22, 158,  40, 209,  58,  62,  89 },  // left = d135
    {  44,  31,  29, 147,  46, 158,  56, 102, 198 },  // left = d117
    {  35,  19,  12, 135,  87, 209,  41,  45, 167 },  // left = d153
    {  55,  25,  21, 118,  95, 215,  38,  39,  66 },  // left = d207
    {  51,  38,  25, 113,  58, 164,  70,  93,  97 },  // left = d63
    {  47,  54,  34, 146, 108, 203,  72, 103, 151 }   // left = tm
  }, {  // above = d117
    {  64,  19,  37, 156,  66, 138,  49,  95, 133 },  // left = dc
    {  46,  27,  80, 150,  55, 124,  55, 121, 135 },  // left = v
    {  36,  23,  27, 165, 149, 166,  54,  64, 118 },  // left = h
    {  53,  21,  36, 131,  63, 163,  60, 109,  81 },  // left = d45
    {  40,  26,  35, 154,  40, 185,  51,  97, 123 },  // left = d135
    {  35,  19,  34, 179,  19,  97,  48, 129, 124 },  // left = d117
    {  36,  20,  26, 136,  62, 164,  33,  77, 154 },  // left = d153
    {  45,  18,  32, 130,  90, 157,  40,  79,  91 },  // left = d207
    {  45,  26,  28, 129,  45, 129,  49, 147, 123 },  // left = d63
    {  38,  44,  51, 136,  74, 162,  57,  97, 121 }   // left = tm
  }, {  // above = d153
    {  75,  17,  22, 136, 138, 185,  32,  34, 166 },  // left = dc
    {  56,  39,  58, 133, 117, 173,  48,  53, 187 },  // left = v
    {  35,  21,  12, 161, 212, 207,  20,  23, 145 },  // left = h
    {  56,  29,  19, 117, 109, 181,  55,  68, 112 },  // left = d45
    {  47,  29,  17, 153,  64, 220,  59,  51, 114 },  // left = d135
    {  46,  16,  24, 136,  76, 147,  41,  64, 172 },  // left = d117
    {  34,  17,  11, 108, 152, 187,  13,  15, 209 },  // left = d153
    {  51,  24,  14, 115, 133, 209,  32,  26, 104 },  // left = d207
    {  55,  30,  18, 122,  79, 179,  44,  88, 116 },  // left = d63
    {  37,  49,  25, 129, 168, 164,  41,  54, 148 }   // left = tm
  }, {  // above = d207
    {  82,  22,  32, 127, 143, 213,  39,  41,  70 },  // left = dc
    {  62,  44,  61, 123, 105, 189,  48,  57,  64 },  // left = v
    {  47,  25,  17, 175, 222, 220,  24,  30,  86 },  // left = h
    {  68,  36,  17, 106, 102, 206,  59,  74,  74 },  // left = d45
    {  57,  39,  23, 151,  68, 216,  55,  63,  58 },  // left = d135
    {  49,  30,  35, 141,  70, 168,  82,  40, 115 },  // left = d117
    {  51,  25,  15, 136, 129, 202,  38,  35, 139 },  // left = d153
    {  68,  26,  16, 111, 141, 215,  29,  28,  28 },  // left = d207
    {  59,  39,  19, 114,  75, 180,  77, 104,  42 },  // left = d63
    {  40,  61,  26, 126, 152, 206,  61,  59,  93 }   // left = tm
  }, {  // above = d63
    {  78,  23,  39, 111, 117, 170,  74, 124,  94 },  // left = dc
    {  48,  34,  86, 101,  92, 146,  78, 179, 134 },  // left = v
    {  47,  22,  24, 138, 187, 178,  68,  69,  59 },  // left = h
    {  56,  25,  33, 105, 112, 187,  95, 177, 129 },  // left = d45
    {  48,  31,  27, 114,  63, 183,  82, 116,  56 },  // left = d135
    {  43,  28,  37, 121,  63, 123,  61, 192, 169 },  // left = d117
    {  42,  17,  24, 109,  97, 177,  56,  76, 122 },  // left = d153
    {  58,  18,  28, 105, 139, 182,  70,  92,  63 },  // left = d207
    {  46,  23,  32,  74,  86, 150,  67, 183,  88 },  // left = d63
    {  36,  38,  48,  92, 122, 165,  88, 137,  91 }   // left = tm
  }, {  // above = tm
    {  65,  70,  60, 155, 159, 199,  61,  60,  81 },  // left = dc
    {  44,  78, 115, 132, 119, 173,  71, 112,  93 },  // left = v
    {  39,  38,  21, 184, 227, 206,  42,  32,  64 },  // left = h
    {  58,  47,  36, 124, 137, 193,  80,  82,  78 },  // left = d45
    {  49,  50,  35, 144,  95, 205,  63,  78,  59 },  // left = d135
    {  41,  53,  52, 148,  71, 142,  65, 128,  51 },  // left = d117
    {  40,  36,  28, 143, 143, 202,  40,  55, 137 },  // left = d153
    {  52,  34,  29, 129, 183, 227,  42,  35,  43 },  // left = d207
    {  42,  44,  44, 104, 105, 164,  64, 130,  80 },  // left = d63
    {  43,  81,  53, 140, 169, 204,  68,  84,  72 }   // left = tm
  }
};

static const vpx_prob default_if_y_probs[BLOCK_SIZE_GROUPS][INTRA_MODES - 1] = {
  {  65,  32,  18, 144, 162, 194,  41,  51,  98 },  // block_size < 8x8
  { 132,  68,  18, 165, 217, 196,  45,  40,  78 },  // block_size < 16x16
  { 173,  80,  19, 176, 240, 193,  64,  35,  46 },  // block_size < 32x32
  { 221, 135,  38, 194, 248, 121,  96,  85,  29 }   // block_size >= 32x32
};

static const vpx_prob default_uv_probs[INTRA_MODES][INTRA_MODES - 1] = {
  { 120,   7,  76, 176, 208, 126,  28,  54, 103 },  // y = dc
  {  48,  12, 154, 155, 139,  90,  34, 117, 119 },  // y = v
  {  67,   6,  25, 204, 243, 158,  13,  21,  96 },  // y = h
  {  97,   5,  44, 131, 176, 139,  48,  68,  97 },  // y = d45
  {  83,   5,  42, 156, 111, 152,  26,  49, 152 },  // y = d135
  {  80,   5,  58, 178,  74,  83,  33,  62, 145 },  // y = d117
  {  86,   5,  32, 154, 192, 168,  14,  22, 163 },  // y = d153
  {  85,   5,  32, 156, 216, 148,  19,  29,  73 },  // y = d207
  {  77,   7,  64, 116, 132, 122,  37, 126, 120 },  // y = d63
  { 101,  21, 107, 181, 192, 103,  19,  67, 125 }   // y = tm
};

#if CONFIG_EXT_PARTITION_TYPES
static const vpx_prob default_partition_probs[PARTITION_CONTEXTS]
                                             [EXT_PARTITION_TYPES - 1] = {
  // 8x8 -> 4x4
  { 199, 122, 141, 128, 128, 128, 128 },  // a/l both not split
  { 147,  63, 159, 128, 128, 128, 128 },  // a split, l not split
  { 148, 133, 118, 128, 128, 128, 128 },  // l split, a not split
  { 121, 104, 114, 128, 128, 128, 128 },  // a/l both split
  // 16x16 -> 8x8
  { 174,  73,  87, 128, 128, 128, 128 },  // a/l both not split
  {  92,  41,  83, 128, 128, 128, 128 },  // a split, l not split
  {  82,  99,  50, 128, 128, 128, 128 },  // l split, a not split
  {  53,  39,  39, 128, 128, 128, 128 },  // a/l both split
  // 32x32 -> 16x16
  { 177,  58,  59, 128, 128, 128, 128 },  // a/l both not split
  {  68,  26,  63, 128, 128, 128, 128 },  // a split, l not split
  {  52,  79,  25, 128, 128, 128, 128 },  // l split, a not split
  {  17,  14,  12, 128, 128, 128, 128 },  // a/l both split
  // 64x64 -> 32x32
  { 222,  34,  30, 128, 128, 128, 128 },  // a/l both not split
  {  72,  16,  44, 128, 128, 128, 128 },  // a split, l not split
  {  58,  32,  12, 128, 128, 128, 128 },  // l split, a not split
  {  10,   7,   6, 128, 128, 128, 128 },  // a/l both split
#if CONFIG_EXT_PARTITION
  // 128x128 -> 64x64
  { 222,  34,  30, 128, 128, 128, 128 },  // a/l both not split
  {  72,  16,  44, 128, 128, 128, 128 },  // a split, l not split
  {  58,  32,  12, 128, 128, 128, 128 },  // l split, a not split
  {  10,   7,   6, 128, 128, 128, 128 },  // a/l both split
#endif  // CONFIG_EXT_PARTITION
};
#else
static const vpx_prob default_partition_probs[PARTITION_CONTEXTS]
                                             [PARTITION_TYPES - 1] = {
  // 8x8 -> 4x4
  { 199, 122, 141 },  // a/l both not split
  { 147,  63, 159 },  // a split, l not split
  { 148, 133, 118 },  // l split, a not split
  { 121, 104, 114 },  // a/l both split
  // 16x16 -> 8x8
  { 174,  73,  87 },  // a/l both not split
  {  92,  41,  83 },  // a split, l not split
  {  82,  99,  50 },  // l split, a not split
  {  53,  39,  39 },  // a/l both split
  // 32x32 -> 16x16
  { 177,  58,  59 },  // a/l both not split
  {  68,  26,  63 },  // a split, l not split
  {  52,  79,  25 },  // l split, a not split
  {  17,  14,  12 },  // a/l both split
  // 64x64 -> 32x32
  { 222,  34,  30 },  // a/l both not split
  {  72,  16,  44 },  // a split, l not split
  {  58,  32,  12 },  // l split, a not split
  {  10,   7,   6 },  // a/l both split
#if CONFIG_EXT_PARTITION
  // 128x128 -> 64x64
  { 222,  34,  30 },  // a/l both not split
  {  72,  16,  44 },  // a split, l not split
  {  58,  32,  12 },  // l split, a not split
  {  10,   7,   6 },  // a/l both split
#endif  // CONFIG_EXT_PARTITION
};
#endif  // CONFIG_EXT_PARTITION_TYPES

#if CONFIG_REF_MV
static const vpx_prob default_newmv_prob[NEWMV_MODE_CONTEXTS] = {
    200, 180, 150, 150, 110, 70, 60,
};

static const vpx_prob default_zeromv_prob[ZEROMV_MODE_CONTEXTS] = {
    192, 64,
};

static const vpx_prob default_refmv_prob[REFMV_MODE_CONTEXTS] = {
    220, 220, 200, 200, 180, 128, 30, 220, 30,
};

static const vpx_prob default_drl_prob[DRL_MODE_CONTEXTS] = {
    128, 160, 180, 128, 160
};

#if CONFIG_EXT_INTER
static const vpx_prob default_new2mv_prob = 180;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV

static const vpx_prob default_inter_mode_probs[INTER_MODE_CONTEXTS]
                                              [INTER_MODES - 1] = {
#if CONFIG_EXT_INTER
  // TODO(zoeliu): To adjust the initial default probs
  {2,       173,   34,   173},  // 0 = both zero mv
  {7,       145,   85,   145},  // 1 = one zero mv + one a predicted mv
  {7,       166,   63,   166},  // 2 = two predicted mvs
  {7,       94,    66,   128},  // 3 = one predicted/zero and one new mv
  {8,       64,    46,   128},  // 4 = two new mvs
  {17,      81,    31,   128},  // 5 = one intra neighbour + x
  {25,      29,    30,    96},  // 6 = two intra neighbours
#else
  {2,       173,   34},  // 0 = both zero mv
  {7,       145,   85},  // 1 = one zero mv + one a predicted mv
  {7,       166,   63},  // 2 = two predicted mvs
  {7,       94,    66},  // 3 = one predicted/zero and one new mv
  {8,       64,    46},  // 4 = two new mvs
  {17,      81,    31},  // 5 = one intra neighbour + x
  {25,      29,    30},  // 6 = two intra neighbours
#endif  // CONFIG_EXT_INTER
};

#if CONFIG_EXT_INTER
static const vpx_prob default_inter_compound_mode_probs
                      [INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES - 1] = {
  { 2, 173,  68, 192, 64, 192, 128, 180, 180},   // 0 = both zero mv
  { 7, 145, 160, 192, 64, 192, 128, 180, 180},   // 1 = 1 zero + 1 predicted
  { 7, 166, 126, 192, 64, 192, 128, 180, 180},   // 2 = two predicted mvs
  { 7,  94, 132, 192, 64, 192, 128, 180, 180},   // 3 = 1 pred/zero, 1 new
  { 8,  64,  64, 192, 64, 192, 128, 180, 180},   // 4 = two new mvs
  {17,  81,  52, 192, 64, 192, 128, 180, 180},   // 5 = one intra neighbour
  {25,  29,  50, 192, 64, 192, 128, 180, 180},   // 6 = two intra neighbours
};

static const vpx_prob default_interintra_prob[BLOCK_SIZE_GROUPS] = {
  208, 208, 208, 208,
};

static const vpx_prob
    default_interintra_mode_prob[BLOCK_SIZE_GROUPS][INTERINTRA_MODES - 1] = {
  {  65,  32,  18, 144, 162, 194,  41,  51,  98 },  // block_size < 8x8
  { 132,  68,  18, 165, 217, 196,  45,  40,  78 },  // block_size < 16x16
  { 173,  80,  19, 176, 240, 193,  64,  35,  46 },  // block_size < 32x32
  { 221, 135,  38, 194, 248, 121,  96,  85,  29 }   // block_size >= 32x32
};

static const vpx_prob default_wedge_interintra_prob[BLOCK_SIZES] = {
  208, 208, 208, 208, 208, 208, 216, 216, 216, 224, 224, 224, 240,
#if CONFIG_EXT_PARTITION
  208, 208, 208
#endif  // CONFIG_EXT_PARTITION
};

static const vpx_prob default_wedge_interinter_prob[BLOCK_SIZES] = {
  208, 208, 208, 208, 208, 208, 216, 216, 216, 224, 224, 224, 240,
#if CONFIG_EXT_PARTITION
  255, 255, 255
#endif  // CONFIG_EXT_PARTITION
};
#endif  // CONFIG_EXT_INTER

// Change this section appropriately once warped motion is supported
#if CONFIG_OBMC && !CONFIG_WARPED_MOTION
const vpx_tree_index vp10_motvar_tree[TREE_SIZE(MOTION_VARIATIONS)] = {
  -SIMPLE_TRANSLATION, -OBMC_CAUSAL
};
static
const vpx_prob default_motvar_prob[BLOCK_SIZES][MOTION_VARIATIONS - 1] = {
  {255},
  {255}, {255}, {151},
  {153}, {144}, {178},
  {165}, {160}, {207},
  {195}, {168}, {244},
#if CONFIG_EXT_PARTITION
  {252}, {252}, {252},
#endif  // CONFIG_EXT_PARTITION
};

#elif !CONFIG_OBMC && CONFIG_WARPED_MOTION

const vpx_tree_index vp10_motvar_tree[TREE_SIZE(MOTION_VARIATIONS)] = {
  -SIMPLE_TRANSLATION, -WARPED_CAUSAL
};
static
const vpx_prob default_motvar_prob[BLOCK_SIZES][MOTION_VARIATIONS - 1] = {
  {255},
  {255}, {255}, {151},
  {153}, {144}, {178},
  {165}, {160}, {207},
  {195}, {168}, {244},
#if CONFIG_EXT_PARTITION
  {252}, {252}, {252},
#endif  // CONFIG_EXT_PARTITION
};

#elif CONFIG_OBMC && CONFIG_WARPED_MOTION

const vpx_tree_index vp10_motvar_tree[TREE_SIZE(MOTION_VARIATIONS)] = {
  -SIMPLE_TRANSLATION, 2,
  -OBMC_CAUSAL, -WARPED_CAUSAL,
};
static
const vpx_prob default_motvar_prob[BLOCK_SIZES][MOTION_VARIATIONS - 1] = {
  {255, 200},
  {255, 200}, {255, 200}, {151, 200},
  {153, 200}, {144, 200}, {178, 200},
  {165, 200}, {160, 200}, {207, 200},
  {195, 200}, {168, 200}, {244, 200},
#if CONFIG_EXT_PARTITION
  {252, 200}, {252, 200}, {252, 200},
#endif  // CONFIG_EXT_PARTITION
};
#endif  // CONFIG_OBMC || !CONFIG_WARPED_MOTION

/* Array indices are identical to previously-existing INTRAMODECONTEXTNODES. */
const vpx_tree_index vp10_intra_mode_tree[TREE_SIZE(INTRA_MODES)] = {
  -DC_PRED, 2,                      /* 0 = DC_NODE */
  -TM_PRED, 4,                      /* 1 = TM_NODE */
  -V_PRED, 6,                       /* 2 = V_NODE */
  8, 12,                            /* 3 = COM_NODE */
  -H_PRED, 10,                      /* 4 = H_NODE */
  -D135_PRED, -D117_PRED,           /* 5 = D135_NODE */
  -D45_PRED, 14,                    /* 6 = D45_NODE */
  -D63_PRED, 16,                    /* 7 = D63_NODE */
  -D153_PRED, -D207_PRED            /* 8 = D153_NODE */
};

const vpx_tree_index vp10_inter_mode_tree[TREE_SIZE(INTER_MODES)] = {
  -INTER_OFFSET(ZEROMV), 2,
  -INTER_OFFSET(NEARESTMV), 4,
#if CONFIG_EXT_INTER
  -INTER_OFFSET(NEARMV), 6,
  -INTER_OFFSET(NEWMV), -INTER_OFFSET(NEWFROMNEARMV)
#else
  -INTER_OFFSET(NEARMV), -INTER_OFFSET(NEWMV)
#endif  // CONFIG_EXT_INTER
};

#if CONFIG_EXT_INTER
const vpx_tree_index vp10_interintra_mode_tree[TREE_SIZE(INTERINTRA_MODES)] = {
  -II_DC_PRED, 2,                   /* 0 = II_DC_NODE     */
  -II_TM_PRED, 4,                   /* 1 = II_TM_NODE     */
  -II_V_PRED, 6,                    /* 2 = II_V_NODE      */
  8, 12,                            /* 3 = II_COM_NODE    */
  -II_H_PRED, 10,                   /* 4 = II_H_NODE      */
  -II_D135_PRED, -II_D117_PRED,     /* 5 = II_D135_NODE   */
  -II_D45_PRED, 14,                 /* 6 = II_D45_NODE    */
  -II_D63_PRED, 16,                 /* 7 = II_D63_NODE    */
  -II_D153_PRED, -II_D207_PRED      /* 8 = II_D153_NODE   */
};

const vpx_tree_index vp10_inter_compound_mode_tree
      [TREE_SIZE(INTER_COMPOUND_MODES)] = {
  -INTER_COMPOUND_OFFSET(ZERO_ZEROMV), 2,
  -INTER_COMPOUND_OFFSET(NEAREST_NEARESTMV), 4,
  6, -INTER_COMPOUND_OFFSET(NEW_NEWMV),
  8, 12,
  -INTER_COMPOUND_OFFSET(NEAR_NEARMV), 10,
  -INTER_COMPOUND_OFFSET(NEAREST_NEARMV),
      -INTER_COMPOUND_OFFSET(NEAR_NEARESTMV),
  14, 16,
  -INTER_COMPOUND_OFFSET(NEAREST_NEWMV), -INTER_COMPOUND_OFFSET(NEW_NEARESTMV),
  -INTER_COMPOUND_OFFSET(NEAR_NEWMV), -INTER_COMPOUND_OFFSET(NEW_NEARMV)
};
#endif  // CONFIG_EXT_INTER

const vpx_tree_index vp10_partition_tree[TREE_SIZE(PARTITION_TYPES)] = {
  -PARTITION_NONE, 2,
  -PARTITION_HORZ, 4,
  -PARTITION_VERT, -PARTITION_SPLIT
};

#if CONFIG_EXT_PARTITION_TYPES
const vpx_tree_index vp10_ext_partition_tree[TREE_SIZE(EXT_PARTITION_TYPES)] = {
  -PARTITION_NONE, 2,
  6, 4,
  8, -PARTITION_SPLIT,
  -PARTITION_HORZ, 10,
  -PARTITION_VERT, 12,
  -PARTITION_HORZ_A, -PARTITION_HORZ_B,
  -PARTITION_VERT_A, -PARTITION_VERT_B
};
#endif  // CONFIG_EXT_PARTITION_TYPES

static const vpx_prob default_intra_inter_p[INTRA_INTER_CONTEXTS] = {
  9, 102, 187, 225
};

static const vpx_prob default_comp_inter_p[COMP_INTER_CONTEXTS] = {
  239, 183, 119,  96,  41
};


#if CONFIG_EXT_REFS
static const vpx_prob default_comp_ref_p[REF_CONTEXTS][FWD_REFS - 1] = {
  // TODO(zoeliu): To adjust the initial prob values.
  {  33,  16,  16 },
  {  77,  74,  74 },
  { 142, 142, 142 },
  { 172, 170, 170 },
  { 238, 247, 247 }
};
static const vpx_prob default_comp_bwdref_p[REF_CONTEXTS][BWD_REFS - 1] = {
  { 16 }, { 74 }, { 142 }, { 170 }, { 247 }
};
#else
static const vpx_prob default_comp_ref_p[REF_CONTEXTS][COMP_REFS - 1] = {
  { 50 }, { 126 }, { 123 }, { 221 }, { 226 }
};
#endif  // CONFIG_EXT_REFS

static const vpx_prob default_single_ref_p[REF_CONTEXTS][SINGLE_REFS - 1] = {
#if CONFIG_EXT_REFS
  {  33,  16,  16,  16,  16 },
  {  77,  74,  74,  74,  74 },
  { 142, 142, 142, 142, 142 },
  { 172, 170, 170, 170, 170 },
  { 238, 247, 247, 247, 247 }
#else
  {  33,  16 },
  {  77,  74 },
  { 142, 142 },
  { 172, 170 },
  { 238, 247 }
#endif  // CONFIG_EXT_REFS
};

const vpx_tree_index vp10_palette_size_tree[TREE_SIZE(PALETTE_SIZES)] = {
    -TWO_COLORS, 2,
    -THREE_COLORS, 4,
    -FOUR_COLORS, 6,
    -FIVE_COLORS, 8,
    -SIX_COLORS, 10,
    -SEVEN_COLORS, -EIGHT_COLORS,
};

// TODO(huisu): tune these probs
const vpx_prob
vp10_default_palette_y_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1] = {
    {  96,  89, 100,  64,  77, 130},
    {  22,  15,  44,  16,  34,  82},
    {  30,  19,  57,  18,  38,  86},
    {  94,  36, 104,  23,  43,  92},
    { 116,  76, 107,  46,  65, 105},
    { 112,  82,  94,  40,  70, 112},
    { 147, 124, 123,  58,  69, 103},
    { 180, 113, 136,  49,  45, 114},
    { 107,  70,  87,  49, 154, 156},
    {  98, 105, 142,  63,  64, 152},
#if CONFIG_EXT_PARTITION
    {  98, 105, 142,  63,  64, 152},
    {  98, 105, 142,  63,  64, 152},
    {  98, 105, 142,  63,  64, 152},
#endif  // CONFIG_EXT_PARTITION
};

const vpx_prob
vp10_default_palette_uv_size_prob[PALETTE_BLOCK_SIZES][PALETTE_SIZES - 1] = {
    { 160, 196, 228, 213, 175, 230},
    {  87, 148, 208, 141, 166, 163},
    {  72, 151, 204, 139, 155, 161},
    {  78, 135, 171, 104, 120, 173},
    {  59,  92, 131,  78,  92, 142},
    {  75, 118, 149,  84,  90, 128},
    {  89,  87,  92,  66,  66, 128},
    {  67,  53,  54,  55,  66,  93},
    { 120, 130,  83, 171,  75, 214},
    {  72,  55,  66,  68,  79, 107},
#if CONFIG_EXT_PARTITION
    {  72,  55,  66,  68,  79, 107},
    {  72,  55,  66,  68,  79, 107},
    {  72,  55,  66,  68,  79, 107},
#endif  // CONFIG_EXT_PARTITION
};

const vpx_prob
vp10_default_palette_y_mode_prob[PALETTE_BLOCK_SIZES][PALETTE_Y_MODE_CONTEXTS]
                                                      = {
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
#if CONFIG_EXT_PARTITION
    { 240,  180,  100, },
    { 240,  180,  100, },
    { 240,  180,  100, },
#endif  // CONFIG_EXT_PARTITION
};


const vpx_prob vp10_default_palette_uv_mode_prob[2] = {
    253, 229
};

const vpx_tree_index
vp10_palette_color_tree[PALETTE_MAX_SIZE - 1][TREE_SIZE(PALETTE_COLORS)] = {
    {  // 2 colors
        -PALETTE_COLOR_ONE, -PALETTE_COLOR_TWO,
    },
    {  // 3 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, -PALETTE_COLOR_THREE,
    },
    {  // 4 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, 4,
        -PALETTE_COLOR_THREE, -PALETTE_COLOR_FOUR,
    },
    {  // 5 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, 4,
        -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, -PALETTE_COLOR_FIVE,
    },
    {  // 6 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, 4,
        -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8,
        -PALETTE_COLOR_FIVE, -PALETTE_COLOR_SIX,
    },
    {  // 7 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, 4,
        -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8,
        -PALETTE_COLOR_FIVE, 10,
        -PALETTE_COLOR_SIX, -PALETTE_COLOR_SEVEN,
    },
    {  // 8 colors
        -PALETTE_COLOR_ONE, 2,
        -PALETTE_COLOR_TWO, 4,
        -PALETTE_COLOR_THREE, 6,
        -PALETTE_COLOR_FOUR, 8,
        -PALETTE_COLOR_FIVE, 10,
        -PALETTE_COLOR_SIX, 12,
        -PALETTE_COLOR_SEVEN, -PALETTE_COLOR_EIGHT,
    },
};

const vpx_prob vp10_default_palette_y_color_prob
[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1] = {
    {  // 2 colors
        { 230, 255, 128, 128, 128, 128, 128 },
        { 214, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 240, 255, 128, 128, 128, 128, 128 },
        {  73, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 130, 255, 128, 128, 128, 128, 128 },
        { 227, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 188, 255, 128, 128, 128, 128, 128 },
        {  75, 255, 128, 128, 128, 128, 128 },
        { 250, 255, 128, 128, 128, 128, 128 },
        { 223, 255, 128, 128, 128, 128, 128 },
        { 252, 255, 128, 128, 128, 128, 128 },
    }, {  // 3 colors
        { 229, 137, 255, 128, 128, 128, 128 },
        { 197, 120, 255, 128, 128, 128, 128 },
        { 107, 195, 255, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        {  27, 151, 255, 128, 128, 128, 128 },
        { 230, 130, 255, 128, 128, 128, 128 },
        {  37, 230, 255, 128, 128, 128, 128 },
        {  67, 221, 255, 128, 128, 128, 128 },
        { 124, 230, 255, 128, 128, 128, 128 },
        { 195, 109, 255, 128, 128, 128, 128 },
        {  99, 122, 255, 128, 128, 128, 128 },
        { 205, 208, 255, 128, 128, 128, 128 },
        {  40, 235, 255, 128, 128, 128, 128 },
        { 251, 132, 255, 128, 128, 128, 128 },
        { 237, 186, 255, 128, 128, 128, 128 },
        { 253, 112, 255, 128, 128, 128, 128 },
    }, {  // 4 colors
        { 195,  87, 128, 255, 128, 128, 128 },
        { 143, 100, 123, 255, 128, 128, 128 },
        {  94, 124, 119, 255, 128, 128, 128 },
        {  77,  91, 130, 255, 128, 128, 128 },
        {  39, 114, 178, 255, 128, 128, 128 },
        { 222,  94, 125, 255, 128, 128, 128 },
        {  44, 203, 132, 255, 128, 128, 128 },
        {  68, 175, 122, 255, 128, 128, 128 },
        { 110, 187, 124, 255, 128, 128, 128 },
        { 152,  91, 128, 255, 128, 128, 128 },
        {  70, 109, 181, 255, 128, 128, 128 },
        { 133, 113, 164, 255, 128, 128, 128 },
        {  47, 205, 133, 255, 128, 128, 128 },
        { 247,  94, 136, 255, 128, 128, 128 },
        { 205, 122, 146, 255, 128, 128, 128 },
        { 251, 100, 141, 255, 128, 128, 128 },
    }, {  // 5 colors
        { 195,  65,  84, 125, 255, 128, 128 },
        { 150,  76,  84, 121, 255, 128, 128 },
        {  94, 110,  81, 117, 255, 128, 128 },
        {  79,  85,  91, 139, 255, 128, 128 },
        {  26, 102, 139, 127, 255, 128, 128 },
        { 220,  73,  91, 119, 255, 128, 128 },
        {  38, 203,  86, 127, 255, 128, 128 },
        {  61, 186,  72, 124, 255, 128, 128 },
        { 132, 199,  84, 128, 255, 128, 128 },
        { 172,  52,  62, 120, 255, 128, 128 },
        { 102,  89, 121, 122, 255, 128, 128 },
        { 182,  48,  69, 186, 255, 128, 128 },
        {  36, 206,  87, 126, 255, 128, 128 },
        { 249,  55,  67, 122, 255, 128, 128 },
        { 218,  88,  75, 122, 255, 128, 128 },
        { 253,  64,  80, 119, 255, 128, 128 },
    }, {  // 6 colors
        { 182,  54,  64,  75, 118, 255, 128 },
        { 126,  67,  70,  76, 116, 255, 128 },
        {  79,  92,  67,  85, 120, 255, 128 },
        {  63,  61,  81, 118, 132, 255, 128 },
        {  21,  80, 105,  83, 119, 255, 128 },
        { 215,  72,  74,  74, 111, 255, 128 },
        {  50, 176,  63,  79, 120, 255, 128 },
        {  72, 148,  66,  77, 120, 255, 128 },
        { 105, 177,  57,  78, 130, 255, 128 },
        { 150,  66,  66,  80, 127, 255, 128 },
        {  81,  76, 109,  85, 116, 255, 128 },
        { 113,  81,  62,  96, 148, 255, 128 },
        {  54, 179,  69,  82, 121, 255, 128 },
        { 244,  47,  48,  67, 118, 255, 128 },
        { 198,  83,  53,  65, 121, 255, 128 },
        { 250,  42,  51,  69, 110, 255, 128 },
    }, {  // 7 colors
        { 182,  45,  54,  62,  74, 113, 255 },
        { 124,  63,  57,  62,  77, 114, 255 },
        {  77,  80,  56,  66,  76, 117, 255 },
        {  63,  57,  69,  98,  85, 131, 255 },
        {  19,  81,  98,  63,  80, 116, 255 },
        { 215,  56,  60,  63,  68, 105, 255 },
        {  50, 174,  50,  60,  79, 118, 255 },
        {  68, 151,  50,  58,  73, 117, 255 },
        { 104, 182,  53,  57,  79, 127, 255 },
        { 156,  50,  51,  63,  77, 111, 255 },
        {  88,  67,  97,  59,  82, 120, 255 },
        { 114,  81,  46,  65, 103, 132, 255 },
        {  55, 166,  57,  66,  82, 120, 255 },
        { 245,  34,  38,  43,  63, 114, 255 },
        { 203,  68,  45,  47,  60, 118, 255 },
        { 250,  35,  37,  47,  66, 110, 255 },
    }, {  // 8 colors
        { 180,  43,  46,  50,  56,  69, 109 },
        { 116,  53,  51,  49,  57,  73, 115 },
        {  79,  70,  49,  50,  59,  74, 117 },
        {  60,  54,  57,  70,  62,  83, 129 },
        {  20,  73,  85,  52,  66,  81, 119 },
        { 213,  56,  52,  49,  53,  62, 104 },
        {  48, 161,  41,  45,  56,  77, 116 },
        {  68, 139,  40,  47,  54,  71, 116 },
        { 123, 166,  42,  43,  52,  76, 130 },
        { 153,  44,  44,  47,  54,  79, 129 },
        {  87,  64,  83,  49,  60,  75, 127 },
        { 131,  68,  43,  48,  73,  96, 130 },
        {  55, 152,  45,  51,  64,  77, 113 },
        { 243,  30,  28,  33,  41,  65, 114 },
        { 202,  56,  35,  36,  42,  63, 123 },
        { 249,  31,  29,  32,  45,  68, 111 },
    }
};

const vpx_prob vp10_default_palette_uv_color_prob
[PALETTE_MAX_SIZE - 1][PALETTE_COLOR_CONTEXTS][PALETTE_COLORS - 1] = {
    {  // 2 colors
        { 228, 255, 128, 128, 128, 128, 128 },
        { 195, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 228, 255, 128, 128, 128, 128, 128 },
        {  71, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 129, 255, 128, 128, 128, 128, 128 },
        { 206, 255, 128, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        { 136, 255, 128, 128, 128, 128, 128 },
        {  98, 255, 128, 128, 128, 128, 128 },
        { 236, 255, 128, 128, 128, 128, 128 },
        { 222, 255, 128, 128, 128, 128, 128 },
        { 249, 255, 128, 128, 128, 128, 128 },
    }, {  // 3 colors
        { 198, 136, 255, 128, 128, 128, 128 },
        { 178, 105, 255, 128, 128, 128, 128 },
        { 100, 206, 255, 128, 128, 128, 128 },
        { 128, 128, 128, 128, 128, 128, 128 },
        {  12, 136, 255, 128, 128, 128, 128 },
        { 219, 134, 255, 128, 128, 128, 128 },
        {  50, 198, 255, 128, 128, 128, 128 },
        {  61, 231, 255, 128, 128, 128, 128 },
        { 110, 209, 255, 128, 128, 128, 128 },
        { 173, 106, 255, 128, 128, 128, 128 },
        { 145, 166, 255, 128, 128, 128, 128 },
        { 156, 175, 255, 128, 128, 128, 128 },
        {  69, 183, 255, 128, 128, 128, 128 },
        { 241, 163, 255, 128, 128, 128, 128 },
        { 224, 160, 255, 128, 128, 128, 128 },
        { 246, 154, 255, 128, 128, 128, 128 },
    }, {  // 4 colors
        { 173,  88, 143, 255, 128, 128, 128 },
        { 146,  81, 127, 255, 128, 128, 128 },
        {  84, 134, 102, 255, 128, 128, 128 },
        {  69, 138, 140, 255, 128, 128, 128 },
        {  31, 103, 200, 255, 128, 128, 128 },
        { 217, 101, 139, 255, 128, 128, 128 },
        {  51, 174, 121, 255, 128, 128, 128 },
        {  64, 177, 109, 255, 128, 128, 128 },
        {  96, 179, 145, 255, 128, 128, 128 },
        { 164,  77, 114, 255, 128, 128, 128 },
        {  87,  94, 156, 255, 128, 128, 128 },
        { 105,  57, 173, 255, 128, 128, 128 },
        {  63, 158, 137, 255, 128, 128, 128 },
        { 236, 102, 156, 255, 128, 128, 128 },
        { 197, 115, 153, 255, 128, 128, 128 },
        { 245, 106, 154, 255, 128, 128, 128 },
    }, {  // 5 colors
        { 179,  64,  97, 129, 255, 128, 128 },
        { 137,  56,  88, 125, 255, 128, 128 },
        {  82, 107,  61, 118, 255, 128, 128 },
        {  59, 113,  86, 115, 255, 128, 128 },
        {  23,  88, 118, 130, 255, 128, 128 },
        { 213,  66,  90, 125, 255, 128, 128 },
        {  37, 181, 103, 121, 255, 128, 128 },
        {  47, 188,  61, 131, 255, 128, 128 },
        { 104, 185, 103, 144, 255, 128, 128 },
        { 163,  39,  76, 112, 255, 128, 128 },
        {  94,  74, 131, 126, 255, 128, 128 },
        { 142,  42, 103, 163, 255, 128, 128 },
        {  53, 162,  99, 149, 255, 128, 128 },
        { 239,  54,  84, 108, 255, 128, 128 },
        { 203,  84, 110, 147, 255, 128, 128 },
        { 248,  70, 105, 151, 255, 128, 128 },
    }, {  // 6 colors
        { 189,  50,  67,  90, 130, 255, 128 },
        { 114,  50,  55,  90, 123, 255, 128 },
        {  66,  76,  54,  82, 128, 255, 128 },
        {  43,  69,  69,  80, 129, 255, 128 },
        {  22,  59,  87,  88, 141, 255, 128 },
        { 203,  49,  68,  87, 122, 255, 128 },
        {  43, 157,  74, 104, 146, 255, 128 },
        {  54, 138,  51,  95, 138, 255, 128 },
        {  82, 171,  58, 102, 146, 255, 128 },
        { 129,  38,  59,  64, 168, 255, 128 },
        {  56,  67, 119,  92, 112, 255, 128 },
        {  96,  62,  53, 132,  82, 255, 128 },
        {  60, 147,  77, 108, 145, 255, 128 },
        { 238,  76,  73,  93, 148, 255, 128 },
        { 189,  86,  73, 103, 157, 255, 128 },
        { 246,  62,  75,  83, 167, 255, 128 },
    }, {  // 7 colors
        { 179,  42,  51,  73,  99, 134, 255 },
        { 119,  52,  52,  61,  64, 114, 255 },
        {  53,  77,  35,  65,  71, 131, 255 },
        {  38,  70,  51,  68,  89, 144, 255 },
        {  23,  65, 128,  73,  97, 131, 255 },
        { 210,  47,  52,  63,  81, 143, 255 },
        {  42, 159,  57,  68,  98, 143, 255 },
        {  49, 153,  45,  82,  93, 143, 255 },
        {  81, 169,  52,  72, 113, 151, 255 },
        { 136,  46,  35,  56,  75,  96, 255 },
        {  57,  84, 109,  47, 107, 131, 255 },
        { 128,  78,  57,  36, 128,  85, 255 },
        {  54, 149,  68,  77,  94, 153, 255 },
        { 243,  58,  50,  71,  81, 167, 255 },
        { 189,  92,  64,  70, 121, 173, 255 },
        { 248,  35,  38,  51,  82, 201, 255 },
    }, {  // 8 colors
        { 201,  40,  36,  42,  64,  92, 123 },
        { 116,  43,  33,  43,  73, 102, 128 },
        {  46,  77,  37,  69,  62,  78, 150 },
        {  40,  65,  52,  50,  76,  89, 133 },
        {  28,  48,  91,  17,  64,  77, 133 },
        { 218,  43,  43,  37,  56,  72, 163 },
        {  41, 155,  44,  83,  82, 129, 180 },
        {  44, 141,  29,  55,  64,  89, 147 },
        {  92, 166,  48,  45,  59, 126, 179 },
        { 169,  35,  49,  41,  36,  99, 139 },
        {  55,  77,  77,  56,  60,  75, 156 },
        { 155,  81,  51,  64,  57, 182, 255 },
        {  60, 134,  49,  49,  93, 128, 174 },
        { 244,  98,  51,  46,  22,  73, 238 },
        { 189,  70,  40,  87,  93,  79, 201 },
        { 248,  54,  49,  40,  29,  42, 227 },
    }
};

static const int palette_color_context_lookup[PALETTE_COLOR_CONTEXTS] = {
    // (3, 0, 0, 0), (3, 2, 0, 0), (3, 3, 2, 0), (3, 3, 2, 2),
    3993,  4235,  4378,  4380,
    // (4, 3, 3, 0), (5, 0, 0, 0), (5, 3, 0, 0), (5, 3, 2, 0),
    5720,  6655,  7018,  7040,
    // (5, 5, 0, 0), (6, 2, 0, 0), (6, 2, 2, 0), (6, 4, 0, 0),
    7260,  8228,  8250,  8470,
    // (7, 3, 0, 0), (8, 0, 0, 0), (8, 2, 0, 0), (10, 0, 0, 0)
    9680, 10648, 10890, 13310
};

const vpx_tree_index vp10_tx_size_tree[TX_SIZES - 1][TREE_SIZE(TX_SIZES)] = {
    {  // Max tx_size is 8X8
        -TX_4X4, -TX_8X8,
    },
    {  // Max tx_size is 16X16
        -TX_4X4, 2,
        -TX_8X8, -TX_16X16,
    },
    {  // Max tx_size is 32X32
        -TX_4X4, 2,
        -TX_8X8, 4,
        -TX_16X16, -TX_32X32,
    },
};

static const vpx_prob
default_tx_size_prob[TX_SIZES - 1][TX_SIZE_CONTEXTS][TX_SIZES - 1] = {
    {  // Max tx_size is 8X8
        { 100, }, { 66, },
    },
    {  // Max tx_size is 16X16
        { 20, 152, }, { 15, 101, },
    },
    {  // Max tx_size is 32X32
        { 3, 136, 37 }, { 5, 52,  13 },
    },
};

int vp10_get_palette_color_context(const uint8_t *color_map, int cols,
                                   int r, int c, int n, int *color_order) {
  int i, j, max, max_idx, temp;
  int scores[PALETTE_MAX_SIZE + 10];
  int weights[4] = {3, 2, 3, 2};
  int color_ctx = 0;
  int color_neighbors[4];

  assert(n <= PALETTE_MAX_SIZE);

  if (c - 1 >= 0)
    color_neighbors[0] = color_map[r * cols + c - 1];
  else
    color_neighbors[0] = -1;
  if (c - 1 >= 0 && r - 1 >= 0)
    color_neighbors[1] = color_map[(r - 1) * cols + c - 1];
  else
    color_neighbors[1] = -1;
  if (r - 1 >= 0)
    color_neighbors[2] = color_map[(r - 1) * cols + c];
  else
    color_neighbors[2] = -1;
  if (r - 1 >= 0 && c + 1 <= cols - 1)
    color_neighbors[3] = color_map[(r - 1) * cols + c + 1];
  else
    color_neighbors[3] = -1;

  for (i = 0; i < PALETTE_MAX_SIZE; ++i)
    color_order[i] = i;
  memset(scores, 0, PALETTE_MAX_SIZE * sizeof(scores[0]));
  for (i = 0; i < 4; ++i) {
    if (color_neighbors[i] >= 0)
      scores[color_neighbors[i]] += weights[i];
  }

  for (i = 0; i < 4; ++i) {
    max = scores[i];
    max_idx = i;
    j = i + 1;
    while (j < n) {
      if (scores[j] > max) {
        max = scores[j];
        max_idx = j;
      }
      ++j;
    }

    if (max_idx != i) {
      temp = scores[i];
      scores[i] = scores[max_idx];
      scores[max_idx] = temp;

      temp = color_order[i];
      color_order[i] = color_order[max_idx];
      color_order[max_idx] = temp;
    }
  }

  for (i = 0; i < 4; ++i)
    color_ctx = color_ctx * 11 + scores[i];

  for (i = 0; i < PALETTE_COLOR_CONTEXTS; ++i)
    if (color_ctx == palette_color_context_lookup[i]) {
      color_ctx = i;
      break;
    }

  if (color_ctx >= PALETTE_COLOR_CONTEXTS)
    color_ctx = 0;

  return color_ctx;
}

#if CONFIG_VAR_TX
static const vpx_prob default_txfm_partition_probs[TXFM_PARTITION_CONTEXTS] = {
    192, 128, 64, 192, 128, 64, 192, 128, 64,
};
#endif

static const vpx_prob default_skip_probs[SKIP_CONTEXTS] = {
  192, 128, 64
};

#if CONFIG_EXT_INTERP
static const vpx_prob default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                                    [SWITCHABLE_FILTERS - 1] = {
#if CONFIG_DUAL_FILTER
  { 235, 192, 128, 128},
  { 36, 243, 208, 128},
  { 34, 16, 128, 128},
  { 36, 243, 48, 128},
  { 34, 16, 128, 128},
  { 149, 160, 128, 128},

  { 235, 192, 128, 128},
  { 36, 243, 208, 128},
  { 34, 16, 128, 128},
  { 36, 243, 48, 128},
  { 34, 16, 128, 128},
  { 149, 160, 128, 128},

  { 235, 192, 128, 128},
  { 36, 243, 208, 128},
  { 34, 16, 128, 128},
  { 36, 243, 48, 128},
  { 34, 16, 128, 128},
  { 149, 160, 128, 128},

  { 235, 192, 128, 128},
  { 36, 243, 208, 128},
  { 34, 16, 128, 128},
  { 36, 243, 48, 128},
  { 34, 16, 128, 128},
  { 149, 160, 128, 128},
#else
  { 235, 192, 128, 128},
  { 36, 243, 208, 128},
  { 34, 16, 128, 128},
  { 36, 243, 48, 128},
  { 34, 16, 128, 128},
  { 149, 160, 128, 128},
#endif
};
#else  // CONFIG_EXT_INTERP
#if CONFIG_DUAL_FILTER
static const vpx_prob default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                                    [SWITCHABLE_FILTERS - 1] = {
  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 149, 144, },

  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 10, 3, },

  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 149, 144, },

  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 10, 3, },
};
#else
static const vpx_prob default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]
                                                    [SWITCHABLE_FILTERS - 1] = {
  { 235, 162, },
  { 36, 255, },
  { 34, 3, },
  { 149, 144, },
};
#endif
#endif  // CONFIG_EXT_INTERP

#if CONFIG_EXT_TX
const vpx_tree_index vp10_ext_tx_inter_tree[EXT_TX_SETS_INTER]
                                           [TREE_SIZE(TX_TYPES)] = {
  { // ToDo(yaowu): remove used entry 0.
    0
  }, {
    -IDTX, 2,
    4, 14,
    6, 8,
    -V_DCT, -H_DCT,
    10, 12,
    -V_ADST, -H_ADST,
    -V_FLIPADST, -H_FLIPADST,
    -DCT_DCT, 16,
    18, 24,
    20, 22,
    -ADST_DCT, -DCT_ADST,
    -FLIPADST_DCT, -DCT_FLIPADST,
    26, 28,
    -ADST_ADST, -FLIPADST_FLIPADST,
    -ADST_FLIPADST, -FLIPADST_ADST
  }, {
    -IDTX, 2,
    4, 6,
    -V_DCT, -H_DCT,
    -DCT_DCT, 8,
    10, 16,
    12, 14,
    -ADST_DCT, -DCT_ADST,
    -FLIPADST_DCT, -DCT_FLIPADST,
    18, 20,
    -ADST_ADST, -FLIPADST_FLIPADST,
    -ADST_FLIPADST, -FLIPADST_ADST
  }, {
    -IDTX, -DCT_DCT,
  }
};

const vpx_tree_index vp10_ext_tx_intra_tree[EXT_TX_SETS_INTRA]
                                           [TREE_SIZE(TX_TYPES)] = {
  {  // ToDo(yaowu): remove unused entry 0.
    0
  }, {
    -IDTX, 2,
    -DCT_DCT, 4,
    6, 8,
    -V_DCT, -H_DCT,
    -ADST_ADST, 10,
    -ADST_DCT, -DCT_ADST,
  }, {
    -IDTX, 2,
    -DCT_DCT, 4,
    -ADST_ADST, 6,
    -ADST_DCT, -DCT_ADST,
  }
};

static const vpx_prob
default_inter_ext_tx_prob[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES - 1] = {
  { // ToDo(yaowu): remove unused entry 0.
    { 0 },
    { 0 },
    { 0 },
#if EXT_TX_SIZES == 4
    { 0 },
#endif
  }, {
    { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128, 128},
    { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128, 128},
    { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128, 128},
#if EXT_TX_SIZES == 4
    { 10, 24, 30, 128, 128, 128, 128, 112, 160, 128, 128, 128, 128, 128, 128},
#endif
  }, {
    { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
    { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
    { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
#if EXT_TX_SIZES == 4
    { 10, 30, 128, 112, 160, 128, 128, 128, 128, 128, 128 },
#endif
  }, {
    { 12, },
    { 12, },
    { 12, },
#if EXT_TX_SIZES == 4
    { 12, },
#endif
  }
};

static const vpx_prob
default_intra_ext_tx_prob[EXT_TX_SETS_INTRA][EXT_TX_SIZES]
                         [INTRA_MODES][TX_TYPES - 1] = {
  { // ToDo(yaowu): remove unused entry 0.
    {
      { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    }, {
      { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    }, {
      { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
#if EXT_TX_SIZES == 4
    }, {
      { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
#endif
    },
  }, {
    {
      {   8, 224,  32, 128, 64, 128, },
      {  10,  32,  32, 128, 16, 192, },
      {  10,  32,  32, 128, 16,  64, },
      {   9, 200,  32, 128, 64, 128, },
      {   8,   8,  32, 128, 224, 128, },
      {  10,  32,  32, 128, 16, 192, },
      {  10,  32,  32, 128, 16,  64, },
      {  10,  23,  32, 128, 80, 176, },
      {  10,  23,  32, 128, 80, 176, },
      {  10,  32,  32, 128, 16,  64, },
    }, {
      {   8, 224, 32, 128,  64, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {   9, 200, 32, 128,  64, 128, },
      {   8,   8, 32, 128, 224, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  32, 32, 128,  16,  64, },
    }, {
      {   8, 224, 32, 128,  64, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {   9, 200, 32, 128,  64, 128, },
      {   8,   8, 32, 128, 224, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  32, 32, 128,  16,  64, },
#if EXT_TX_SIZES == 4
    }, {
      {   8, 224, 32, 128,  64, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {   9, 200, 32, 128,  64, 128, },
      {   8,   8, 32, 128, 224, 128, },
      {  10,  32, 32, 128,  16, 192, },
      {  10,  32, 32, 128,  16,  64, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  23, 32, 128,  80, 176, },
      {  10,  32, 32, 128,  16,  64, },
#endif
    },
  }, {
    {
      {   8, 224,  64, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {   9, 200,  64, 128, },
      {   8,   8, 224, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {  10,  23,  80, 176, },
      {  10,  23,  80, 176, },
      {  10,  32,  16,  64, },
    }, {
      {   8, 224,  64, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {   9, 200,  64, 128, },
      {   8,   8, 224, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {  10,  23,  80, 176, },
      {  10,  23,  80, 176, },
      {  10,  32,  16,  64, },
    }, {
      {   8, 224,  64, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {   9, 200,  64, 128, },
      {   8,   8, 224, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {  10,  23,  80, 176, },
      {  10,  23,  80, 176, },
      {  10,  32,  16,  64, },
#if EXT_TX_SIZES == 4
    }, {
      {   8, 224,  64, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {   9, 200,  64, 128, },
      {   8,   8, 224, 128, },
      {  10,  32,  16, 192, },
      {  10,  32,  16,  64, },
      {  10,  23,  80, 176, },
      {  10,  23,  80, 176, },
      {  10,  32,  16,  64, },
#endif
    },
  },
};

#else

const vpx_tree_index vp10_ext_tx_tree[TREE_SIZE(TX_TYPES)] = {
  -DCT_DCT, 2,
  -ADST_ADST, 4,
  -ADST_DCT, -DCT_ADST
};

static const vpx_prob default_intra_ext_tx_prob[EXT_TX_SIZES]
                                               [TX_TYPES][TX_TYPES - 1] = {
  {{240, 85, 128}, {4, 1, 248}, {4, 1, 8}, {4, 248, 128}},
  {{244, 85, 128}, {8, 2, 248}, {8, 2, 8}, {8, 248, 128}},
  {{248, 85, 128}, {16, 4, 248}, {16, 4, 8}, {16, 248, 128}},
};

static const vpx_prob default_inter_ext_tx_prob[EXT_TX_SIZES]
                                               [TX_TYPES - 1] = {
  {160, 85, 128},
  {176, 85, 128},
  {192, 85, 128},
};
#endif  // CONFIG_EXT_TX

#if CONFIG_EXT_INTRA
static const vpx_prob
default_intra_filter_probs[INTRA_FILTERS + 1][INTRA_FILTERS - 1] = {
  { 98,  63,  60,  },
  { 98,  82,  80,  },
  { 94,  65, 103,  },
  { 49,  25,  24,  },
  { 72,  38,  50,  },
};
static const vpx_prob default_ext_intra_probs[2] = {230, 230};

const vpx_tree_index vp10_intra_filter_tree[TREE_SIZE(INTRA_FILTERS)] = {
  -INTRA_FILTER_LINEAR, 2,
  -INTRA_FILTER_8TAP, 4,
  -INTRA_FILTER_8TAP_SHARP, -INTRA_FILTER_8TAP_SMOOTH,
};
#endif  // CONFIG_EXT_INTRA

#if CONFIG_SUPERTX
static const vpx_prob default_supertx_prob[PARTITION_SUPERTX_CONTEXTS]
                                          [TX_SIZES] = {
  { 1, 160, 160, 170 },
  { 1, 200, 200, 210 },
};
#endif  // CONFIG_SUPERTX

// FIXME(someone) need real defaults here
static const struct segmentation_probs default_seg_probs = {
  { 128, 128, 128, 128, 128, 128, 128 },
  { 128, 128, 128 },
};

static void init_mode_probs(FRAME_CONTEXT *fc) {
  vp10_copy(fc->uv_mode_prob, default_uv_probs);
  vp10_copy(fc->y_mode_prob, default_if_y_probs);
  vp10_copy(fc->switchable_interp_prob, default_switchable_interp_prob);
  vp10_copy(fc->partition_prob, default_partition_probs);
  vp10_copy(fc->intra_inter_prob, default_intra_inter_p);
  vp10_copy(fc->comp_inter_prob, default_comp_inter_p);
  vp10_copy(fc->comp_ref_prob, default_comp_ref_p);
#if CONFIG_EXT_REFS
  vp10_copy(fc->comp_bwdref_prob, default_comp_bwdref_p);
#endif  // CONFIG_EXT_REFS
  vp10_copy(fc->single_ref_prob, default_single_ref_p);
  vp10_copy(fc->tx_size_probs, default_tx_size_prob);
#if CONFIG_VAR_TX
  vp10_copy(fc->txfm_partition_prob, default_txfm_partition_probs);
#endif
  vp10_copy(fc->skip_probs, default_skip_probs);
#if CONFIG_REF_MV
  vp10_copy(fc->newmv_prob, default_newmv_prob);
  vp10_copy(fc->zeromv_prob, default_zeromv_prob);
  vp10_copy(fc->refmv_prob, default_refmv_prob);
  vp10_copy(fc->drl_prob, default_drl_prob);
#if CONFIG_EXT_INTER
  fc->new2mv_prob = default_new2mv_prob;
#endif  // CONFIG_EXT_INTER
#endif  // CONFIG_REF_MV
  vp10_copy(fc->inter_mode_probs, default_inter_mode_probs);
#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  vp10_copy(fc->motvar_prob, default_motvar_prob);
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION
#if CONFIG_EXT_INTER
  vp10_copy(fc->inter_compound_mode_probs, default_inter_compound_mode_probs);
  vp10_copy(fc->interintra_prob, default_interintra_prob);
  vp10_copy(fc->interintra_mode_prob, default_interintra_mode_prob);
  vp10_copy(fc->wedge_interintra_prob, default_wedge_interintra_prob);
  vp10_copy(fc->wedge_interinter_prob, default_wedge_interinter_prob);
#endif  // CONFIG_EXT_INTER
#if CONFIG_SUPERTX
  vp10_copy(fc->supertx_prob, default_supertx_prob);
#endif  // CONFIG_SUPERTX
  vp10_copy(fc->seg.tree_probs, default_seg_probs.tree_probs);
  vp10_copy(fc->seg.pred_probs, default_seg_probs.pred_probs);
#if CONFIG_EXT_INTRA
  vp10_copy(fc->ext_intra_probs, default_ext_intra_probs);
  vp10_copy(fc->intra_filter_probs, default_intra_filter_probs);
#endif  // CONFIG_EXT_INTRA
  vp10_copy(fc->inter_ext_tx_prob, default_inter_ext_tx_prob);
  vp10_copy(fc->intra_ext_tx_prob, default_intra_ext_tx_prob);
}

#if CONFIG_EXT_INTERP
const vpx_tree_index vp10_switchable_interp_tree
[TREE_SIZE(SWITCHABLE_FILTERS)] = {
  -EIGHTTAP_REGULAR, 2,
  4, 6,
  -EIGHTTAP_SMOOTH, -EIGHTTAP_SMOOTH2,
  -MULTITAP_SHARP, -MULTITAP_SHARP2,
};
#else
const vpx_tree_index vp10_switchable_interp_tree
[TREE_SIZE(SWITCHABLE_FILTERS)] = {
  -EIGHTTAP_REGULAR, 2,
  -EIGHTTAP_SMOOTH, -MULTITAP_SHARP
};
#endif  // CONFIG_EXT_INTERP

void vp10_adapt_inter_frame_probs(VP10_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = cm->fc;
  const FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];
  const FRAME_COUNTS *counts = &cm->counts;

  for (i = 0; i < INTRA_INTER_CONTEXTS; i++)
    fc->intra_inter_prob[i] = vp10_mode_mv_merge_probs(
        pre_fc->intra_inter_prob[i], counts->intra_inter[i]);
  for (i = 0; i < COMP_INTER_CONTEXTS; i++)
    fc->comp_inter_prob[i] = vp10_mode_mv_merge_probs(
        pre_fc->comp_inter_prob[i], counts->comp_inter[i]);

#if CONFIG_EXT_REFS
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (FWD_REFS - 1); j++)
      fc->comp_ref_prob[i][j] = mode_mv_merge_probs(
          pre_fc->comp_ref_prob[i][j], counts->comp_ref[i][j]);
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (BWD_REFS - 1); j++)
      fc->comp_bwdref_prob[i][j] = mode_mv_merge_probs(
          pre_fc->comp_bwdref_prob[i][j], counts->comp_bwdref[i][j]);
#else
  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (COMP_REFS - 1); j++)
      fc->comp_ref_prob[i][j] = mode_mv_merge_probs(
          pre_fc->comp_ref_prob[i][j], counts->comp_ref[i][j]);
#endif  // CONFIG_EXT_REFS

  for (i = 0; i < REF_CONTEXTS; i++)
    for (j = 0; j < (SINGLE_REFS - 1); j++)
      fc->single_ref_prob[i][j] = vp10_mode_mv_merge_probs(
          pre_fc->single_ref_prob[i][j], counts->single_ref[i][j]);

#if CONFIG_REF_MV
  for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i)
    fc->newmv_prob[i] = vp10_mode_mv_merge_probs(pre_fc->newmv_prob[i],
                                                 counts->newmv_mode[i]);
  for (i = 0; i < ZEROMV_MODE_CONTEXTS; ++i)
    fc->zeromv_prob[i] = vp10_mode_mv_merge_probs(pre_fc->zeromv_prob[i],
                                                  counts->zeromv_mode[i]);
  for (i = 0; i < REFMV_MODE_CONTEXTS; ++i)
    fc->refmv_prob[i] = vp10_mode_mv_merge_probs(pre_fc->refmv_prob[i],
                                                 counts->refmv_mode[i]);

  for (i = 0; i < DRL_MODE_CONTEXTS; ++i)
    fc->drl_prob[i] = vp10_mode_mv_merge_probs(pre_fc->drl_prob[i],
                                               counts->drl_mode[i]);
#if CONFIG_EXT_INTER
  fc->new2mv_prob = vp10_mode_mv_merge_probs(pre_fc->new2mv_prob,
                                             counts->new2mv_mode);
#endif  // CONFIG_EXT_INTER
#else
  for (i = 0; i < INTER_MODE_CONTEXTS; i++)
    vpx_tree_merge_probs(vp10_inter_mode_tree, pre_fc->inter_mode_probs[i],
                         counts->inter_mode[i], fc->inter_mode_probs[i]);
#endif

#if CONFIG_OBMC || CONFIG_WARPED_MOTION
  for (i = BLOCK_8X8; i < BLOCK_SIZES; ++i)
    vpx_tree_merge_probs(vp10_motvar_tree, pre_fc->motvar_prob[i],
                         counts->motvar[i], fc->motvar_prob[i]);
#endif  // CONFIG_OBMC || CONFIG_WARPED_MOTION

#if CONFIG_SUPERTX
  for (i = 0; i < PARTITION_SUPERTX_CONTEXTS; ++i) {
    int j;
    for (j = 1; j < TX_SIZES; ++j) {
      fc->supertx_prob[i][j] = vp10_mode_mv_merge_probs(
          pre_fc->supertx_prob[i][j], counts->supertx[i][j]);
    }
  }
#endif  // CONFIG_SUPERTX

#if CONFIG_EXT_INTER
  for (i = 0; i < INTER_MODE_CONTEXTS; i++)
    vpx_tree_merge_probs(vp10_inter_compound_mode_tree,
                         pre_fc->inter_compound_mode_probs[i],
                         counts->inter_compound_mode[i],
                         fc->inter_compound_mode_probs[i]);
  for (i = 0; i < BLOCK_SIZE_GROUPS; ++i) {
    if (is_interintra_allowed_bsize_group(i))
      fc->interintra_prob[i] = vp10_mode_mv_merge_probs(
          pre_fc->interintra_prob[i], counts->interintra[i]);
  }
  for (i = 0; i < BLOCK_SIZE_GROUPS; i++) {
    vpx_tree_merge_probs(
        vp10_interintra_mode_tree, pre_fc->interintra_mode_prob[i],
        counts->interintra_mode[i], fc->interintra_mode_prob[i]);
  }
  for (i = 0; i < BLOCK_SIZES; ++i) {
    if (is_interintra_allowed_bsize(i) && is_interintra_wedge_used(i))
      fc->wedge_interintra_prob[i] = vp10_mode_mv_merge_probs(
          pre_fc->wedge_interintra_prob[i], counts->wedge_interintra[i]);
  }
  for (i = 0; i < BLOCK_SIZES; ++i) {
    if (is_interinter_wedge_used(i))
      fc->wedge_interinter_prob[i] = vp10_mode_mv_merge_probs(
          pre_fc->wedge_interinter_prob[i], counts->wedge_interinter[i]);
  }
#endif  // CONFIG_EXT_INTER

  for (i = 0; i < BLOCK_SIZE_GROUPS; i++)
    vpx_tree_merge_probs(vp10_intra_mode_tree, pre_fc->y_mode_prob[i],
                counts->y_mode[i], fc->y_mode_prob[i]);

  if (cm->interp_filter == SWITCHABLE) {
    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; i++)
      vpx_tree_merge_probs(vp10_switchable_interp_tree,
                           pre_fc->switchable_interp_prob[i],
                           counts->switchable_interp[i],
                           fc->switchable_interp_prob[i]);
  }
}

void vp10_adapt_intra_frame_probs(VP10_COMMON *cm) {
  int i, j;
  FRAME_CONTEXT *fc = cm->fc;
  const FRAME_CONTEXT *pre_fc = &cm->frame_contexts[cm->frame_context_idx];
  const FRAME_COUNTS *counts = &cm->counts;

  if (cm->tx_mode == TX_MODE_SELECT) {
    for (i = 0; i < TX_SIZES - 1; ++i) {
      for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
        vpx_tree_merge_probs(vp10_tx_size_tree[i],
                             pre_fc->tx_size_probs[i][j],
                             counts->tx_size[i][j],
                             fc->tx_size_probs[i][j]);
    }
  }

#if CONFIG_VAR_TX
  if (cm->tx_mode == TX_MODE_SELECT)
    for (i = 0; i < TXFM_PARTITION_CONTEXTS; ++i)
      fc->txfm_partition_prob[i] =
          vp10_mode_mv_merge_probs(pre_fc->txfm_partition_prob[i],
                              counts->txfm_partition[i]);
#endif

  for (i = 0; i < SKIP_CONTEXTS; ++i)
    fc->skip_probs[i] = vp10_mode_mv_merge_probs(
        pre_fc->skip_probs[i], counts->skip[i]);

#if CONFIG_EXT_TX
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    int s;
    for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
      if (use_inter_ext_tx_for_txsize[s][i]) {
        vpx_tree_merge_probs(vp10_ext_tx_inter_tree[s],
                             pre_fc->inter_ext_tx_prob[s][i],
                             counts->inter_ext_tx[s][i],
                             fc->inter_ext_tx_prob[s][i]);
      }
    }
    for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
      if (use_intra_ext_tx_for_txsize[s][i]) {
        int j;
        for (j = 0; j < INTRA_MODES; ++j)
          vpx_tree_merge_probs(vp10_ext_tx_intra_tree[s],
                               pre_fc->intra_ext_tx_prob[s][i][j],
                               counts->intra_ext_tx[s][i][j],
                               fc->intra_ext_tx_prob[s][i][j]);
      }
    }
  }
#else
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    for (j = 0; j < TX_TYPES; ++j)
      vpx_tree_merge_probs(vp10_ext_tx_tree,
                           pre_fc->intra_ext_tx_prob[i][j],
                           counts->intra_ext_tx[i][j],
                           fc->intra_ext_tx_prob[i][j]);
  }
  for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
    vpx_tree_merge_probs(vp10_ext_tx_tree,
                         pre_fc->inter_ext_tx_prob[i],
                         counts->inter_ext_tx[i],
                         fc->inter_ext_tx_prob[i]);
  }
#endif  // CONFIG_EXT_TX

  if (cm->seg.temporal_update) {
    for (i = 0; i < PREDICTION_PROBS; i++)
      fc->seg.pred_probs[i] = vp10_mode_mv_merge_probs(
          pre_fc->seg.pred_probs[i], counts->seg.pred[i]);

    vpx_tree_merge_probs(vp10_segment_tree, pre_fc->seg.tree_probs,
                         counts->seg.tree_mispred, fc->seg.tree_probs);
  } else {
    vpx_tree_merge_probs(vp10_segment_tree, pre_fc->seg.tree_probs,
                         counts->seg.tree_total, fc->seg.tree_probs);
  }

  for (i = 0; i < INTRA_MODES; ++i)
    vpx_tree_merge_probs(vp10_intra_mode_tree, pre_fc->uv_mode_prob[i],
                         counts->uv_mode[i], fc->uv_mode_prob[i]);

#if CONFIG_EXT_PARTITION_TYPES
  vpx_tree_merge_probs(vp10_partition_tree, pre_fc->partition_prob[0],
                       counts->partition[0], fc->partition_prob[0]);
  for (i = 1; i < PARTITION_CONTEXTS; i++)
    vpx_tree_merge_probs(vp10_ext_partition_tree, pre_fc->partition_prob[i],
                         counts->partition[i], fc->partition_prob[i]);
#else
  for (i = 0; i < PARTITION_CONTEXTS; i++)
    vpx_tree_merge_probs(vp10_partition_tree, pre_fc->partition_prob[i],
                         counts->partition[i], fc->partition_prob[i]);
#endif  // CONFIG_EXT_PARTITION_TYPES

#if CONFIG_EXT_INTRA
  for (i = 0; i < PLANE_TYPES; ++i) {
    fc->ext_intra_probs[i] = vp10_mode_mv_merge_probs(
              pre_fc->ext_intra_probs[i], counts->ext_intra[i]);
  }

  for (i = 0; i < INTRA_FILTERS + 1; ++i)
    vpx_tree_merge_probs(vp10_intra_filter_tree, pre_fc->intra_filter_probs[i],
                         counts->intra_filter[i], fc->intra_filter_probs[i]);
#endif  // CONFIG_EXT_INTRA
}

static void set_default_lf_deltas(struct loopfilter *lf) {
  lf->mode_ref_delta_enabled = 1;
  lf->mode_ref_delta_update = 1;

  lf->ref_deltas[INTRA_FRAME] = 1;
  lf->ref_deltas[LAST_FRAME] = 0;
#if CONFIG_EXT_REFS
  lf->ref_deltas[LAST2_FRAME] = lf->ref_deltas[LAST_FRAME];
  lf->ref_deltas[LAST3_FRAME] = lf->ref_deltas[LAST_FRAME];
  lf->ref_deltas[BWDREF_FRAME] = lf->ref_deltas[LAST_FRAME];
#endif  // CONFIG_EXT_REFS
  lf->ref_deltas[GOLDEN_FRAME] = -1;
  lf->ref_deltas[ALTREF_FRAME] = -1;

  lf->mode_deltas[0] = 0;
  lf->mode_deltas[1] = 0;
}

void vp10_setup_past_independence(VP10_COMMON *cm) {
  // Reset the segment feature data to the default stats:
  // Features disabled, 0, with delta coding (Default state).
  struct loopfilter *const lf = &cm->lf;

  int i;
  vp10_clearall_segfeatures(&cm->seg);
  cm->seg.abs_delta = SEGMENT_DELTADATA;

  if (cm->last_frame_seg_map && !cm->frame_parallel_decode)
    memset(cm->last_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  if (cm->current_frame_seg_map)
    memset(cm->current_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  // Reset the mode ref deltas for loop filter
  vp10_zero(lf->last_ref_deltas);
  vp10_zero(lf->last_mode_deltas);
  set_default_lf_deltas(lf);

  // To force update of the sharpness
  lf->last_sharpness_level = -1;
#if CONFIG_LOOP_RESTORATION
  cm->rst_info.restoration_level = -1;
#endif  // CONFIG_LOOP_RESTORATION

  vp10_default_coef_probs(cm);
  init_mode_probs(cm->fc);
  vp10_init_mv_probs(cm);
  cm->fc->initialized = 1;

  if (cm->frame_type == KEY_FRAME || cm->error_resilient_mode ||
      cm->reset_frame_context == RESET_FRAME_CONTEXT_ALL) {
    // Reset all frame contexts.
    for (i = 0; i < FRAME_CONTEXTS; ++i)
      cm->frame_contexts[i] = *cm->fc;
  } else if (cm->reset_frame_context == RESET_FRAME_CONTEXT_CURRENT) {
    // Reset only the frame context specified in the frame header.
    cm->frame_contexts[cm->frame_context_idx] = *cm->fc;
  }

  // prev_mip will only be allocated in encoder.
  if (frame_is_intra_only(cm) && cm->prev_mip && !cm->frame_parallel_decode)
    memset(cm->prev_mip, 0,
           cm->mi_stride * (cm->mi_rows + 1) * sizeof(*cm->prev_mip));

  cm->frame_context_idx = 0;
}
