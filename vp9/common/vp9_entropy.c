/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <stdio.h>

#include "vp9/common/vp9_entropy.h"
#include "string.h"
#include "vp9/common/vp9_blockd.h"
#include "vp9/common/vp9_onyxc_int.h"
#include "vp9/common/vp9_entropymode.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_coefupdateprobs.h"

const int vp9_i8x8_block[4] = {0, 2, 8, 10};

DECLARE_ALIGNED(16, const uint8_t, vp9_norm[256]) = {
  0, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// Unified coefficient band structure used by all block sizes
DECLARE_ALIGNED(16, const int, vp9_coef_bands[32]) = {
  0, 1, 2, 3, 3, 3, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 5,
  5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5
};
DECLARE_ALIGNED(16, const int, vp9_coef_bands4x4[16]) = {
  0, 1, 2, 3, 3, 3, 4, 4,
  4, 4, 5, 5, 5, 5, 5, 5
};

DECLARE_ALIGNED(16, const uint8_t, vp9_pt_energy_class[MAX_ENTROPY_TOKENS]) = {
  0, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5
};

DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_4x4[16]) = {
  0,  1,  4,  8,
  5,  2,  3,  6,
  9, 12, 13, 10,
  7, 11, 14, 15,
};

DECLARE_ALIGNED(16, const int, vp9_col_scan_4x4[16]) = {
  0, 4,  8, 12,
  1, 5,  9, 13,
  2, 6, 10, 14,
  3, 7, 11, 15
};

DECLARE_ALIGNED(16, const int, vp9_row_scan_4x4[16]) = {
  0,   1,  2,  3,
  4,   5,  6,  7,
  8,   9, 10, 11,
  12, 13, 14, 15
};

DECLARE_ALIGNED(64, const int, vp9_default_zig_zag1d_8x8[64]) = {
  0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
  12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
};

DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_16x16[256]) = {
  0,   1,  16,  32,  17,   2,   3,  18,
  33,  48,  64,  49,  34,  19,   4,   5,
  20,  35,  50,  65,  80,  96,  81,  66,
  51,  36,  21,   6,   7,  22,  37,  52,
  67,  82,  97, 112, 128, 113,  98,  83,
  68,  53,  38,  23,   8,   9,  24,  39,
  54,  69,  84,  99, 114, 129, 144, 160,
  145, 130, 115, 100,  85,  70,  55,  40,
  25,  10,  11,  26,  41,  56,  71,  86,
  101, 116, 131, 146, 161, 176, 192, 177,
  162, 147, 132, 117, 102,  87,  72,  57,
  42,  27,  12,  13,  28,  43,  58, 73,
  88, 103, 118, 133, 148, 163, 178, 193,
  208, 224, 209, 194, 179, 164, 149, 134,
  119, 104,  89,  74,  59,  44,  29,  14,
  15,  30, 45,  60,  75,  90, 105, 120,
  135, 150, 165, 180, 195, 210, 225, 240,
  241, 226, 211, 196, 181, 166, 151, 136,
  121, 106,  91,  76,  61,  46,  31,  47,
  62,  77, 92, 107, 122, 137, 152, 167,
  182, 197, 212, 227, 242, 243, 228, 213,
  198, 183, 168, 153, 138, 123, 108, 93,
  78,  63,  79,  94, 109, 124, 139, 154,
  169, 184, 199, 214, 229, 244, 245, 230,
  215, 200, 185, 170, 155, 140, 125, 110,
  95, 111, 126, 141, 156, 171, 186, 201,
  216, 231, 246, 247, 232, 217, 202, 187,
  172, 157, 142, 127, 143, 158, 173, 188,
  203, 218, 233, 248, 249, 234, 219, 204,
  189, 174, 159, 175, 190, 205, 220, 235,
  250, 251, 236, 221, 206, 191, 207, 222,
  237, 252, 253, 238, 223, 239, 254, 255,
};

DECLARE_ALIGNED(16, const int, vp9_default_zig_zag1d_32x32[1024]) = {
    0,    1,   32,   64,   33,    2,    3,   34,   65,   96,  128,   97,   66,   35,    4,    5,   36,   67,   98,  129,  160,  192,  161,  130,   99,   68,   37,    6,    7,   38,   69,  100,
  131,  162,  193,  224,  256,  225,  194,  163,  132,  101,   70,   39,    8,    9,   40,   71,  102,  133,  164,  195,  226,  257,  288,  320,  289,  258,  227,  196,  165,  134,  103,   72,
   41,   10,   11,   42,   73,  104,  135,  166,  197,  228,  259,  290,  321,  352,  384,  353,  322,  291,  260,  229,  198,  167,  136,  105,   74,   43,   12,   13,   44,   75,  106,  137,
  168,  199,  230,  261,  292,  323,  354,  385,  416,  448,  417,  386,  355,  324,  293,  262,  231,  200,  169,  138,  107,   76,   45,   14,   15,   46,   77,  108,  139,  170,  201,  232,
  263,  294,  325,  356,  387,  418,  449,  480,  512,  481,  450,  419,  388,  357,  326,  295,  264,  233,  202,  171,  140,  109,   78,   47,   16,   17,   48,   79,  110,  141,  172,  203,
  234,  265,  296,  327,  358,  389,  420,  451,  482,  513,  544,  576,  545,  514,  483,  452,  421,  390,  359,  328,  297,  266,  235,  204,  173,  142,  111,   80,   49,   18,   19,   50,
   81,  112,  143,  174,  205,  236,  267,  298,  329,  360,  391,  422,  453,  484,  515,  546,  577,  608,  640,  609,  578,  547,  516,  485,  454,  423,  392,  361,  330,  299,  268,  237,
  206,  175,  144,  113,   82,   51,   20,   21,   52,   83,  114,  145,  176,  207,  238,  269,  300,  331,  362,  393,  424,  455,  486,  517,  548,  579,  610,  641,  672,  704,  673,  642,
  611,  580,  549,  518,  487,  456,  425,  394,  363,  332,  301,  270,  239,  208,  177,  146,  115,   84,   53,   22,   23,   54,   85,  116,  147,  178,  209,  240,  271,  302,  333,  364,
  395,  426,  457,  488,  519,  550,  581,  612,  643,  674,  705,  736,  768,  737,  706,  675,  644,  613,  582,  551,  520,  489,  458,  427,  396,  365,  334,  303,  272,  241,  210,  179,
  148,  117,   86,   55,   24,   25,   56,   87,  118,  149,  180,  211,  242,  273,  304,  335,  366,  397,  428,  459,  490,  521,  552,  583,  614,  645,  676,  707,  738,  769,  800,  832,
  801,  770,  739,  708,  677,  646,  615,  584,  553,  522,  491,  460,  429,  398,  367,  336,  305,  274,  243,  212,  181,  150,  119,   88,   57,   26,   27,   58,   89,  120,  151,  182,
  213,  244,  275,  306,  337,  368,  399,  430,  461,  492,  523,  554,  585,  616,  647,  678,  709,  740,  771,  802,  833,  864,  896,  865,  834,  803,  772,  741,  710,  679,  648,  617,
  586,  555,  524,  493,  462,  431,  400,  369,  338,  307,  276,  245,  214,  183,  152,  121,   90,   59,   28,   29,   60,   91,  122,  153,  184,  215,  246,  277,  308,  339,  370,  401,
  432,  463,  494,  525,  556,  587,  618,  649,  680,  711,  742,  773,  804,  835,  866,  897,  928,  960,  929,  898,  867,  836,  805,  774,  743,  712,  681,  650,  619,  588,  557,  526,
  495,  464,  433,  402,  371,  340,  309,  278,  247,  216,  185,  154,  123,   92,   61,   30,   31,   62,   93,  124,  155,  186,  217,  248,  279,  310,  341,  372,  403,  434,  465,  496,
  527,  558,  589,  620,  651,  682,  713,  744,  775,  806,  837,  868,  899,  930,  961,  992,  993,  962,  931,  900,  869,  838,  807,  776,  745,  714,  683,  652,  621,  590,  559,  528,
  497,  466,  435,  404,  373,  342,  311,  280,  249,  218,  187,  156,  125,   94,   63,   95,  126,  157,  188,  219,  250,  281,  312,  343,  374,  405,  436,  467,  498,  529,  560,  591,
  622,  653,  684,  715,  746,  777,  808,  839,  870,  901,  932,  963,  994,  995,  964,  933,  902,  871,  840,  809,  778,  747,  716,  685,  654,  623,  592,  561,  530,  499,  468,  437,
  406,  375,  344,  313,  282,  251,  220,  189,  158,  127,  159,  190,  221,  252,  283,  314,  345,  376,  407,  438,  469,  500,  531,  562,  593,  624,  655,  686,  717,  748,  779,  810,
  841,  872,  903,  934,  965,  996,  997,  966,  935,  904,  873,  842,  811,  780,  749,  718,  687,  656,  625,  594,  563,  532,  501,  470,  439,  408,  377,  346,  315,  284,  253,  222,
  191,  223,  254,  285,  316,  347,  378,  409,  440,  471,  502,  533,  564,  595,  626,  657,  688,  719,  750,  781,  812,  843,  874,  905,  936,  967,  998,  999,  968,  937,  906,  875,
  844,  813,  782,  751,  720,  689,  658,  627,  596,  565,  534,  503,  472,  441,  410,  379,  348,  317,  286,  255,  287,  318,  349,  380,  411,  442,  473,  504,  535,  566,  597,  628,
  659,  690,  721,  752,  783,  814,  845,  876,  907,  938,  969, 1000, 1001,  970,  939,  908,  877,  846,  815,  784,  753,  722,  691,  660,  629,  598,  567,  536,  505,  474,  443,  412,
  381,  350,  319,  351,  382,  413,  444,  475,  506,  537,  568,  599,  630,  661,  692,  723,  754,  785,  816,  847,  878,  909,  940,  971, 1002, 1003,  972,  941,  910,  879,  848,  817,
  786,  755,  724,  693,  662,  631,  600,  569,  538,  507,  476,  445,  414,  383,  415,  446,  477,  508,  539,  570,  601,  632,  663,  694,  725,  756,  787,  818,  849,  880,  911,  942,
  973, 1004, 1005,  974,  943,  912,  881,  850,  819,  788,  757,  726,  695,  664,  633,  602,  571,  540,  509,  478,  447,  479,  510,  541,  572,  603,  634,  665,  696,  727,  758,  789,
  820,  851,  882,  913,  944,  975, 1006, 1007,  976,  945,  914,  883,  852,  821,  790,  759,  728,  697,  666,  635,  604,  573,  542,  511,  543,  574,  605,  636,  667,  698,  729,  760,
  791,  822,  853,  884,  915,  946,  977, 1008, 1009,  978,  947,  916,  885,  854,  823,  792,  761,  730,  699,  668,  637,  606,  575,  607,  638,  669,  700,  731,  762,  793,  824,  855,
  886,  917,  948,  979, 1010, 1011,  980,  949,  918,  887,  856,  825,  794,  763,  732,  701,  670,  639,  671,  702,  733,  764,  795,  826,  857,  888,  919,  950,  981, 1012, 1013,  982,
  951,  920,  889,  858,  827,  796,  765,  734,  703,  735,  766,  797,  828,  859,  890,  921,  952,  983, 1014, 1015,  984,  953,  922,  891,  860,  829,  798,  767,  799,  830,  861,  892,
  923,  954,  985, 1016, 1017,  986,  955,  924,  893,  862,  831,  863,  894,  925,  956,  987, 1018, 1019,  988,  957,  926,  895,  927,  958,  989, 1020, 1021,  990,  959,  991, 1022, 1023,
};

/* Array indices are identical to previously-existing CONTEXT_NODE indices */

const vp9_tree_index vp9_coef_tree[ 22] =     /* corresponding _CONTEXT_NODEs */
{
  -DCT_EOB_TOKEN, 2,                             /* 0 = EOB */
  -ZERO_TOKEN, 4,                               /* 1 = ZERO */
  -ONE_TOKEN, 6,                               /* 2 = ONE */
  8, 12,                                      /* 3 = LOW_VAL */
  -TWO_TOKEN, 10,                            /* 4 = TWO */
  -THREE_TOKEN, -FOUR_TOKEN,                /* 5 = THREE */
  14, 16,                                    /* 6 = HIGH_LOW */
  -DCT_VAL_CATEGORY1, -DCT_VAL_CATEGORY2,   /* 7 = CAT_ONE */
  18, 20,                                   /* 8 = CAT_THREEFOUR */
  -DCT_VAL_CATEGORY3, -DCT_VAL_CATEGORY4,  /* 9 = CAT_THREE */
  -DCT_VAL_CATEGORY5, -DCT_VAL_CATEGORY6   /* 10 = CAT_FIVE */
};

struct vp9_token_struct vp9_coef_encodings[MAX_ENTROPY_TOKENS];

/* Trees for extra bits.  Probabilities are constant and
   do not depend on previously encoded bits */

static const vp9_prob Pcat1[] = { 159};
static const vp9_prob Pcat2[] = { 165, 145};
static const vp9_prob Pcat3[] = { 173, 148, 140};
static const vp9_prob Pcat4[] = { 176, 155, 140, 135};
static const vp9_prob Pcat5[] = { 180, 157, 141, 134, 130};
static const vp9_prob Pcat6[] = {
  254, 254, 254, 252, 249, 243, 230, 196, 177, 153, 140, 133, 130, 129
};

#if CONFIG_CODE_NONZEROCOUNT
const vp9_tree_index vp9_nzc4x4_tree[2 * NZC4X4_NODES] = {
  -NZC_0, 2,
  4, 6,
  -NZC_1, -NZC_2,
  -NZC_3TO4, 8,
  -NZC_5TO8, -NZC_9TO16,
};
struct vp9_token_struct vp9_nzc4x4_encodings[NZC4X4_TOKENS];

const vp9_tree_index vp9_nzc8x8_tree[2 * NZC8X8_NODES] = {
  -NZC_0, 2,
  4, 6,
  -NZC_1, -NZC_2,
  8, 10,
  -NZC_3TO4, -NZC_5TO8,
  -NZC_9TO16, 12,
  -NZC_17TO32, -NZC_33TO64,
};
struct vp9_token_struct vp9_nzc8x8_encodings[NZC8X8_TOKENS];

const vp9_tree_index vp9_nzc16x16_tree[2 * NZC16X16_NODES] = {
  -NZC_0, 2,
  4, 6,
  -NZC_1, -NZC_2,
  8, 10,
  -NZC_3TO4, -NZC_5TO8,
  12, 14,
  -NZC_9TO16, -NZC_17TO32,
  -NZC_33TO64, 16,
  -NZC_65TO128, -NZC_129TO256,
};
struct vp9_token_struct vp9_nzc16x16_encodings[NZC16X16_TOKENS];

const vp9_tree_index vp9_nzc32x32_tree[2 * NZC32X32_NODES] = {
  -NZC_0, 2,
  4, 6,
  -NZC_1, -NZC_2,
  8, 10,
  -NZC_3TO4, -NZC_5TO8,
  12, 14,
  -NZC_9TO16, -NZC_17TO32,
  16, 18,
  -NZC_33TO64, -NZC_65TO128,
  -NZC_129TO256, 20,
  -NZC_257TO512, -NZC_513TO1024,
};
struct vp9_token_struct vp9_nzc32x32_encodings[NZC32X32_TOKENS];

const vp9_prob Pcat_nzc[MAX_NZC_CONTEXTS]
                       [NZC_TOKENS_EXTRA][NZC_BITS_EXTRA] = { {
    // Bit probabilities are in least to most significance order
    {176,   0,   0,   0,   0,   0,   0,   0,   0},   // 3 - 4
    {164, 192,   0,   0,   0,   0,   0,   0,   0},   // 5 - 8
    {154, 184, 208,   0,   0,   0,   0,   0,   0},   // 9 - 16
    {144, 176, 200, 216,   0,   0,   0,   0,   0},   // 17 - 32
    {140, 172, 192, 208, 224,   0,   0,   0,   0},   // 33 - 64
    {136, 168, 188, 200, 220, 232,   0,   0,   0},   // 65 - 128
    {132, 164, 184, 196, 216, 228, 240,   0,   0},   // 129 - 256
    {130, 162, 178, 194, 212, 226, 240, 248,   0},   // 257 - 512
    {128, 160, 176, 192, 208, 224, 240, 248, 254},   // 513 - 1024
  }, {
    {168,   0,   0,   0,   0,   0,   0,   0,   0},   // 3 - 4
    {152, 184,   0,   0,   0,   0,   0,   0,   0},   // 5 - 8
    {152, 184, 208,   0,   0,   0,   0,   0,   0},   // 9 - 16
    {144, 176, 200, 216,   0,   0,   0,   0,   0},   // 17 - 32
    {140, 172, 192, 208, 224,   0,   0,   0,   0},   // 33 - 64
    {136, 168, 188, 200, 220, 232,   0,   0,   0},   // 65 - 128
    {132, 164, 184, 196, 216, 228, 240,   0,   0},   // 129 - 256
    {130, 162, 178, 194, 212, 226, 240, 248,   0},   // 257 - 512
    {128, 160, 176, 192, 208, 224, 240, 248, 254},   // 513 - 1024
  }, {
    {160,   0,   0,   0,   0,   0,   0,   0,   0},   // 3 - 4
    {152, 176,   0,   0,   0,   0,   0,   0,   0},   // 5 - 8
    {150, 184, 208,   0,   0,   0,   0,   0,   0},   // 9 - 16
    {144, 176, 200, 216,   0,   0,   0,   0,   0},   // 17 - 32
    {140, 172, 192, 208, 224,   0,   0,   0,   0},   // 33 - 64
    {136, 168, 188, 200, 220, 232,   0,   0,   0},   // 65 - 128
    {132, 164, 184, 196, 216, 228, 240,   0,   0},   // 129 - 256
    {130, 162, 178, 194, 212, 226, 240, 248,   0},   // 257 - 512
    {128, 160, 176, 192, 208, 224, 240, 248, 254},   // 513 - 1024
  },
};

#endif  // CONFIG_CODE_NONZEROCOUNT

static vp9_tree_index cat1[2], cat2[4], cat3[6], cat4[8], cat5[10], cat6[28];

static void init_bit_tree(vp9_tree_index *p, int n) {
  int i = 0;

  while (++i < n) {
    p[0] = p[1] = i << 1;
    p += 2;
  }

  p[0] = p[1] = 0;
}

static void init_bit_trees() {
  init_bit_tree(cat1, 1);
  init_bit_tree(cat2, 2);
  init_bit_tree(cat3, 3);
  init_bit_tree(cat4, 4);
  init_bit_tree(cat5, 5);
  init_bit_tree(cat6, 14);
}

vp9_extra_bit_struct vp9_extra_bits[12] = {
  { 0, 0, 0, 0},
  { 0, 0, 0, 1},
  { 0, 0, 0, 2},
  { 0, 0, 0, 3},
  { 0, 0, 0, 4},
  { cat1, Pcat1, 1, 5},
  { cat2, Pcat2, 2, 7},
  { cat3, Pcat3, 3, 11},
  { cat4, Pcat4, 4, 19},
  { cat5, Pcat5, 5, 35},
  { cat6, Pcat6, 14, 67},
  { 0, 0, 0, 0}
};

#include "vp9/common/vp9_default_coef_probs.h"

// This function updates and then returns n AC coefficient context
// This is currently a placeholder function to allow experimentation
// using various context models based on the energy earlier tokens
// within the current block.
//
// For now it just returns the previously used context.
int vp9_get_coef_context(int * recent_energy, int token) {
  // int token_energy;
  // int av_energy;

  /*token_energy = ((token != DCT_EOB_TOKEN) ? token : 0);
  if (!token_energy) {
    if (!(*recent_energy)) {
      av_energy = 0;
    } else {
      av_energy = 1;
    }
  } else {
    av_energy = ((token_energy + *recent_energy + 1) >> 1) + 1;
    if (av_energy > DCT_VAL_CATEGORY6)
      av_energy = DCT_VAL_CATEGORY6;
  }
  *recent_energy = token_energy;*/

  return vp9_pt_energy_class[token];
};

void vp9_default_coef_probs(VP9_COMMON *pc) {
#if CONFIG_CODE_NONZEROCOUNT
  int h, g;
#endif
  vpx_memcpy(pc->fc.coef_probs_4x4, default_coef_probs_4x4,
             sizeof(pc->fc.coef_probs_4x4));
  vpx_memcpy(pc->fc.coef_probs_8x8, default_coef_probs_8x8,
             sizeof(pc->fc.coef_probs_8x8));
  vpx_memcpy(pc->fc.coef_probs_16x16, default_coef_probs_16x16,
             sizeof(pc->fc.coef_probs_16x16));
  vpx_memcpy(pc->fc.coef_probs_32x32, default_coef_probs_32x32,
             sizeof(pc->fc.coef_probs_32x32));
#if CONFIG_CODE_NONZEROCOUNT
  for (h = 0; h < MAX_NZC_CONTEXTS; ++h) {
    for (g = 0; g < REF_TYPES; ++g) {
      int i;
      unsigned int branch_ct4x4[NZC4X4_NODES][2];
      unsigned int branch_ct8x8[NZC8X8_NODES][2];
      unsigned int branch_ct16x16[NZC16X16_NODES][2];
      unsigned int branch_ct32x32[NZC32X32_NODES][2];
      for (i = 0; i < BLOCK_TYPES; ++i) {
        vp9_tree_probs_from_distribution(
          NZC4X4_TOKENS, vp9_nzc4x4_encodings, vp9_nzc4x4_tree,
          pc->fc.nzc_probs_4x4[h][g][i], branch_ct4x4,
          default_nzc4x4_counts[h][g][i]);
      }
      for (i = 0; i < BLOCK_TYPES; ++i) {
        vp9_tree_probs_from_distribution(
          NZC8X8_TOKENS, vp9_nzc8x8_encodings, vp9_nzc8x8_tree,
          pc->fc.nzc_probs_8x8[h][g][i], branch_ct8x8,
          default_nzc8x8_counts[h][g][i]);
      }
      for (i = 0; i < BLOCK_TYPES; ++i) {
        vp9_tree_probs_from_distribution(
          NZC16X16_TOKENS, vp9_nzc16x16_encodings, vp9_nzc16x16_tree,
          pc->fc.nzc_probs_16x16[h][g][i], branch_ct16x16,
          default_nzc16x16_counts[h][g][i]);
      }
      for (i = 0; i < BLOCK_TYPES; ++i) {
        vp9_tree_probs_from_distribution(
          NZC32X32_TOKENS, vp9_nzc32x32_encodings, vp9_nzc32x32_tree,
          pc->fc.nzc_probs_32x32[h][g][i], branch_ct32x32,
          default_nzc32x32_counts[h][g][i]);
      }
    }
  }
#endif  // CONFIG_CODE_NONZEROCOUNTyy
}

void vp9_coef_tree_initialize() {
  init_bit_trees();
  vp9_tokens_from_tree(vp9_coef_encodings, vp9_coef_tree);
#if CONFIG_CODE_NONZEROCOUNT
  vp9_tokens_from_tree(vp9_nzc4x4_encodings, vp9_nzc4x4_tree);
  vp9_tokens_from_tree(vp9_nzc8x8_encodings, vp9_nzc8x8_tree);
  vp9_tokens_from_tree(vp9_nzc16x16_encodings, vp9_nzc16x16_tree);
  vp9_tokens_from_tree(vp9_nzc32x32_encodings, vp9_nzc32x32_tree);
#endif
}

#if CONFIG_CODE_NONZEROCOUNT

#define mb_in_cur_tile(cm, mb_row, mb_col)      \
    ((mb_col) >= (cm)->cur_tile_mb_col_start && \
     (mb_col) <= (cm)->cur_tile_mb_col_end   && \
     (mb_row) >= 0)

#define choose_nzc_context(nzc_exp, t2, t1)     \
    ((nzc_exp) >= ((t2) << 6) ? 2 : (nzc_exp) >= ((t1) << 6) ? 1 : 0)

#define NZC_T2_32X32    32
#define NZC_T1_32X32     8
#define NZC_T2_16X16    16
#define NZC_T1_16X16     4
#define NZC_T2_8X8       8
#define NZC_T1_8X8       2
#define NZC_T2_4X4       4
#define NZC_T1_4X4       1

// Transforms a mb16 block index to a sb64 block index
static inline int mb16_to_sb64_index(int mb_row, int mb_col, int block) {
  int r = (mb_row & 3);
  int c = (mb_col & 3);
  int b;
  if (block < 16) {  // Y
    int ib = block >> 2;
    int jb = block & 3;
    ib += r * 4;
    jb += c * 4;
    b = ib * 16 + jb;
    assert(b < 256);
    return b;
  } else {  // UV
    int base = block - (block & 3);
    int ib = (block - base) >> 1;
    int jb = (block - base) & 1;
    ib += r * 2;
    jb += c * 2;
    b = base * 16 + ib * 8 + jb;
    assert(b >= 256 && b < 384);
    return b;
  }
}

// Transforms a mb16 block index to a sb32 block index
static inline int mb16_to_sb32_index(int mb_row, int mb_col, int block) {
  int r = (mb_row & 1);
  int c = (mb_col & 1);
  int b;
  if (block < 16) {  // Y
    int ib = block >> 2;
    int jb = block & 3;
    ib += r * 4;
    jb += c * 4;
    b = ib * 8 + jb;
    assert(b < 64);
    return b;
  } else {  // UV
    int base = block - (block & 3);
    int ib = (block - base) >> 1;
    int jb = (block - base) & 1;
    ib += r * 2;
    jb += c * 2;
    b = base * 4 + ib * 4 + jb;
    assert(b >= 64 && b < 96);
    return b;
  }
}

static inline int block_to_txfm_index(int block, TX_SIZE tx_size, int s) {
  // s is the log of the number of 4x4 blocks in each row/col of larger block
  int b, ib, jb, nb;
  ib = block >> s;
  jb = block - (ib << s);
  ib >>= tx_size;
  jb >>= tx_size;
  nb = 1 << (s - tx_size);
  b = (ib * nb + jb) << (2 * tx_size);
  return b;
}

/* BEGIN - Helper functions to get the y nzcs */
static unsigned int get_nzc_4x4_y_sb64(MB_MODE_INFO *mi, int block) {
  int b;
  assert(block < 256);
  b = block_to_txfm_index(block, mi->txfm_size, 4);
  assert(b < 256);
  return mi->nzcs[b] << (6 - 2 * mi->txfm_size);
}

static unsigned int get_nzc_4x4_y_sb32(MB_MODE_INFO *mi, int block) {
  int b;
  assert(block < 64);
  b = block_to_txfm_index(block, mi->txfm_size, 3);
  assert(b < 64);
  return mi->nzcs[b] << (6 - 2 * mi->txfm_size);
}

static unsigned int get_nzc_4x4_y_mb16(MB_MODE_INFO *mi, int block) {
  int b;
  assert(block < 16);
  b = block_to_txfm_index(block, mi->txfm_size, 2);
  assert(b < 16);
  return mi->nzcs[b] << (6 - 2 * mi->txfm_size);
}
/* END - Helper functions to get the y nzcs */

/* Function to get y nzc where block index is in mb16 terms */
static unsigned int get_nzc_4x4_y(VP9_COMMON *cm, MODE_INFO *m,
                                  int mb_row, int mb_col, int block) {
  // NOTE: All values returned are at 64 times the true value at 4x4 scale
  MB_MODE_INFO *const mi = &m->mbmi;
  const int mis = cm->mode_info_stride;
  if (mi->mb_skip_coeff || !mb_in_cur_tile(cm, mb_row, mb_col))
    return 0;
  if (mi->sb_type == BLOCK_SIZE_SB64X64) {
    int r = mb_row & 3;
    int c = mb_col & 3;
    m -= c + r * mis;
    if (m->mbmi.mb_skip_coeff || !mb_in_cur_tile(cm, mb_row - r, mb_col - c))
      return 0;
    else
      return get_nzc_4x4_y_sb64(
          &m->mbmi, mb16_to_sb64_index(mb_row, mb_col, block));
  } else if (mi->sb_type == BLOCK_SIZE_SB32X32) {
    int r = mb_row & 1;
    int c = mb_col & 1;
    m -= c + r * mis;
    if (m->mbmi.mb_skip_coeff || !mb_in_cur_tile(cm, mb_row - r, mb_col - c))
      return 0;
    else
      return get_nzc_4x4_y_sb32(
          &m->mbmi, mb16_to_sb32_index(mb_row, mb_col, block));
  } else {
    if (m->mbmi.mb_skip_coeff || !mb_in_cur_tile(cm, mb_row, mb_col))
      return 0;
    return get_nzc_4x4_y_mb16(mi, block);
  }
}

/* BEGIN - Helper functions to get the uv nzcs */
static unsigned int get_nzc_4x4_uv_sb64(MB_MODE_INFO *mi, int block) {
  int b;
  int base, uvtxfm_size;
  assert(block >= 256 && block < 384);
  uvtxfm_size = mi->txfm_size;
  base = 256 + (block & 64);
  block -= base;
  b = base + block_to_txfm_index(block, uvtxfm_size, 3);
  assert(b >= 256 && b < 384);
  return mi->nzcs[b] << (6 - 2 * uvtxfm_size);
}

static unsigned int get_nzc_4x4_uv_sb32(MB_MODE_INFO *mi, int block) {
  int b;
  int base, uvtxfm_size;
  assert(block >= 64 && block < 96);
  if (mi->txfm_size == TX_32X32)
    uvtxfm_size = TX_16X16;
  else
    uvtxfm_size = mi->txfm_size;
  base = 64 + (block & 16);
  block -= base;
  b = base + block_to_txfm_index(block, uvtxfm_size, 2);
  assert(b >= 64 && b < 96);
  return mi->nzcs[b] << (6 - 2 * uvtxfm_size);
}

static unsigned int get_nzc_4x4_uv_mb16(MB_MODE_INFO *mi, int block) {
  int b;
  int base, uvtxfm_size;
  assert(block >= 16 && block < 24);
  if (mi->txfm_size == TX_8X8 &&
      (mi->mode == SPLITMV || mi->mode == I8X8_PRED))
    uvtxfm_size = TX_4X4;
  else if (mi->txfm_size == TX_16X16)
    uvtxfm_size = TX_8X8;
  else
    uvtxfm_size = mi->txfm_size;
  base = 16 + (block & 4);
  block -= base;
  b = base + block_to_txfm_index(block, uvtxfm_size, 1);
  assert(b >= 16 && b < 24);
  return mi->nzcs[b] << (6 - 2 * uvtxfm_size);
}
/* END - Helper functions to get the uv nzcs */

/* Function to get uv nzc where block index is in mb16 terms */
static unsigned int get_nzc_4x4_uv(VP9_COMMON *cm, MODE_INFO *m,
                                   int mb_row, int mb_col, int block) {
  // NOTE: All values returned are at 64 times the true value at 4x4 scale
  MB_MODE_INFO *const mi = &m->mbmi;
  const int mis = cm->mode_info_stride;
  if (mi->mb_skip_coeff || !mb_in_cur_tile(cm, mb_row, mb_col))
    return 0;
  if (mi->sb_type == BLOCK_SIZE_SB64X64) {
    int r = mb_row & 3;
    int c = mb_col & 3;
    m -= c + r * mis;
    if (m->mbmi.mb_skip_coeff || !mb_in_cur_tile(cm, mb_row - r, mb_col - c))
      return 0;
    else
      return get_nzc_4x4_uv_sb64(
          &m->mbmi, mb16_to_sb64_index(mb_row, mb_col, block));
  } else if (mi->sb_type == BLOCK_SIZE_SB32X32) {
    int r = mb_row & 1;
    int c = mb_col & 1;
    m -= c + r * mis;
    if (m->mbmi.mb_skip_coeff || !mb_in_cur_tile(cm, mb_row - r, mb_col - c))
      return 0;
    else
    return get_nzc_4x4_uv_sb32(
        &m->mbmi, mb16_to_sb32_index(mb_row, mb_col, block));
  } else {
    return get_nzc_4x4_uv_mb16(mi, block);
  }
}

int vp9_get_nzc_context_y_sb64(VP9_COMMON *cm, MODE_INFO *cur,
                               int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  assert(block < 256);
  switch (txfm_size) {
    case TX_32X32:
      assert((block & 63) == 0);
      if (block < 128) {
        int o = (block >> 6) * 2;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 12) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 13) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 14) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 15) +
            get_nzc_4x4_y(cm, cur - mis + o + 1,
                          mb_row - 1, mb_col + o + 1, 12) +
            get_nzc_4x4_y(cm, cur - mis + o + 1,
                          mb_row - 1, mb_col + o + 1, 13) +
            get_nzc_4x4_y(cm, cur - mis + o + 1,
                          mb_row - 1, mb_col + o + 1, 14) +
            get_nzc_4x4_y(cm, cur - mis + o + 1,
                          mb_row - 1, mb_col + o + 1, 15);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 128] << 3;
      }
      if ((block & 127) == 0) {
        int o = (block >> 7) * 2;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 3) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 7) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 11) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 15) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis + mis,
                          mb_row + o + 1, mb_col - 1, 3) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis + mis,
                          mb_row + o + 1, mb_col - 1, 7) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis + mis,
                          mb_row + o + 1, mb_col - 1, 11) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis + mis,
                          mb_row + o + 1, mb_col - 1, 15);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 64] << 3;
      }
      nzc_exp <<= 2;
      // Note nzc_exp is 64 times the average value expected at 32x32 scale
      return choose_nzc_context(nzc_exp, NZC_T2_32X32, NZC_T1_32X32);
      break;

    case TX_16X16:
      assert((block & 15) == 0);
      if (block < 64) {
        int o = block >> 4;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 12) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 13) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 14) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 15);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 64] << 4;
      }
      if ((block & 63) == 0) {
        int o = block >> 6;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 3) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 7) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 11) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 15);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 16] << 4;
      }
      nzc_exp <<= 1;
      // Note nzc_exp is 64 times the average value expected at 16x16 scale
      return choose_nzc_context(nzc_exp, NZC_T2_16X16, NZC_T1_16X16);
      break;

    case TX_8X8:
      assert((block & 3) == 0);
      if (block < 32) {
        int o = block >> 3;
        int p = ((block >> 2) & 1) ? 14 : 12;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, p) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, p + 1);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 32] << 5;
      }
      if ((block & 31) == 0) {
        int o = block >> 6;
        int p = ((block >> 5) & 1) ? 11 : 3;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, p) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, p + 4);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 4] << 5;
      }
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);
      break;

    case TX_4X4:
      if (block < 16) {
        int o = block >> 2;
        int p = block & 3;
        nzc_exp = get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o,
                                12 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 16] << 6);
      }
      if ((block & 15) == 0) {
        int o = block >> 6;
        int p = (block >> 4) & 3;
        nzc_exp += get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                                 3 + 4 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);
      break;

    default:
      return 0;
  }
}

int vp9_get_nzc_context_y_sb32(VP9_COMMON *cm, MODE_INFO *cur,
                               int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  assert(block < 64);
  switch (txfm_size) {
    case TX_32X32:
      assert(block == 0);
      nzc_exp =
          (get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 12) +
           get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 13) +
           get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 14) +
           get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 15) +
           get_nzc_4x4_y(cm, cur - mis + 1, mb_row - 1, mb_col + 1, 12) +
           get_nzc_4x4_y(cm, cur - mis + 1, mb_row - 1, mb_col + 1, 13) +
           get_nzc_4x4_y(cm, cur - mis + 1, mb_row - 1, mb_col + 1, 14) +
           get_nzc_4x4_y(cm, cur - mis + 1, mb_row - 1, mb_col + 1, 15) +
           get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 3) +
           get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 7) +
           get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 11) +
           get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 15) +
           get_nzc_4x4_y(cm, cur - 1 + mis, mb_row + 1, mb_col - 1, 3) +
           get_nzc_4x4_y(cm, cur - 1 + mis, mb_row + 1, mb_col - 1, 7) +
           get_nzc_4x4_y(cm, cur - 1 + mis, mb_row + 1, mb_col - 1, 11) +
           get_nzc_4x4_y(cm, cur - 1 + mis, mb_row + 1, mb_col - 1, 15)) << 2;
      // Note nzc_exp is 64 times the average value expected at 32x32 scale
      return choose_nzc_context(nzc_exp, NZC_T2_32X32, NZC_T1_32X32);
      break;

    case TX_16X16:
      assert((block & 15) == 0);
      if (block < 32) {
        int o = (block >> 4) & 1;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 12) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 13) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 14) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, 15);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 32] << 4;
      }
      if ((block & 31) == 0) {
        int o = block >> 5;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 3) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 7) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 11) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, 15);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 16] << 4;
      }
      nzc_exp <<= 1;
      // Note nzc_exp is 64 times the average value expected at 16x16 scale
      return choose_nzc_context(nzc_exp, NZC_T2_16X16, NZC_T1_16X16);
      break;

    case TX_8X8:
      assert((block & 3) == 0);
      if (block < 16) {
        int o = block >> 3;
        int p = ((block >> 2) & 1) ? 14 : 12;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, p) +
            get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o, p + 1);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 16] << 5;
      }
      if ((block & 15) == 0) {
        int o = block >> 5;
        int p = ((block >> 4) & 1) ? 11 : 3;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, p) +
            get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1, p + 4);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 4] << 5;
      }
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);
      break;

    case TX_4X4:
      if (block < 8) {
        int o = block >> 2;
        int p = block & 3;
        nzc_exp = get_nzc_4x4_y(cm, cur - mis + o, mb_row - 1, mb_col + o,
                                12 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 8] << 6);
      }
      if ((block & 7) == 0) {
        int o = block >> 5;
        int p = (block >> 3) & 3;
        nzc_exp += get_nzc_4x4_y(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                                 3 + 4 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);
      break;

    default:
      return 0;
      break;
  }
}

int vp9_get_nzc_context_y_mb16(VP9_COMMON *cm, MODE_INFO *cur,
                               int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  assert(block < 16);
  switch (txfm_size) {
    case TX_16X16:
      assert(block == 0);
      nzc_exp =
          get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 12) +
          get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 13) +
          get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 14) +
          get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, 15) +
          get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 3) +
          get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 7) +
          get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 11) +
          get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, 15);
      nzc_exp <<= 1;
      // Note nzc_exp is 64 times the average value expected at 16x16 scale
      return choose_nzc_context(nzc_exp, NZC_T2_16X16, NZC_T1_16X16);

    case TX_8X8:
      assert((block & 3) == 0);
      if (block < 8) {
        int p = ((block >> 2) & 1) ? 14 : 12;
        nzc_exp =
            get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, p) +
            get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col, p + 1);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 8] << 5;
      }
      if ((block & 7) == 0) {
        int p = ((block >> 3) & 1) ? 11 : 3;
        nzc_exp +=
            get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, p) +
            get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1, p + 4);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 4] << 5;
      }
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);

    case TX_4X4:
      if (block < 4) {
        int p = block & 3;
        nzc_exp = get_nzc_4x4_y(cm, cur - mis, mb_row - 1, mb_col,
                                12 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 4] << 6);
      }
      if ((block & 3) == 0) {
        int p = (block >> 2) & 3;
        nzc_exp += get_nzc_4x4_y(cm, cur - 1, mb_row, mb_col - 1,
                                 3 + 4 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);

    default:
      return 0;
      break;
  }
}

int vp9_get_nzc_context_uv_sb64(VP9_COMMON *cm, MODE_INFO *cur,
                                int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  const int base = block - (block & 63);
  const int boff = (block & 63);
  const int base_mb16 = base >> 4;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  TX_SIZE txfm_size_uv;

  assert(block >= 256 && block < 384);
  txfm_size_uv = txfm_size;

  switch (txfm_size_uv) {
    case TX_32X32:
      assert(block == 256 || block == 320);
      nzc_exp =
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - mis + 1, mb_row - 1, mb_col + 1,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis + 1, mb_row - 1, mb_col + 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - mis + 2, mb_row - 1, mb_col + 2,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis + 2, mb_row - 1, mb_col + 2,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - mis + 3, mb_row - 1, mb_col + 3,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis + 3, mb_row - 1, mb_col + 3,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1, mb_row, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1, mb_row, mb_col - 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row + 1, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row + 1, mb_col - 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1 + 2 * mis, mb_row + 2, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1 + 2 * mis, mb_row + 2, mb_col - 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1 + 3 * mis, mb_row + 3, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1 + 3 * mis, mb_row + 3, mb_col - 1,
                         base_mb16 + 3);
      nzc_exp <<= 2;
      // Note nzc_exp is 64 times the average value expected at 32x32 scale
      return choose_nzc_context(nzc_exp, NZC_T2_32X32, NZC_T1_32X32);

    case TX_16X16:
      // uv txfm_size 16x16
      assert((block & 15) == 0);
      if (boff < 32) {
        int o = (boff >> 4) & 1;
        nzc_exp =
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 2) +
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 3) +
            get_nzc_4x4_uv(cm, cur - mis + o + 1, mb_row - 1, mb_col + o + 1,
                           base_mb16 + 2) +
            get_nzc_4x4_uv(cm, cur - mis + o + 1, mb_row - 1, mb_col + o + 1,
                           base_mb16 + 3);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 32] << 4;
      }
      if ((boff & 31) == 0) {
        int o = boff >> 5;
        nzc_exp +=
            get_nzc_4x4_uv(cm, cur - 1 + o * mis,
                           mb_row + o, mb_col - 1, base_mb16 + 1) +
            get_nzc_4x4_uv(cm, cur - 1 + o * mis,
                           mb_row + o, mb_col - 1, base_mb16 + 3) +
            get_nzc_4x4_uv(cm, cur - 1 + o * mis + mis,
                           mb_row + o + 1, mb_col - 1, base_mb16 + 1) +
            get_nzc_4x4_uv(cm, cur - 1 + o * mis + mis,
                           mb_row + o + 1, mb_col - 1, base_mb16 + 3);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 16] << 4;
      }
      nzc_exp <<= 1;
      // Note nzc_exp is 64 times the average value expected at 16x16 scale
      return choose_nzc_context(nzc_exp, NZC_T2_16X16, NZC_T1_16X16);

    case TX_8X8:
      assert((block & 3) == 0);
      if (boff < 16) {
        int o = boff >> 2;
        nzc_exp =
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 2) +
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 3);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 16] << 5;
      }
      if ((boff & 15) == 0) {
        int o = boff >> 4;
        nzc_exp +=
            get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                           base_mb16 + 1) +
            get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                           base_mb16 + 3);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 4] << 5;
      }
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);

    case TX_4X4:
      if (boff < 8) {
        int o = boff >> 1;
        int p = boff & 1;
        nzc_exp = get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                                 base_mb16 + 2 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 8] << 6);
      }
      if ((boff & 7) == 0) {
        int o = boff >> 4;
        int p = (boff >> 3) & 1;
        nzc_exp += get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                                  base_mb16 + 1 + 2 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);

    default:
      return 0;
  }
}

int vp9_get_nzc_context_uv_sb32(VP9_COMMON *cm, MODE_INFO *cur,
                                int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  const int base = block - (block & 15);
  const int boff = (block & 15);
  const int base_mb16 = base >> 2;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  TX_SIZE txfm_size_uv;

  assert(block >= 64 && block < 96);
  if (txfm_size == TX_32X32)
    txfm_size_uv = TX_16X16;
  else
    txfm_size_uv = txfm_size;

  switch (txfm_size_uv) {
    case TX_16X16:
      // uv txfm_size 16x16
      assert(block == 64 || block == 80);
      nzc_exp =
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - mis + 1, mb_row - 1, mb_col + 1,
                         base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis + 1, mb_row - 1, mb_col + 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row, mb_col - 1,
                         base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row + 1, mb_col - 1,
                         base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1 + mis, mb_row + 1, mb_col - 1,
                         base_mb16 + 3);
      nzc_exp <<= 1;
      // Note nzc_exp is 64 times the average value expected at 16x16 scale
      return choose_nzc_context(nzc_exp, NZC_T2_16X16, NZC_T1_16X16);
      break;

    case TX_8X8:
      assert((block & 3) == 0);
      if (boff < 8) {
        int o = boff >> 2;
        nzc_exp =
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 2) +
            get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                           base_mb16 + 3);
      } else {
        nzc_exp = cur->mbmi.nzcs[block - 8] << 5;
      }
      if ((boff & 7) == 0) {
        int o = boff >> 3;
        nzc_exp +=
            get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                           base_mb16 + 1) +
            get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                           base_mb16 + 3);
      } else {
        nzc_exp += cur->mbmi.nzcs[block - 4] << 5;
      }
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);

    case TX_4X4:
      if (boff < 4) {
        int o = boff >> 1;
        int p = boff & 1;
        nzc_exp = get_nzc_4x4_uv(cm, cur - mis + o, mb_row - 1, mb_col + o,
                                 base_mb16 + 2 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 4] << 6);
      }
      if ((boff & 3) == 0) {
        int o = boff >> 3;
        int p = (boff >> 2) & 1;
        nzc_exp += get_nzc_4x4_uv(cm, cur - 1 + o * mis, mb_row + o, mb_col - 1,
                                  base_mb16 + 1 + 2 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);

    default:
      return 0;
  }
}

int vp9_get_nzc_context_uv_mb16(VP9_COMMON *cm, MODE_INFO *cur,
                                int mb_row, int mb_col, int block) {
  // returns an index in [0, MAX_NZC_CONTEXTS - 1] to reflect how busy
  // neighboring blocks are
  int mis = cm->mode_info_stride;
  int nzc_exp = 0;
  const int base = block - (block & 3);
  const int boff = (block & 3);
  const int base_mb16 = base;
  TX_SIZE txfm_size = cur->mbmi.txfm_size;
  TX_SIZE txfm_size_uv;

  assert(block >= 16 && block < 24);
  if (txfm_size == TX_16X16)
    txfm_size_uv = TX_8X8;
  else if (txfm_size == TX_8X8 &&
           (cur->mbmi.mode == I8X8_PRED || cur->mbmi.mode == SPLITMV))
    txfm_size_uv = TX_4X4;
  else
    txfm_size_uv = txfm_size;

  switch (txfm_size_uv) {
    case TX_8X8:
      assert((block & 3) == 0);
      nzc_exp =
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col, base_mb16 + 2) +
          get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col, base_mb16 + 3) +
          get_nzc_4x4_uv(cm, cur - 1, mb_row, mb_col - 1, base_mb16 + 1) +
          get_nzc_4x4_uv(cm, cur - 1, mb_row, mb_col - 1, base_mb16 + 3);
      // Note nzc_exp is 64 times the average value expected at 8x8 scale
      return choose_nzc_context(nzc_exp, NZC_T2_8X8, NZC_T1_8X8);

    case TX_4X4:
      if (boff < 2) {
        int p = boff & 1;
        nzc_exp = get_nzc_4x4_uv(cm, cur - mis, mb_row - 1, mb_col,
                                 base_mb16 + 2 + p);
      } else {
        nzc_exp = (cur->mbmi.nzcs[block - 2] << 6);
      }
      if ((boff & 1) == 0) {
        int p = (boff >> 1) & 1;
        nzc_exp += get_nzc_4x4_uv(cm, cur - 1, mb_row, mb_col - 1,
                                  base_mb16 + 1 + 2 * p);
      } else {
        nzc_exp += (cur->mbmi.nzcs[block - 1] << 6);
      }
      nzc_exp >>= 1;
      // Note nzc_exp is 64 times the average value expected at 4x4 scale
      return choose_nzc_context(nzc_exp, NZC_T2_4X4, NZC_T1_4X4);

    default:
      return 0;
  }
}

int vp9_get_nzc_context(VP9_COMMON *cm, MACROBLOCKD *xd, int block) {
  if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB64X64) {
    assert(block < 384);
    if (block < 256)
      return vp9_get_nzc_context_y_sb64(cm, xd->mode_info_context,
                                        get_mb_row(xd), get_mb_col(xd), block);
    else
      return vp9_get_nzc_context_uv_sb64(cm, xd->mode_info_context,
                                         get_mb_row(xd), get_mb_col(xd), block);
  } else if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB32X32) {
    assert(block < 96);
    if (block < 64)
      return vp9_get_nzc_context_y_sb32(cm, xd->mode_info_context,
                                        get_mb_row(xd), get_mb_col(xd), block);
    else
      return vp9_get_nzc_context_uv_sb32(cm, xd->mode_info_context,
                                         get_mb_row(xd), get_mb_col(xd), block);
  } else {
    assert(block < 64);
    if (block < 16)
      return vp9_get_nzc_context_y_mb16(cm, xd->mode_info_context,
                                        get_mb_row(xd), get_mb_col(xd), block);
    else
      return vp9_get_nzc_context_uv_mb16(cm, xd->mode_info_context,
                                         get_mb_row(xd), get_mb_col(xd), block);
  }
}

static void update_nzc(VP9_COMMON *cm,
                       uint16_t nzc,
                       int nzc_context,
                       TX_SIZE tx_size,
                       int ref,
                       int type) {
  int c;
  c = codenzc(nzc);
  if (tx_size == TX_32X32)
    cm->fc.nzc_counts_32x32[nzc_context][ref][type][c]++;
  else if (tx_size == TX_16X16)
    cm->fc.nzc_counts_16x16[nzc_context][ref][type][c]++;
  else if (tx_size == TX_8X8)
    cm->fc.nzc_counts_8x8[nzc_context][ref][type][c]++;
  else if (tx_size == TX_4X4)
    cm->fc.nzc_counts_4x4[nzc_context][ref][type][c]++;
  else
    assert(0);
  // TODO(debargha): Handle extra bits later if needed
}

static void update_nzcs_sb64(VP9_COMMON *cm,
                             MACROBLOCKD *xd,
                             int mb_row,
                             int mb_col) {
  MODE_INFO *m = xd->mode_info_context;
  MB_MODE_INFO *const mi = &m->mbmi;
  int j, nzc_context;
  const int ref = m->mbmi.ref_frame != INTRA_FRAME;

  assert(mb_col == get_mb_col(xd));
  assert(mb_row == get_mb_row(xd));

  if (mi->mb_skip_coeff)
    return;

  switch (mi->txfm_size) {
    case TX_32X32:
      for (j = 0; j < 256; j += 64) {
        nzc_context = vp9_get_nzc_context_y_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_32X32, ref, 0);
      }
      for (j = 256; j < 384; j += 64) {
        nzc_context = vp9_get_nzc_context_uv_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_32X32, ref, 1);
      }
      break;

    case TX_16X16:
      for (j = 0; j < 256; j += 16) {
        nzc_context = vp9_get_nzc_context_y_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 0);
      }
      for (j = 256; j < 384; j += 16) {
        nzc_context = vp9_get_nzc_context_uv_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 1);
      }
      break;

    case TX_8X8:
      for (j = 0; j < 256; j += 4) {
        nzc_context = vp9_get_nzc_context_y_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 0);
      }
      for (j = 256; j < 384; j += 4) {
        nzc_context = vp9_get_nzc_context_uv_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 1);
      }
      break;

    case TX_4X4:
      for (j = 0; j < 256; ++j) {
        nzc_context = vp9_get_nzc_context_y_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 0);
      }
      for (j = 256; j < 384; ++j) {
        nzc_context = vp9_get_nzc_context_uv_sb64(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 1);
      }
      break;

    default:
      break;
  }
}

static void update_nzcs_sb32(VP9_COMMON *cm,
                            MACROBLOCKD *xd,
                            int mb_row,
                            int mb_col) {
  MODE_INFO *m = xd->mode_info_context;
  MB_MODE_INFO *const mi = &m->mbmi;
  int j, nzc_context;
  const int ref = m->mbmi.ref_frame != INTRA_FRAME;

  assert(mb_col == get_mb_col(xd));
  assert(mb_row == get_mb_row(xd));

  if (mi->mb_skip_coeff)
    return;

  switch (mi->txfm_size) {
    case TX_32X32:
      for (j = 0; j < 64; j += 64) {
        nzc_context = vp9_get_nzc_context_y_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_32X32, ref, 0);
      }
      for (j = 64; j < 96; j += 16) {
        nzc_context = vp9_get_nzc_context_uv_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 1);
      }
      break;

    case TX_16X16:
      for (j = 0; j < 64; j += 16) {
        nzc_context = vp9_get_nzc_context_y_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 0);
      }
      for (j = 64; j < 96; j += 16) {
        nzc_context = vp9_get_nzc_context_uv_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 1);
      }
      break;

    case TX_8X8:
      for (j = 0; j < 64; j += 4) {
        nzc_context = vp9_get_nzc_context_y_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 0);
      }
      for (j = 64; j < 96; j += 4) {
        nzc_context = vp9_get_nzc_context_uv_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 1);
      }
      break;

    case TX_4X4:
      for (j = 0; j < 64; ++j) {
        nzc_context = vp9_get_nzc_context_y_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 0);
      }
      for (j = 64; j < 96; ++j) {
        nzc_context = vp9_get_nzc_context_uv_sb32(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 1);
      }
      break;

    default:
      break;
  }
}

static void update_nzcs_mb16(VP9_COMMON *cm,
                             MACROBLOCKD *xd,
                             int mb_row,
                             int mb_col) {
  MODE_INFO *m = xd->mode_info_context;
  MB_MODE_INFO *const mi = &m->mbmi;
  int j, nzc_context;
  const int ref = m->mbmi.ref_frame != INTRA_FRAME;

  assert(mb_col == get_mb_col(xd));
  assert(mb_row == get_mb_row(xd));

  if (mi->mb_skip_coeff)
    return;

  switch (mi->txfm_size) {
    case TX_16X16:
      for (j = 0; j < 16; j += 16) {
        nzc_context = vp9_get_nzc_context_y_mb16(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_16X16, ref, 0);
      }
      for (j = 16; j < 24; j += 4) {
        nzc_context = vp9_get_nzc_context_uv_mb16(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 1);
      }
      break;

    case TX_8X8:
      for (j = 0; j < 16; j += 4) {
        nzc_context = vp9_get_nzc_context_y_mb16(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 0);
      }
      if (mi->mode == I8X8_PRED || mi->mode == SPLITMV) {
        for (j = 16; j < 24; ++j) {
          nzc_context = vp9_get_nzc_context_uv_mb16(cm, m, mb_row, mb_col, j);
          update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 1);
        }
      } else {
        for (j = 16; j < 24; j += 4) {
          nzc_context = vp9_get_nzc_context_uv_mb16(cm, m, mb_row, mb_col, j);
          update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_8X8, ref, 1);
        }
      }
      break;

    case TX_4X4:
      for (j = 0; j < 16; ++j) {
        nzc_context = vp9_get_nzc_context_y_mb16(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 0);
      }
      for (j = 16; j < 24; ++j) {
        nzc_context = vp9_get_nzc_context_uv_mb16(cm, m, mb_row, mb_col, j);
        update_nzc(cm, m->mbmi.nzcs[j], nzc_context, TX_4X4, ref, 1);
      }
      break;

    default:
      break;
  }
}

void vp9_update_nzc_counts(VP9_COMMON *cm,
                           MACROBLOCKD *xd,
                           int mb_row,
                           int mb_col) {
  if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB64X64)
    update_nzcs_sb64(cm, xd, mb_row, mb_col);
  else if (xd->mode_info_context->mbmi.sb_type == BLOCK_SIZE_SB32X32)
    update_nzcs_sb32(cm, xd, mb_row, mb_col);
  else
    update_nzcs_mb16(cm, xd, mb_row, mb_col);
}
#endif  // CONFIG_CODE_NONZEROCOUNT

// #define COEF_COUNT_TESTING

#define COEF_COUNT_SAT 24
#define COEF_MAX_UPDATE_FACTOR 112
#define COEF_COUNT_SAT_KEY 24
#define COEF_MAX_UPDATE_FACTOR_KEY 112
#define COEF_COUNT_SAT_AFTER_KEY 24
#define COEF_MAX_UPDATE_FACTOR_AFTER_KEY 128

static void adapt_coef_probs(vp9_coeff_probs *dst_coef_probs,
                             vp9_coeff_probs *pre_coef_probs,
                             int block_types, vp9_coeff_count *coef_counts,
                             int count_sat, int update_factor) {
  int t, i, j, k, l, count;
  unsigned int branch_ct[ENTROPY_NODES][2];
  vp9_prob coef_probs[ENTROPY_NODES];
  int factor;

  for (i = 0; i < block_types; ++i)
    for (j = 0; j < REF_TYPES; ++j)
      for (k = 0; k < COEF_BANDS; ++k)
        for (l = 0; l < PREV_COEF_CONTEXTS; ++l) {
          if (l >= 3 && k == 0)
            continue;
          vp9_tree_probs_from_distribution(MAX_ENTROPY_TOKENS,
                                           vp9_coef_encodings, vp9_coef_tree,
                                           coef_probs, branch_ct,
                                           coef_counts[i][j][k][l]);
          for (t = 0; t < ENTROPY_NODES; ++t) {
            count = branch_ct[t][0] + branch_ct[t][1];
            count = count > count_sat ? count_sat : count;
            factor = (update_factor * count / count_sat);
            dst_coef_probs[i][j][k][l][t] =
                weighted_prob(pre_coef_probs[i][j][k][l][t],
                              coef_probs[t], factor);
          }
        }
}

void vp9_adapt_coef_probs(VP9_COMMON *cm) {
  int count_sat;
  int update_factor; /* denominator 256 */

  // printf("Frame type: %d\n", cm->frame_type);
  if (cm->frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_KEY;
    count_sat = COEF_COUNT_SAT_KEY;
  } else if (cm->last_frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_AFTER_KEY;  /* adapt quickly */
    count_sat = COEF_COUNT_SAT_AFTER_KEY;
  } else {
    update_factor = COEF_MAX_UPDATE_FACTOR;
    count_sat = COEF_COUNT_SAT;
  }

  adapt_coef_probs(cm->fc.coef_probs_4x4, cm->fc.pre_coef_probs_4x4,
                   BLOCK_TYPES, cm->fc.coef_counts_4x4,
                   count_sat, update_factor);
  adapt_coef_probs(cm->fc.coef_probs_8x8, cm->fc.pre_coef_probs_8x8,
                   BLOCK_TYPES, cm->fc.coef_counts_8x8,
                   count_sat, update_factor);
  adapt_coef_probs(cm->fc.coef_probs_16x16, cm->fc.pre_coef_probs_16x16,
                   BLOCK_TYPES, cm->fc.coef_counts_16x16,
                   count_sat, update_factor);
  adapt_coef_probs(cm->fc.coef_probs_32x32, cm->fc.pre_coef_probs_32x32,
                   BLOCK_TYPES, cm->fc.coef_counts_32x32,
                   count_sat, update_factor);
}

#if CONFIG_CODE_NONZEROCOUNT
static void adapt_nzc_probs(VP9_COMMON *cm,
                            int block_size,
                            int count_sat,
                            int update_factor) {
  int c, r, b, n;
  int count, factor;
  unsigned int nzc_branch_ct[NZC32X32_NODES][2];
  vp9_prob nzc_probs[NZC32X32_NODES];
  int tokens, nodes;
  const vp9_tree_index *nzc_tree;
  const struct vp9_token_struct *nzc_encodings;
  vp9_prob *dst_nzc_probs;
  vp9_prob *pre_nzc_probs;
  unsigned int *nzc_counts;

  if (block_size == 32) {
    tokens = NZC32X32_TOKENS;
    nzc_tree = vp9_nzc32x32_tree;
    nzc_encodings = vp9_nzc32x32_encodings;
    dst_nzc_probs = cm->fc.nzc_probs_32x32[0][0][0];
    pre_nzc_probs = cm->fc.pre_nzc_probs_32x32[0][0][0];
    nzc_counts = cm->fc.nzc_counts_32x32[0][0][0];
  } else if (block_size == 16) {
    tokens = NZC16X16_TOKENS;
    nzc_tree = vp9_nzc16x16_tree;
    nzc_encodings = vp9_nzc16x16_encodings;
    dst_nzc_probs = cm->fc.nzc_probs_16x16[0][0][0];
    pre_nzc_probs = cm->fc.pre_nzc_probs_16x16[0][0][0];
    nzc_counts = cm->fc.nzc_counts_16x16[0][0][0];
  } else if (block_size == 8) {
    tokens = NZC8X8_TOKENS;
    nzc_tree = vp9_nzc8x8_tree;
    nzc_encodings = vp9_nzc8x8_encodings;
    dst_nzc_probs = cm->fc.nzc_probs_8x8[0][0][0];
    pre_nzc_probs = cm->fc.pre_nzc_probs_8x8[0][0][0];
    nzc_counts = cm->fc.nzc_counts_8x8[0][0][0];
  } else {
    nzc_tree = vp9_nzc4x4_tree;
    nzc_encodings = vp9_nzc4x4_encodings;
    tokens = NZC4X4_TOKENS;
    dst_nzc_probs = cm->fc.nzc_probs_4x4[0][0][0];
    pre_nzc_probs = cm->fc.pre_nzc_probs_4x4[0][0][0];
    nzc_counts = cm->fc.nzc_counts_4x4[0][0][0];
  }
  nodes = tokens - 1;
  for (c = 0; c < MAX_NZC_CONTEXTS; ++c)
    for (r = 0; r < REF_TYPES; ++r)
      for (b = 0; b < BLOCK_TYPES; ++b) {
        int offset = c * REF_TYPES * BLOCK_TYPES + r * BLOCK_TYPES + b;
        int offset_nodes = offset * nodes;
        int offset_tokens = offset * tokens;
        vp9_tree_probs_from_distribution(tokens,
                                         nzc_encodings, nzc_tree,
                                         nzc_probs, nzc_branch_ct,
                                         nzc_counts + offset_tokens);
        for (n = 0; n < nodes; ++n) {
          count = nzc_branch_ct[n][0] + nzc_branch_ct[n][1];
          count = count > count_sat ? count_sat : count;
          factor = (update_factor * count / count_sat);
          dst_nzc_probs[offset_nodes + n] =
              weighted_prob(pre_nzc_probs[offset_nodes + n],
                            nzc_probs[n], factor);
        }
      }
}

// #define NZC_COUNT_TESTING
void vp9_adapt_nzc_probs(VP9_COMMON *cm) {
  int count_sat;
  int update_factor; /* denominator 256 */
#ifdef NZC_COUNT_TESTING
  int c, r, b, t;
  printf("\n");
  for (c = 0; c < MAX_NZC_CONTEXTS; ++c)
    for (r = 0; r < REF_TYPES; ++r) {
      for (b = 0; b < BLOCK_TYPES; ++b) {
        printf("    {");
        for (t = 0; t < NZC4X4_TOKENS; ++t) {
          printf(" %d,", cm->fc.nzc_counts_4x4[c][r][b][t]);
        }
        printf("}\n");
      }
      printf("\n");
    }
#endif

  if (cm->frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_KEY;
    count_sat = COEF_COUNT_SAT_KEY;
  } else if (cm->last_frame_type == KEY_FRAME) {
    update_factor = COEF_MAX_UPDATE_FACTOR_AFTER_KEY;  /* adapt quickly */
    count_sat = COEF_COUNT_SAT_AFTER_KEY;
  } else {
    update_factor = COEF_MAX_UPDATE_FACTOR;
    count_sat = COEF_COUNT_SAT;
  }

  adapt_nzc_probs(cm, 4, count_sat, update_factor);
  adapt_nzc_probs(cm, 8, count_sat, update_factor);
  adapt_nzc_probs(cm, 16, count_sat, update_factor);
  adapt_nzc_probs(cm, 32, count_sat, update_factor);
}
#endif  // CONFIG_CODE_NONZEROCOUNT
