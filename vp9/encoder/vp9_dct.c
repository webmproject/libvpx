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
#include "vpx_ports/config.h"
#include "vp9/common/vp9_systemdependent.h"

#include "vp9/common/vp9_blockd.h"

// TODO: these transforms can be converted into integer forms to reduce
//       the complexity
static const float dct_4[16] = {
  0.500000000000000,  0.500000000000000,  0.500000000000000,  0.500000000000000,
  0.653281482438188,  0.270598050073099, -0.270598050073099, -0.653281482438188,
  0.500000000000000, -0.500000000000000, -0.500000000000000,  0.500000000000000,
  0.270598050073099, -0.653281482438188,  0.653281482438188, -0.270598050073099
};

static const float adst_4[16] = {
  0.228013428883779,  0.428525073124360,  0.577350269189626,  0.656538502008139,
  0.577350269189626,  0.577350269189626,  0.000000000000000, -0.577350269189626,
  0.656538502008139, -0.228013428883779, -0.577350269189626,  0.428525073124359,
  0.428525073124360, -0.656538502008139,  0.577350269189626, -0.228013428883779
};

static const float dct_8[64] = {
  0.353553390593274,   0.353553390593274,   0.353553390593274,   0.353553390593274,
  0.353553390593274,   0.353553390593274,   0.353553390593274,   0.353553390593274,
  0.490392640201615,   0.415734806151273,   0.277785116509801,   0.097545161008064,
 -0.097545161008064,  -0.277785116509801,  -0.415734806151273,  -0.490392640201615,
  0.461939766255643,   0.191341716182545,  -0.191341716182545,  -0.461939766255643,
 -0.461939766255643,  -0.191341716182545,   0.191341716182545,   0.461939766255643,
  0.415734806151273,  -0.097545161008064,  -0.490392640201615,  -0.277785116509801,
  0.277785116509801,   0.490392640201615,   0.097545161008064,  -0.415734806151273,
  0.353553390593274,  -0.353553390593274,  -0.353553390593274,   0.353553390593274,
  0.353553390593274,  -0.353553390593274,  -0.353553390593274,   0.353553390593274,
  0.277785116509801,  -0.490392640201615,   0.097545161008064,   0.415734806151273,
 -0.415734806151273,  -0.097545161008064,   0.490392640201615,  -0.277785116509801,
  0.191341716182545,  -0.461939766255643,   0.461939766255643,  -0.191341716182545,
 -0.191341716182545,   0.461939766255643,  -0.461939766255643,   0.191341716182545,
  0.097545161008064,  -0.277785116509801,   0.415734806151273,  -0.490392640201615,
  0.490392640201615,  -0.415734806151273,   0.277785116509801,  -0.097545161008064
};

static const float adst_8[64] = {
  0.089131608307533,   0.175227946595735,   0.255357107325376,   0.326790388032145,
  0.387095214016349,   0.434217976756762,   0.466553967085785,   0.483002021635509,
  0.255357107325376,   0.434217976756762,   0.483002021635509,   0.387095214016349,
  0.175227946595735,  -0.089131608307533,  -0.326790388032145,  -0.466553967085785,
  0.387095214016349,   0.466553967085785,   0.175227946595735,  -0.255357107325376,
 -0.483002021635509,  -0.326790388032145,   0.089131608307533,   0.434217976756762,
  0.466553967085785,   0.255357107325376,  -0.326790388032145,  -0.434217976756762,
  0.089131608307533,   0.483002021635509,   0.175227946595735,  -0.387095214016348,
  0.483002021635509,  -0.089131608307533,  -0.466553967085785,   0.175227946595735,
  0.434217976756762,  -0.255357107325376,  -0.387095214016348,   0.326790388032145,
  0.434217976756762,  -0.387095214016348,  -0.089131608307533,   0.466553967085786,
 -0.326790388032145,  -0.175227946595735,   0.483002021635509,  -0.255357107325375,
  0.326790388032145,  -0.483002021635509,   0.387095214016349,  -0.089131608307534,
 -0.255357107325377,   0.466553967085785,  -0.434217976756762,   0.175227946595736,
  0.175227946595735,  -0.326790388032145,   0.434217976756762,  -0.483002021635509,
  0.466553967085785,  -0.387095214016348,   0.255357107325376,  -0.089131608307532
};

/* Converted the transforms to integers. */
static const int16_t dct_i4[16] = {
  16384,  16384,  16384,  16384,
  21407,   8867,  -8867, -21407,
  16384, -16384, -16384,  16384,
   8867, -21407,  21407,  -8867
};

static const int16_t adst_i4[16] = {
   7472,  14042,  18919,  21513,
  18919,  18919,      0, -18919,
  21513,  -7472, -18919,  14042,
  14042, -21513,  18919,  -7472
};

static const int16_t dct_i8[64] = {
   11585,  11585,  11585,  11585,
   11585,  11585,  11585,  11585,
   16069,  13623,   9102,   3196,
   -3196,  -9102, -13623, -16069,
   15137,   6270,  -6270, -15137,
  -15137,  -6270,   6270,  15137,
   13623,  -3196, -16069,  -9102,
    9102,  16069,   3196, -13623,
   11585, -11585, -11585,  11585,
   11585, -11585, -11585,  11585,
    9102, -16069,   3196,  13623,
  -13623,  -3196,  16069,  -9102,
    6270, -15137,  15137,  -6270,
   -6270,  15137, -15137,   6270,
    3196,  -9102,  13623, -16069,
   16069, -13623,   9102,  -3196
};

static const int16_t adst_i8[64] = {
    2921,   5742,   8368,  10708,
   12684,  14228,  15288,  15827,
    8368,  14228,  15827,  12684,
    5742,  -2921, -10708, -15288,
   12684,  15288,   5742,  -8368,
  -15827, -10708,   2921,  14228,
   15288,   8368, -10708, -14228,
    2921,  15827,   5742, -12684,
   15827,  -2921, -15288,   5742,
   14228,  -8368, -12684,  10708,
   14228, -12684,  -2921,  15288,
  -10708,  -5742,  15827,  -8368,
   10708, -15827,  12684,  -2921,
   -8368,  15288, -14228,   5742,
    5742, -10708,  14228, -15827,
   15288, -12684,   8368,  -2921
};

static const float dct_16[256] = {
  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,
  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,  0.250000,
  0.351851,  0.338330,  0.311806,  0.273300,  0.224292,  0.166664,  0.102631,  0.034654,
 -0.034654, -0.102631, -0.166664, -0.224292, -0.273300, -0.311806, -0.338330, -0.351851,
  0.346760,  0.293969,  0.196424,  0.068975, -0.068975, -0.196424, -0.293969, -0.346760,
 -0.346760, -0.293969, -0.196424, -0.068975,  0.068975,  0.196424,  0.293969,  0.346760,
  0.338330,  0.224292,  0.034654, -0.166664, -0.311806, -0.351851, -0.273300, -0.102631,
  0.102631,  0.273300,  0.351851,  0.311806,  0.166664, -0.034654, -0.224292, -0.338330,
  0.326641,  0.135299, -0.135299, -0.326641, -0.326641, -0.135299,  0.135299,  0.326641,
  0.326641,  0.135299, -0.135299, -0.326641, -0.326641, -0.135299,  0.135299,  0.326641,
  0.311806,  0.034654, -0.273300, -0.338330, -0.102631,  0.224292,  0.351851,  0.166664,
 -0.166664, -0.351851, -0.224292,  0.102631,  0.338330,  0.273300, -0.034654, -0.311806,
  0.293969, -0.068975, -0.346760, -0.196424,  0.196424,  0.346760,  0.068975, -0.293969,
 -0.293969,  0.068975,  0.346760,  0.196424, -0.196424, -0.346760, -0.068975,  0.293969,
  0.273300, -0.166664, -0.338330,  0.034654,  0.351851,  0.102631, -0.311806, -0.224292,
  0.224292,  0.311806, -0.102631, -0.351851, -0.034654,  0.338330,  0.166664, -0.273300,
  0.250000, -0.250000, -0.250000,  0.250000,  0.250000, -0.250000, -0.250000,  0.250000,
  0.250000, -0.250000, -0.250000,  0.250000,  0.250000, -0.250000, -0.250000,  0.250000,
  0.224292, -0.311806, -0.102631,  0.351851, -0.034654, -0.338330,  0.166664,  0.273300,
 -0.273300, -0.166664,  0.338330,  0.034654, -0.351851,  0.102631,  0.311806, -0.224292,
  0.196424, -0.346760,  0.068975,  0.293969, -0.293969, -0.068975,  0.346760, -0.196424,
 -0.196424,  0.346760, -0.068975, -0.293969,  0.293969,  0.068975, -0.346760,  0.196424,
  0.166664, -0.351851,  0.224292,  0.102631, -0.338330,  0.273300,  0.034654, -0.311806,
  0.311806, -0.034654, -0.273300,  0.338330, -0.102631, -0.224292,  0.351851, -0.166664,
  0.135299, -0.326641,  0.326641, -0.135299, -0.135299,  0.326641, -0.326641,  0.135299,
  0.135299, -0.326641,  0.326641, -0.135299, -0.135299,  0.326641, -0.326641,  0.135299,
  0.102631, -0.273300,  0.351851, -0.311806,  0.166664,  0.034654, -0.224292,  0.338330,
 -0.338330,  0.224292, -0.034654, -0.166664,  0.311806, -0.351851,  0.273300, -0.102631,
  0.068975, -0.196424,  0.293969, -0.346760,  0.346760, -0.293969,  0.196424, -0.068975,
 -0.068975,  0.196424, -0.293969,  0.346760, -0.346760,  0.293969, -0.196424,  0.068975,
  0.034654, -0.102631,  0.166664, -0.224292,  0.273300, -0.311806,  0.338330, -0.351851,
  0.351851, -0.338330,  0.311806, -0.273300,  0.224292, -0.166664,  0.102631, -0.034654
};

static const float adst_16[256] = {
  0.033094,  0.065889,  0.098087,  0.129396,  0.159534,  0.188227,  0.215215,  0.240255,
  0.263118,  0.283599,  0.301511,  0.316693,  0.329007,  0.338341,  0.344612,  0.347761,
  0.098087,  0.188227,  0.263118,  0.316693,  0.344612,  0.344612,  0.316693,  0.263118,
  0.188227,  0.098087,  0.000000, -0.098087, -0.188227, -0.263118, -0.316693, -0.344612,
  0.159534,  0.283599,  0.344612,  0.329007,  0.240255,  0.098087, -0.065889, -0.215215,
 -0.316693, -0.347761, -0.301511, -0.188227, -0.033094,  0.129396,  0.263118,  0.338341,
  0.215215,  0.338341,  0.316693,  0.159534, -0.065889, -0.263118, -0.347761, -0.283599,
 -0.098087,  0.129396,  0.301511,  0.344612,  0.240255,  0.033094, -0.188227, -0.329007,
  0.263118,  0.344612,  0.188227, -0.098087, -0.316693, -0.316693, -0.098087,  0.188227,
  0.344612,  0.263118,  0.000000, -0.263118, -0.344612, -0.188227,  0.098087,  0.316693,
  0.301511,  0.301511,  0.000000, -0.301511, -0.301511, -0.000000,  0.301511,  0.301511,
  0.000000, -0.301511, -0.301511, -0.000000,  0.301511,  0.301511,  0.000000, -0.301511,
  0.329007,  0.215215, -0.188227, -0.338341, -0.033094,  0.316693,  0.240255, -0.159534,
 -0.344612, -0.065889,  0.301511,  0.263118, -0.129396, -0.347761, -0.098087,  0.283599,
  0.344612,  0.098087, -0.316693, -0.188227,  0.263118,  0.263118, -0.188227, -0.316693,
  0.098087,  0.344612,  0.000000, -0.344612, -0.098087,  0.316693,  0.188227, -0.263118,
  0.347761, -0.033094, -0.344612,  0.065889,  0.338341, -0.098087, -0.329007,  0.129396,
  0.316693, -0.159534, -0.301511,  0.188227,  0.283599, -0.215215, -0.263118,  0.240255,
  0.338341, -0.159534, -0.263118,  0.283599,  0.129396, -0.344612,  0.033094,  0.329007,
 -0.188227, -0.240255,  0.301511,  0.098087, -0.347761,  0.065889,  0.316693, -0.215215,
  0.316693, -0.263118, -0.098087,  0.344612, -0.188227, -0.188227,  0.344612, -0.098087,
 -0.263118,  0.316693,  0.000000, -0.316693,  0.263118,  0.098087, -0.344612,  0.188227,
  0.283599, -0.329007,  0.098087,  0.215215, -0.347761,  0.188227,  0.129396, -0.338341,
  0.263118,  0.033094, -0.301511,  0.316693, -0.065889, -0.240255,  0.344612, -0.159534,
  0.240255, -0.347761,  0.263118, -0.033094, -0.215215,  0.344612, -0.283599,  0.065889,
  0.188227, -0.338341,  0.301511, -0.098087, -0.159534,  0.329007, -0.316693,  0.129396,
  0.188227, -0.316693,  0.344612, -0.263118,  0.098087,  0.098087, -0.263118,  0.344612,
 -0.316693,  0.188227,  0.000000, -0.188227,  0.316693, -0.344612,  0.263118, -0.098087,
  0.129396, -0.240255,  0.316693, -0.347761,  0.329007, -0.263118,  0.159534, -0.033094,
 -0.098087,  0.215215, -0.301511,  0.344612, -0.338341,  0.283599, -0.188227,  0.065889,
  0.065889, -0.129396,  0.188227, -0.240255,  0.283599, -0.316693,  0.338341, -0.347761,
  0.344612, -0.329007,  0.301511, -0.263118,  0.215215, -0.159534,  0.098087, -0.033094
};

/* Converted the transforms to integers. */
static const int16_t dct_i16[256] = {
    8192,   8192,   8192,   8192,   8192,   8192,   8192,   8192,
    8192,   8192,   8192,   8192,   8192,   8192,   8192,   8192,
   11529,  11086,  10217,   8955,   7350,   5461,   3363,   1136,
   -1136,  -3363,  -5461,  -7350,  -8955, -10217, -11086, -11529,
   11363,   9633,   6436,   2260,  -2260,  -6436,  -9633, -11363,
  -11363,  -9633,  -6436,  -2260,   2260,   6436,   9633,  11363,
   11086,   7350,   1136,  -5461, -10217, -11529,  -8955,  -3363,
    3363,   8955,  11529,  10217,   5461,  -1136,  -7350, -11086,
   10703,   4433,  -4433, -10703, -10703,  -4433,   4433,  10703,
   10703,   4433,  -4433, -10703, -10703,  -4433,   4433,  10703,
   10217,   1136,  -8955, -11086,  -3363,   7350,  11529,   5461,
   -5461, -11529,  -7350,   3363,  11086,   8955,  -1136, -10217,
    9633,  -2260, -11363,  -6436,   6436,  11363,   2260,  -9633,
   -9633,   2260,  11363,   6436,  -6436, -11363,  -2260,   9633,
    8955,  -5461, -11086,   1136,  11529,   3363, -10217,  -7350,
    7350,  10217,  -3363, -11529,  -1136,  11086,   5461,  -8955,
    8192,  -8192,  -8192,   8192,   8192,  -8192,  -8192,   8192,
    8192,  -8192,  -8192,   8192,   8192,  -8192,  -8192,   8192,
    7350, -10217,  -3363,  11529,  -1136, -11086,   5461,   8955,
   -8955,  -5461,  11086,   1136, -11529,   3363,  10217,  -7350,
    6436, -11363,   2260,   9633,  -9633,  -2260,  11363,  -6436,
   -6436,  11363,  -2260,  -9633,   9633,   2260, -11363,   6436,
    5461, -11529,   7350,   3363, -11086,   8955,   1136, -10217,
   10217,  -1136,  -8955,  11086,  -3363,  -7350,  11529,  -5461,
    4433, -10703,  10703,  -4433,  -4433,  10703, -10703,   4433,
    4433, -10703,  10703,  -4433,  -4433,  10703, -10703,   4433,
    3363,  -8955,  11529, -10217,   5461,   1136,  -7350,  11086,
  -11086,   7350,  -1136,  -5461,  10217, -11529,   8955,  -3363,
    2260,  -6436,   9633, -11363,  11363,  -9633,   6436,  -2260,
   -2260,   6436,  -9633,  11363, -11363,   9633,  -6436,   2260,
    1136,  -3363,   5461,  -7350,   8955, -10217,  11086, -11529,
   11529, -11086,  10217,  -8955,   7350,  -5461,   3363,  -1136
};

static const int16_t adst_i16[256] = {
    1084,   2159,   3214,   4240,   5228,   6168,   7052,   7873,
    8622,   9293,   9880,  10377,  10781,  11087,  11292,  11395,
    3214,   6168,   8622,  10377,  11292,  11292,  10377,   8622,
    6168,   3214,      0,  -3214,  -6168,  -8622, -10377, -11292,
    5228,   9293,  11292,  10781,   7873,   3214,  -2159,  -7052,
  -10377, -11395,  -9880,  -6168,  -1084,   4240,   8622,  11087,
    7052,  11087,  10377,   5228,  -2159,  -8622, -11395,  -9293,
   -3214,   4240,   9880,  11292,   7873,   1084,  -6168, -10781,
    8622,  11292,   6168,  -3214, -10377, -10377,  -3214,   6168,
   11292,   8622,      0,  -8622, -11292,  -6168,   3214,  10377,
    9880,   9880,      0,  -9880,  -9880,      0,   9880,   9880,
       0,  -9880,  -9880,      0,   9880,   9880,      0,  -9880,
   10781,   7052,  -6168, -11087,  -1084,  10377,   7873,  -5228,
  -11292,  -2159,   9880,   8622,  -4240, -11395,  -3214,   9293,
   11292,   3214, -10377,  -6168,   8622,   8622,  -6168, -10377,
    3214,  11292,      0, -11292,  -3214,  10377,   6168,  -8622,
   11395,  -1084, -11292,   2159,  11087,  -3214, -10781,   4240,
   10377,  -5228,  -9880,   6168,   9293,  -7052,  -8622,   7873,
   11087,  -5228,  -8622,   9293,   4240, -11292,   1084,  10781,
   -6168,  -7873,   9880,   3214, -11395,   2159,  10377,  -7052,
   10377,  -8622,  -3214,  11292,  -6168,  -6168,  11292,  -3214,
   -8622,  10377,      0, -10377,   8622,   3214, -11292,   6168,
    9293, -10781,   3214,   7052, -11395,   6168,   4240, -11087,
    8622,   1084,  -9880,  10377,  -2159,  -7873,  11292,  -5228,
    7873, -11395,   8622,  -1084,  -7052,  11292,  -9293,   2159,
    6168, -11087,   9880,  -3214,  -5228,  10781, -10377,   4240,
    6168, -10377,  11292,  -8622,   3214,   3214,  -8622,  11292,
  -10377,   6168,      0,  -6168,  10377, -11292,   8622,  -3214,
    4240,  -7873,  10377, -11395,  10781,  -8622,   5228,  -1084,
   -3214,   7052,  -9880,  11292, -11087,   9293,  -6168,   2159,
    2159,  -4240,   6168,  -7873,   9293, -10377,  11087, -11395,
   11292, -10781,   9880,  -8622,   7052,  -5228,   3214,  -1084
};

static const int xC1S7 = 16069;
static const int xC2S6 = 15137;
static const int xC3S5 = 13623;
static const int xC4S4 = 11585;
static const int xC5S3 =  9102;
static const int xC6S2 =  6270;
static const int xC7S1 =  3196;

#define SHIFT_BITS 14
#define DOROUND(X) X += (1<<(SHIFT_BITS-1));

#define FINAL_SHIFT 3
#define FINAL_ROUNDING (1<<(FINAL_SHIFT -1))
#define IN_SHIFT (FINAL_SHIFT+1)


void vp9_short_fdct8x8_c(short *InputData, short *OutputData, int pitch) {
  int loop;
  int short_pitch = pitch >> 1;
  int is07, is12, is34, is56;
  int is0734, is1256;
  int id07, id12, id34, id56;
  int irot_input_x, irot_input_y;
  int icommon_product1;      // Re-used product  (c4s4 * (s12 - s56))
  int icommon_product2;      // Re-used product  (c4s4 * (d12 + d56))
  int temp1, temp2;          // intermediate variable for computation

  int  InterData[64];
  int  *ip = InterData;
  short *op = OutputData;

  for (loop = 0; loop < 8; loop++) {
    // Pre calculate some common sums and differences.
    is07 = (InputData[0] + InputData[7]) << IN_SHIFT;
    is12 = (InputData[1] + InputData[2]) << IN_SHIFT;
    is34 = (InputData[3] + InputData[4]) << IN_SHIFT;
    is56 = (InputData[5] + InputData[6]) << IN_SHIFT;
    id07 = (InputData[0] - InputData[7]) << IN_SHIFT;
    id12 = (InputData[1] - InputData[2]) << IN_SHIFT;
    id34 = (InputData[3] - InputData[4]) << IN_SHIFT;
    id56 = (InputData[5] - InputData[6]) << IN_SHIFT;

    is0734 = is07 + is34;
    is1256 = is12 + is56;

    // Pre-Calculate some common product terms.
    icommon_product1 = xC4S4 * (is12 - is56);
    DOROUND(icommon_product1)
    icommon_product1 >>= SHIFT_BITS;

    icommon_product2 = xC4S4 * (id12 + id56);
    DOROUND(icommon_product2)
    icommon_product2 >>= SHIFT_BITS;


    ip[0] = (xC4S4 * (is0734 + is1256));
    DOROUND(ip[0]);
    ip[0] >>= SHIFT_BITS;

    ip[4] = (xC4S4 * (is0734 - is1256));
    DOROUND(ip[4]);
    ip[4] >>= SHIFT_BITS;

    // Define inputs to rotation for outputs 2 and 6
    irot_input_x = id12 - id56;
    irot_input_y = is07 - is34;

    // Apply rotation for outputs 2 and 6.
    temp1 = xC6S2 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[2] = temp1 + temp2;

    temp1 = xC6S2 * irot_input_y;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_x;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[6] = temp1 - temp2;

    // Define inputs to rotation for outputs 1 and 7
    irot_input_x = icommon_product1 + id07;
    irot_input_y = -(id34 + icommon_product2);

    // Apply rotation for outputs 1 and 7.
    temp1 = xC1S7 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC7S1 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[1] = temp1 - temp2;

    temp1 = xC7S1 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC1S7 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[7] = temp1 + temp2;

    // Define inputs to rotation for outputs 3 and 5
    irot_input_x = id07 - icommon_product1;
    irot_input_y = id34 - icommon_product2;

    // Apply rotation for outputs 3 and 5.
    temp1 = xC3S5 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC5S3 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[3] = temp1 - temp2;


    temp1 = xC5S3 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC3S5 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    ip[5] = temp1 + temp2;

    // Increment data pointer for next row
    InputData += short_pitch;
    ip += 8;
  }

  // Performed DCT on rows, now transform the columns
  ip = InterData;
  for (loop = 0; loop < 8; loop++) {
    // Pre calculate some common sums and differences.
    is07 = ip[0 * 8] + ip[7 * 8];
    is12 = ip[1 * 8] + ip[2 * 8];
    is34 = ip[3 * 8] + ip[4 * 8];
    is56 = ip[5 * 8] + ip[6 * 8];

    id07 = ip[0 * 8] - ip[7 * 8];
    id12 = ip[1 * 8] - ip[2 * 8];
    id34 = ip[3 * 8] - ip[4 * 8];
    id56 = ip[5 * 8] - ip[6 * 8];

    is0734 = is07 + is34;
    is1256 = is12 + is56;

    // Pre-Calculate some common product terms
    icommon_product1 = xC4S4 * (is12 - is56);
    icommon_product2 = xC4S4 * (id12 + id56);
    DOROUND(icommon_product1)
    DOROUND(icommon_product2)
    icommon_product1 >>= SHIFT_BITS;
    icommon_product2 >>= SHIFT_BITS;


    temp1 = xC4S4 * (is0734 + is1256);
    temp2 = xC4S4 * (is0734 - is1256);
    DOROUND(temp1);
    DOROUND(temp2);
    temp1 >>= SHIFT_BITS;

    temp2 >>= SHIFT_BITS;
    op[0 * 8] = (temp1 + FINAL_ROUNDING) >> FINAL_SHIFT;
    op[4 * 8] = (temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 2 and 6
    irot_input_x = id12 - id56;
    irot_input_y = is07 - is34;

    // Apply rotation for outputs 2 and 6.
    temp1 = xC6S2 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[2 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    temp1 = xC6S2 * irot_input_y;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC2S6 * irot_input_x;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[6 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 1 and 7
    irot_input_x = icommon_product1 + id07;
    irot_input_y = -(id34 + icommon_product2);

    // Apply rotation for outputs 1 and 7.
    temp1 = xC1S7 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC7S1 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[1 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    temp1 = xC7S1 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC1S7 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[7 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Define inputs to rotation for outputs 3 and 5
    irot_input_x = id07 - icommon_product1;
    irot_input_y = id34 - icommon_product2;

    // Apply rotation for outputs 3 and 5.
    temp1 = xC3S5 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC5S3 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[3 * 8] = (temp1 - temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;


    temp1 = xC5S3 * irot_input_x;
    DOROUND(temp1);
    temp1 >>= SHIFT_BITS;
    temp2 = xC3S5 * irot_input_y;
    DOROUND(temp2);
    temp2 >>= SHIFT_BITS;
    op[5 * 8] = (temp1 + temp2 + FINAL_ROUNDING) >> FINAL_SHIFT;

    // Increment data pointer for next column.
    ip++;
    op++;
  }
}

void vp9_short_fhaar2x2_c(short *input, short *output, int pitch) {
  /* [1 1; 1 -1] orthogonal transform */
  /* use position: 0,1, 4, 8 */
  int i;
  short *ip1 = input;
  short *op1 = output;
  for (i = 0; i < 16; i++) {
    op1[i] = 0;
  }

  op1[0] = (ip1[0] + ip1[1] + ip1[4] + ip1[8] + 1) >> 1;
  op1[1] = (ip1[0] - ip1[1] + ip1[4] - ip1[8]) >> 1;
  op1[4] = (ip1[0] + ip1[1] - ip1[4] - ip1[8]) >> 1;
  op1[8] = (ip1[0] - ip1[1] - ip1[4] + ip1[8]) >> 1;
}

/* For test */
#define TEST_INT 1
#if TEST_INT
#define vp9_fht_int_c vp9_fht_c
#else
#define vp9_fht_float_c vp9_fht_c
#endif

void vp9_fht_float_c(const int16_t *input, int pitch, int16_t *output,
               TX_TYPE tx_type, int tx_dim) {
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
  {
    int i, j, k;
    float bufa[256], bufb[256];  // buffers are for floating-point test purpose
                                 // the implementation could be simplified in
                                 // conjunction with integer transform
    const int16_t *ip = input;
    int16_t *op = output;

    float *pfa = &bufa[0];
    float *pfb = &bufb[0];

    // pointers to vertical and horizontal transforms
    const float *ptv, *pth;

    assert(tx_type != DCT_DCT);
    // load and convert residual array into floating-point
    for (j = 0; j < tx_dim; j++) {
      for (i = 0; i < tx_dim; i++) {
        pfa[i] = (float)ip[i];
      }
      pfa += tx_dim;
      ip  += pitch / 2;
    }

    // vertical transformation
    pfa = &bufa[0];
    pfb = &bufb[0];

    switch (tx_type) {
      case ADST_ADST :
      case ADST_DCT  :
        ptv = (tx_dim == 4) ? &adst_4[0] :
                              ((tx_dim == 8) ? &adst_8[0] : &adst_16[0]);
        break;

      default :
        ptv = (tx_dim == 4) ? &dct_4[0] :
                              ((tx_dim == 8) ? &dct_8[0] : &dct_16[0]);
        break;
    }

    for (j = 0; j < tx_dim; j++) {
      for (i = 0; i < tx_dim; i++) {
        pfb[i] = 0;
        for (k = 0; k < tx_dim; k++) {
          pfb[i] += ptv[k] * pfa[(k * tx_dim)];
        }
        pfa += 1;
      }
      pfb += tx_dim;
      ptv += tx_dim;
      pfa = &bufa[0];
    }

    // horizontal transformation
    pfa = &bufa[0];
    pfb = &bufb[0];

    switch (tx_type) {
      case ADST_ADST :
      case  DCT_ADST :
        pth = (tx_dim == 4) ? &adst_4[0] :
                              ((tx_dim == 8) ? &adst_8[0] : &adst_16[0]);
        break;

      default :
        pth = (tx_dim == 4) ? &dct_4[0] :
                              ((tx_dim == 8) ? &dct_8[0] : &dct_16[0]);
        break;
    }

    for (j = 0; j < tx_dim; j++) {
      for (i = 0; i < tx_dim; i++) {
        pfa[i] = 0;
        for (k = 0; k < tx_dim; k++) {
          pfa[i] += pfb[k] * pth[k];
        }
        pth += tx_dim;
      }

      pfa += tx_dim;
      pfb += tx_dim;
      // pth -= tx_dim * tx_dim;

      switch (tx_type) {
        case ADST_ADST :
        case  DCT_ADST :
          pth = (tx_dim == 4) ? &adst_4[0] :
                                ((tx_dim == 8) ? &adst_8[0] : &adst_16[0]);
          break;

        default :
          pth = (tx_dim == 4) ? &dct_4[0] :
                                ((tx_dim == 8) ? &dct_8[0] : &dct_16[0]);
          break;
      }
    }

    // convert to short integer format and load BLOCKD buffer
    op = output;
    pfa = &bufa[0];

    for (j = 0; j < tx_dim; j++) {
      for (i = 0; i < tx_dim; i++) {
        op[i] = (pfa[i] > 0 ) ? (int16_t)( 8 * pfa[i] + 0.49) :
                                     -(int16_t)(- 8 * pfa[i] + 0.49);
      }
      op  += tx_dim;
      pfa += tx_dim;
    }
  }
  vp9_clear_system_state();  // Make it simd safe : __asm emms;
}

/* Converted the transforms to integer form. */
#define VERTICAL_SHIFT 11
#define VERTICAL_ROUNDING ((1 << (VERTICAL_SHIFT - 1)) - 1)
#define HORIZONTAL_SHIFT 16
#define HORIZONTAL_ROUNDING ((1 << (HORIZONTAL_SHIFT - 1)) - 1)
void vp9_fht_int_c(const int16_t *input, int pitch, int16_t *output,
                   TX_TYPE tx_type, int tx_dim) {
  int i, j, k;
  int16_t imbuf[256];

  const int16_t *ip = input;
  int16_t *op = output;
  int16_t *im = &imbuf[0];

  /* pointers to vertical and horizontal transforms. */
  const int16_t *ptv = NULL, *pth = NULL;

  switch (tx_type) {
    case ADST_ADST :
      ptv = pth = (tx_dim == 4) ? &adst_i4[0]
                                  : ((tx_dim == 8) ? &adst_i8[0]
                                                     : &adst_i16[0]);
      break;
    case ADST_DCT  :
      ptv = (tx_dim == 4) ? &adst_i4[0]
                            : ((tx_dim == 8) ? &adst_i8[0] : &adst_i16[0]);
      pth = (tx_dim == 4) ? &dct_i4[0]
                            : ((tx_dim == 8) ? &dct_i8[0] : &dct_i16[0]);
      break;
    case  DCT_ADST :
      ptv = (tx_dim == 4) ? &dct_i4[0]
                            : ((tx_dim == 8) ? &dct_i8[0] : &dct_i16[0]);
      pth = (tx_dim == 4) ? &adst_i4[0]
                            : ((tx_dim == 8) ? &adst_i8[0] : &adst_i16[0]);
      break;
    case  DCT_DCT :
      ptv = pth = (tx_dim == 4) ? &dct_i4[0]
                                  : ((tx_dim == 8) ? &dct_i8[0] : &dct_i16[0]);
      break;
    default:
      assert(0);
      break;
  }

  /* vertical transformation */
  for (j = 0; j < tx_dim; j++) {
    for (i = 0; i < tx_dim; i++) {
      int temp = 0;

      for (k = 0; k < tx_dim; k++) {
        temp += ptv[k] * ip[(k * (pitch >> 1))];
      }

      im[i] = (int16_t)((temp + VERTICAL_ROUNDING) >> VERTICAL_SHIFT);
      ip++;
    }
    im += tx_dim;  // 16
    ptv += tx_dim;
    ip = input;
  }

  /* horizontal transformation */
  im = &imbuf[0];

  for (j = 0; j < tx_dim; j++) {
    const int16_t *pthc = pth;

    for (i = 0; i < tx_dim; i++) {
      int temp = 0;

      for (k = 0; k < tx_dim; k++) {
        temp += im[k] * pthc[k];
      }

      op[i] = (int16_t)((temp + HORIZONTAL_ROUNDING) >> HORIZONTAL_SHIFT);
      pthc += tx_dim;
    }

    im += tx_dim;  // 16
    op += tx_dim;
  }
}

void vp9_short_fdct4x4_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;

  for (i = 0; i < 4; i++) {
    a1 = ((ip[0] + ip[3]) << 5);
    b1 = ((ip[1] + ip[2]) << 5);
    c1 = ((ip[1] - ip[2]) << 5);
    d1 = ((ip[0] - ip[3]) << 5);

    op[0] = a1 + b1;
    op[2] = a1 - b1;

    op[1] = (c1 * 2217 + d1 * 5352 +  14500) >> 12;
    op[3] = (d1 * 2217 - c1 * 5352 +   7500) >> 12;

    ip += pitch / 2;
    op += 4;

  }
  ip = output;
  op = output;
  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[12];
    b1 = ip[4] + ip[8];
    c1 = ip[4] - ip[8];
    d1 = ip[0] - ip[12];

    op[0]  = (a1 + b1 + 7) >> 4;
    op[8]  = (a1 - b1 + 7) >> 4;

    op[4]  = ((c1 * 2217 + d1 * 5352 +  12000) >> 16) + (d1 != 0);
    op[12] = (d1 * 2217 - c1 * 5352 +  51000) >> 16;

    ip++;
    op++;
  }
}

void vp9_short_fdct8x4_c(short *input, short *output, int pitch)
{
    vp9_short_fdct4x4_c(input,   output,    pitch);
    vp9_short_fdct4x4_c(input + 4, output + 16, pitch);
}

void vp9_short_walsh4x4_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0 * pitch_short] + ip[3 * pitch_short];
    b1 = ip[1 * pitch_short] + ip[2 * pitch_short];
    c1 = ip[1 * pitch_short] - ip[2 * pitch_short];
    d1 = ip[0 * pitch_short] - ip[3 * pitch_short];

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = (a1 + b1 + 1) >> 1;
    op[1] = (c1 + d1) >> 1;
    op[2] = (a1 - b1) >> 1;
    op[3] = (d1 - c1) >> 1;

    ip += 4;
    op += 4;
  }
}

#if CONFIG_LOSSLESS
void vp9_short_walsh4x4_lossless_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = (ip[0 * pitch_short] + ip[3 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    b1 = (ip[1 * pitch_short] + ip[2 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    c1 = (ip[1 * pitch_short] - ip[2 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;
    d1 = (ip[0 * pitch_short] - ip[3 * pitch_short]) >> Y2_WHT_UPSCALE_FACTOR;

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = ((a1 + b1 + 1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[1] = ((c1 + d1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[2] = ((a1 - b1) >> 1) << Y2_WHT_UPSCALE_FACTOR;
    op[3] = ((d1 - c1) >> 1) << Y2_WHT_UPSCALE_FACTOR;

    ip += 4;
    op += 4;
  }
}

void vp9_short_walsh4x4_x8_c(short *input, short *output, int pitch) {
  int i;
  int a1, b1, c1, d1;
  short *ip = input;
  short *op = output;
  int pitch_short = pitch >> 1;

  for (i = 0; i < 4; i++) {
    a1 = ip[0 * pitch_short] + ip[3 * pitch_short];
    b1 = ip[1 * pitch_short] + ip[2 * pitch_short];
    c1 = ip[1 * pitch_short] - ip[2 * pitch_short];
    d1 = ip[0 * pitch_short] - ip[3 * pitch_short];

    op[0] = (a1 + b1 + 1) >> 1;
    op[4] = (c1 + d1) >> 1;
    op[8] = (a1 - b1) >> 1;
    op[12] = (d1 - c1) >> 1;

    ip++;
    op++;
  }
  ip = output;
  op = output;

  for (i = 0; i < 4; i++) {
    a1 = ip[0] + ip[3];
    b1 = ip[1] + ip[2];
    c1 = ip[1] - ip[2];
    d1 = ip[0] - ip[3];

    op[0] = ((a1 + b1 + 1) >> 1) << WHT_UPSCALE_FACTOR;
    op[1] = ((c1 + d1) >> 1) << WHT_UPSCALE_FACTOR;
    op[2] = ((a1 - b1) >> 1) << WHT_UPSCALE_FACTOR;
    op[3] = ((d1 - c1) >> 1) << WHT_UPSCALE_FACTOR;

    ip += 4;
    op += 4;
  }
}

void vp9_short_walsh8x4_x8_c(short *input, short *output, int pitch) {
  vp9_short_walsh4x4_x8_c(input,   output,    pitch);
  vp9_short_walsh4x4_x8_c(input + 4, output + 16, pitch);
}
#endif

#define TEST_INT_16x16_DCT 1
#if !TEST_INT_16x16_DCT
static const double C1 = 0.995184726672197;
static const double C2 = 0.98078528040323;
static const double C3 = 0.956940335732209;
static const double C4 = 0.923879532511287;
static const double C5 = 0.881921264348355;
static const double C6 = 0.831469612302545;
static const double C7 = 0.773010453362737;
static const double C8 = 0.707106781186548;
static const double C9 = 0.634393284163646;
static const double C10 = 0.555570233019602;
static const double C11 = 0.471396736825998;
static const double C12 = 0.38268343236509;
static const double C13 = 0.290284677254462;
static const double C14 = 0.195090322016128;
static const double C15 = 0.098017140329561;

static void dct16x16_1d(double input[16], double output[16]) {
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
  {
    double step[16];
    double intermediate[16];
    double temp1, temp2;

    // step 1
    step[ 0] = input[0] + input[15];
    step[ 1] = input[1] + input[14];
    step[ 2] = input[2] + input[13];
    step[ 3] = input[3] + input[12];
    step[ 4] = input[4] + input[11];
    step[ 5] = input[5] + input[10];
    step[ 6] = input[6] + input[ 9];
    step[ 7] = input[7] + input[ 8];
    step[ 8] = input[7] - input[ 8];
    step[ 9] = input[6] - input[ 9];
    step[10] = input[5] - input[10];
    step[11] = input[4] - input[11];
    step[12] = input[3] - input[12];
    step[13] = input[2] - input[13];
    step[14] = input[1] - input[14];
    step[15] = input[0] - input[15];

    // step 2
    output[0] = step[0] + step[7];
    output[1] = step[1] + step[6];
    output[2] = step[2] + step[5];
    output[3] = step[3] + step[4];
    output[4] = step[3] - step[4];
    output[5] = step[2] - step[5];
    output[6] = step[1] - step[6];
    output[7] = step[0] - step[7];

    temp1 = step[ 8]*C7;
    temp2 = step[15]*C9;
    output[ 8] = temp1 + temp2;

    temp1 = step[ 9]*C11;
    temp2 = step[14]*C5;
    output[ 9] = temp1 - temp2;

    temp1 = step[10]*C3;
    temp2 = step[13]*C13;
    output[10] = temp1 + temp2;

    temp1 = step[11]*C15;
    temp2 = step[12]*C1;
    output[11] = temp1 - temp2;

    temp1 = step[11]*C1;
    temp2 = step[12]*C15;
    output[12] = temp2 + temp1;

    temp1 = step[10]*C13;
    temp2 = step[13]*C3;
    output[13] = temp2 - temp1;

    temp1 = step[ 9]*C5;
    temp2 = step[14]*C11;
    output[14] = temp2 + temp1;

    temp1 = step[ 8]*C9;
    temp2 = step[15]*C7;
    output[15] = temp2 - temp1;

    // step 3
    step[ 0] = output[0] + output[3];
    step[ 1] = output[1] + output[2];
    step[ 2] = output[1] - output[2];
    step[ 3] = output[0] - output[3];

    temp1 = output[4]*C14;
    temp2 = output[7]*C2;
    step[ 4] = temp1 + temp2;

    temp1 = output[5]*C10;
    temp2 = output[6]*C6;
    step[ 5] = temp1 + temp2;

    temp1 = output[5]*C6;
    temp2 = output[6]*C10;
    step[ 6] = temp2 - temp1;

    temp1 = output[4]*C2;
    temp2 = output[7]*C14;
    step[ 7] = temp2 - temp1;

    step[ 8] = output[ 8] + output[11];
    step[ 9] = output[ 9] + output[10];
    step[10] = output[ 9] - output[10];
    step[11] = output[ 8] - output[11];

    step[12] = output[12] + output[15];
    step[13] = output[13] + output[14];
    step[14] = output[13] - output[14];
    step[15] = output[12] - output[15];

    // step 4
    output[ 0] = (step[ 0] + step[ 1]);
    output[ 8] = (step[ 0] - step[ 1]);

    temp1 = step[2]*C12;
    temp2 = step[3]*C4;
    temp1 = temp1 + temp2;
    output[ 4] = 2*(temp1*C8);

    temp1 = step[2]*C4;
    temp2 = step[3]*C12;
    temp1 = temp2 - temp1;
    output[12] = 2*(temp1*C8);

    output[ 2] = 2*((step[4] + step[ 5])*C8);
    output[14] = 2*((step[7] - step[ 6])*C8);

    temp1 = step[4] - step[5];
    temp2 = step[6] + step[7];
    output[ 6] = (temp1 + temp2);
    output[10] = (temp1 - temp2);

    intermediate[8] = step[8] + step[14];
    intermediate[9] = step[9] + step[15];

    temp1 = intermediate[8]*C12;
    temp2 = intermediate[9]*C4;
    temp1 = temp1 - temp2;
    output[3] = 2*(temp1*C8);

    temp1 = intermediate[8]*C4;
    temp2 = intermediate[9]*C12;
    temp1 = temp2 + temp1;
    output[13] = 2*(temp1*C8);

    output[ 9] = 2*((step[10] + step[11])*C8);

    intermediate[11] = step[10] - step[11];
    intermediate[12] = step[12] + step[13];
    intermediate[13] = step[12] - step[13];
    intermediate[14] = step[ 8] - step[14];
    intermediate[15] = step[ 9] - step[15];

    output[15] = (intermediate[11] + intermediate[12]);
    output[ 1] = -(intermediate[11] - intermediate[12]);

    output[ 7] = 2*(intermediate[13]*C8);

    temp1 = intermediate[14]*C12;
    temp2 = intermediate[15]*C4;
    temp1 = temp1 - temp2;
    output[11] = -2*(temp1*C8);

    temp1 = intermediate[14]*C4;
    temp2 = intermediate[15]*C12;
    temp1 = temp2 + temp1;
    output[ 5] = 2*(temp1*C8);
  }
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}

void vp9_short_fdct16x16_c(short *input, short *out, int pitch) {
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
  {
    int shortpitch = pitch >> 1;
    int i, j;
    double output[256];
    // First transform columns
    for (i = 0; i < 16; i++) {
        double temp_in[16], temp_out[16];
        for (j = 0; j < 16; j++)
            temp_in[j] = input[j*shortpitch + i];
        dct16x16_1d(temp_in, temp_out);
        for (j = 0; j < 16; j++)
            output[j*16 + i] = temp_out[j];
    }
    // Then transform rows
    for (i = 0; i < 16; ++i) {
        double temp_in[16], temp_out[16];
        for (j = 0; j < 16; ++j)
            temp_in[j] = output[j + i*16];
        dct16x16_1d(temp_in, temp_out);
        for (j = 0; j < 16; ++j)
            output[j + i*16] = temp_out[j];
    }
    // Scale by some magic number
    for (i = 0; i < 256; i++)
        out[i] = (short)round(output[i]/2);
  }
  vp9_clear_system_state(); // Make it simd safe : __asm emms;
}

#else
static const int16_t C1 = 16305;
static const int16_t C2 = 16069;
static const int16_t C3 = 15679;
static const int16_t C4 = 15137;
static const int16_t C5 = 14449;
static const int16_t C6 = 13623;
static const int16_t C7 = 12665;
static const int16_t C8 = 11585;
static const int16_t C9 = 10394;
static const int16_t C10 = 9102;
static const int16_t C11 = 7723;
static const int16_t C12 = 6270;
static const int16_t C13 = 4756;
static const int16_t C14 = 3196;
static const int16_t C15 = 1606;

#define RIGHT_SHIFT 14
#define ROUNDING (1 << (RIGHT_SHIFT - 1))

static void dct16x16_1d(int16_t input[16], int16_t output[16],
                        int last_shift_bits) {
    int16_t step[16];
    int intermediate[16];
    int temp1, temp2;
    int final_shift = RIGHT_SHIFT;
    int final_rounding = ROUNDING;
    int output_shift = 0;
    int output_rounding = 0;

    final_shift += last_shift_bits;
    if (final_shift > 0)
    final_rounding = 1 << (final_shift - 1);

    output_shift += last_shift_bits;
    if (output_shift > 0)
      output_rounding = 1 << (output_shift - 1);

    // step 1
    step[ 0] = input[0] + input[15];
    step[ 1] = input[1] + input[14];
    step[ 2] = input[2] + input[13];
    step[ 3] = input[3] + input[12];
    step[ 4] = input[4] + input[11];
    step[ 5] = input[5] + input[10];
    step[ 6] = input[6] + input[ 9];
    step[ 7] = input[7] + input[ 8];
    step[ 8] = input[7] - input[ 8];
    step[ 9] = input[6] - input[ 9];
    step[10] = input[5] - input[10];
    step[11] = input[4] - input[11];
    step[12] = input[3] - input[12];
    step[13] = input[2] - input[13];
    step[14] = input[1] - input[14];
    step[15] = input[0] - input[15];

    // step 2
    output[0] = step[0] + step[7];
    output[1] = step[1] + step[6];
    output[2] = step[2] + step[5];
    output[3] = step[3] + step[4];
    output[4] = step[3] - step[4];
    output[5] = step[2] - step[5];
    output[6] = step[1] - step[6];
    output[7] = step[0] - step[7];

    temp1 = step[ 8] * C7;
    temp2 = step[15] * C9;
    output[ 8] = (temp1 + temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 9] * C11;
    temp2 = step[14] * C5;
    output[ 9] = (temp1 - temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[10] * C3;
    temp2 = step[13] * C13;
    output[10] = (temp1 + temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[11] * C15;
    temp2 = step[12] * C1;
    output[11] = (temp1 - temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[11] * C1;
    temp2 = step[12] * C15;
    output[12] = (temp2 + temp1 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[10] * C13;
    temp2 = step[13] * C3;
    output[13] = (temp2 - temp1 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 9] * C5;
    temp2 = step[14] * C11;
    output[14] = (temp2 + temp1 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[ 8] * C9;
    temp2 = step[15] * C7;
    output[15] = (temp2 - temp1 + ROUNDING) >> RIGHT_SHIFT;

    // step 3
    step[ 0] = output[0] + output[3];
    step[ 1] = output[1] + output[2];
    step[ 2] = output[1] - output[2];
    step[ 3] = output[0] - output[3];

    temp1 = output[4] * C14;
    temp2 = output[7] * C2;
    step[ 4] = (temp1 + temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[5] * C10;
    temp2 = output[6] * C6;
    step[ 5] = (temp1 + temp2 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[5] * C6;
    temp2 = output[6] * C10;
    step[ 6] = (temp2 - temp1 + ROUNDING) >> RIGHT_SHIFT;

    temp1 = output[4] * C2;
    temp2 = output[7] * C14;
    step[ 7] = (temp2 - temp1 + ROUNDING) >> RIGHT_SHIFT;

    step[ 8] = output[ 8] + output[11];
    step[ 9] = output[ 9] + output[10];
    step[10] = output[ 9] - output[10];
    step[11] = output[ 8] - output[11];

    step[12] = output[12] + output[15];
    step[13] = output[13] + output[14];
    step[14] = output[13] - output[14];
    step[15] = output[12] - output[15];

    // step 4
    output[ 0] = (step[ 0] + step[ 1] + output_rounding) >> output_shift;
    output[ 8] = (step[ 0] - step[ 1] + output_rounding) >> output_shift;

    temp1 = step[2] * C12;
    temp2 = step[3] * C4;
    temp1 = (temp1 + temp2 + final_rounding) >> final_shift;
    output[ 4] = (2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;

    temp1 = step[2] * C4;
    temp2 = step[3] * C12;
    temp1 = (temp2 - temp1 + final_rounding) >> final_shift;
    output[12] = (2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;

    output[ 2] = (2 * ((step[4] + step[ 5]) * C8) + final_rounding)
        >> final_shift;
    output[14] = (2 * ((step[7] - step[ 6]) * C8) + final_rounding)
        >> final_shift;

    temp1 = step[4] - step[5];
    temp2 = step[6] + step[7];
    output[ 6] = (temp1 + temp2 + output_rounding) >> output_shift;
    output[10] = (temp1 - temp2 + output_rounding) >> output_shift;

    intermediate[8] = step[8] + step[14];
    intermediate[9] = step[9] + step[15];

    temp1 = intermediate[8] * C12;
    temp2 = intermediate[9] * C4;
    temp1 = (temp1 - temp2 + final_rounding) >> final_shift;
    output[3] = (2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;

    temp1 = intermediate[8] * C4;
    temp2 = intermediate[9] * C12;
    temp1 = (temp2 + temp1 + final_rounding) >> final_shift;
    output[13] = (2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;

    output[ 9] = (2 * ((step[10] + step[11]) * C8) + final_rounding)
        >> final_shift;

    intermediate[11] = step[10] - step[11];
    intermediate[12] = step[12] + step[13];
    intermediate[13] = step[12] - step[13];
    intermediate[14] = step[ 8] - step[14];
    intermediate[15] = step[ 9] - step[15];

    output[15] = (intermediate[11] + intermediate[12] + output_rounding)
        >> output_shift;
    output[ 1] = -(intermediate[11] - intermediate[12] + output_rounding)
        >> output_shift;

    output[ 7] = (2 * (intermediate[13] * C8) + final_rounding) >> final_shift;

    temp1 = intermediate[14] * C12;
    temp2 = intermediate[15] * C4;
    temp1 = (temp1 - temp2 + final_rounding) >> final_shift;
    output[11] = (-2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;

    temp1 = intermediate[14] * C4;
    temp2 = intermediate[15] * C12;
    temp1 = (temp2 + temp1 + final_rounding) >> final_shift;
    output[ 5] = (2 * (temp1 * C8) + ROUNDING) >> RIGHT_SHIFT;
}

void vp9_short_fdct16x16_c(int16_t *input, int16_t *out, int pitch) {
    int shortpitch = pitch >> 1;
    int i, j;
    int16_t output[256];
    int16_t *outptr = &output[0];

    // First transform columns
    for (i = 0; i < 16; i++) {
        int16_t temp_in[16];
        int16_t temp_out[16];
        for (j = 0; j < 16; j++)
            temp_in[j] = input[j * shortpitch + i];
        dct16x16_1d(temp_in, temp_out, 0);
        for (j = 0; j < 16; j++)
            output[j * 16 + i] = temp_out[j];
    }

    // Then transform rows
    for (i = 0; i < 16; ++i) {
        dct16x16_1d(outptr, out, 1);
        outptr += 16;
        out += 16;
    }
}
#undef RIGHT_SHIFT
#undef ROUNDING
#endif
