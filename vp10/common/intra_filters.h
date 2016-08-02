/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_COMMON_INTRA_FILTERS_H_
#define VP10_COMMON_INTRA_FILTERS_H_

#define FILTER_INTRA_PREC_BITS (10)

static int filter_intra_taps_4[TX_SIZES][INTRA_MODES][4] = {
    {
        {735, 881, -537, -54},
        {1005, 519, -488, -11},
        {383, 990, -343, -6},
        {442, 805, -542, 319},
        {658, 616, -133, -116},
        {875, 442, -141, -151},
        {386, 741, -23, -80},
        {390, 1027, -446, 51},
        {679, 606, -523, 262},
        {903, 922, -778, -23},
    },
    {
        {648, 803, -444, 16},
        {972, 620, -576, 7},
        {561, 967, -499, -5},
        {585, 762, -468, 144},
        {596, 619, -182, -9},
        {895, 459, -176, -153},
        {557, 722, -126, -129},
        {601, 839, -523, 105},
        {562, 709, -499, 251},
        {803, 872, -695, 43},
    },
    {
        {423, 728, -347, 111},
        {963, 685, -665, 23},
        {281, 1024, -480, 216},
        {640, 596, -437, 78},
        {429, 669, -259, 99},
        {740, 646, -415, 23},
        {568, 771, -346, 40},
        {404, 833, -486, 209},
        {398, 712, -423, 307},
        {939, 935, -887, 17},
    },
    {
        {477, 737, -393, 150},
        {881, 630, -546, 67},
        {506, 984, -443, -20},
        {114, 459, -270, 528},
        {433, 528, 14, 3},
        {837, 470, -301, -30},
        {181, 777, 89, -107},
        {-29, 716, -232, 259},
        {589, 646, -495, 255},
        {740, 884, -728, 77},
    },
};

#endif  // VP10_COMMON_INTRA_FILTERS_H_
