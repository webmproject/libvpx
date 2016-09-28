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

#ifndef AV1_COMMON_X86_HIGHBD_CONVOLVE_FILTER_SSE4_H_
#define AV1_COMMON_X86_HIGHBD_CONVOLVE_FILTER_SSE4_H_

#include "./aom_config.h"
#include "av1/common/filter.h"

#if USE_TEMPORALFILTER_12TAP
DECLARE_ALIGNED(16, extern const int16_t,
                av1_sub_pel_filters_temporalfilter_12_highbd_ver_signal_dir[15]
                                                                           [6]
                                                                           [8]);
#endif

#if CONFIG_EXT_INTERP
DECLARE_ALIGNED(16, extern const int16_t,
                av1_sub_pel_filters_12sharp_highbd_ver_signal_dir[15][6][8]);
DECLARE_ALIGNED(16, extern const int16_t,
                av1_sub_pel_filters_10sharp_highbd_ver_signal_dir[15][6][8]);
#endif

#endif  // AV1_COMMON_X86_HIGHBD_CONVOLVE_FILTER_SSE4_H_
