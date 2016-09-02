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
#ifndef AV1_COMMON_CLPF_H_
#define AV1_COMMON_CLPF_H_

#include "av1/common/reconinter.h"

// Configuration
#define CLPF_ALLOW_PIXEL_PARALLELISM \
  1  // 1 = SIMD friendly (adds a buffer requirement)
#define CLPF_ALLOW_BLOCK_PARALLELISM \
  0  // 1 = MT friendly (degrades quality slighty)
#define CLPF_FILTER_ALL_PLANES \
  0  // 1 = filter both luma and chroma, 0 = filter only luma

void av1_clpf_frame(const YV12_BUFFER_CONFIG *frame, const AV1_COMMON *cm,
                    MACROBLOCKD *xd);

#endif
