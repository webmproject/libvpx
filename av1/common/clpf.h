/*
Copyright (c) 2016, Cisco Systems
(Replace with proper AOM header)
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
