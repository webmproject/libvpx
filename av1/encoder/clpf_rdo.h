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

#ifndef AV1_ENCODER_CLPF_H_
#define AV1_ENCODER_CLPF_H_

#include "av1/common/reconinter.h"

int av1_clpf_decision(int k, int l, const YV12_BUFFER_CONFIG *rec,
                      const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                      int block_size, int w, int h, unsigned int strength,
                      unsigned int fb_size_log2, uint8_t *res);

void av1_clpf_test_frame(const YV12_BUFFER_CONFIG *rec,
                         const YV12_BUFFER_CONFIG *org, const AV1_COMMON *cm,
                         int *best_strength, int *best_bs);

#endif
