/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_MBGRAPH_H_
#define AV1_ENCODER_MBGRAPH_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  struct {
    int err;
    union {
      int_mv mv;
      PREDICTION_MODE mode;
    } m;
  } ref[TOTAL_REFS_PER_FRAME];
} MBGRAPH_MB_STATS;

typedef struct { MBGRAPH_MB_STATS *mb_stats; } MBGRAPH_FRAME_STATS;

struct AV1_COMP;

void av1_update_mbgraph_stats(struct AV1_COMP *cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_MBGRAPH_H_
