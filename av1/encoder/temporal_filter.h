/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_TEMPORAL_FILTER_H_
#define AV1_ENCODER_TEMPORAL_FILTER_H_

#ifdef __cplusplus
extern "C" {
#endif

void av1_temporal_filter(AV1_COMP *cpi, int distance);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_TEMPORAL_FILTER_H_
