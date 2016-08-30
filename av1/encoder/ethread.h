/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_ETHREAD_H_
#define AV1_ENCODER_ETHREAD_H_

#ifdef __cplusplus
extern "C" {
#endif

struct AV1_COMP;
struct ThreadData;

typedef struct EncWorkerData {
  struct AV1_COMP *cpi;
  struct ThreadData *td;
  int start;
} EncWorkerData;

void av1_encode_tiles_mt(struct AV1_COMP *cpi);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_ETHREAD_H_
