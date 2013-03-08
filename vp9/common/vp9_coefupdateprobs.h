/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_VP9_COEFUPDATEPROBS_H_
#define VP9_COMMON_VP9_COEFUPDATEPROBS_H_

/* Update probabilities for the nodes in the token entropy tree.
   Generated file included by vp9_entropy.c */
#define COEF_UPDATE_PROB 252
#define COEF_UPDATE_PROB_8X8 252
#define COEF_UPDATE_PROB_16X16 252

#if CONFIG_CODE_NONZEROCOUNT
#define NZC_UPDATE_PROB_4X4     252
#define NZC_UPDATE_PROB_8X8     252
#define NZC_UPDATE_PROB_16X16   252
#define NZC_UPDATE_PROB_32X32   252
#define NZC_UPDATE_PROB_PCAT    252
#endif

#endif  // VP9_COMMON_VP9_COEFUPDATEPROBS_H__
