/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AV1_ENCODER_SUBEXP_H_
#define AV1_ENCODER_SUBEXP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "aom_dsp/prob.h"

struct aom_writer;

void av1_write_prob_diff_update(struct aom_writer *w, aom_prob newp,
                                aom_prob oldp);

void av1_cond_prob_diff_update(struct aom_writer *w, aom_prob *oldp,
                               const unsigned int ct[2]);

int av1_prob_diff_update_savings_search(const unsigned int *ct, aom_prob oldp,
                                        aom_prob *bestp, aom_prob upd);

int av1_prob_diff_update_savings_search_model(const unsigned int *ct,
                                              const aom_prob *oldp,
                                              aom_prob *bestp, aom_prob upd,
                                              int stepsize);
int av1_cond_prob_diff_update_savings(aom_prob *oldp, const unsigned int ct[2]);

#if CONFIG_ENTROPY
int av1_prob_update_search_subframe(unsigned int ct[][2], aom_prob oldp,
                                    aom_prob *bestp, aom_prob upd, int n);
int av1_prob_update_search_model_subframe(
    unsigned int ct[ENTROPY_NODES][COEF_PROBS_BUFS][2], const aom_prob *oldp,
    aom_prob *bestp, aom_prob upd, int stepsize, int n);
#endif  // CONFIG_ENTROPY

//
// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
void aom_write_primitive_symmetric(aom_writer *w, int word,
                                   unsigned int mag_bits);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_SUBEXP_H_
