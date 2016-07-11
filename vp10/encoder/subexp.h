/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP10_ENCODER_SUBEXP_H_
#define VP10_ENCODER_SUBEXP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "vpx_dsp/prob.h"

struct vp10_writer;

void vp10_write_prob_diff_update(struct vp10_writer *w,
                                vpx_prob newp, vpx_prob oldp);

void vp10_cond_prob_diff_update(struct vp10_writer *w, vpx_prob *oldp,
                               const unsigned int ct[2]);

int vp10_prob_diff_update_savings_search(const unsigned int *ct,
                                        vpx_prob oldp, vpx_prob *bestp,
                                        vpx_prob upd);


int vp10_prob_diff_update_savings_search_model(const unsigned int *ct,
                                              const vpx_prob *oldp,
                                              vpx_prob *bestp,
                                              vpx_prob upd,
                                              int stepsize);
int vp10_cond_prob_diff_update_savings(vpx_prob *oldp,
                                       const unsigned int ct[2]);

#if CONFIG_ENTROPY
int vp10_prob_update_search_subframe(unsigned int ct[][2],
                                     vpx_prob oldp, vpx_prob *bestp,
                                     vpx_prob upd, int n);
int vp10_prob_update_search_model_subframe(unsigned int ct[ENTROPY_NODES]
                                                          [COEF_PROBS_BUFS][2],
                                           const vpx_prob *oldp,
                                           vpx_prob *bestp, vpx_prob upd,
                                           int stepsize, int n);
#endif  // CONFIG_ENTROPY

//
// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
void vp10_write_primitive_symmetric(vp10_writer *w, int word,
                                   unsigned int mag_bits);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_SUBEXP_H_
