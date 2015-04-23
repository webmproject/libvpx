/*
 *  Copyright (c) 2013 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_ENCODER_VP9_SUBEXP_H_
#define VP9_ENCODER_VP9_SUBEXP_H_

#ifdef __cplusplus
extern "C" {
#endif

void vp9_write_prob_diff_update(vp9_writer *w,
                                vp9_prob newp, vp9_prob oldp);

void vp9_cond_prob_diff_update(vp9_writer *w, vp9_prob *oldp,
                               unsigned int *ct);

int vp9_cond_prob_diff_update_savings(vp9_prob *oldp,
                                      const unsigned int ct[2]);

int vp9_prob_diff_update_savings_search(const unsigned int *ct,
                                        vp9_prob oldp, vp9_prob *bestp,
                                        vp9_prob upd);

int vp9_prob_diff_update_savings_search_model(const unsigned int *ct,
                                              const vp9_prob *oldp,
                                              vp9_prob *bestp,
                                              vp9_prob upd);

// num_values is the number of values word can take
void vp9_write_primitive_uniform(vp9_writer *w, int word,
                                 unsigned int num_values);

// k is the parameter of the subexponential code
void vp9_write_primitive_subexp(vp9_writer *w, int word,
                                unsigned int k);
//
// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
void vp9_write_primitive_symmetric(vp9_writer *w, int word,
                                   unsigned int mag_bits);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP9_ENCODER_VP9_SUBEXP_H_
