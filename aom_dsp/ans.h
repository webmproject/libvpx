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

#ifndef AOM_DSP_ANS_H_
#define AOM_DSP_ANS_H_
// Constants, types and utilities for Asymmetric Numeral Systems
// http://arxiv.org/abs/1311.2540v2

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/prob.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

typedef uint8_t AnsP8;
#define ANS_P8_PRECISION 256u
#define ANS_P8_SHIFT 8
typedef uint16_t AnsP10;
#define ANS_P10_PRECISION 1024u
#define RANS_PROB_BITS 10

#define RANS_PRECISION ANS_P10_PRECISION

#define L_BASE (ANS_P10_PRECISION * 4)  // L_BASE % precision must be 0
#define IO_BASE 256
// Range I = { L_BASE, L_BASE + 1, ..., L_BASE * IO_BASE - 1 }

// This is now just a boring cdf. It starts with an explicit zero.
// TODO(aconverse): Remove starting zero.
typedef uint16_t rans_lut[16];

void aom_rans_build_cdf_from_pdf(const AnsP10 token_probs[], rans_lut cdf_tab);

void aom_rans_merge_prob8_pdf(AnsP10 *const out_pdf, const AnsP8 node_prob,
                              const AnsP10 *const src_pdf, int in_syms);
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
#endif  // AOM_DSP_ANS_H_
