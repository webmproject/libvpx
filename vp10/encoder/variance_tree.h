/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_VARIANCE_TREE_H_
#define VP10_ENCODER_VARIANCE_TREE_H_

#include <assert.h>

#include "./vpx_config.h"

#include "vpx/vpx_integer.h"

#include "vp10/common/enums.h"

#ifdef __cplusplus
extern "C" {
#endif

struct VP10Common;
struct ThreadData;

typedef struct {
  int64_t sum_square_error;
  int64_t sum_error;
  int log2_count;
  int variance;
} var;

typedef struct {
  var none;
  var horz[2];
  var vert[2];
} partition_variance;

typedef struct VAR_TREE {
  int force_split;
  partition_variance variances;
  struct VAR_TREE *split[4];
  BLOCK_SIZE bsize;
  const uint8_t *src;
  const uint8_t *ref;
  int src_stride;
  int ref_stride;
  int width;
  int height;
#if CONFIG_VPX_HIGHBITDEPTH
  int highbd;
#endif  // CONFIG_VPX_HIGHBITDEPTH
} VAR_TREE;

void vp10_setup_var_tree(struct VP10Common *cm, struct ThreadData *td);
void vp10_free_var_tree(struct ThreadData *td);

// Set variance values given sum square error, sum error, count.
static INLINE void fill_variance(int64_t s2, int64_t s, int c, var *v) {
  v->sum_square_error = s2;
  v->sum_error = s;
  v->log2_count = c;
  v->variance = (int)(256 * (v->sum_square_error -
      ((v->sum_error * v->sum_error) >> v->log2_count)) >> v->log2_count);
}

static INLINE void sum_2_variances(const var *a, const var *b, var *r) {
  assert(a->log2_count == b->log2_count);
  fill_variance(a->sum_square_error + b->sum_square_error,
                a->sum_error + b->sum_error, a->log2_count + 1, r);
}

static INLINE void fill_variance_node(VAR_TREE *vt) {
  sum_2_variances(&vt->split[0]->variances.none,
                  &vt->split[1]->variances.none,
                  &vt->variances.horz[0]);
  sum_2_variances(&vt->split[2]->variances.none,
                  &vt->split[3]->variances.none,
                  &vt->variances.horz[1]);
  sum_2_variances(&vt->split[0]->variances.none,
                  &vt->split[2]->variances.none,
                  &vt->variances.vert[0]);
  sum_2_variances(&vt->split[1]->variances.none,
                  &vt->split[3]->variances.none,
                  &vt->variances.vert[1]);
  sum_2_variances(&vt->variances.vert[0],
                  &vt->variances.vert[1],
                  &vt->variances.none);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif /* VP10_ENCODER_VARIANCE_TREE_H_ */
