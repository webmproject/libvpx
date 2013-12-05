/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <assert.h>

#include "./vpx_config.h"
#include "vp9/common/vp9_treecoder.h"

static unsigned int convert_distribution(unsigned int i, vp9_tree tree,
                                         unsigned int branch_ct[][2],
                                         const unsigned int num_events[]) {
  unsigned int left, right;

  if (tree[i] <= 0)
    left = num_events[-tree[i]];
  else
    left = convert_distribution(tree[i], tree, branch_ct, num_events);

  if (tree[i + 1] <= 0)
    right = num_events[-tree[i + 1]];
  else
    right = convert_distribution(tree[i + 1], tree, branch_ct, num_events);

  branch_ct[i >> 1][0] = left;
  branch_ct[i >> 1][1] = right;
  return left + right;
}

void vp9_tree_probs_from_distribution(vp9_tree tree,
                                      unsigned int branch_ct[/* n-1 */][2],
                                      const unsigned int num_events[/* n */]) {
  convert_distribution(0, tree, branch_ct, num_events);
}


