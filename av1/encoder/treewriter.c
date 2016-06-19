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

#include "av1/encoder/treewriter.h"

static void tree2tok(struct av1_token *tokens, const aom_tree_index *tree,
                     int i, int v, int l) {
  v += v;
  ++l;

  do {
    const aom_tree_index j = tree[i++];
    if (j <= 0) {
      tokens[-j].value = v;
      tokens[-j].len = l;
    } else {
      tree2tok(tokens, tree, j, v, l);
    }
  } while (++v & 1);
}

void av1_tokens_from_tree(struct av1_token *tokens,
                          const aom_tree_index *tree) {
  tree2tok(tokens, tree, 0, 0, 0);
}

/* This code assumes that tree contains as unique leaf nodes the integer values
    0 to len - 1 and produces the forward and inverse mapping tables in ind[]
    and inv[] respectively. */
void av1_indices_from_tree(int *ind, int *inv, int len,
                           const aom_tree_index *tree) {
  int i;
  int index;
  for (i = index = 0; i < TREE_SIZE(len); i++) {
    const aom_tree_index j = tree[i];
    if (j <= 0) {
      inv[index] = -j;
      ind[-j] = index++;
    }
  }
}

static unsigned int convert_distribution(unsigned int i, aom_tree tree,
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

void av1_tree_probs_from_distribution(aom_tree tree,
                                      unsigned int branch_ct[/* n-1 */][2],
                                      const unsigned int num_events[/* n */]) {
  convert_distribution(0, tree, branch_ct, num_events);
}
