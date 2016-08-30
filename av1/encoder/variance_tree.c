/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "av1/encoder/variance_tree.h"
#include "av1/encoder/encoder.h"

void av1_setup_var_tree(struct AV1Common *cm, ThreadData *td) {
  int i, j;
#if CONFIG_EXT_PARTITION
  const int leaf_nodes = 1024;
  const int tree_nodes = 1024 + 256 + 64 + 16 + 4 + 1;
#else
  const int leaf_nodes = 256;
  const int tree_nodes = 256 + 64 + 16 + 4 + 1;
#endif  // CONFIG_EXT_PARTITION
  int index = 0;
  VAR_TREE *this_var;
  int nodes;

  aom_free(td->var_tree);
  CHECK_MEM_ERROR(cm, td->var_tree,
                  aom_calloc(tree_nodes, sizeof(*td->var_tree)));

  this_var = &td->var_tree[0];

  // Sets up all the leaf nodes in the tree.
  for (index = 0; index < leaf_nodes; ++index) {
    VAR_TREE *const leaf = &td->var_tree[index];
    leaf->split[0] = NULL;
  }

  // Each node has 4 leaf nodes, fill in the child pointers
  // from leafs to the root.
  for (nodes = leaf_nodes >> 2; nodes > 0; nodes >>= 2) {
    for (i = 0; i < nodes; ++i, ++index) {
      VAR_TREE *const node = &td->var_tree[index];
      for (j = 0; j < 4; j++) node->split[j] = this_var++;
    }
  }

  // Set up the root node for the largest superblock size
  i = MAX_MIB_SIZE_LOG2 - MIN_MIB_SIZE_LOG2;
  td->var_root[i] = &td->var_tree[tree_nodes - 1];
  // Set up the root nodes for the rest of the possible superblock sizes
  while (--i >= 0) {
    td->var_root[i] = td->var_root[i + 1]->split[0];
  }
}

void av1_free_var_tree(ThreadData *td) {
  aom_free(td->var_tree);
  td->var_tree = NULL;
}
