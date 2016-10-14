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

#include "./aom_config.h"

#if CONFIG_DAALA_EC
#include <string.h>
#endif

#include "aom_dsp/prob.h"

#if CONFIG_DAALA_EC
#include "aom_dsp/entcode.h"
#endif

const uint8_t aom_norm[256] = {
  0, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

static unsigned int tree_merge_probs_impl(unsigned int i,
                                          const aom_tree_index *tree,
                                          const aom_prob *pre_probs,
                                          const unsigned int *counts,
                                          aom_prob *probs) {
  const int l = tree[i];
  const unsigned int left_count =
      (l <= 0) ? counts[-l]
               : tree_merge_probs_impl(l, tree, pre_probs, counts, probs);
  const int r = tree[i + 1];
  const unsigned int right_count =
      (r <= 0) ? counts[-r]
               : tree_merge_probs_impl(r, tree, pre_probs, counts, probs);
  const unsigned int ct[2] = { left_count, right_count };
  probs[i >> 1] = mode_mv_merge_probs(pre_probs[i >> 1], ct);
  return left_count + right_count;
}

void aom_tree_merge_probs(const aom_tree_index *tree, const aom_prob *pre_probs,
                          const unsigned int *counts, aom_prob *probs) {
  tree_merge_probs_impl(0, tree, pre_probs, counts, probs);
}

#if CONFIG_DAALA_EC
typedef struct tree_node tree_node;

struct tree_node {
  aom_tree_index index;
  uint8_t probs[16];
  uint8_t prob;
  int path;
  int len;
  int l;
  int r;
  uint16_t pdf;
};

/* Compute the probability of this node in Q23 */
static uint32_t tree_node_prob(tree_node n, int i) {
  uint32_t prob;
  /* 1.0 in Q23 */
  prob = 16777216;
  for (; i < n.len; i++) {
    prob = prob * n.probs[i] >> 8;
  }
  return prob;
}

static int tree_node_cmp(tree_node a, tree_node b) {
  int i;
  uint32_t pa;
  uint32_t pb;
  for (i = 0; i < OD_MINI(a.len, b.len) && a.probs[i] == b.probs[i]; i++) {
  }
  pa = tree_node_prob(a, i);
  pb = tree_node_prob(b, i);
  return pa > pb ? 1 : pa < pb ? -1 : 0;
}

/* Given a Q15 probability for symbol subtree rooted at tree[n], this function
    computes the probability of each symbol (defined as a node that has no
    children). */
static uint16_t tree_node_compute_probs(tree_node *tree, int n, uint16_t pdf) {
  if (tree[n].l == 0) {
    /* This prevents probability computations in Q15 that underflow from
        producing a symbol that has zero probability. */
    if (pdf == 0) pdf = 1;
    tree[n].pdf = pdf;
    return pdf;
  } else {
    /* We process the smaller probability first,  */
    if (tree[n].prob < 128) {
      uint16_t lp;
      uint16_t rp;
      lp = (((uint32_t)pdf) * tree[n].prob + 128) >> 8;
      lp = tree_node_compute_probs(tree, tree[n].l, lp);
      rp = tree_node_compute_probs(tree, tree[n].r, lp > pdf ? 0 : pdf - lp);
      return lp + rp;
    } else {
      uint16_t rp;
      uint16_t lp;
      rp = (((uint32_t)pdf) * (256 - tree[n].prob) + 128) >> 8;
      rp = tree_node_compute_probs(tree, tree[n].r, rp);
      lp = tree_node_compute_probs(tree, tree[n].l, rp > pdf ? 0 : pdf - rp);
      return lp + rp;
    }
  }
}

static int tree_node_extract(tree_node *tree, int n, int symb, uint16_t *pdf,
                             aom_tree_index *index, int *path, int *len) {
  if (tree[n].l == 0) {
    pdf[symb] = tree[n].pdf;
    if (index != NULL) index[symb] = tree[n].index;
    if (path != NULL) path[symb] = tree[n].path;
    if (len != NULL) len[symb] = tree[n].len;
    return symb + 1;
  } else {
    symb = tree_node_extract(tree, tree[n].l, symb, pdf, index, path, len);
    return tree_node_extract(tree, tree[n].r, symb, pdf, index, path, len);
  }
}

int tree_to_cdf(const aom_tree_index *tree, const aom_prob *probs,
                aom_tree_index root, uint16_t *cdf, aom_tree_index *index,
                int *path, int *len) {
  tree_node symb[2 * 16 - 1];
  int nodes;
  int next[16];
  int size;
  int nsymbs;
  int i;
  /* Create the root node with probability 1 in Q15. */
  symb[0].index = root;
  symb[0].path = 0;
  symb[0].len = 0;
  symb[0].l = symb[0].r = 0;
  nodes = 1;
  next[0] = 0;
  size = 1;
  nsymbs = 1;
  while (size > 0 && nsymbs < 16) {
    int m;
    tree_node n;
    aom_tree_index j;
    uint8_t prob;
    m = 0;
    /* Find the internal node with the largest probability. */
    for (i = 1; i < size; i++) {
      if (tree_node_cmp(symb[next[i]], symb[next[m]]) > 0) m = i;
    }
    i = next[m];
    memmove(&next[m], &next[m + 1], sizeof(*next) * (size - (m + 1)));
    size--;
    /* Split this symbol into two symbols */
    n = symb[i];
    j = n.index;
    prob = probs[j >> 1];
    /* Left */
    n.index = tree[j];
    n.path <<= 1;
    n.len++;
    n.probs[n.len - 1] = prob;
    symb[nodes] = n;
    if (n.index > 0) {
      next[size++] = nodes;
    }
    /* Right */
    n.index = tree[j + 1];
    n.path += 1;
    n.probs[n.len - 1] = 256 - prob;
    symb[nodes + 1] = n;
    if (n.index > 0) {
      next[size++] = nodes + 1;
    }
    symb[i].prob = prob;
    symb[i].l = nodes;
    symb[i].r = nodes + 1;
    nodes += 2;
    nsymbs++;
  }
  /* Compute the probabilities of each symbol in Q15 */
  tree_node_compute_probs(symb, 0, 32768);
  /* Extract the cdf, index, path and length */
  tree_node_extract(symb, 0, 0, cdf, index, path, len);
  /* Convert to CDF */
  for (i = 1; i < nsymbs; i++) {
    cdf[i] = cdf[i - 1] + cdf[i];
  }
  return nsymbs;
}
#endif
