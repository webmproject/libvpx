/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VP9_COMMON_VP9_TREECODER_H_
#define VP9_COMMON_VP9_TREECODER_H_

typedef unsigned char vp9_prob;

#define vp9_prob_half ( (vp9_prob) 128)

typedef signed char vp9_tree_index;
struct bool_coder_spec;

typedef struct bool_coder_spec bool_coder_spec;
typedef struct bool_writer bool_writer;
typedef struct bool_reader bool_reader;

typedef const bool_coder_spec c_bool_coder_spec;
typedef const bool_writer c_bool_writer;
typedef const bool_reader c_bool_reader;



# define vp9_complement( x) (255 - x)


/* We build coding trees compactly in arrays.
   Each node of the tree is a pair of vp9_tree_indices.
   Array index often references a corresponding probability table.
   Index <= 0 means done encoding/decoding and value = -Index,
   Index > 0 means need another bit, specification at index.
   Nonnegative indices are always even;  processing begins at node 0. */

typedef const vp9_tree_index vp9_tree[], *vp9_tree_p;


typedef const struct vp9_token_struct {
  int value;
  int Len;
} vp9_token;

/* Construct encoding array from tree. */

void vp9_tokens_from_tree(struct vp9_token_struct *, vp9_tree);
void vp9_tokens_from_tree_offset(struct vp9_token_struct *, vp9_tree,
                                 int offset);


/* Convert array of token occurrence counts into a table of probabilities
   for the associated binary encoding tree.  Also writes count of branches
   taken for each node on the tree; this facilitiates decisions as to
   probability updates. */

void vp9_tree_probs_from_distribution(
  int n,                      /* n = size of alphabet */
  vp9_token tok               [ /* n */ ],
  vp9_tree tree,
  vp9_prob probs          [ /* n-1 */ ],
  unsigned int branch_ct       [ /* n-1 */ ] [2],
  const unsigned int num_events[ /* n */ ],
  unsigned int Pfactor,
  int Round
);

static __inline int clip_prob(int p) {
  if (p > 255)
    return 255;
  else if (p < 1)
    return 1;
  return p;
}

vp9_prob vp9_bin_prob_from_distribution(const unsigned int counts[2]);

#endif
