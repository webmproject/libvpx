/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "mvref_common.h"

#if CONFIG_NEW_MVREF

#define MVREF_NEIGHBOURS 8
static int mv_ref_search[MVREF_NEIGHBOURS][2] =
  { {0,-1},{-1,0},{-1,-1},{0,-2},{-2,0},{-1,-2},{-2,-1},{-2,-2} };
static int ref_distance_weight[MVREF_NEIGHBOURS] =
  { 3,3,2,1,1,1,1,1 };
  //{ 4,4,2,1,1,1,1,1 };

// clamp_mv
#define MV_BORDER (16 << 3) // Allow 16 pels in 1/8th pel units
static void clamp_mv(const MACROBLOCKD *xd, int_mv *mv) {

  if (mv->as_mv.col < (xd->mb_to_left_edge - MV_BORDER))
    mv->as_mv.col = xd->mb_to_left_edge - MV_BORDER;
  else if (mv->as_mv.col > xd->mb_to_right_edge + MV_BORDER)
    mv->as_mv.col = xd->mb_to_right_edge + MV_BORDER;

  if (mv->as_mv.row < (xd->mb_to_top_edge - MV_BORDER))
    mv->as_mv.row = xd->mb_to_top_edge - MV_BORDER;
  else if (mv->as_mv.row > xd->mb_to_bottom_edge + MV_BORDER)
    mv->as_mv.row = xd->mb_to_bottom_edge + MV_BORDER;
}

// Code for selecting / building and entropy coding a motion vector reference
// Returns a seperation value for two vectors.
// This is taken as the sum of the abs x and y difference.
unsigned int mv_distance(int_mv *mv1, int_mv *mv2) {
  return (abs(mv1->as_mv.row - mv2->as_mv.row) +
          abs(mv1->as_mv.col - mv2->as_mv.col));
}

// Gets a best matching candidate refenence motion vector
// from the given mode info structure (if available)
int get_candidate_mvref(
  const MODE_INFO *candidate_mi,
  MV_REFERENCE_FRAME ref_frame,
  MV_REFERENCE_FRAME *candidate_ref_frame,
  int_mv *candidate_mv
) {

  int ret_val = FALSE;

  if (ref_frame == candidate_mi->mbmi.ref_frame) {
    candidate_mv->as_int = candidate_mi->mbmi.mv[FIRST_REF].as_int;
    *candidate_ref_frame = ref_frame;
    ret_val = TRUE;

  } else if (ref_frame == candidate_mi->mbmi.second_ref_frame) {
    candidate_mv->as_int = candidate_mi->mbmi.mv[SECOND_REF].as_int;
    *candidate_ref_frame = ref_frame;
    ret_val = TRUE;

  } else if (candidate_mi->mbmi.ref_frame != INTRA_FRAME) {
    candidate_mv->as_int = candidate_mi->mbmi.mv[FIRST_REF].as_int;
    *candidate_ref_frame = candidate_mi->mbmi.ref_frame;
    ret_val = TRUE;

  } else if (candidate_mi->mbmi.second_ref_frame != INTRA_FRAME) {
    candidate_mv->as_int = candidate_mi->mbmi.mv[SECOND_REF].as_int;
    *candidate_ref_frame = candidate_mi->mbmi.second_ref_frame;
    ret_val = TRUE;
  }

  return ret_val;
}

// Performs mv adjustment based on reference frame and clamps the MV
// if it goes off the edge of the buffer.
void scale_mv(
  MACROBLOCKD *xd,
  MV_REFERENCE_FRAME this_ref_frame,
  MV_REFERENCE_FRAME candidate_ref_frame,
  int_mv *candidate_mv,
  int *ref_sign_bias
) {

  if (candidate_ref_frame != this_ref_frame) {

    //int frame_distances[MAX_REF_FRAMES];
    //int last_distance = 1;
    //int gf_distance = xd->frames_since_golden;
    //int arf_distance = xd->frames_till_alt_ref_frame;

    // Sign inversion where appropriate.
    if (ref_sign_bias[candidate_ref_frame] != ref_sign_bias[this_ref_frame]) {
      candidate_mv->as_mv.row = -candidate_mv->as_mv.row;
      candidate_mv->as_mv.col = -candidate_mv->as_mv.col;
    }

    // Scale based on frame distance if the reference frames not the same.
    /*frame_distances[INTRA_FRAME] = 1;   // should never be used
    frame_distances[LAST_FRAME] = 1;
    frame_distances[GOLDEN_FRAME] =
      (xd->frames_since_golden) ? xd->frames_since_golden : 1;
    frame_distances[ALTREF_FRAME] =
      (xd->frames_till_alt_ref_frame) ? xd->frames_till_alt_ref_frame : 1;

    if (frame_distances[this_ref_frame] &&
        frame_distances[candidate_ref_frame]) {
      candidate_mv->as_mv.row =
        (short)(((int)(candidate_mv->as_mv.row) *
                 frame_distances[this_ref_frame]) /
                frame_distances[candidate_ref_frame]);

      candidate_mv->as_mv.col =
        (short)(((int)(candidate_mv->as_mv.col) *
                 frame_distances[this_ref_frame]) /
                frame_distances[candidate_ref_frame]);
    }
    */
  }

  // Clamp the MV so it does not point out of the frame buffer
  clamp_mv(xd, candidate_mv);
}

// Adds a new candidate reference vector to the list if indeed it is new.
// If it is not new then the score of the existing candidate that it matches
// is increased and the list is resorted.
void addmv_and_shuffle(
  int_mv *mv_list,
  int *mv_scores,
  int *index,
  int_mv candidate_mv,
  int weight
) {

  int i = *index;
  int duplicate_found = FALSE;

  // Check for duplicates. If there is one increment its score.
  while (i > 0) {
    i--;
    if (candidate_mv.as_int == mv_list[i].as_int) {
      duplicate_found = TRUE;
      mv_scores[i] += weight;
      break;
    }
  }

  // If no duplicate was found add the new vector and give it a weight
  if (!duplicate_found) {
    mv_list[*index].as_int = candidate_mv.as_int;
    mv_scores[*index] = weight;
    i = *index;
    (*index)++;
  }

  // Reshuffle the list so that highest scoring mvs at the top.
  while (i > 0) {
    if (mv_scores[i] > mv_scores[i-1]) {
      int tmp_score = mv_scores[i-1];
      int_mv tmp_mv = mv_list[i-1];

      mv_scores[i-1] = mv_scores[i];
      mv_list[i-1] = mv_list[i];
      mv_scores[i] = tmp_score;
      mv_list[i] = tmp_mv;
      i--;
    } else
      break;
  }
}


// Measure the distance of each reference candidate from the actual
// residual vector and return the nearest
unsigned int pick_best_mv_ref( int_mv target_mv,
                               int_mv * mv_ref_list,
                               int_mv * best_ref ) {

  int i;
  int best_index = 0;
  unsigned int distance, distance2;

  distance = mv_distance(&target_mv, &mv_ref_list[0]);

  for (i = 1; i < MAX_MV_REFS; ++i ) {
    distance2 =
      mv_distance(&target_mv, &mv_ref_list[i]);
    if (distance2 < distance) {
      distance = distance2;
      best_index = i;
    }
  }

  (*best_ref).as_int = mv_ref_list[best_index].as_int;

  return best_index;
}

// This function searches the neighbourhood of a given MB/SB and populates a
// list of candidate reference vectors.
//
void find_mv_refs(
  MACROBLOCKD *xd,
  MODE_INFO *here,
  MODE_INFO *lf_here,
  MV_REFERENCE_FRAME ref_frame,
  int_mv *mv_ref_list,
  int *ref_sign_bias
) {

  int i;
  MODE_INFO *candidate_mi;
  int_mv candidate_mvs[MAX_MV_REFS];
  int_mv c_refmv;
  MV_REFERENCE_FRAME c_ref_frame;
  int candidate_scores[MAX_MV_REFS];
  int index = 0;
  int ref_weight = 0;
  int valid_mv_ref;

  // Blank the reference vector lists and other local structures.
  vpx_memset(mv_ref_list, 0, sizeof(int_mv) * MAX_MV_REFS);
  vpx_memset(candidate_mvs, 0, sizeof(int_mv) * MAX_MV_REFS);
  vpx_memset(candidate_scores, 0, sizeof(candidate_scores));

  // Populate a list with candidate reference vectors from the
  // spatial neighbours.
  for (i = 0; i < 2; ++i) {
    if (((mv_ref_search[i][0] << 7) >= xd->mb_to_left_edge) &&
        ((mv_ref_search[i][1] << 7) >= xd->mb_to_top_edge)) {

      candidate_mi = here + mv_ref_search[i][0] +
                     (mv_ref_search[i][1] * xd->mode_info_stride);

      valid_mv_ref = get_candidate_mvref(candidate_mi, ref_frame,
                                         &c_ref_frame, &c_refmv);

      if (valid_mv_ref) {
          scale_mv(xd, ref_frame, c_ref_frame, &c_refmv, ref_sign_bias );
          ref_weight = ref_distance_weight[i] +
                       ((c_ref_frame == ref_frame) << 3);

          addmv_and_shuffle(candidate_mvs, candidate_scores,
                            &index, c_refmv, ref_weight);
      }
    }
  }

  // Look at the corresponding vector in the last frame
  candidate_mi = lf_here;
  valid_mv_ref = get_candidate_mvref(candidate_mi, ref_frame,
                                     &c_ref_frame, &c_refmv);
  if (valid_mv_ref) {
    scale_mv(xd, ref_frame, c_ref_frame, &c_refmv, ref_sign_bias );
    ref_weight = 2 + ((c_ref_frame == ref_frame) << 3);
    addmv_and_shuffle(candidate_mvs, candidate_scores,
                      &index, c_refmv, ref_weight);
  }

  // Populate a list with candidate reference vectors from the
  // spatial neighbours.
  for (i = 2; i < MVREF_NEIGHBOURS; ++i) {
    if (((mv_ref_search[i][0] << 7) >= xd->mb_to_left_edge) &&
        ((mv_ref_search[i][1] << 7) >= xd->mb_to_top_edge)) {

      candidate_mi = here + mv_ref_search[i][0] +
                     (mv_ref_search[i][1] * xd->mode_info_stride);

      valid_mv_ref = get_candidate_mvref(candidate_mi, ref_frame,
                                         &c_ref_frame, &c_refmv);

      if (valid_mv_ref) {
          scale_mv(xd, ref_frame, c_ref_frame, &c_refmv, ref_sign_bias );
          ref_weight = ref_distance_weight[i] +
                       ((c_ref_frame == ref_frame) << 3);

          addmv_and_shuffle(candidate_mvs, candidate_scores,
                            &index, c_refmv, ref_weight);
      }
    }
  }

  // 0,0 is always a valid reference.
  for (i = 0; i < index; ++i)
    if (candidate_mvs[i].as_int == 0)
      break;
  if (i == index) {
    c_refmv.as_int = 0;
    addmv_and_shuffle(candidate_mvs, candidate_scores,
                      &index, c_refmv, 1);
  }

  // Copy over the candidate list.
  vpx_memcpy(mv_ref_list, candidate_mvs, sizeof(candidate_mvs));
}

#endif
