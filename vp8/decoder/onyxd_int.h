/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_VP8D_INT_H
#define __INC_VP8D_INT_H
#include "vpx_ports/config.h"
#include "vp8/common/onyxd.h"
#include "treereader.h"
#include "vp8/common/onyxc_int.h"
#include "dequantize.h"

// #define DEC_DEBUG

typedef struct {
  int ithread;
  void *ptr1;
  void *ptr2;
} DECODETHREAD_DATA;

typedef struct {
  MACROBLOCKD  mbd;
  int mb_row;
  int current_mb_col;
  short *coef_ptr;
} MB_ROW_DEC;

typedef struct {
  int64_t time_stamp;
  int size;
} DATARATE;

typedef struct {
  int const *scan;
  int const *scan_8x8;
  UINT8 const *ptr_block2leftabove;
  vp8_tree_index const *vp8_coef_tree_ptr;
  unsigned char *norm_ptr;
  UINT8 *ptr_coef_bands_x;
  UINT8 *ptr_coef_bands_x_8x8;

  ENTROPY_CONTEXT_PLANES *A;
  ENTROPY_CONTEXT_PLANES *L;

  INT16 *qcoeff_start_ptr;
  BOOL_DECODER *current_bc;

  vp8_prob const *coef_probs[BLOCK_TYPES];
  vp8_prob const *coef_probs_8x8[BLOCK_TYPES_8X8];
#if CONFIG_TX16X16
  vp8_prob const *coef_probs_16X16[BLOCK_TYPES_16X16];
#endif

  UINT8 eob[25];

} DETOK;

typedef struct VP8Decompressor {
  DECLARE_ALIGNED(16, MACROBLOCKD, mb);

  DECLARE_ALIGNED(16, VP8_COMMON, common);

  vp8_reader bc, bc2;

  VP8D_CONFIG oxcf;


  const unsigned char *Source;
  unsigned int   source_sz;

  vp8_reader *mbc;
  int64_t last_time_stamp;
  int   ready_for_new_data;

  DATARATE dr[16];

  DETOK detoken;

#if CONFIG_RUNTIME_CPU_DETECT
  vp8_dequant_rtcd_vtable_t        dequant;
#endif

  vp8_prob prob_skip_false;

  int decoded_key_frame;

  int interleaved_decoding;

} VP8D_COMP;

int vp8_decode_frame(VP8D_COMP *cpi);
void vp8_dmachine_specific_config(VP8D_COMP *pbi);


#if CONFIG_DEBUG
#define CHECK_MEM_ERROR(lval,expr) do {\
    lval = (expr); \
    if(!lval) \
      vpx_internal_error(&pbi->common.error, VPX_CODEC_MEM_ERROR,\
                         "Failed to allocate "#lval" at %s:%d", \
                         __FILE__,__LINE__);\
  } while(0)
#else
#define CHECK_MEM_ERROR(lval,expr) do {\
    lval = (expr); \
    if(!lval) \
      vpx_internal_error(&pbi->common.error, VPX_CODEC_MEM_ERROR,\
                         "Failed to allocate "#lval);\
  } while(0)
#endif

#endif
