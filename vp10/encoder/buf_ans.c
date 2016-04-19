/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <string.h>

#include "vp10/common/common.h"
#include "vp10/encoder/buf_ans.h"
#include "vp10/encoder/encoder.h"
#include "vpx_mem/vpx_mem.h"

void vp10_buf_ans_alloc(struct BufAnsCoder *c, struct VP10Common *cm,
                        int size_hint) {
  c->cm = cm;
  c->size = size_hint;
  CHECK_MEM_ERROR(cm, c->buf, vpx_malloc(c->size * sizeof(*c->buf)));
  // Initialize to overfull to trigger the assert in write.
  c->offset = c->size + 1;
}

void vp10_buf_ans_free(struct BufAnsCoder *c) {
  vpx_free(c->buf);
  c->buf = NULL;
  c->size = 0;
}

void vp10_buf_ans_grow(struct BufAnsCoder *c) {
  struct buffered_ans_symbol *new_buf = NULL;
  int new_size = c->size * 2;
  CHECK_MEM_ERROR(c->cm, new_buf, vpx_malloc(new_size * sizeof(*new_buf)));
  memcpy(new_buf, c->buf, c->size * sizeof(*c->buf));
  vpx_free(c->buf);
  c->buf = new_buf;
  c->size = new_size;
}
