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

#include "aom_dsp/buf_ans.h"
#include "aom_mem/aom_mem.h"
#include "aom/internal/aom_codec_internal.h"

void aom_buf_ans_alloc(struct BufAnsCoder *c,
                       struct aom_internal_error_info *error, int size_hint) {
  c->error = error;
  c->size = size_hint;
  AOM_CHECK_MEM_ERROR(error, c->buf, aom_malloc(c->size * sizeof(*c->buf)));
  // Initialize to overfull to trigger the assert in write.
  c->offset = c->size + 1;
}

void aom_buf_ans_free(struct BufAnsCoder *c) {
  aom_free(c->buf);
  c->buf = NULL;
  c->size = 0;
}

void aom_buf_ans_grow(struct BufAnsCoder *c) {
  struct buffered_ans_symbol *new_buf = NULL;
  int new_size = c->size * 2;
  AOM_CHECK_MEM_ERROR(c->error, new_buf,
                      aom_malloc(new_size * sizeof(*new_buf)));
  memcpy(new_buf, c->buf, c->size * sizeof(*c->buf));
  aom_free(c->buf);
  c->buf = new_buf;
  c->size = new_size;
}
