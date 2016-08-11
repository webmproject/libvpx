/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_util/debug_util.h"
#include <assert.h>
#include <stdio.h>
#if CONFIG_BITSTREAM_DEBUG
#define QUEUE_MAX_SIZE 2000000
static int result_queue[QUEUE_MAX_SIZE];
static int prob_queue[QUEUE_MAX_SIZE];
static int queue_r = 0;
static int queue_w = 0;
static int queue_prev_w = -1;
static int skip_r = 0;
static int skip_w = 0;

void bitstream_queue_set_skip_write(int skip) { skip_w = skip; }

void bitstream_queue_set_skip_read(int skip) { skip_r = skip; }

void bitstream_queue_record_write(void) { queue_prev_w = queue_w; }

void bitstream_queue_reset_write(void) { queue_w = queue_prev_w; }

int bitstream_queue_get_write(void) { return queue_w; }

int bitstream_queue_get_read(void) { return queue_r; }

void bitstream_queue_pop(int* result, int* prob) {
  if (!skip_r) {
    if (queue_w == queue_r) {
      printf("buffer underflow queue_w %d queue_r %d\n", queue_w, queue_r);
      assert(0);
    }
    *result = result_queue[queue_r];
    *prob = prob_queue[queue_r];
    queue_r = (queue_r + 1) % QUEUE_MAX_SIZE;
  }
}

void bitstream_queue_push(int result, int prob) {
  if (!skip_w) {
    result_queue[queue_w] = result;
    prob_queue[queue_w] = prob;
    queue_w = (queue_w + 1) % QUEUE_MAX_SIZE;
    if (queue_w == queue_r) {
      printf("buffer overflow queue_w %d queue_r %d\n", queue_w, queue_r);
      assert(0);
    }
  }
}
#endif  // CONFIG_BITSTREAM_DEBUG
