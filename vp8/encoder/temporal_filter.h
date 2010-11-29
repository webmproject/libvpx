/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef __INC_VP8_TEMPORAL_FILTER_H
#define __INC_VP8_TEMPORAL_FILTER_H

#define prototype_filter(sym)\
    void (sym) \
    ( \
     unsigned char *frame1, \
     unsigned int stride, \
     unsigned char *frame2, \
     unsigned int block_size, \
     int strength, \
     int filter_weight, \
     unsigned int *accumulator, \
     unsigned int *count \
    )

#ifndef vp8_temporal_filter
#define vp8_temporal_filter vp8_apply_temporal_filter_c
#endif
extern prototype_filter(vp8_temporal_filter);

typedef struct
{
    prototype_filter(*filter);
} vp8_temporal_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define TEMPORAL_INVOKE(ctx,fn) (ctx)->fn
#else
#define TEMPORAL_INVOKE(ctx,fn) vp8_temporal_##fn
#endif

#endif // __INC_VP8_TEMPORAL_FILTER_H
