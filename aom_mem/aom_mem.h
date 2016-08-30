/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef AOM_MEM_AOM_MEM_H_
#define AOM_MEM_AOM_MEM_H_

#include "aom_config.h"
#if defined(__uClinux__)
#include <lddk.h>
#endif

#include <stdlib.h>
#include <stddef.h>

#if defined(__cplusplus)
extern "C" {
#endif

void *aom_memalign(size_t align, size_t size);
void *aom_malloc(size_t size);
void *aom_calloc(size_t num, size_t size);
void *aom_realloc(void *memblk, size_t size);
void aom_free(void *memblk);

#if CONFIG_AOM_HIGHBITDEPTH
void *aom_memset16(void *dest, int val, size_t length);
#endif

#include <string.h>

#ifdef AOM_MEM_PLTFRM
#include AOM_MEM_PLTFRM
#endif

#if defined(__cplusplus)
}
#endif

#endif  // AOM_MEM_AOM_MEM_H_
