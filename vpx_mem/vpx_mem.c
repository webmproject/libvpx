/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vpx_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/vpx_mem_intrnl.h"
#include "vpx/vpx_integer.h"

static INLINE size_t *get_malloc_address_location(void *const mem) {
  return ((size_t *)mem) - 1;
}

static INLINE size_t get_aligned_malloc_size(size_t size, size_t align) {
  return size + align - 1 + ADDRESS_STORAGE_SIZE;
}

static INLINE void set_actual_malloc_address(void *const mem,
                                             const void *const malloc_addr) {
  size_t *const malloc_addr_location = get_malloc_address_location(mem);
  *malloc_addr_location = (size_t)malloc_addr;
}

static INLINE void *get_actual_malloc_address(void *const mem) {
  size_t *const malloc_addr_location = get_malloc_address_location(mem);
  return (void *)(*malloc_addr_location);
}

void *vpx_memalign(size_t align, size_t size) {
  void *x = NULL;
  const size_t aligned_size = get_aligned_malloc_size(size, align);
  void *const addr = malloc(aligned_size);
  if (addr) {
    x = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE, (int)align);
    set_actual_malloc_address(x, addr);
  }
  return x;
}

void *vpx_malloc(size_t size) { return vpx_memalign(DEFAULT_ALIGNMENT, size); }

void *vpx_calloc(size_t num, size_t size) {
  const size_t total_size = num * size;
  void *const x = vpx_malloc(total_size);
  if (x) memset(x, 0, total_size);
  return x;
}

void *vpx_realloc(void *memblk, size_t size) {
  void *new_addr = NULL;

  /*
  The realloc() function changes the size of the object pointed to by
  ptr to the size specified by size, and returns a pointer to the
  possibly moved block. The contents are unchanged up to the lesser
  of the new and old sizes. If ptr is null, realloc() behaves like
  malloc() for the specified size. If size is zero (0) and ptr is
  not a null pointer, the object pointed to is freed.
  */
  if (!memblk)
    new_addr = vpx_malloc(size);
  else if (!size)
    vpx_free(memblk);
  else {
    void *addr = get_actual_malloc_address(memblk);
    const size_t aligned_size =
        get_aligned_malloc_size(size, DEFAULT_ALIGNMENT);
    memblk = NULL;
    addr = realloc(addr, aligned_size);
    if (addr) {
      new_addr = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE,
                            DEFAULT_ALIGNMENT);
      set_actual_malloc_address(new_addr, addr);
    }
  }

  return new_addr;
}

void vpx_free(void *memblk) {
  if (memblk) {
    void *addr = get_actual_malloc_address(memblk);
    free(addr);
  }
}

#if CONFIG_VP9_HIGHBITDEPTH
void *vpx_memset16(void *dest, int val, size_t length) {
  size_t i;
  uint16_t *dest16 = (uint16_t *)dest;
  for (i = 0; i < length; i++) *dest16++ = val;
  return dest;
}
#endif  // CONFIG_VP9_HIGHBITDEPTH
