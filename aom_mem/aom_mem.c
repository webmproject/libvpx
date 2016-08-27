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

#define __AOM_MEM_C__

#include "aom_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/aom_mem_intrnl.h"
#include "aom/aom_integer.h"

static size_t GetAlignedMallocSize(size_t size, size_t align) {
  return size + align - 1 + ADDRESS_STORAGE_SIZE;
}

static size_t *GetMallocAddressLocation(void *const mem) {
  return ((size_t *)mem) - 1;
}

static void SetActualMallocAddress(void *const mem,
                                   const void *const malloc_addr) {
  size_t *const malloc_addr_location = GetMallocAddressLocation(mem);
  *malloc_addr_location = (size_t)malloc_addr;
}

static void *GetActualMallocAddress(void *const mem) {
  const size_t *const malloc_addr_location = GetMallocAddressLocation(mem);
  return (void *)(*malloc_addr_location);
}

void *aom_memalign(size_t align, size_t size) {
  void *x = NULL;
  const size_t aligned_size = GetAlignedMallocSize(size, align);
  void *const addr = malloc(aligned_size);
  if (addr) {
    x = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE, align);
    SetActualMallocAddress(x, addr);
  }
  return x;
}

void *aom_malloc(size_t size) { return aom_memalign(DEFAULT_ALIGNMENT, size); }

void *aom_calloc(size_t num, size_t size) {
  const size_t total_size = num * size;
  void *const x = aom_malloc(total_size);
  if (x) memset(x, 0, total_size);
  return x;
}

void *aom_realloc(void *memblk, size_t size) {
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
    new_addr = aom_malloc(size);
  else if (!size)
    aom_free(memblk);
  else {
    void *addr = GetActualMallocAddress(memblk);
    const size_t aligned_size = GetAlignedMallocSize(size, DEFAULT_ALIGNMENT);
    memblk = NULL;
    addr = realloc(addr, aligned_size);
    if (addr) {
      new_addr = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE,
                            DEFAULT_ALIGNMENT);
      SetActualMallocAddress(new_addr, addr);
    }
  }

  return new_addr;
}

void aom_free(void *memblk) {
  if (memblk) {
    void *addr = GetActualMallocAddress(memblk);
    free(addr);
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
void *aom_memset16(void *dest, int val, size_t length) {
  size_t i;
  uint16_t *dest16 = (uint16_t *)dest;
  for (i = 0; i < length; i++) *dest16++ = val;
  return dest;
}
#endif  // CONFIG_AOM_HIGHBITDEPTH
