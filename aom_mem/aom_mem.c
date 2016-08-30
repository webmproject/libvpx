/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#define __AOM_MEM_C__

#include "aom_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/aom_mem_intrnl.h"
#include "aom/aom_integer.h"

void *aom_memalign(size_t align, size_t size) {
  void *addr, *x = NULL;

  addr = malloc(size + align - 1 + ADDRESS_STORAGE_SIZE);

  if (addr) {
    x = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE, (int)align);
    /* save the actual malloc address */
    ((size_t *)x)[-1] = (size_t)addr;
  }

  return x;
}

void *aom_malloc(size_t size) { return aom_memalign(DEFAULT_ALIGNMENT, size); }

void *aom_calloc(size_t num, size_t size) {
  void *x;

  x = aom_memalign(DEFAULT_ALIGNMENT, num * size);

  if (x) memset(x, 0, num * size);

  return x;
}

void *aom_realloc(void *memblk, size_t size) {
  void *addr, *new_addr = NULL;
  int align = DEFAULT_ALIGNMENT;

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
    addr = (void *)(((size_t *)memblk)[-1]);
    memblk = NULL;

    new_addr = realloc(addr, size + align + ADDRESS_STORAGE_SIZE);

    if (new_addr) {
      addr = new_addr;
      new_addr =
          (void *)(((size_t)((unsigned char *)new_addr + ADDRESS_STORAGE_SIZE) +
                    (align - 1)) &
                   (size_t)-align);
      /* save the actual malloc address */
      ((size_t *)new_addr)[-1] = (size_t)addr;
    }
  }

  return new_addr;
}

void aom_free(void *memblk) {
  if (memblk) {
    void *addr = (void *)(((size_t *)memblk)[-1]);
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
