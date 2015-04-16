/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#define __VPX_MEM_C__

#include "vpx_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/vpx_mem_intrnl.h"
#include "vpx/vpx_integer.h"

#if CONFIG_MEM_TRACKER
#ifndef VPX_NO_GLOBALS
static unsigned long g_alloc_count = 0;
#else
#include "vpx_global_handling.h"
#define g_alloc_count vpxglobalm(vpxmem,g_alloc_count)
#endif
#endif

#if USE_GLOBAL_FUNCTION_POINTERS
struct GLOBAL_FUNC_POINTERS {
  g_malloc_func g_malloc;
  g_calloc_func g_calloc;
  g_realloc_func g_realloc;
  g_free_func g_free;
  g_memcpy_func g_memcpy;
  g_memset_func g_memset;
  g_memmove_func g_memmove;
} *g_func = NULL;

# define VPX_MALLOC_L  g_func->g_malloc
# define VPX_REALLOC_L g_func->g_realloc
# define VPX_FREE_L    g_func->g_free
# define VPX_MEMCPY_L  g_func->g_memcpy
# define VPX_MEMSET_L  g_func->g_memset
# define VPX_MEMMOVE_L g_func->g_memmove
#else
# define VPX_MALLOC_L  malloc
# define VPX_REALLOC_L realloc
# define VPX_FREE_L    free
# define VPX_MEMCPY_L  memcpy
# define VPX_MEMSET_L  memset
# define VPX_MEMMOVE_L memmove
#endif /* USE_GLOBAL_FUNCTION_POINTERS */

unsigned int vpx_mem_get_version() {
  unsigned int ver = ((unsigned int)(unsigned char)VPX_MEM_VERSION_CHIEF << 24 |
                      (unsigned int)(unsigned char)VPX_MEM_VERSION_MAJOR << 16 |
                      (unsigned int)(unsigned char)VPX_MEM_VERSION_MINOR << 8  |
                      (unsigned int)(unsigned char)VPX_MEM_VERSION_PATCH);
  return ver;
}

void *vpx_memalign(size_t align, size_t size) {
  void *addr,
       * x = NULL;

  addr = VPX_MALLOC_L(size + align - 1 + ADDRESS_STORAGE_SIZE);

  if (addr) {
    x = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE, (int)align);
    /* save the actual malloc address */
    ((size_t *)x)[-1] = (size_t)addr;
  }

  return x;
}

void *vpx_malloc(size_t size) {
  return vpx_memalign(DEFAULT_ALIGNMENT, size);
}

void *vpx_calloc(size_t num, size_t size) {
  void *x;

  x = vpx_memalign(DEFAULT_ALIGNMENT, num * size);

  if (x)
    VPX_MEMSET_L(x, 0, num * size);

  return x;
}

void *vpx_realloc(void *memblk, size_t size) {
  void *addr,
       * new_addr = NULL;
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
    new_addr = vpx_malloc(size);
  else if (!size)
    vpx_free(memblk);
  else {
    addr   = (void *)(((size_t *)memblk)[-1]);
    memblk = NULL;

    new_addr = VPX_REALLOC_L(addr, size + align + ADDRESS_STORAGE_SIZE);

    if (new_addr) {
      addr = new_addr;
      new_addr = (void *)(((size_t)
                           ((unsigned char *)new_addr + ADDRESS_STORAGE_SIZE) + (align - 1)) &
                          (size_t) - align);
      /* save the actual malloc address */
      ((size_t *)new_addr)[-1] = (size_t)addr;
    }
  }

  return new_addr;
}

void vpx_free(void *memblk) {
  if (memblk) {
    void *addr = (void *)(((size_t *)memblk)[-1]);
    VPX_FREE_L(addr);
  }
}

#if CONFIG_MEM_TRACKER
void *xvpx_memalign(size_t align, size_t size, char *file, int line) {
#if TRY_BOUNDS_CHECK
  unsigned char *x_bounds;
#endif

  void *x;

  if (g_alloc_count == 0) {
#if TRY_BOUNDS_CHECK
    int i_rv = vpx_memory_tracker_init(BOUNDS_CHECK_PAD_SIZE, BOUNDS_CHECK_VALUE);
#else
    int i_rv = vpx_memory_tracker_init(0, 0);
#endif

    if (i_rv < 0) {
      _P(printf("ERROR xvpx_malloc MEM_TRACK_USAGE error vpx_memory_tracker_init().\n");)
    }
  }

#if TRY_BOUNDS_CHECK
  {
    int i;
    unsigned int tempme = BOUNDS_CHECK_VALUE;

    x_bounds = vpx_memalign(align, size + (BOUNDS_CHECK_PAD_SIZE * 2));

    if (x_bounds) {
      /*we're aligning the address twice here but to keep things
        consistent we want to have the padding come before the stored
        address so no matter what free function gets called we will
        attempt to free the correct address*/
      x_bounds = (unsigned char *)(((size_t *)x_bounds)[-1]);
      x = align_addr(x_bounds + BOUNDS_CHECK_PAD_SIZE + ADDRESS_STORAGE_SIZE,
                     (int)align);
      /* save the actual malloc address */
      ((size_t *)x)[-1] = (size_t)x_bounds;

      for (i = 0; i < BOUNDS_CHECK_PAD_SIZE; i += sizeof(unsigned int)) {
        VPX_MEMCPY_L(x_bounds + i, &tempme, sizeof(unsigned int));
        VPX_MEMCPY_L((unsigned char *)x + size + i,
                     &tempme, sizeof(unsigned int));
      }
    } else
      x = NULL;
  }
#else
  x = vpx_memalign(align, size);
#endif /*TRY_BOUNDS_CHECK*/

  g_alloc_count++;

  vpx_memory_tracker_add((size_t)x, (unsigned int)size, file, line, 1);

  return x;
}

void *xvpx_malloc(size_t size, char *file, int line) {
  return xvpx_memalign(DEFAULT_ALIGNMENT, size, file, line);
}

void *xvpx_calloc(size_t num, size_t size, char *file, int line) {
  void *x = xvpx_memalign(DEFAULT_ALIGNMENT, num * size, file, line);

  if (x)
    VPX_MEMSET_L(x, 0, num * size);

  return x;
}

void *xvpx_realloc(void *memblk, size_t size, char *file, int line) {
  struct mem_block *p = NULL;
  int orig_size = 0,
      orig_line = 0;
  char *orig_file = NULL;

#if TRY_BOUNDS_CHECK
  unsigned char *x_bounds = memblk ?
                            (unsigned char *)(((size_t *)memblk)[-1]) :
                            NULL;
#endif

  void *x;

  if (g_alloc_count == 0) {
#if TRY_BOUNDS_CHECK

    if (!vpx_memory_tracker_init(BOUNDS_CHECK_PAD_SIZE, BOUNDS_CHECK_VALUE))
#else
    if (!vpx_memory_tracker_init(0, 0))
#endif
    {
      _P(printf("ERROR xvpx_malloc MEM_TRACK_USAGE error vpx_memory_tracker_init().\n");)
    }
  }

  if ((p = vpx_memory_tracker_find((size_t)memblk))) {
    orig_size = p->size;
    orig_file = p->file;
    orig_line = p->line;
  }

#if TRY_BOUNDS_CHECK_ON_FREE
  vpx_memory_tracker_check_integrity(file, line);
#endif

  /* have to do this regardless of success, because
   * the memory that does get realloc'd may change
   * the bounds values of this block
   */
  vpx_memory_tracker_remove((size_t)memblk);

#if TRY_BOUNDS_CHECK
  {
    int i;
    unsigned int tempme = BOUNDS_CHECK_VALUE;

    x_bounds = vpx_realloc(memblk, size + (BOUNDS_CHECK_PAD_SIZE * 2));

    if (x_bounds) {
      x_bounds = (unsigned char *)(((size_t *)x_bounds)[-1]);
      x = align_addr(x_bounds + BOUNDS_CHECK_PAD_SIZE + ADDRESS_STORAGE_SIZE,
                     (int)DEFAULT_ALIGNMENT);
      /* save the actual malloc address */
      ((size_t *)x)[-1] = (size_t)x_bounds;

      for (i = 0; i < BOUNDS_CHECK_PAD_SIZE; i += sizeof(unsigned int)) {
        VPX_MEMCPY_L(x_bounds + i, &tempme, sizeof(unsigned int));
        VPX_MEMCPY_L((unsigned char *)x + size + i,
                     &tempme, sizeof(unsigned int));
      }
    } else
      x = NULL;
  }
#else
  x = vpx_realloc(memblk, size);
#endif /*TRY_BOUNDS_CHECK*/

  if (!memblk) ++g_alloc_count;

  if (x)
    vpx_memory_tracker_add((size_t)x, (unsigned int)size, file, line, 1);
  else
    vpx_memory_tracker_add((size_t)memblk, orig_size, orig_file, orig_line, 1);

  return x;
}

void xvpx_free(void *p_address, char *file, int line) {
#if TRY_BOUNDS_CHECK
  unsigned char *p_bounds_address = (unsigned char *)p_address;
  /*p_bounds_address -= BOUNDS_CHECK_PAD_SIZE;*/
#endif

#if !TRY_BOUNDS_CHECK_ON_FREE
  (void)file;
  (void)line;
#endif

  if (p_address) {
#if TRY_BOUNDS_CHECK_ON_FREE
    vpx_memory_tracker_check_integrity(file, line);
#endif

    /* if the addr isn't found in the list, assume it was allocated via
     * vpx_ calls not xvpx_, therefore it does not contain any padding
     */
    if (vpx_memory_tracker_remove((size_t)p_address) == -2) {
      p_bounds_address = p_address;
      _P(fprintf(stderr, "[vpx_mem][xvpx_free] addr: %p not found in"
                 " list; freed from file:%s"
                 " line:%d\n", p_address, file, line));
    } else
      --g_alloc_count;

#if TRY_BOUNDS_CHECK
    vpx_free(p_bounds_address);
#else
    vpx_free(p_address);
#endif

    if (!g_alloc_count)
      vpx_memory_tracker_destroy();
  }
}

#endif /*CONFIG_MEM_TRACKER*/

void *vpx_memcpy(void *dest, const void *source, size_t length) {
  return VPX_MEMCPY_L(dest, source, length);
}

void *vpx_memset(void *dest, int val, size_t length) {
  return VPX_MEMSET_L(dest, val, length);
}

#if CONFIG_VP9 && CONFIG_VP9_HIGHBITDEPTH
void *vpx_memset16(void *dest, int val, size_t length) {
  int i;
  void *orig = dest;
  uint16_t *dest16 = dest;
  for (i = 0; i < length; i++)
    *dest16++ = val;
  return orig;
}
#endif  // CONFIG_VP9 && CONFIG_VP9_HIGHBITDEPTH

void *vpx_memmove(void *dest, const void *src, size_t count) {
  return VPX_MEMMOVE_L(dest, src, count);
}

#if USE_GLOBAL_FUNCTION_POINTERS
# if CONFIG_MEM_TRACKER
extern int vpx_memory_tracker_set_functions(g_malloc_func g_malloc_l
, g_calloc_func g_calloc_l
, g_realloc_func g_realloc_l
, g_free_func g_free_l
, g_memcpy_func g_memcpy_l
, g_memset_func g_memset_l
, g_memmove_func g_memmove_l);
# endif
#endif /*USE_GLOBAL_FUNCTION_POINTERS*/
int vpx_mem_set_functions(g_malloc_func g_malloc_l
, g_calloc_func g_calloc_l
, g_realloc_func g_realloc_l
, g_free_func g_free_l
, g_memcpy_func g_memcpy_l
, g_memset_func g_memset_l
, g_memmove_func g_memmove_l) {
#if USE_GLOBAL_FUNCTION_POINTERS

  /* If use global functions is turned on then the
  application must set the global functions before
  it does anything else or vpx_mem will have
  unpredictable results. */
  if (!g_func) {
    g_func = (struct GLOBAL_FUNC_POINTERS *)
             g_malloc_l(sizeof(struct GLOBAL_FUNC_POINTERS));

    if (!g_func) {
      return -1;
    }
  }

#if CONFIG_MEM_TRACKER
  {
    int rv = 0;
    rv = vpx_memory_tracker_set_functions(g_malloc_l
, g_calloc_l
, g_realloc_l
, g_free_l
, g_memcpy_l
, g_memset_l
, g_memmove_l);

    if (rv < 0) {
      return rv;
    }
  }
#endif

  g_func->g_malloc  = g_malloc_l;
  g_func->g_calloc  = g_calloc_l;
  g_func->g_realloc = g_realloc_l;
  g_func->g_free    = g_free_l;
  g_func->g_memcpy  = g_memcpy_l;
  g_func->g_memset  = g_memset_l;
  g_func->g_memmove = g_memmove_l;

  return 0;
#else
  (void)g_malloc_l;
  (void)g_calloc_l;
  (void)g_realloc_l;
  (void)g_free_l;
  (void)g_memcpy_l;
  (void)g_memset_l;
  (void)g_memmove_l;
  return -1;
#endif
}

int vpx_mem_unset_functions() {
#if USE_GLOBAL_FUNCTION_POINTERS

  if (g_func) {
    g_free_func temp_free = g_func->g_free;
    temp_free(g_func);
    g_func = NULL;
  }

#endif
  return 0;
}
