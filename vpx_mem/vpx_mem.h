/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VPX_MEM_VPX_MEM_H_
#define VPX_MEM_VPX_MEM_H_

#include "vpx_config.h"
#if defined(__uClinux__)
# include <lddk.h>
#endif

#ifndef REPLACE_BUILTIN_FUNCTIONS
# define REPLACE_BUILTIN_FUNCTIONS 0  /* replace builtin functions with their
vpx_ equivalents */
#endif

#include <stdlib.h>
#include <stddef.h>

#if defined(__cplusplus)
extern "C" {
#endif

  void *vpx_memalign(size_t align, size_t size);
  void *vpx_malloc(size_t size);
  void *vpx_calloc(size_t num, size_t size);
  void *vpx_realloc(void *memblk, size_t size);
  void vpx_free(void *memblk);

  void *vpx_memcpy(void *dest, const void *src, size_t length);
  void *vpx_memset(void *dest, int val, size_t length);
#if CONFIG_VP9 && CONFIG_VP9_HIGHBITDEPTH
  void *vpx_memset16(void *dest, int val, size_t length);
#endif
  void *vpx_memmove(void *dest, const void *src, size_t count);

  /* some defines for backward compatibility */
#define DMEM_GENERAL 0

// (*)<

#if REPLACE_BUILTIN_FUNCTIONS
# ifndef __VPX_MEM_C__
#  define memalign vpx_memalign
#  define malloc   vpx_malloc
#  define calloc   vpx_calloc
#  define realloc  vpx_realloc
#  define free     vpx_free
#  define memcpy   vpx_memcpy
#  define memmove  vpx_memmove
#  define memset   vpx_memset
# endif
#endif

#ifndef __VPX_MEM_C__
# include <string.h>
# define vpx_memcpy  memcpy
# define vpx_memset  memset
# define vpx_memmove memmove
#endif

#ifdef VPX_MEM_PLTFRM
# include VPX_MEM_PLTFRM
#endif

#if defined(__cplusplus)
}
#endif

#endif  // VPX_MEM_VPX_MEM_H_
