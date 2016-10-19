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
#ifndef AV1_COMMON_ODINTRIN_H_
#define AV1_COMMON_ODINTRIN_H_

#include "aom/aom_integer.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_ports/bitops.h"
#include "av1/common/enums.h"

#ifdef __cplusplus
extern "C" {
#endif

/*Smallest blocks are 4x4*/
#define OD_LOG_BSIZE0 (2)
/*There are 5 block sizes total (4x4, 8x8, 16x16, 32x32 and 64x64).*/
#define OD_NBSIZES (5)
/*The log of the maximum length of the side of a block.*/
#define OD_LOG_BSIZE_MAX (OD_LOG_BSIZE0 + OD_NBSIZES - 1)
/*The maximum length of the side of a block.*/
#define OD_BSIZE_MAX (1 << OD_LOG_BSIZE_MAX)

typedef int od_coeff;

typedef int16_t od_dering_in;

#define OD_DIVU_DMAX (1024)

extern uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2];

/*Translate unsigned division by small divisors into multiplications.*/
#define OD_DIVU_SMALL(_x, _d)                                     \
  ((uint32_t)((OD_DIVU_SMALL_CONSTS[(_d)-1][0] * (uint64_t)(_x) + \
               OD_DIVU_SMALL_CONSTS[(_d)-1][1]) >>                \
              32) >>                                              \
   (OD_ILOG(_d) - 1))

#define OD_DIVU(_x, _d) \
  (((_d) < OD_DIVU_DMAX) ? (OD_DIVU_SMALL((_x), (_d))) : ((_x) / (_d)))

#define OD_MINI AOMMIN
#define OD_MAXI AOMMAX
#define OD_CLAMPI(min, val, max) clamp((val), (min), (max))

#define OD_CLZ0 (1)
#define OD_CLZ(x) (-get_msb(x))
#define OD_ILOG_NZ(x) (OD_CLZ0 - OD_CLZ(x))
/*Note that __builtin_clz is not defined when x == 0, according to the gcc
   documentation (and that of the x86 BSR instruction that implements it), so
   we have to special-case it.
  We define a special version of the macro to use when x can be zero.*/
#define OD_ILOG(x) ((x) ? OD_ILOG_NZ(x) : 0)

#define OD_LOG2 AOMLOG2

/*Enable special features for gcc and compatible compilers.*/
#if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#define OD_GNUC_PREREQ(maj, min, pat)                                \
  ((__GNUC__ << 16) + (__GNUC_MINOR__ << 8) + __GNUC_PATCHLEVEL__ >= \
   ((maj) << 16) + ((min) << 8) + pat)  // NOLINT
#else
#define OD_GNUC_PREREQ(maj, min, pat) (0)
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
#define OD_WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))
#else
#define OD_WARN_UNUSED_RESULT
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
#define OD_ARG_NONNULL(x) __attribute__((__nonnull__(x)))
#else
#define OD_ARG_NONNULL(x)
#endif

#if defined(OD_ENABLE_ASSERTIONS)
#if OD_GNUC_PREREQ(2, 5, 0)
__attribute__((noreturn))
#endif
void od_fatal_impl(const char *_str, const char *_file, int _line);

#define OD_FATAL(_str) (od_fatal_impl(_str, __FILE__, __LINE__))

#define OD_ASSERT(_cond)                     \
  do {                                       \
    if (!(_cond)) {                          \
      OD_FATAL("assertion failed: " #_cond); \
    }                                        \
  } while (0)

#define OD_ASSERT2(_cond, _message)                        \
  do {                                                     \
    if (!(_cond)) {                                        \
      OD_FATAL("assertion failed: " #_cond "\n" _message); \
    }                                                      \
  } while (0)

#define OD_ALWAYS_TRUE(_cond) OD_ASSERT(_cond)

#else
#define OD_ASSERT(_cond)
#define OD_ASSERT2(_cond, _message)
#define OD_ALWAYS_TRUE(_cond) ((void)(_cond))
#endif

/** Copy n elements of memory from src to dst. The 0* term provides
    compile-time type checking  */
#if !defined(OVERRIDE_OD_COPY)
#define OD_COPY(dst, src, n) \
  (memcpy((dst), (src), sizeof(*(dst)) * (n) + 0 * ((dst) - (src))))
#endif

/** Copy n elements of memory from src to dst, allowing overlapping regions.
    The 0* term provides compile-time type checking */
#if !defined(OVERRIDE_OD_MOVE)
#define OD_MOVE(dst, src, n) \
  (memmove((dst), (src), sizeof(*(dst)) * (n) + 0 * ((dst) - (src))))
#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_COMMON_ODINTRIN_H_
