/*
 * Copyright (c) 2001-2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#if !defined(_entcode_H)
#define _entcode_H (1)
#include <limits.h>
#include <stddef.h>
#include "av1/common/odintrin.h"

/*Set this flag 1 to enable a "reduced overhead" version of the entropy coder.
  This uses a partition function that more accurately follows the input
   probability estimates at the expense of some additional CPU cost (though
   still an order of magnitude less than a full division).

  In classic arithmetic coding, the partition function maps a value x in the
   range [0, ft] to a value in y in [0, r] with 0 < ft <= r via
    y = x*r/ft.
  Any deviation from this value increases coding inefficiency.

  To avoid divisions, we require ft <= r < 2*ft (enforcing it by shifting up
   ft if necessary), and replace that function with
    y = x + OD_MINI(x, r - ft).
  This counts values of x smaller than r - ft double compared to values larger
   than r - ft, which over-estimates the probability of symbols at the start of
   the alphabet, and under-estimates the probability of symbols at the end of
   the alphabet.
  The overall coding inefficiency assuming accurate probability models and
   independent symbols is in the 1% range, which is similar to that of CABAC.

  To reduce overhead even further, we split this into two cases:
  1) r - ft > ft - (r - ft).
     That is, we have more values of x that are double-counted than
      single-counted.
     In this case, we still double-count the first 2*r - 3*ft values of x, but
      after that we alternate between single-counting and double-counting for
      the rest.
  2) r - ft < ft - (r - ft).
     That is, we have more values of x that are single-counted than
      double-counted.
     In this case, we alternate between single-counting and double-counting for
      the first 2*(r - ft) values of x, and single-count the rest.
  For two equiprobable symbols in different places in the alphabet, this
   reduces the maximum ratio of over-estimation to under-estimation from 2:1
   for the previous partition function to either 4:3 or 3:2 (for each of the
   two cases above, respectively), assuming symbol probabilities significantly
   greater than 1/32768.
  That reduces the worst-case per-symbol overhead from 1 bit to 0.58 bits.

  The resulting function is
    e = OD_MAXI(2*r - 3*ft, 0);
    y = x + OD_MINI(x, e) + OD_MINI(OD_MAXI(x - e, 0) >> 1, r - ft).
  Here, e is a value that is greater than 0 in case 1, and 0 in case 2.
  This function is about 3 times as expensive to evaluate as the high-overhead
   version, but still an order of magnitude cheaper than a division, since it
   is composed only of very simple operations.
  Because we want to fit in 16-bit registers and must use unsigned values to do
   so, we use saturating subtraction to enforce the maximums with 0.

  Enabling this reduces the measured overhead in ectest from 0.805% to 0.621%
   (vs. 0.022% for the division-based partition function with r much greater
   than ft).
  It improves performance on ntt-short-1 by about 0.3%.*/
#define OD_EC_REDUCED_OVERHEAD (1)

/*OPT: od_ec_window must be at least 32 bits, but if you have fast arithmetic
   on a larger type, you can speed up the decoder by using it here.*/
typedef uint32_t od_ec_window;

#define OD_EC_WINDOW_SIZE ((int)sizeof(od_ec_window) * CHAR_BIT)

/*Unsigned subtraction with unsigned saturation.
  This implementation of the macro is intentionally chosen to increase the
   number of common subexpressions in the reduced-overhead partition function.
  This matters for C code, but it would not for hardware with a saturating
   subtraction instruction.*/
#define OD_SUBSATU(a, b) ((a)-OD_MINI(a, b))

/*The number of bits to use for the range-coded part of unsigned integers.*/
#define OD_EC_UINT_BITS (4)

/*The resolution of fractional-precision bit usage measurements, i.e.,
   3 => 1/8th bits.*/
#define OD_BITRES (3)

extern const uint16_t OD_UNIFORM_CDFS_Q15[135];

/*Returns a Q15 CDF for a uniform probability distribution of the given size.
  n: The size of the distribution.
     This must be at least 2, and no more than 16.*/
#define OD_UNIFORM_CDF_Q15(n) (OD_UNIFORM_CDFS_Q15 + ((n) * ((n)-1) >> 1) - 1)

/*See entcode.c for further documentation.*/

OD_WARN_UNUSED_RESULT uint32_t od_ec_tell_frac(uint32_t nbits_total,
                                               uint32_t rng);

#endif
