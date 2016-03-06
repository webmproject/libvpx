/*Daala video codec
Copyright (c) 2001-2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "entcode.h"

/*CDFs for uniform probability distributions of small sizes (2 through 16,
   inclusive).*/
const uint16_t OD_UNIFORM_CDFS_Q15[135] = {
  16384, 32768, 10923, 21845, 32768, 8192,  16384, 24576, 32768, 6554,  13107,
  19661, 26214, 32768, 5461,  10923, 16384, 21845, 27307, 32768, 4681,  9362,
  14043, 18725, 23406, 28087, 32768, 4096,  8192,  12288, 16384, 20480, 24576,
  28672, 32768, 3641,  7282,  10923, 14564, 18204, 21845, 25486, 29127, 32768,
  3277,  6554,  9830,  13107, 16384, 19661, 22938, 26214, 29491, 32768, 2979,
  5958,  8937,  11916, 14895, 17873, 20852, 23831, 26810, 29789, 32768, 2731,
  5461,  8192,  10923, 13653, 16384, 19115, 21845, 24576, 27307, 30037, 32768,
  2521,  5041,  7562,  10082, 12603, 15124, 17644, 20165, 22686, 25206, 27727,
  30247, 32768, 2341,  4681,  7022,  9362,  11703, 14043, 16384, 18725, 21065,
  23406, 25746, 28087, 30427, 32768, 2185,  4369,  6554,  8738,  10923, 13107,
  15292, 17476, 19661, 21845, 24030, 26214, 28399, 30583, 32768, 2048,  4096,
  6144,  8192,  10240, 12288, 14336, 16384, 18432, 20480, 22528, 24576, 26624,
  28672, 30720, 32768
};

/*Given the current total integer number of bits used and the current value of
   rng, computes the fraction number of bits used to OD_BITRES precision.
  This is used by od_ec_enc_tell_frac() and od_ec_dec_tell_frac().
  nbits_total: The number of whole bits currently used, i.e., the value
                returned by od_ec_enc_tell() or od_ec_dec_tell().
  rng: The current value of rng from either the encoder or decoder state.
  Return: The number of bits scaled by 2**OD_BITRES.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
uint32_t od_ec_tell_frac(uint32_t nbits_total, uint32_t rng) {
  uint32_t nbits;
  int l;
  int i;
  /*To handle the non-integral number of bits still left in the encoder/decoder
     state, we compute the worst-case number of bits of val that must be
     encoded to ensure that the value is inside the range for any possible
     subsequent bits.
    The computation here is independent of val itself (the decoder does not
     even track that value), even though the real number of bits used after
     od_ec_enc_done() may be 1 smaller if rng is a power of two and the
     corresponding trailing bits of val are all zeros.
    If we did try to track that special case, then coding a value with a
     probability of 1/(1 << n) might sometimes appear to use more than n bits.
    This may help explain the surprising result that a newly initialized
     encoder or decoder claims to have used 1 bit.*/
  nbits = nbits_total << OD_BITRES;
  l = 0;
  for (i = OD_BITRES; i-- > 0;) {
    int b;
    rng = rng * rng >> 15;
    b = (int)(rng >> 16);
    l = l << 1 | b;
    rng >>= b;
  }
  return nbits - l;
}
