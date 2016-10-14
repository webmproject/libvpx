/*Daala video codec
Copyright (c) 2001-2013 Daala project contributors.  All rights reserved.

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

#if !defined(_entdec_H)
#define _entdec_H (1)
#include <limits.h>
#include "aom_dsp/entcode.h"
#if OD_ACCOUNTING
#include "./accounting.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct od_ec_dec od_ec_dec;

#if OD_ACCOUNTING
#define OD_ACC_STR , char *acc_str
#define od_ec_decode_bool(dec, fz, ft, str) od_ec_decode_bool_(dec, fz, ft, str)
#define od_ec_decode_bool_q15(dec, fz, str) od_ec_decode_bool_q15_(dec, fz, str)
#define od_ec_decode_cdf(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_(dec, cdf, nsyms, str)
#define od_ec_decode_cdf_q15(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_q15_(dec, cdf, nsyms, str)
#define od_ec_decode_cdf_unscaled(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_unscaled_(dec, cdf, nsyms, str)
#define od_ec_decode_cdf_unscaled_dyadic(dec, cdf, nsyms, ftb, str) \
  od_ec_decode_cdf_unscaled_dyadic_(dec, cdf, nsyms, ftb, str)
#define od_ec_dec_uint(dec, ft, str) od_ec_dec_uint_(dec, ft, str)
#define od_ec_dec_bits(dec, ftb, str) od_ec_dec_bits_(dec, ftb, str)
#else
#define OD_ACC_STR
#define od_ec_decode_bool(dec, fz, ft, str) od_ec_decode_bool_(dec, fz, ft)
#define od_ec_decode_bool_q15(dec, fz, str) od_ec_decode_bool_q15_(dec, fz)
#define od_ec_decode_cdf(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_(dec, cdf, nsyms)
#define od_ec_decode_cdf_q15(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_q15_(dec, cdf, nsyms)
#define od_ec_decode_cdf_unscaled(dec, cdf, nsyms, str) \
  od_ec_decode_cdf_unscaled_(dec, cdf, nsyms)
#define od_ec_decode_cdf_unscaled_dyadic(dec, cdf, nsyms, ftb, str) \
  od_ec_decode_cdf_unscaled_dyadic_(dec, cdf, nsyms, ftb)
#define od_ec_dec_uint(dec, ft, str) od_ec_dec_uint_(dec, ft)
#define od_ec_dec_bits(dec, ftb, str) od_ec_dec_bits_(dec, ftb)
#endif

/*The entropy decoder context.*/
struct od_ec_dec {
  /*The start of the current input buffer.*/
  const unsigned char *buf;
  /*The read pointer for the raw bits.*/
  const unsigned char *eptr;
  /*Bits that will be read from/written at the end.*/
  od_ec_window end_window;
  /*Number of valid bits in end_window.*/
  int nend_bits;
  /*An offset used to keep track of tell after reaching the end of the stream.
    This is constant throughout most of the decoding process, but becomes
     important once we hit the end of the buffer and stop incrementing pointers
     (and instead pretend cnt/nend_bits have lots of bits).*/
  int32_t tell_offs;
  /*The end of the current input buffer.*/
  const unsigned char *end;
  /*The read pointer for the entropy-coded bits.*/
  const unsigned char *bptr;
  /*The difference between the coded value and the low end of the current
     range.*/
  od_ec_window dif;
  /*The number of values in the current range.*/
  uint16_t rng;
  /*The number of bits of data in the current value.*/
  int16_t cnt;
  /*Nonzero if an error occurred.*/
  int error;
#if OD_ACCOUNTING
  od_accounting_internal *acct;
#endif
};

/*See entdec.c for further documentation.*/

void od_ec_dec_init(od_ec_dec *dec, const unsigned char *buf, uint32_t storage)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT int od_ec_decode_bool_(od_ec_dec *dec, unsigned fz,
                                             unsigned ft OD_ACC_STR)
    OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT int od_ec_decode_bool_q15_(od_ec_dec *dec,
                                                 unsigned fz OD_ACC_STR)
    OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT int od_ec_decode_cdf_(od_ec_dec *dec, const uint16_t *cdf,
                                            int nsyms OD_ACC_STR)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);
OD_WARN_UNUSED_RESULT int od_ec_decode_cdf_q15_(od_ec_dec *dec,
                                                const uint16_t *cdf,
                                                int nsyms OD_ACC_STR)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);
OD_WARN_UNUSED_RESULT int od_ec_decode_cdf_unscaled_(od_ec_dec *dec,
                                                     const uint16_t *cdf,
                                                     int nsyms OD_ACC_STR)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);
OD_WARN_UNUSED_RESULT int od_ec_decode_cdf_unscaled_dyadic_(
    od_ec_dec *dec, const uint16_t *cdf, int nsyms, unsigned _ftb OD_ACC_STR)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT uint32_t od_ec_dec_uint_(od_ec_dec *dec,
                                               uint32_t ft OD_ACC_STR)
    OD_ARG_NONNULL(1);

OD_WARN_UNUSED_RESULT uint32_t od_ec_dec_bits_(od_ec_dec *dec,
                                               unsigned ftb OD_ACC_STR)
    OD_ARG_NONNULL(1);

OD_WARN_UNUSED_RESULT int od_ec_dec_tell(od_ec_dec *dec) OD_ARG_NONNULL(1);
uint32_t od_ec_dec_tell_frac(od_ec_dec *dec) OD_ARG_NONNULL(1);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
