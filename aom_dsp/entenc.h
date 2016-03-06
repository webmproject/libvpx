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

#if !defined(_entenc_H)
#define _entenc_H (1)
#include <stddef.h>
#include "entcode.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct od_ec_enc od_ec_enc;

#define OD_MEASURE_EC_OVERHEAD (0)

/*The entropy encoder context.*/
struct od_ec_enc {
  /*Buffered output.
    This contains only the raw bits until the final call to od_ec_enc_done(),
     where all the arithmetic-coded data gets prepended to it.*/
  unsigned char *buf;
  /*The size of the buffer.*/
  uint32_t storage;
  /*The offset at which the last byte containing raw bits was written.*/
  uint32_t end_offs;
  /*Bits that will be read from/written at the end.*/
  od_ec_window end_window;
  /*Number of valid bits in end_window.*/
  int nend_bits;
  /*A buffer for output bytes with their associated carry flags.*/
  uint16_t *precarry_buf;
  /*The size of the pre-carry buffer.*/
  uint32_t precarry_storage;
  /*The offset at which the next entropy-coded byte will be written.*/
  uint32_t offs;
  /*The low end of the current range.*/
  od_ec_window low;
  /*The number of values in the current range.*/
  uint16_t rng;
  /*The number of bits of data in the current value.*/
  int16_t cnt;
  /*Nonzero if an error occurred.*/
  int error;
#if OD_MEASURE_EC_OVERHEAD
  double entropy;
  int nb_symbols;
#endif
};

/*See entenc.c for further documentation.*/

void od_ec_enc_init(od_ec_enc *enc, uint32_t size) OD_ARG_NONNULL(1);
void od_ec_enc_reset(od_ec_enc *enc) OD_ARG_NONNULL(1);
void od_ec_enc_clear(od_ec_enc *enc) OD_ARG_NONNULL(1);

void od_ec_encode_bool(od_ec_enc *enc, int val, unsigned fz, unsigned _ft)
    OD_ARG_NONNULL(1);
void od_ec_encode_bool_q15(od_ec_enc *enc, int val, unsigned fz_q15)
    OD_ARG_NONNULL(1);
void od_ec_encode_cdf(od_ec_enc *enc, int s, const uint16_t *cdf, int nsyms)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(3);
void od_ec_encode_cdf_q15(od_ec_enc *enc, int s, const uint16_t *cdf, int nsyms)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(3);
void od_ec_encode_cdf_unscaled(od_ec_enc *enc, int s, const uint16_t *cdf,
                               int nsyms) OD_ARG_NONNULL(1) OD_ARG_NONNULL(3);
void od_ec_encode_cdf_unscaled_dyadic(od_ec_enc *enc, int s,
                                      const uint16_t *cdf, int nsyms,
                                      unsigned ftb) OD_ARG_NONNULL(1)
    OD_ARG_NONNULL(3);

void od_ec_enc_uint(od_ec_enc *enc, uint32_t fl, uint32_t ft) OD_ARG_NONNULL(1);

void od_ec_enc_bits(od_ec_enc *enc, uint32_t fl, unsigned ftb)
    OD_ARG_NONNULL(1);

void od_ec_enc_patch_initial_bits(od_ec_enc *enc, unsigned val, int nbits)
    OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT unsigned char *od_ec_enc_done(od_ec_enc *enc,
                                                    uint32_t *nbytes)
    OD_ARG_NONNULL(1) OD_ARG_NONNULL(2);

OD_WARN_UNUSED_RESULT int od_ec_enc_tell(od_ec_enc *enc) OD_ARG_NONNULL(1);
OD_WARN_UNUSED_RESULT uint32_t od_ec_enc_tell_frac(od_ec_enc *enc)
    OD_ARG_NONNULL(1);

void od_ec_enc_checkpoint(od_ec_enc *dst, const od_ec_enc *src);
void od_ec_enc_rollback(od_ec_enc *dst, const od_ec_enc *src);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
