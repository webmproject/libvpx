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

#ifdef HAVE_CONFIG_H
#include "./config.h"
#endif

#include "aom_dsp/entdec.h"

/*A range decoder.
  This is an entropy decoder based upon \cite{Mar79}, which is itself a
   rediscovery of the FIFO arithmetic code introduced by \cite{Pas76}.
  It is very similar to arithmetic encoding, except that encoding is done with
   digits in any base, instead of with bits, and so it is faster when using
   larger bases (i.e.: a byte).
  The author claims an average waste of $\frac{1}{2}\log_b(2b)$ bits, where $b$
   is the base, longer than the theoretical optimum, but to my knowledge there
   is no published justification for this claim.
  This only seems true when using near-infinite precision arithmetic so that
   the process is carried out with no rounding errors.

  An excellent description of implementation details is available at
   http://www.arturocampos.com/ac_range.html
  A recent work \cite{MNW98} which proposes several changes to arithmetic
   encoding for efficiency actually re-discovers many of the principles
   behind range encoding, and presents a good theoretical analysis of them.

  End of stream is handled by writing out the smallest number of bits that
   ensures that the stream will be correctly decoded regardless of the value of
   any subsequent bits.
  od_ec_dec_tell() can be used to determine how many bits were needed to decode
   all the symbols thus far; other data can be packed in the remaining bits of
   the input buffer.
  @PHDTHESIS{Pas76,
    author="Richard Clark Pasco",
    title="Source coding algorithms for fast data compression",
    school="Dept. of Electrical Engineering, Stanford University",
    address="Stanford, CA",
    month=May,
    year=1976,
    URL="http://www.richpasco.org/scaffdc.pdf"
  }
  @INPROCEEDINGS{Mar79,
   author="Martin, G.N.N.",
   title="Range encoding: an algorithm for removing redundancy from a digitised
    message",
   booktitle="Video & Data Recording Conference",
   year=1979,
   address="Southampton",
   month=Jul,
   URL="http://www.compressconsult.com/rangecoder/rngcod.pdf.gz"
  }
  @ARTICLE{MNW98,
   author="Alistair Moffat and Radford Neal and Ian H. Witten",
   title="Arithmetic Coding Revisited",
   journal="{ACM} Transactions on Information Systems",
   year=1998,
   volume=16,
   number=3,
   pages="256--294",
   month=Jul,
   URL="http://researchcommons.waikato.ac.nz/bitstream/handle/10289/78/content.pdf"
  }*/

/*This is meant to be a large, positive constant that can still be efficiently
   loaded as an immediate (on platforms like ARM, for example).
  Even relatively modest values like 100 would work fine.*/
#define OD_EC_LOTS_OF_BITS (0x4000)

static void od_ec_dec_refill(od_ec_dec *dec) {
  int s;
  od_ec_window dif;
  int16_t cnt;
  const unsigned char *bptr;
  const unsigned char *end;
  dif = dec->dif;
  cnt = dec->cnt;
  bptr = dec->bptr;
  end = dec->end;
  s = OD_EC_WINDOW_SIZE - 9 - (cnt + 15);
  for (; s >= 0 && bptr < end; s -= 8, bptr++) {
    OD_ASSERT(s <= OD_EC_WINDOW_SIZE - 8);
    dif |= (od_ec_window)bptr[0] << s;
    cnt += 8;
  }
  if (bptr >= end) {
    dec->tell_offs += OD_EC_LOTS_OF_BITS - cnt;
    cnt = OD_EC_LOTS_OF_BITS;
  }
  dec->dif = dif;
  dec->cnt = cnt;
  dec->bptr = bptr;
}

/*Takes updated dif and range values, renormalizes them so that
   32768 <= rng < 65536 (reading more bytes from the stream into dif if
   necessary), and stores them back in the decoder context.
  dif: The new value of dif.
  rng: The new value of the range.
  ret: The value to return.
  Return: ret.
          This allows the compiler to jump to this function via a tail-call.*/
static int od_ec_dec_normalize(od_ec_dec *dec, od_ec_window dif, unsigned rng,
                               int ret) {
  int d;
  OD_ASSERT(rng <= 65535U);
  d = 16 - OD_ILOG_NZ(rng);
  dec->cnt -= d;
  dec->dif = dif << d;
  dec->rng = rng << d;
  if (dec->cnt < 0) od_ec_dec_refill(dec);
  return ret;
}

/*Initializes the decoder.
  buf: The input buffer to use.
  Return: 0 on success, or a negative value on error.*/
void od_ec_dec_init(od_ec_dec *dec, const unsigned char *buf,
                    uint32_t storage) {
  dec->buf = buf;
  dec->eptr = buf + storage;
  dec->end_window = 0;
  dec->nend_bits = 0;
  dec->tell_offs = 10 - (OD_EC_WINDOW_SIZE - 8);
  dec->end = buf + storage;
  dec->bptr = buf;
  dec->dif = 0;
  dec->rng = 0x8000;
  dec->cnt = -15;
  dec->error = 0;
  od_ec_dec_refill(dec);
}

/*Decode a bit that has an fz/ft probability of being a zero.
  fz: The probability that the bit is zero, scaled by _ft.
  ft: The total probability.
      This must be at least 16384 and no more than 32768.
  Return: The value decoded (0 or 1).*/
int od_ec_decode_bool(od_ec_dec *dec, unsigned fz, unsigned ft) {
  od_ec_window dif;
  od_ec_window vw;
  unsigned r;
  int s;
  unsigned v;
  int ret;
  OD_ASSERT(0 < fz);
  OD_ASSERT(fz < ft);
  OD_ASSERT(16384 <= ft);
  OD_ASSERT(ft <= 32768U);
  dif = dec->dif;
  r = dec->rng;
  OD_ASSERT(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
  OD_ASSERT(ft <= r);
  s = r - ft >= ft;
  ft <<= s;
  fz <<= s;
  OD_ASSERT(r - ft < ft);
#if OD_EC_REDUCED_OVERHEAD
  {
    unsigned d;
    unsigned e;
    d = r - ft;
    e = OD_SUBSATU(2 * d, ft);
    v = fz + OD_MINI(fz, e) + OD_MINI(OD_SUBSATU(fz, e) >> 1, d);
  }
#else
  v = fz + OD_MINI(fz, r - ft);
#endif
  vw = (od_ec_window)v << (OD_EC_WINDOW_SIZE - 16);
  ret = dif >= vw;
  if (ret) dif -= vw;
  r = ret ? r - v : v;
  return od_ec_dec_normalize(dec, dif, r, ret);
}

/*Decode a bit that has an fz probability of being a zero in Q15.
  This is a simpler, lower overhead version of od_ec_decode_bool() for use when
   ft == 32768.
  To be decoded properly by this function, symbols cannot have been encoded by
   od_ec_encode(), but must have been encoded with one of the equivalent _q15()
   or _dyadic() functions instead.
  fz: The probability that the bit is zero, scaled by 32768.
  Return: The value decoded (0 or 1).*/
int od_ec_decode_bool_q15(od_ec_dec *dec, unsigned fz) {
  od_ec_window dif;
  od_ec_window vw;
  unsigned r;
  unsigned r_new;
  unsigned v;
  int ret;
  OD_ASSERT(0 < fz);
  OD_ASSERT(fz < 32768U);
  dif = dec->dif;
  r = dec->rng;
  OD_ASSERT(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
  OD_ASSERT(32768U <= r);
  v = fz * (uint32_t)r >> 15;
  vw = (od_ec_window)v << (OD_EC_WINDOW_SIZE - 16);
  ret = 0;
  r_new = v;
  if (dif >= vw) {
    r_new = r - v;
    dif -= vw;
    ret = 1;
  }
  return od_ec_dec_normalize(dec, dif, r_new, ret);
}

/*Decodes a symbol given a cumulative distribution function (CDF) table.
  cdf: The CDF, such that symbol s falls in the range
        [s > 0 ? cdf[s - 1] : 0, cdf[s]).
       The values must be monotonically non-increasing, and cdf[nsyms - 1]
        must be at least 16384, and no more than 32768.
  nsyms: The number of symbols in the alphabet.
         This should be at most 16.
  Return: The decoded symbol s.*/
int od_ec_decode_cdf(od_ec_dec *dec, const uint16_t *cdf, int nsyms) {
  od_ec_window dif;
  unsigned r;
  unsigned c;
  unsigned d;
#if OD_EC_REDUCED_OVERHEAD
  unsigned e;
#endif
  int s;
  unsigned u;
  unsigned v;
  unsigned q;
  unsigned fl;
  unsigned fh;
  unsigned ft;
  int ret;
  dif = dec->dif;
  r = dec->rng;
  OD_ASSERT(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
  OD_ASSERT(nsyms > 0);
  ft = cdf[nsyms - 1];
  OD_ASSERT(16384 <= ft);
  OD_ASSERT(ft <= 32768U);
  OD_ASSERT(ft <= r);
  s = r - ft >= ft;
  ft <<= s;
  d = r - ft;
  OD_ASSERT(d < ft);
  c = (unsigned)(dif >> (OD_EC_WINDOW_SIZE - 16));
  q = OD_MAXI((int)(c >> 1), (int)(c - d));
#if OD_EC_REDUCED_OVERHEAD
  e = OD_SUBSATU(2 * d, ft);
  /*The correctness of this inverse partition function is not obvious, but it
     was checked exhaustively for all possible values of r, ft, and c.
    TODO: It should be possible to optimize this better than the compiler,
     given that we do not care about the accuracy of negative results (as we
     will not use them).
    It would also be nice to get rid of the 32-bit dividend, as it requires a
     32x32->64 bit multiply to invert.*/
  q = OD_MAXI((int)q, (int)((2 * (int32_t)c + 1 - (int32_t)e) / 3));
#endif
  q >>= s;
  OD_ASSERT(q<ft>> s);
  fl = 0;
  ret = 0;
  for (fh = cdf[ret]; fh <= q; fh = cdf[++ret]) fl = fh;
  OD_ASSERT(fh <= ft >> s);
  fl <<= s;
  fh <<= s;
#if OD_EC_REDUCED_OVERHEAD
  u = fl + OD_MINI(fl, e) + OD_MINI(OD_SUBSATU(fl, e) >> 1, d);
  v = fh + OD_MINI(fh, e) + OD_MINI(OD_SUBSATU(fh, e) >> 1, d);
#else
  u = fl + OD_MINI(fl, d);
  v = fh + OD_MINI(fh, d);
#endif
  r = v - u;
  dif -= (od_ec_window)u << (OD_EC_WINDOW_SIZE - 16);
  return od_ec_dec_normalize(dec, dif, r, ret);
}

/*Decodes a symbol given a cumulative distribution function (CDF) table.
  cdf: The CDF, such that symbol s falls in the range
        [s > 0 ? cdf[s - 1] : 0, cdf[s]).
       The values must be monotonically non-increasing, and cdf[nsyms - 1]
       must be at least 2, and no more than 32768.
  nsyms: The number of symbols in the alphabet.
         This should be at most 16.
  Return: The decoded symbol s.*/
int od_ec_decode_cdf_unscaled(od_ec_dec *dec, const uint16_t *cdf, int nsyms) {
  od_ec_window dif;
  unsigned r;
  unsigned c;
  unsigned d;
#if OD_EC_REDUCED_OVERHEAD
  unsigned e;
#endif
  int s;
  unsigned u;
  unsigned v;
  unsigned q;
  unsigned fl;
  unsigned fh;
  unsigned ft;
  int ret;
  dif = dec->dif;
  r = dec->rng;
  OD_ASSERT(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
  OD_ASSERT(nsyms > 0);
  ft = cdf[nsyms - 1];
  OD_ASSERT(2 <= ft);
  OD_ASSERT(ft <= 32768U);
  s = 15 - OD_ILOG_NZ(ft - 1);
  ft <<= s;
  OD_ASSERT(ft <= r);
  if (r - ft >= ft) {
    ft <<= 1;
    s++;
  }
  d = r - ft;
  OD_ASSERT(d < ft);
  c = (unsigned)(dif >> (OD_EC_WINDOW_SIZE - 16));
  q = OD_MAXI((int)(c >> 1), (int)(c - d));
#if OD_EC_REDUCED_OVERHEAD
  e = OD_SUBSATU(2 * d, ft);
  /*TODO: See TODO above.*/
  q = OD_MAXI((int)q, (int)((2 * (int32_t)c + 1 - (int32_t)e) / 3));
#endif
  q >>= s;
  OD_ASSERT(q<ft>> s);
  fl = 0;
  ret = 0;
  for (fh = cdf[ret]; fh <= q; fh = cdf[++ret]) fl = fh;
  OD_ASSERT(fh <= ft >> s);
  fl <<= s;
  fh <<= s;
#if OD_EC_REDUCED_OVERHEAD
  u = fl + OD_MINI(fl, e) + OD_MINI(OD_SUBSATU(fl, e) >> 1, d);
  v = fh + OD_MINI(fh, e) + OD_MINI(OD_SUBSATU(fh, e) >> 1, d);
#else
  u = fl + OD_MINI(fl, d);
  v = fh + OD_MINI(fh, d);
#endif
  r = v - u;
  dif -= (od_ec_window)u << (OD_EC_WINDOW_SIZE - 16);
  return od_ec_dec_normalize(dec, dif, r, ret);
}

/*Decodes a symbol given a cumulative distribution function (CDF) table that
   sums to a power of two.
  This is a simpler, lower overhead version of od_ec_decode_cdf() for use when
   cdf[nsyms - 1] is a power of two.
  To be decoded properly by this function, symbols cannot have been encoded by
   od_ec_encode(), but must have been encoded with one of the equivalent _q15()
   functions instead.
  cdf: The CDF, such that symbol s falls in the range
        [s > 0 ? cdf[s - 1] : 0, cdf[s]).
       The values must be monotonically non-increasing, and cdf[nsyms - 1]
       must be exactly 1 << ftb.
  nsyms: The number of symbols in the alphabet.
         This should be at most 16.
  ftb: The number of bits of precision in the cumulative distribution.
       This must be no more than 15.
  Return: The decoded symbol s.*/
int od_ec_decode_cdf_unscaled_dyadic(od_ec_dec *dec, const uint16_t *cdf,
                                     int nsyms, unsigned ftb) {
  od_ec_window dif;
  unsigned r;
  unsigned c;
  unsigned u;
  unsigned v;
  int ret;
  (void)nsyms;
  dif = dec->dif;
  r = dec->rng;
  OD_ASSERT(dif >> (OD_EC_WINDOW_SIZE - 16) < r);
  OD_ASSERT(ftb <= 15);
  OD_ASSERT(cdf[nsyms - 1] == 1U << ftb);
  OD_ASSERT(32768U <= r);
  c = (unsigned)(dif >> (OD_EC_WINDOW_SIZE - 16));
  v = 0;
  ret = -1;
  do {
    u = v;
    v = cdf[++ret] * (uint32_t)r >> ftb;
  } while (v <= c);
  OD_ASSERT(v <= r);
  r = v - u;
  dif -= (od_ec_window)u << (OD_EC_WINDOW_SIZE - 16);
  return od_ec_dec_normalize(dec, dif, r, ret);
}

/*Decodes a symbol given a cumulative distribution function (CDF) table in Q15.
  This is a simpler, lower overhead version of od_ec_decode_cdf() for use when
   cdf[nsyms - 1] == 32768.
  To be decoded properly by this function, symbols cannot have been encoded by
   od_ec_encode(), but must have been encoded with one of the equivalent _q15()
   or dyadic() functions instead.
  cdf: The CDF, such that symbol s falls in the range
        [s > 0 ? cdf[s - 1] : 0, cdf[s]).
       The values must be monotonically non-increasing, and cdf[nsyms - 1]
        must be 32768.
  nsyms: The number of symbols in the alphabet.
         This should be at most 16.
  Return: The decoded symbol s.*/
int od_ec_decode_cdf_q15(od_ec_dec *dec, const uint16_t *cdf, int nsyms) {
  return od_ec_decode_cdf_unscaled_dyadic(dec, cdf, nsyms, 15);
}

/*Extracts a raw unsigned integer with a non-power-of-2 range from the stream.
  The integer must have been encoded with od_ec_enc_uint().
  ft: The number of integers that can be decoded (one more than the max).
      This must be at least 2, and no more than 2**29.
  Return: The decoded bits.*/
uint32_t od_ec_dec_uint(od_ec_dec *dec, uint32_t ft) {
  OD_ASSERT(ft >= 2);
  OD_ASSERT(ft <= (uint32_t)1 << (25 + OD_EC_UINT_BITS));
  if (ft > 1U << OD_EC_UINT_BITS) {
    uint32_t t;
    int ft1;
    int ftb;
    ft--;
    ftb = OD_ILOG_NZ(ft) - OD_EC_UINT_BITS;
    ft1 = (int)(ft >> ftb) + 1;
    t = od_ec_decode_cdf_q15(dec, OD_UNIFORM_CDF_Q15(ft1), ft1);
    t = t << ftb | od_ec_dec_bits(dec, ftb, "");
    if (t <= ft) return t;
    dec->error = 1;
    return ft;
  }
  return od_ec_decode_cdf_q15(dec, OD_UNIFORM_CDF_Q15(ft), (int)ft);
}

/*Extracts a sequence of raw bits from the stream.
  The bits must have been encoded with od_ec_enc_bits().
  ftb: The number of bits to extract.
       This must be between 0 and 25, inclusive.
  Return: The decoded bits.*/
uint32_t od_ec_dec_bits_(od_ec_dec *dec, unsigned ftb) {
  od_ec_window window;
  int available;
  uint32_t ret;
  OD_ASSERT(ftb <= 25);
  window = dec->end_window;
  available = dec->nend_bits;
  if ((unsigned)available < ftb) {
    const unsigned char *buf;
    const unsigned char *eptr;
    buf = dec->buf;
    eptr = dec->eptr;
    OD_ASSERT(available <= OD_EC_WINDOW_SIZE - 8);
    do {
      if (eptr <= buf) {
        dec->tell_offs += OD_EC_LOTS_OF_BITS - available;
        available = OD_EC_LOTS_OF_BITS;
        break;
      }
      window |= (od_ec_window) * --eptr << available;
      available += 8;
    } while (available <= OD_EC_WINDOW_SIZE - 8);
    dec->eptr = eptr;
  }
  ret = (uint32_t)window & (((uint32_t)1 << ftb) - 1);
  window >>= ftb;
  available -= ftb;
  dec->end_window = window;
  dec->nend_bits = available;
  return ret;
}

/*Returns the number of bits "used" by the decoded symbols so far.
  This same number can be computed in either the encoder or the decoder, and is
   suitable for making coding decisions.
  Return: The number of bits.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
int od_ec_dec_tell(const od_ec_dec *dec) {
  return ((dec->end - dec->eptr) + (dec->bptr - dec->buf)) * 8 - dec->cnt -
         dec->nend_bits + dec->tell_offs;
}

/*Returns the number of bits "used" by the decoded symbols so far.
  This same number can be computed in either the encoder or the decoder, and is
   suitable for making coding decisions.
  Return: The number of bits scaled by 2**OD_BITRES.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
uint32_t od_ec_dec_tell_frac(const od_ec_dec *dec) {
  return od_ec_tell_frac(od_ec_dec_tell(dec), dec->rng);
}
