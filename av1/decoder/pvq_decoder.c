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

/* clang-format off */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include "./aom_config.h"
#include "aom_dsp/entcode.h"
#include "aom_dsp/entdec.h"
#include "av1/common/odintrin.h"
#include "av1/common/partition.h"
#include "av1/common/pvq_state.h"
#include "av1/decoder/decint.h"
#include "av1/decoder/pvq_decoder.h"

static void od_decode_pvq_codeword(od_ec_dec *ec, od_pvq_codeword_ctx *ctx,
 od_coeff *y, int n, int k) {
  int i;
  od_decode_band_pvq_splits(ec, ctx, y, n, k, 0);
  for (i = 0; i < n; i++) {
    if (y[i] && od_ec_dec_bits(ec, 1, "pvq:sign")) y[i] = -y[i];
  }
}

/** Inverse of neg_interleave; decodes the interleaved gain.
 *
 * @param [in]      x      quantized/interleaved gain to decode
 * @param [in]      ref    quantized gain of the reference
 * @return                 original quantized gain value
 */
static int neg_deinterleave(int x, int ref) {
  if (x < 2*ref-1) {
    if (x & 1) return ref - 1 - (x >> 1);
    else return ref + (x >> 1);
  }
  else return x+1;
}

/** Synthesizes one parition of coefficient values from a PVQ-encoded
 * vector.
 *
 * @param [out]     xcoeff  output coefficient partition (x in math doc)
 * @param [in]      ypulse  PVQ-encoded values (y in math doc); in the noref
 *                          case, this vector has n entries, in the
 *                          reference case it contains n-1 entries
 *                          (the m-th entry is not included)
 * @param [in]      ref     reference vector (prediction)
 * @param [in]      n       number of elements in this partition
 * @param [in]      gr      gain of the reference vector (prediction)
 * @param [in]      noref   indicates presence or lack of prediction
 * @param [in]      g       decoded quantized vector gain
 * @param [in]      theta   decoded theta (prediction error)
 * @param [in]      qm      QM with magnitude compensation
 * @param [in]      qm_inv  Inverse of QM with magnitude compensation
 */
static void pvq_synthesis(od_coeff *xcoeff, od_coeff *ypulse, od_val16 *r16,
 int n, od_val32 gr, int noref, od_val32 g, od_val32 theta, const int16_t *qm_inv,
 int shift) {
  int s;
  int m;
  /* Sign of the Householder reflection vector */
  s = 0;
  /* Direction of the Householder reflection vector */
  m = noref ? 0 : od_compute_householder(r16, n, gr, &s, shift);
  od_pvq_synthesis_partial(xcoeff, ypulse, r16, n, noref, g, theta, m, s,
   qm_inv);
}

typedef struct {
  od_coeff *ref;
  int nb_coeffs;
  int allow_flip;
} cfl_ctx;

/** Decodes a single vector of integers (eg, a partition within a
 *  coefficient block) encoded using PVQ
 *
 * @param [in,out] ec      range encoder
 * @param [in]     q0      scale/quantizer
 * @param [in]     n       number of coefficients in partition
 * @param [in,out] model   entropy decoder state
 * @param [in,out] adapt   adaptation context
 * @param [in,out] exg     ExQ16 expectation of decoded gain value
 * @param [in,out] ext     ExQ16 expectation of decoded theta value
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [out]    noref   boolean indicating absence of reference
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     robust  stream is robust to error in the reference
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [in]     pli     plane index
 * @param [in]     cdf_ctx selects which cdf context to use
 * @param [in,out] skip_rest whether to skip further bands in each direction
 * @param [in]     band    index of the band being decoded
 * @param [in]     band    index of the band being decoded
 * @param [out]    skip    skip flag with range [0,1]
 * @param [in]     qm      QM with magnitude compensation
 * @param [in]     qm_inv  Inverse of QM with magnitude compensation
 */
static void pvq_decode_partition(od_ec_dec *ec,
                                 int q0,
                                 int n,
                                 generic_encoder model[3],
                                 od_adapt_ctx *adapt,
                                 int *exg,
                                 int *ext,
                                 od_coeff *ref,
                                 od_coeff *out,
                                 int *noref,
                                 od_val16 beta,
                                 int robust,
                                 int is_keyframe,
                                 int pli,
                                 int cdf_ctx,
                                 cfl_ctx *cfl,
                                 int has_skip,
                                 int *skip_rest,
                                 int band,
                                 int *skip,
                                 const int16_t *qm,
                                 const int16_t *qm_inv) {
  int k;
  od_val32 qcg;
  int max_theta;
  int itheta;
  od_val32 theta;
  od_val32 gr;
  od_val32 gain_offset;
  od_coeff y[MAXN];
  int qg;
  int nodesync;
  int id;
  int i;
  od_val16 ref16[MAXN];
  int rshift;
  theta = 0;
  gr = 0;
  gain_offset = 0;
  /* We always use the robust bitstream for keyframes to avoid having
     PVQ and entropy decoding depending on each other, hurting parallelism. */
  nodesync = robust || is_keyframe;
  /* Skip is per-direction. For band=0, we can use any of the flags. */
  if (skip_rest[(band + 2) % 3]) {
    qg = 0;
    if (is_keyframe) {
      itheta = -1;
      *noref = 1;
    }
    else {
      itheta = 0;
      *noref = 0;
    }
  }
  else {
    /* Jointly decode gain, itheta and noref for small values. Then we handle
       larger gain. We need to wait for itheta because in the !nodesync case
       it depends on max_theta, which depends on the gain. */
    id = od_decode_cdf_adapt(ec, &adapt->pvq.pvq_gaintheta_cdf[cdf_ctx][0],
     8 + 7*has_skip, adapt->pvq.pvq_gaintheta_increment,
     "pvq:gaintheta");
    if (!is_keyframe && id >= 10) id++;
    if (is_keyframe && id >= 8) id++;
    if (id >= 8) {
      id -= 8;
      skip_rest[0] = skip_rest[1] = skip_rest[2] = 1;
    }
    qg = id & 1;
    itheta = (id >> 1) - 1;
    *noref = (itheta == -1);
  }
  /* The CfL flip bit is only decoded on the first band that has noref=0. */
  if (cfl->allow_flip && !*noref) {
    int flip;
    flip = od_ec_dec_bits(ec, 1, "cfl:flip");
    if (flip) {
      for (i = 0; i < cfl->nb_coeffs; i++) cfl->ref[i] = -cfl->ref[i];
    }
    cfl->allow_flip = 0;
  }
  if (qg > 0) {
    int tmp;
    tmp = *exg;
    qg = 1 + generic_decode(ec, &model[!*noref], -1, &tmp, 2, "pvq:gain");
    OD_IIR_DIADIC(*exg, qg << 16, 2);
  }
  *skip = 0;
#if defined(OD_FLOAT_PVQ)
  rshift = 0;
#else
  /* Shift needed to make the reference fit in 15 bits, so that the Householder
     vector can fit in 16 bits. */
  rshift = OD_MAXI(0, od_vector_log_mag(ref, n) - 14);
#endif
  for (i = 0; i < n; i++) {
#if defined(OD_FLOAT_PVQ)
    ref16[i] = ref[i]*(double)qm[i]*OD_QM_SCALE_1;
#else
    ref16[i] = OD_SHR_ROUND(ref[i]*qm[i], OD_QM_SHIFT + rshift);
#endif
  }
  if(!*noref){
    /* we have a reference; compute its gain */
    od_val32 cgr;
    int icgr;
    int cfl_enabled;
    cfl_enabled = pli != 0 && is_keyframe && !OD_DISABLE_CFL;
    cgr = od_pvq_compute_gain(ref16, n, q0, &gr, beta, rshift);
    if (cfl_enabled) cgr = OD_CGAIN_SCALE;
#if defined(OD_FLOAT_PVQ)
    icgr = (int)floor(.5 + cgr);
#else
    icgr = OD_SHR_ROUND(cgr, OD_CGAIN_SHIFT);
#endif
    /* quantized gain is interleave encoded when there's a reference;
       deinterleave it now */
    if (is_keyframe) qg = neg_deinterleave(qg, icgr);
    else {
      qg = neg_deinterleave(qg, icgr + 1) - 1;
      if (qg == 0) *skip = (icgr ? OD_PVQ_SKIP_ZERO : OD_PVQ_SKIP_COPY);
    }
    if (qg == icgr && itheta == 0 && !cfl_enabled) *skip = OD_PVQ_SKIP_COPY;
    gain_offset = cgr - OD_SHL(icgr, OD_CGAIN_SHIFT);
    qcg = OD_SHL(qg, OD_CGAIN_SHIFT) + gain_offset;
    /* read and decode first-stage PVQ error theta */
    max_theta = od_pvq_compute_max_theta(qcg, beta);
    if (itheta > 1 && (nodesync || max_theta > 3)) {
      int tmp;
      tmp = *ext;
      itheta = 2 + generic_decode(ec, &model[2], nodesync ? -1 : max_theta - 3,
       &tmp, 2, "pvq:theta");
      OD_IIR_DIADIC(*ext, itheta << 16, 2);
    }
    theta = od_pvq_compute_theta(itheta, max_theta);
  }
  else{
    itheta = 0;
    if (!is_keyframe) qg++;
    qcg = OD_SHL(qg, OD_CGAIN_SHIFT);
    if (qg == 0) *skip = OD_PVQ_SKIP_ZERO;
  }

  k = od_pvq_compute_k(qcg, itheta, theta, *noref, n, beta, nodesync);
  if (k != 0) {
    /* when noref==0, y is actually size n-1 */
    od_decode_pvq_codeword(ec, &adapt->pvq.pvq_codeword_ctx, y, n - !*noref,
     k);
  }
  else {
    OD_CLEAR(y, n);
  }
  if (*skip) {
    if (*skip == OD_PVQ_SKIP_COPY) OD_COPY(out, ref, n);
    else OD_CLEAR(out, n);
  }
  else {
    od_val32 g;
    g = od_gain_expand(qcg, q0, beta);
    pvq_synthesis(out, y, ref16, n, gr, *noref, g, theta, qm_inv, rshift);
  }
  *skip = !!*skip;
}

/** Decodes a coefficient block (except for DC) encoded using PVQ
 *
 * @param [in,out] dec     daala decoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [in]     q0      quantizer
 * @param [in]     pli     plane index
 * @param [in]     bs      log of the block size minus two
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     robust  stream is robust to error in the reference
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [out]    flags   bitmask of the per band skip and noref flags
 * @param [in]     block_skip skip flag for the block (range 0-3)
 * @param [in]     qm      QM with magnitude compensation
 * @param [in]     qm_inv  Inverse of QM with magnitude compensation
 */
void od_pvq_decode(daala_dec_ctx *dec,
                   od_coeff *ref,
                   od_coeff *out,
                   int q0,
                   int pli,
                   int bs,
                   const od_val16 *beta,
                   int robust,
                   int is_keyframe,
                   unsigned int *flags,
                   int block_skip,
                   const int16_t *qm,
                   const int16_t *qm_inv){

  int noref[PVQ_MAX_PARTITIONS];
  int skip[PVQ_MAX_PARTITIONS];
  int *exg;
  int *ext;
  int nb_bands;
  int i;
  const int *off;
  int size[PVQ_MAX_PARTITIONS];
  generic_encoder *model;
  int skip_rest[3] = {0};
  cfl_ctx cfl;
  /* const unsigned char *pvq_qm; */
  /*Default to skip=1 and noref=0 for all bands.*/
  for (i = 0; i < PVQ_MAX_PARTITIONS; i++) {
    noref[i] = 0;
    skip[i] = 1;
  }
  /* TODO(yushin): Enable this for activity masking,
     when pvq_qm_q4 is available in AOM. */
  /*pvq_qm = &dec->state.pvq_qm_q4[pli][0];*/
  exg = &dec->state.adapt.pvq.pvq_exg[pli][bs][0];
  ext = dec->state.adapt.pvq.pvq_ext + bs*PVQ_MAX_PARTITIONS;
  model = dec->state.adapt.pvq.pvq_param_model;
  nb_bands = OD_BAND_OFFSETS[bs][0];
  off = &OD_BAND_OFFSETS[bs][1];
  OD_ASSERT(block_skip < 4);
  out[0] = block_skip & 1;
  if (!(block_skip >> 1)) {
    if (is_keyframe) for (i = 1; i < 1 << (2*bs + 4); i++) out[i] = 0;
    else for (i = 1; i < 1 << (2*bs + 4); i++) out[i] = ref[i];
  }
  else {
    for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
    cfl.ref = ref;
    cfl.nb_coeffs = off[nb_bands];
    cfl.allow_flip = pli != 0 && is_keyframe;
    for (i = 0; i < nb_bands; i++) {
      int q;
      /* TODO(yushin): Enable this for activity masking,
         when pvq_qm_q4 is available in AOM. */
      /*q = OD_MAXI(1, q0*pvq_qm[od_qm_get_index(bs, i + 1)] >> 4);*/
      q = OD_MAXI(1, q0);
      pvq_decode_partition(dec->ec, q, size[i],
       model, &dec->state.adapt, exg + i, ext + i, ref + off[i], out + off[i],
       &noref[i], beta[i], robust, is_keyframe, pli,
       (pli != 0)*OD_NBSIZES*PVQ_MAX_PARTITIONS + bs*PVQ_MAX_PARTITIONS + i,
       &cfl, i == 0 && (i < nb_bands - 1), skip_rest, i, &skip[i],
       qm + off[i], qm_inv + off[i]);
      if (i == 0 && !skip_rest[0] && bs > 0) {
        int skip_dir;
        int j;
        skip_dir = od_decode_cdf_adapt(dec->ec,
         &dec->state.adapt.pvq.pvq_skip_dir_cdf[(pli != 0) + 2*(bs - 1)][0], 7,
         dec->state.adapt.pvq.pvq_skip_dir_increment, "pvq:skiprest");
        for (j = 0; j < 3; j++) skip_rest[j] = !!(skip_dir & (1 << j));
      }
    }
  }
  *flags = 0;
  for (i = nb_bands - 1; i >= 0; i--) {
    *flags <<= 1;
    *flags |= noref[i]&1;
    *flags <<= 1;
    *flags |= skip[i]&1;
  }
}
