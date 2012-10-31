/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef POSTPROC_H
#define POSTPROC_H

#define prototype_postproc_inplace(sym)\
  void sym (unsigned char *dst, int pitch, int rows, int cols,int flimit)

#define prototype_postproc(sym)\
  void sym (unsigned char *src, unsigned char *dst, int src_pitch,\
            int dst_pitch, int rows, int cols, int flimit)

#define prototype_postproc_addnoise(sym) \
  void sym (unsigned char *s, char *noise, char blackclamp[16],\
            char whiteclamp[16], char bothclamp[16],\
            unsigned int w, unsigned int h, int pitch)

#define prototype_postproc_blend_mb_inner(sym)\
  void sym (unsigned char *y, unsigned char *u, unsigned char *v,\
            int y1, int u1, int v1, int alpha, int stride)

#define prototype_postproc_blend_mb_outer(sym)\
  void sym (unsigned char *y, unsigned char *u, unsigned char *v,\
            int y1, int u1, int v1, int alpha, int stride)

#define prototype_postproc_blend_b(sym)\
  void sym (unsigned char *y, unsigned char *u, unsigned char *v,\
            int y1, int u1, int v1, int alpha, int stride)

#if ARCH_X86 || ARCH_X86_64
#include "x86/postproc_x86.h"
#endif

#ifndef vp9_postproc_down
#define vp9_postproc_down vp9_mbpost_proc_down_c
#endif
extern prototype_postproc_inplace(vp9_postproc_down);

#ifndef vp9_postproc_across
#define vp9_postproc_across vp9_mbpost_proc_across_ip_c
#endif
extern prototype_postproc_inplace(vp9_postproc_across);

#ifndef vp9_postproc_downacross
#define vp9_postproc_downacross vp9_post_proc_down_and_across_c
#endif
extern prototype_postproc(vp9_postproc_downacross);

#ifndef vp9_postproc_addnoise
#define vp9_postproc_addnoise vp9_plane_add_noise_c
#endif
extern prototype_postproc_addnoise(vp9_postproc_addnoise);

#ifndef vp9_postproc_blend_mb_inner
#define vp9_postproc_blend_mb_inner vp9_blend_mb_inner_c
#endif
extern prototype_postproc_blend_mb_inner(vp9_postproc_blend_mb_inner);

#ifndef vp9_postproc_blend_mb_outer
#define vp9_postproc_blend_mb_outer vp9_blend_mb_outer_c
#endif
extern prototype_postproc_blend_mb_outer(vp9_postproc_blend_mb_outer);

#ifndef vp9_postproc_blend_b
#define vp9_postproc_blend_b vp9_blend_b_c
#endif
extern prototype_postproc_blend_b(vp9_postproc_blend_b);

typedef prototype_postproc((*vp9_postproc_fn_t));
typedef prototype_postproc_inplace((*vp9_postproc_inplace_fn_t));
typedef prototype_postproc_addnoise((*vp9_postproc_addnoise_fn_t));
typedef prototype_postproc_blend_mb_inner((*vp9_postproc_blend_mb_inner_fn_t));
typedef prototype_postproc_blend_mb_outer((*vp9_postproc_blend_mb_outer_fn_t));
typedef prototype_postproc_blend_b((*vp9_postproc_blend_b_fn_t));
typedef struct {
  vp9_postproc_inplace_fn_t           down;
  vp9_postproc_inplace_fn_t           across;
  vp9_postproc_fn_t                   downacross;
  vp9_postproc_addnoise_fn_t          addnoise;
  vp9_postproc_blend_mb_inner_fn_t    blend_mb_inner;
  vp9_postproc_blend_mb_outer_fn_t    blend_mb_outer;
  vp9_postproc_blend_b_fn_t           blend_b;
} vp9_postproc_rtcd_vtable_t;

#if CONFIG_RUNTIME_CPU_DETECT
#define POSTPROC_INVOKE(ctx,fn) (ctx)->fn
#else
#define POSTPROC_INVOKE(ctx,fn) vp9_postproc_##fn
#endif

#include "vpx_ports/mem.h"
struct postproc_state {
  int           last_q;
  int           last_noise;
  char          noise[3072];
  DECLARE_ALIGNED(16, char, blackclamp[16]);
  DECLARE_ALIGNED(16, char, whiteclamp[16]);
  DECLARE_ALIGNED(16, char, bothclamp[16]);
};
#include "onyxc_int.h"
#include "ppflags.h"
int vp9_post_proc_frame(struct VP9Common *oci, YV12_BUFFER_CONFIG *dest,
                        vp9_ppflags_t *flags);


void vp9_de_noise(YV12_BUFFER_CONFIG         *source,
                  YV12_BUFFER_CONFIG         *post,
                  int                         q,
                  int                         low_var_thresh,
                  int                         flag,
                  vp9_postproc_rtcd_vtable_t *rtcd);

void vp9_deblock(YV12_BUFFER_CONFIG         *source,
                 YV12_BUFFER_CONFIG         *post,
                 int                         q,
                 int                         low_var_thresh,
                 int                         flag,
                 vp9_postproc_rtcd_vtable_t *rtcd);
#endif
