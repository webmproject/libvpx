/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP9_COMMON_MIPS_MSA_VP9_MACROS_MSA_H_
#define VP9_COMMON_MIPS_MSA_VP9_MACROS_MSA_H_

#include <msa.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

#if HAVE_MSA
/* load macros */
#define LOAD_UB(psrc) *((const v16u8 *)(psrc))
#define LOAD_SB(psrc) *((const v16i8 *)(psrc))
#define LOAD_UH(psrc) *((const v8u16 *)(psrc))
#define LOAD_SH(psrc) *((const v8i16 *)(psrc))
#define LOAD_UW(psrc) *((const v4u32 *)(psrc))
#define LOAD_SW(psrc) *((const v4i32 *)(psrc))
#define LOAD_UD(psrc) *((const v2u64 *)(psrc))
#define LOAD_SD(psrc) *((const v2i64 *)(psrc))

/* store macros */
#define STORE_UB(vec, pdest) *((v16u8 *)(pdest)) = (vec)
#define STORE_SB(vec, pdest) *((v16i8 *)(pdest)) = (vec)
#define STORE_UH(vec, pdest) *((v8u16 *)(pdest)) = (vec)
#define STORE_SH(vec, pdest) *((v8i16 *)(pdest)) = (vec)
#define STORE_UW(vec, pdest) *((v4u32 *)(pdest)) = (vec)
#define STORE_SW(vec, pdest) *((v4i32 *)(pdest)) = (vec)
#define STORE_UD(vec, pdest) *((v2u64 *)(pdest)) = (vec)
#define STORE_SD(vec, pdest) *((v2i64 *)(pdest)) = (vec)

#if (__mips_isa_rev >= 6)
#define LOAD_WORD(psrc) ({                         \
  const uint8_t *src_m = (const uint8_t *)(psrc);  \
  uint32_t val_m;                                  \
                                                   \
  __asm__ __volatile__ (                           \
      "lw  %[val_m],  %[src_m]  \n\t"              \
                                                   \
      : [val_m] "=r" (val_m)                       \
      : [src_m] "m" (*src_m)                       \
  );                                               \
                                                   \
  val_m;                                           \
})

#if (__mips == 64)
#define LOAD_DWORD(psrc) ({                        \
  const uint8_t *src_m = (const uint8_t *)(psrc);  \
  uint64_t val_m = 0;                              \
                                                   \
  __asm__ __volatile__ (                           \
      "ld  %[val_m],  %[src_m]  \n\t"              \
                                                   \
      : [val_m] "=r" (val_m)                       \
      : [src_m] "m" (*src_m)                       \
  );                                               \
                                                   \
  val_m;                                           \
})
#else  // !(__mips == 64)
#define LOAD_DWORD(psrc) ({                                      \
  const uint8_t *src1_m = (const uint8_t *)(psrc);               \
  const uint8_t *src2_m = ((const uint8_t *)(psrc)) + 4;         \
  uint32_t val0_m, val1_m;                                       \
  uint64_t genval_m = 0;                                         \
                                                                 \
  __asm__ __volatile__ (                                         \
      "lw  %[val0_m],  %[src1_m]  \n\t"                          \
                                                                 \
      : [val0_m] "=r" (val0_m)                                   \
      : [src1_m] "m" (*src1_m)                                   \
  );                                                             \
                                                                 \
  __asm__ __volatile__ (                                         \
      "lw  %[val1_m],  %[src2_m]  \n\t"                          \
                                                                 \
      : [val1_m] "=r" (val1_m)                                   \
      : [src2_m] "m" (*src2_m)                                   \
  );                                                             \
                                                                 \
  genval_m = (uint64_t)(val1_m);                                 \
  genval_m = (uint64_t)((genval_m << 32) & 0xFFFFFFFF00000000);  \
  genval_m = (uint64_t)(genval_m | (uint64_t)val0_m);            \
                                                                 \
  genval_m;                                                      \
})
#endif  // (__mips == 64)
#define STORE_WORD_WITH_OFFSET_1(pdst, val) {    \
  uint8_t *dst_ptr_m = ((uint8_t *)(pdst)) + 1;  \
  const uint32_t val_m = (val);                  \
                                                 \
  __asm__ __volatile__ (                         \
      "sw  %[val_m],  %[dst_ptr_m]  \n\t"        \
                                                 \
      : [dst_ptr_m] "=m" (*dst_ptr_m)            \
      : [val_m] "r" (val_m)                      \
  );                                             \
}

#define STORE_WORD(pdst, val) {            \
  uint8_t *dst_ptr_m = (uint8_t *)(pdst);  \
  const uint32_t val_m = (val);            \
                                           \
  __asm__ __volatile__ (                   \
      "sw  %[val_m],  %[dst_ptr_m]  \n\t"  \
                                           \
      : [dst_ptr_m] "=m" (*dst_ptr_m)      \
      : [val_m] "r" (val_m)                \
  );                                       \
}

#define STORE_DWORD(pdst, val) {           \
  uint8_t *dst_ptr_m = (uint8_t *)(pdst);  \
  const uint64_t val_m = (val);            \
                                           \
  __asm__ __volatile__ (                   \
      "sd  %[val_m],  %[dst_ptr_m]  \n\t"  \
                                           \
      : [dst_ptr_m] "=m" (*dst_ptr_m)      \
      : [val_m] "r" (val_m)                \
  );                                       \
}
#else  // !(__mips_isa_rev >= 6)
#define LOAD_WORD(psrc) ({                         \
  const uint8_t *src_m = (const uint8_t *)(psrc);  \
  uint32_t val_m;                                  \
                                                   \
  __asm__ __volatile__ (                           \
      "ulw  %[val_m],  %[src_m]  \n\t"             \
                                                   \
      : [val_m] "=r" (val_m)                       \
      : [src_m] "m" (*src_m)                       \
  );                                               \
                                                   \
  val_m;                                           \
})

#if (__mips == 64)
#define LOAD_DWORD(psrc) ({                        \
  const uint8_t *src_m = (const uint8_t *)(psrc);  \
  uint64_t val_m = 0;                              \
                                                   \
  __asm__ __volatile__ (                           \
      "uld  %[val_m],  %[src_m]  \n\t"             \
                                                   \
      : [val_m] "=r" (val_m)                       \
      : [src_m] "m" (*src_m)                       \
  );                                               \
                                                   \
  val_m;                                           \
})
#else  // !(__mips == 64)
#define LOAD_DWORD(psrc) ({                                      \
  const uint8_t *src1_m = (const uint8_t *)(psrc);               \
  const uint8_t *src2_m = ((const uint8_t *)(psrc)) + 4;         \
  uint32_t val0_m, val1_m;                                       \
  uint64_t genval_m = 0;                                         \
                                                                 \
  __asm__ __volatile__ (                                         \
      "ulw  %[val0_m],  %[src1_m]  \n\t"                         \
                                                                 \
      : [val0_m] "=r" (val0_m)                                   \
      : [src1_m] "m" (*src1_m)                                   \
  );                                                             \
                                                                 \
  __asm__ __volatile__ (                                         \
      "ulw  %[val1_m],  %[src2_m]  \n\t"                         \
                                                                 \
      : [val1_m] "=r" (val1_m)                                   \
      : [src2_m] "m" (*src2_m)                                   \
  );                                                             \
                                                                 \
  genval_m = (uint64_t)(val1_m);                                 \
  genval_m = (uint64_t)((genval_m << 32) & 0xFFFFFFFF00000000);  \
  genval_m = (uint64_t)(genval_m | (uint64_t)val0_m);            \
                                                                 \
  genval_m;                                                      \
})
#endif  // (__mips == 64)

#define STORE_WORD_WITH_OFFSET_1(pdst, val) {    \
  uint8_t *dst_ptr_m = ((uint8_t *)(pdst)) + 1;  \
  const uint32_t val_m = (val);                  \
                                                 \
  __asm__ __volatile__ (                         \
      "usw  %[val_m],  %[dst_ptr_m]  \n\t"       \
                                                 \
      : [dst_ptr_m] "=m" (*dst_ptr_m)            \
      : [val_m] "r" (val_m)                      \
  );                                             \
}

#define STORE_WORD(pdst, val) {             \
  uint8_t *dst_ptr_m = (uint8_t *)(pdst);   \
  const uint32_t val_m = (val);             \
                                            \
  __asm__ __volatile__ (                    \
      "usw  %[val_m],  %[dst_ptr_m]  \n\t"  \
                                            \
      : [dst_ptr_m] "=m" (*dst_ptr_m)       \
      : [val_m] "r" (val_m)                 \
  );                                        \
}

#define STORE_DWORD(pdst, val) {                            \
  uint8_t *dst1_m = (uint8_t *)(pdst);                      \
  uint8_t *dst2_m = ((uint8_t *)(pdst)) + 4;                \
  uint32_t val0_m, val1_m;                                  \
                                                            \
  val0_m = (uint32_t)((val) & 0x00000000FFFFFFFF);          \
  val1_m = (uint32_t)(((val) >> 32) & 0x00000000FFFFFFFF);  \
                                                            \
  __asm__ __volatile__ (                                    \
      "usw  %[val0_m],  %[dst1_m]  \n\t"                    \
      "usw  %[val1_m],  %[dst2_m]  \n\t"                    \
                                                            \
      : [dst1_m] "=m" (*dst1_m), [dst2_m] "=m" (*dst2_m)    \
      : [val0_m] "r" (val0_m), [val1_m] "r" (val1_m)        \
  );                                                        \
}
#endif  // (__mips_isa_rev >= 6)

#define LOAD_2VECS_UB(psrc, stride,   \
                      val0, val1) {   \
  val0 = LOAD_UB(psrc + 0 * stride);  \
  val1 = LOAD_UB(psrc + 1 * stride);  \
}

#define LOAD_4VECS_UB(psrc, stride,              \
                      val0, val1, val2, val3) {  \
  val0 = LOAD_UB(psrc + 0 * stride);             \
  val1 = LOAD_UB(psrc + 1 * stride);             \
  val2 = LOAD_UB(psrc + 2 * stride);             \
  val3 = LOAD_UB(psrc + 3 * stride);             \
}

#define LOAD_4VECS_SB(psrc, stride,              \
                      val0, val1, val2, val3) {  \
  val0 = LOAD_SB(psrc + 0 * stride);             \
  val1 = LOAD_SB(psrc + 1 * stride);             \
  val2 = LOAD_SB(psrc + 2 * stride);             \
  val3 = LOAD_SB(psrc + 3 * stride);             \
}

#define LOAD_5VECS_UB(psrc, stride,                    \
                      out0, out1, out2, out3, out4) {  \
  LOAD_4VECS_UB((psrc), (stride),                      \
                (out0), (out1), (out2), (out3));       \
  out4 = LOAD_UB(psrc + 4 * stride);                   \
}

#define LOAD_5VECS_SB(psrc, stride,                    \
                      out0, out1, out2, out3, out4) {  \
  LOAD_4VECS_SB((psrc), (stride),                      \
                (out0), (out1), (out2), (out3));       \
  out4 = LOAD_SB(psrc + 4 * stride);                   \
}

#define LOAD_7VECS_SB(psrc, stride,            \
                      val0, val1, val2, val3,  \
                      val4, val5, val6) {      \
  val0 = LOAD_SB((psrc) + 0 * (stride));       \
  val1 = LOAD_SB((psrc) + 1 * (stride));       \
  val2 = LOAD_SB((psrc) + 2 * (stride));       \
  val3 = LOAD_SB((psrc) + 3 * (stride));       \
  val4 = LOAD_SB((psrc) + 4 * (stride));       \
  val5 = LOAD_SB((psrc) + 5 * (stride));       \
  val6 = LOAD_SB((psrc) + 6 * (stride));       \
}

#define LOAD_8VECS_UB(psrc, stride,               \
                      out0, out1, out2, out3,     \
                      out4, out5, out6, out7) {   \
  LOAD_4VECS_UB((psrc), (stride),                 \
                (out0), (out1), (out2), (out3));  \
  LOAD_4VECS_UB((psrc + 4 * stride), (stride),    \
                (out4), (out5), (out6), (out7));  \
}

#define LOAD_8VECS_SB(psrc, stride,               \
                      out0, out1, out2, out3,     \
                      out4, out5, out6, out7) {   \
  LOAD_4VECS_SB((psrc), (stride),                 \
                (out0), (out1), (out2), (out3));  \
  LOAD_4VECS_SB((psrc + 4 * stride), (stride),    \
                (out4), (out5), (out6), (out7));  \
}

#define STORE_4VECS_UB(dst_out, pitch,         \
                       in0, in1, in2, in3) {   \
  STORE_UB((in0), (dst_out));                  \
  STORE_UB((in1), ((dst_out) + (pitch)));      \
  STORE_UB((in2), ((dst_out) + 2 * (pitch)));  \
  STORE_UB((in3), ((dst_out) + 3 * (pitch)));  \
}

#define STORE_8VECS_UB(dst_out, pitch_in,               \
                       in0, in1, in2, in3,              \
                       in4, in5, in6, in7) {            \
  STORE_4VECS_UB(dst_out, pitch_in,                     \
                 in0, in1, in2, in3);                   \
  STORE_4VECS_UB((dst_out + 4 * (pitch_in)), pitch_in,  \
                 in4, in5, in6, in7);                   \
}

#define VEC_INSERT_4W_UB(src, src0, src1, src2, src3) {  \
  src = (v16u8)__msa_insert_w((v4i32)(src), 0, (src0));  \
  src = (v16u8)__msa_insert_w((v4i32)(src), 1, (src1));  \
  src = (v16u8)__msa_insert_w((v4i32)(src), 2, (src2));  \
  src = (v16u8)__msa_insert_w((v4i32)(src), 3, (src3));  \
}

#define VEC_INSERT_2DW_UB(src, src0, src1) {             \
  src = (v16u8)__msa_insert_d((v2i64)(src), 0, (src0));  \
  src = (v16u8)__msa_insert_d((v2i64)(src), 1, (src1));  \
}

/* interleave macros */
/* no in-place support */
#define ILV_B_LRLR_UB(in0, in1, in2, in3,                  \
                      out0, out1, out2, out3) {            \
  out0 = (v16u8)__msa_ilvl_b((v16i8)(in1), (v16i8)(in0));  \
  out1 = (v16u8)__msa_ilvr_b((v16i8)(in1), (v16i8)(in0));  \
  out2 = (v16u8)__msa_ilvl_b((v16i8)(in3), (v16i8)(in2));  \
  out3 = (v16u8)__msa_ilvr_b((v16i8)(in3), (v16i8)(in2));  \
}

#define ILVR_B_2VECS_UB(in0_r, in1_r, in0_l, in1_l,            \
                        out0, out1) {                          \
  out0 = (v16u8)__msa_ilvr_b((v16i8)(in0_l), (v16i8)(in0_r));  \
  out1 = (v16u8)__msa_ilvr_b((v16i8)(in1_l), (v16i8)(in1_r));  \
}

#define ILVR_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,     \
                        out0, out1) {                   \
  out0 = __msa_ilvr_b((v16i8)(in0_l), (v16i8)(in0_r));  \
  out1 = __msa_ilvr_b((v16i8)(in1_l), (v16i8)(in1_r));  \
}

#define ILVR_B_4VECS_UB(in0_r, in1_r, in2_r, in3_r,  \
                        in0_l, in1_l, in2_l, in3_l,  \
                        out0, out1, out2, out3) {    \
  ILVR_B_2VECS_UB(in0_r, in1_r, in0_l, in1_l,        \
                  out0, out1);                       \
  ILVR_B_2VECS_UB(in2_r, in3_r, in2_l, in3_l,        \
                  out2, out3);                       \
}

#define ILVR_B_4VECS_SB(in0_r, in1_r, in2_r, in3_r,  \
                        in0_l, in1_l, in2_l, in3_l,  \
                        out0, out1, out2, out3) {    \
  ILVR_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,        \
                  out0, out1);                       \
  ILVR_B_2VECS_SB(in2_r, in3_r, in2_l, in3_l,        \
                  out2, out3);                       \
}

#define ILVR_B_6VECS_SB(in0_r, in1_r, in2_r,   \
                        in3_r, in4_r, in5_r,   \
                        in0_l, in1_l, in2_l,   \
                        in3_l, in4_l, in5_l,   \
                        out0, out1, out2,      \
                        out3, out4, out5) {    \
  ILVR_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,  \
                  out0, out1);                 \
  ILVR_B_2VECS_SB(in2_r, in3_r, in2_l, in3_l,  \
                  out2, out3);                 \
  ILVR_B_2VECS_SB(in4_r, in5_r, in4_l, in5_l,  \
                  out4, out5);                 \
}

#define ILVR_B_8VECS_SB(in0_r, in1_r, in2_r, in3_r,  \
                        in4_r, in5_r, in6_r, in7_r,  \
                        in0_l, in1_l, in2_l, in3_l,  \
                        in4_l, in5_l, in6_l, in7_l,  \
                        out0, out1, out2, out3,      \
                        out4, out5, out6, out7) {    \
  ILVR_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,        \
                  out0, out1);                       \
  ILVR_B_2VECS_SB(in2_r, in3_r, in2_l, in3_l,        \
                  out2, out3);                       \
  ILVR_B_2VECS_SB(in4_r, in5_r, in4_l, in5_l,        \
                  out4, out5);                       \
  ILVR_B_2VECS_SB(in6_r, in7_r, in6_l, in7_l,        \
                  out6, out7);                       \
}

#define ILVL_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,     \
                        out0, out1) {                   \
  out0 = __msa_ilvl_b((v16i8)(in0_l), (v16i8)(in0_r));  \
  out1 = __msa_ilvl_b((v16i8)(in1_l), (v16i8)(in1_r));  \
}

#define ILVL_B_4VECS_SB(in0_r, in1_r, in2_r, in3_r,  \
                        in0_l, in1_l, in2_l, in3_l,  \
                        out0, out1, out2, out3) {    \
  ILVL_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,        \
                  out0, out1);                       \
  ILVL_B_2VECS_SB(in2_r, in3_r, in2_l, in3_l,        \
                  out2, out3);                       \
}

#define ILVL_B_6VECS_SB(in0_r, in1_r, in2_r,   \
                        in3_r, in4_r, in5_r,   \
                        in0_l, in1_l, in2_l,   \
                        in3_l, in4_l, in5_l,   \
                        out0, out1, out2,      \
                        out3, out4, out5) {    \
  ILVL_B_2VECS_SB(in0_r, in1_r, in0_l, in1_l,  \
                  out0, out1);                 \
  ILVL_B_2VECS_SB(in2_r, in3_r, in2_l, in3_l,  \
                  out2, out3);                 \
  ILVL_B_2VECS_SB(in4_r, in5_r, in4_l, in5_l,  \
                  out4, out5);                 \
}

#define ILVR_D_2VECS_SB(out0, in0_l, in0_r,                    \
                        out1, in1_l, in1_r) {                  \
  out0 = (v16i8)__msa_ilvr_d((v2i64)(in0_l), (v2i64)(in0_r));  \
  out1 = (v16i8)__msa_ilvr_d((v2i64)(in1_l), (v2i64)(in1_r));  \
}

#define ILVR_D_3VECS_SB(out0, in0_l, in0_r,                    \
                        out1, in1_l, in1_r,                    \
                        out2, in2_l, in2_r) {                  \
  ILVR_D_2VECS_SB(out0, in0_l, in0_r,                          \
                  out1, in1_l, in1_r);                         \
  out2 = (v16i8)__msa_ilvr_d((v2i64)(in2_l), (v2i64)(in2_r));  \
}

#define ILVR_D_4VECS_SB(out0, in0_l, in0_r,    \
                        out1, in1_l, in1_r,    \
                        out2, in2_l, in2_r,    \
                        out3, in3_l, in3_r) {  \
  ILVR_D_2VECS_SB(out0, in0_l, in0_r,          \
                  out1, in1_l, in1_r);         \
  ILVR_D_2VECS_SB(out2, in2_l, in2_r,          \
                  out3, in3_l, in3_r);         \
}

#define XORI_B_2VECS_UB(val0, val1,               \
                        out0, out1, xor_val) {    \
  out0 = __msa_xori_b((v16u8)(val0), (xor_val));  \
  out1 = __msa_xori_b((v16u8)(val1), (xor_val));  \
}

#define XORI_B_2VECS_SB(val0, val1,                      \
                        out0, out1, xor_val) {           \
  out0 = (v16i8)__msa_xori_b((v16u8)(val0), (xor_val));  \
  out1 = (v16i8)__msa_xori_b((v16u8)(val1), (xor_val));  \
}

#define XORI_B_3VECS_SB(val0, val1, val2,                \
                        out0, out1, out2, xor_val) {     \
  XORI_B_2VECS_SB(val0, val1,  out0, out1, xor_val);     \
  out2 = (v16i8)__msa_xori_b((v16u8)(val2), (xor_val));  \
}

#define XORI_B_4VECS_UB(val0, val1, val2, val3,      \
                        out0, out1, out2, out3,      \
                        xor_val) {                   \
  XORI_B_2VECS_UB(val0, val1, out0, out1, xor_val);  \
  XORI_B_2VECS_UB(val2, val3, out2, out3, xor_val);  \
}

#define XORI_B_4VECS_SB(val0, val1, val2, val3,      \
                        out0, out1, out2, out3,      \
                        xor_val) {                   \
  XORI_B_2VECS_SB(val0, val1, out0, out1, xor_val);  \
  XORI_B_2VECS_SB(val2, val3, out2, out3, xor_val);  \
}

#define XORI_B_7VECS_SB(val0, val1, val2, val3,      \
                        val4, val5, val6,            \
                        out0, out1, out2, out3,      \
                        out4, out5, out6,            \
                        xor_val) {                   \
  XORI_B_4VECS_SB(val0, val1, val2, val3,            \
                  out0, out1, out2, out3, xor_val);  \
  XORI_B_3VECS_SB(val4, val5, val6,                  \
                  out4, out5, out6, xor_val);        \
}

#define SRARI_H_4VECS_UH(val0, val1, val2, val3,                  \
                         out0, out1, out2, out3,                  \
                         shift_right_val) {                       \
  out0 = (v8u16)__msa_srari_h((v8i16)(val0), (shift_right_val));  \
  out1 = (v8u16)__msa_srari_h((v8i16)(val1), (shift_right_val));  \
  out2 = (v8u16)__msa_srari_h((v8i16)(val2), (shift_right_val));  \
  out3 = (v8u16)__msa_srari_h((v8i16)(val3), (shift_right_val));  \
}

#define SRARI_SATURATE_UNSIGNED_H(input, right_shift_val, sat_val) ({  \
  v8u16 out_m;                                                         \
                                                                       \
  out_m = (v8u16)__msa_srari_h((v8i16)(input), (right_shift_val));     \
  out_m = __msa_sat_u_h(out_m, (sat_val));                             \
  out_m;                                                               \
})

#define SRARI_SATURATE_SIGNED_H(input, right_shift_val, sat_val) ({  \
  v8i16 out_m;                                                       \
                                                                     \
  out_m = __msa_srari_h((v8i16)(input), (right_shift_val));          \
  out_m = __msa_sat_s_h(out_m, (sat_val));                           \
  out_m;                                                             \
})

#define PCKEV_2B_XORI128_STORE_4_BYTES_4(in1, in2,        \
                                         pdst, stride) {  \
  uint32_t out0_m, out1_m, out2_m, out3_m;                \
  v16i8 tmp0_m;                                           \
  uint8_t *dst_m = (uint8_t *)(pdst);                     \
                                                          \
  tmp0_m = __msa_pckev_b((v16i8)(in2), (v16i8)(in1));     \
  tmp0_m = (v16i8)__msa_xori_b((v16u8)tmp0_m, 128);       \
                                                          \
  out0_m = __msa_copy_u_w((v4i32)tmp0_m, 0);              \
  out1_m = __msa_copy_u_w((v4i32)tmp0_m, 1);              \
  out2_m = __msa_copy_u_w((v4i32)tmp0_m, 2);              \
  out3_m = __msa_copy_u_w((v4i32)tmp0_m, 3);              \
                                                          \
  STORE_WORD(dst_m, out0_m);                              \
  dst_m += stride;                                        \
  STORE_WORD(dst_m, out1_m);                              \
  dst_m += stride;                                        \
  STORE_WORD(dst_m, out2_m);                              \
  dst_m += stride;                                        \
  STORE_WORD(dst_m, out3_m);                              \
}

#define PCKEV_B_4_XORI128_STORE_8_BYTES_4(in1, in2,        \
                                          in3, in4,        \
                                          pdst, stride) {  \
  uint64_t out0_m, out1_m, out2_m, out3_m;                 \
  v16i8 tmp0_m, tmp1_m;                                    \
  uint8_t *dst_m = (uint8_t *)(pdst);                      \
                                                           \
  tmp0_m = __msa_pckev_b((v16i8)(in2), (v16i8)(in1));      \
  tmp1_m = __msa_pckev_b((v16i8)(in4), (v16i8)(in3));      \
                                                           \
  tmp0_m = (v16i8)__msa_xori_b((v16u8)tmp0_m, 128);        \
  tmp1_m = (v16i8)__msa_xori_b((v16u8)tmp1_m, 128);        \
                                                           \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);               \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);               \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);               \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);               \
                                                           \
  STORE_DWORD(dst_m, out0_m);                              \
  dst_m += stride;                                         \
  STORE_DWORD(dst_m, out1_m);                              \
  dst_m += stride;                                         \
  STORE_DWORD(dst_m, out2_m);                              \
  dst_m += stride;                                         \
  STORE_DWORD(dst_m, out3_m);                              \
}

/* Only for signed vecs */
#define PCKEV_B_XORI128_STORE_VEC(in1, in2, pdest) {  \
  v16i8 tmp_m;                                        \
                                                      \
  tmp_m = __msa_pckev_b((v16i8)(in1), (v16i8)(in2));  \
  tmp_m = (v16i8)__msa_xori_b((v16u8)tmp_m, 128);     \
  STORE_SB(tmp_m, (pdest));                           \
}

/* Only for signed vecs */
#define PCKEV_B_4_XORI128_AVG_STORE_8_BYTES_4(in1, dst0,       \
                                              in2, dst1,       \
                                              in3, dst2,       \
                                              in4, dst3,       \
                                              pdst, stride) {  \
  uint64_t out0_m, out1_m, out2_m, out3_m;                     \
  v16u8 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                        \
  uint8_t *dst_m = (uint8_t *)(pdst);                          \
                                                               \
  tmp0_m = (v16u8)__msa_pckev_b((v16i8)(in2), (v16i8)(in1));   \
  tmp1_m = (v16u8)__msa_pckev_b((v16i8)(in4), (v16i8)(in3));   \
                                                               \
  tmp2_m = (v16u8)__msa_ilvr_d((v2i64)(dst1), (v2i64)(dst0));  \
  tmp3_m = (v16u8)__msa_ilvr_d((v2i64)(dst3), (v2i64)(dst2));  \
                                                               \
  tmp0_m = __msa_xori_b(tmp0_m, 128);                          \
  tmp1_m = __msa_xori_b(tmp1_m, 128);                          \
                                                               \
  tmp0_m = __msa_aver_u_b(tmp0_m, tmp2_m);                     \
  tmp1_m = __msa_aver_u_b(tmp1_m, tmp3_m);                     \
                                                               \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);                   \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);                   \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);                   \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);                   \
                                                               \
  STORE_DWORD(dst_m, out0_m);                                  \
  dst_m += stride;                                             \
  STORE_DWORD(dst_m, out1_m);                                  \
  dst_m += stride;                                             \
  STORE_DWORD(dst_m, out2_m);                                  \
  dst_m += stride;                                             \
  STORE_DWORD(dst_m, out3_m);                                  \
}

/* Only for signed vecs */
#define PCKEV_B_XORI128_AVG_STORE_VEC(in1, in2, dst, pdest) {  \
  v16u8 tmp_m;                                                 \
                                                               \
  tmp_m = (v16u8)__msa_pckev_b((v16i8)(in1), (v16i8)(in2));    \
  tmp_m = __msa_xori_b(tmp_m, 128);                            \
  tmp_m = __msa_aver_u_b(tmp_m, (v16u8)(dst));                 \
  STORE_UB(tmp_m, (pdest));                                    \
}

#define PCKEV_B_STORE_8_BYTES_4(in1, in2, in3, in4,    \
                                pdst, stride) {        \
  uint64_t out0_m, out1_m, out2_m, out3_m;             \
  v16i8 tmp0_m, tmp1_m;                                \
  uint8_t *dst_m = (uint8_t *)(pdst);                  \
                                                       \
  tmp0_m = __msa_pckev_b((v16i8)(in2), (v16i8)(in1));  \
  tmp1_m = __msa_pckev_b((v16i8)(in4), (v16i8)(in3));  \
                                                       \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);           \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);           \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);           \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);           \
                                                       \
  STORE_DWORD(dst_m, out0_m);                          \
  dst_m += stride;                                     \
  STORE_DWORD(dst_m, out1_m);                          \
  dst_m += stride;                                     \
  STORE_DWORD(dst_m, out2_m);                          \
  dst_m += stride;                                     \
  STORE_DWORD(dst_m, out3_m);                          \
}

/* Only for unsigned vecs */
#define PCKEV_B_STORE_VEC(in1, in2, pdest) {          \
  v16i8 tmp_m;                                        \
                                                      \
  tmp_m = __msa_pckev_b((v16i8)(in1), (v16i8)(in2));  \
  STORE_SB(tmp_m, (pdest));                           \
}

#define PCKEV_B_AVG_STORE_8_BYTES_4(in1, dst0, in2, dst1,       \
                                    in3, dst2, in4, dst3,       \
                                    pdst, stride) {             \
  uint64_t out0_m, out1_m, out2_m, out3_m;                      \
  v16u8 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                         \
  uint8_t *dst_m = (uint8_t *)(pdst);                           \
                                                                \
  tmp0_m = (v16u8)__msa_pckev_b((v16i8)(in2), (v16i8)(in1));    \
  tmp1_m = (v16u8)__msa_pckev_b((v16i8)(in4), (v16i8)(in3));    \
                                                                \
  tmp2_m = (v16u8)__msa_pckev_d((v2i64)(dst1), (v2i64)(dst0));  \
  tmp3_m = (v16u8)__msa_pckev_d((v2i64)(dst3), (v2i64)(dst2));  \
                                                                \
  tmp0_m = __msa_aver_u_b(tmp0_m, tmp2_m);                      \
  tmp1_m = __msa_aver_u_b(tmp1_m, tmp3_m);                      \
                                                                \
  out0_m = __msa_copy_u_d((v2i64)tmp0_m, 0);                    \
  out1_m = __msa_copy_u_d((v2i64)tmp0_m, 1);                    \
  out2_m = __msa_copy_u_d((v2i64)tmp1_m, 0);                    \
  out3_m = __msa_copy_u_d((v2i64)tmp1_m, 1);                    \
                                                                \
  STORE_DWORD(dst_m, out0_m);                                   \
  dst_m += stride;                                              \
  STORE_DWORD(dst_m, out1_m);                                   \
  dst_m += stride;                                              \
  STORE_DWORD(dst_m, out2_m);                                   \
  dst_m += stride;                                              \
  STORE_DWORD(dst_m, out3_m);                                   \
}

#define PCKEV_B_AVG_STORE_VEC(in1, in2, dst, pdest) {        \
  v16u8 tmp_m;                                               \
                                                             \
  tmp_m = (v16u8)__msa_pckev_b((v16i8)(in1), (v16i8)(in2));  \
  tmp_m = __msa_aver_u_b(tmp_m, (v16u8)(dst));               \
  STORE_UB(tmp_m, (pdest));                                  \
}
#endif  /* HAVE_MSA */
#endif  /* VP9_COMMON_MIPS_MSA_VP9_MACROS_MSA_H_ */
