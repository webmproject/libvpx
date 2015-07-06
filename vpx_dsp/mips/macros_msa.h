/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_MIPS_MACROS_MSA_H_
#define VPX_DSP_MIPS_MACROS_MSA_H_

#include <msa.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

#define LD_B(RTYPE, psrc) *((const RTYPE *)(psrc))
#define LD_UB(...) LD_B(v16u8, __VA_ARGS__)
#define LD_SB(...) LD_B(v16i8, __VA_ARGS__)

#define LD_H(RTYPE, psrc) *((const RTYPE *)(psrc))
#define LD_UH(...) LD_H(v8u16, __VA_ARGS__)
#define LD_SH(...) LD_H(v8i16, __VA_ARGS__)

#define ST_H(RTYPE, in, pdst) *((RTYPE *)(pdst)) = (in)
#define ST_SH(...) ST_H(v8i16, __VA_ARGS__)

#if (__mips_isa_rev >= 6)
#define LW(psrc) ({                                 \
  const uint8_t *psrc_m = (const uint8_t *)(psrc);  \
  uint32_t val_m;                                   \
                                                    \
  __asm__ __volatile__ (                            \
      "lw  %[val_m],  %[psrc_m]  \n\t"              \
                                                    \
      : [val_m] "=r" (val_m)                        \
      : [psrc_m] "m" (*psrc_m)                      \
  );                                                \
                                                    \
  val_m;                                            \
})

#if (__mips == 64)
#define LD(psrc) ({                                 \
  const uint8_t *psrc_m = (const uint8_t *)(psrc);  \
  uint64_t val_m = 0;                               \
                                                    \
  __asm__ __volatile__ (                            \
      "ld  %[val_m],  %[psrc_m]  \n\t"              \
                                                    \
      : [val_m] "=r" (val_m)                        \
      : [psrc_m] "m" (*psrc_m)                      \
  );                                                \
                                                    \
  val_m;                                            \
})
#else  // !(__mips == 64)
#define LD(psrc) ({                                        \
  const uint8_t *psrc_m = (const uint8_t *)(psrc);         \
  uint32_t val0_m, val1_m;                                 \
  uint64_t val_m = 0;                                      \
                                                           \
  val0_m = LW(psrc_m);                                     \
  val1_m = LW(psrc_m + 4);                                 \
                                                           \
  val_m = (uint64_t)(val1_m);                              \
  val_m = (uint64_t)((val_m << 32) & 0xFFFFFFFF00000000);  \
  val_m = (uint64_t)(val_m | (uint64_t)val0_m);            \
                                                           \
  val_m;                                                   \
})
#endif  // (__mips == 64)

#define SW(val, pdst) {                 \
  uint8_t *pdst_m = (uint8_t *)(pdst);  \
  const uint32_t val_m = (val);         \
                                        \
  __asm__ __volatile__ (                \
      "sw  %[val_m],  %[pdst_m]  \n\t"  \
                                        \
      : [pdst_m] "=m" (*pdst_m)         \
      : [val_m] "r" (val_m)             \
  );                                    \
}

#define SD(val, pdst) {                 \
  uint8_t *pdst_m = (uint8_t *)(pdst);  \
  const uint64_t val_m = (val);         \
                                        \
  __asm__ __volatile__ (                \
      "sd  %[val_m],  %[pdst_m]  \n\t"  \
                                        \
      : [pdst_m] "=m" (*pdst_m)         \
      : [val_m] "r" (val_m)             \
  );                                    \
}
#else  // !(__mips_isa_rev >= 6)
#define LW(psrc) ({                                 \
  const uint8_t *psrc_m = (const uint8_t *)(psrc);  \
  uint32_t val_m;                                   \
                                                    \
  __asm__ __volatile__ (                            \
      "ulw  %[val_m],  %[psrc_m]  \n\t"             \
                                                    \
      : [val_m] "=r" (val_m)                        \
      : [psrc_m] "m" (*psrc_m)                      \
  );                                                \
                                                    \
  val_m;                                            \
})

#define SW(val, pdst) {                  \
  uint8_t *pdst_m = (uint8_t *)(pdst);   \
  const uint32_t val_m = (val);          \
                                         \
  __asm__ __volatile__ (                 \
      "usw  %[val_m],  %[pdst_m]  \n\t"  \
                                         \
      : [pdst_m] "=m" (*pdst_m)          \
      : [val_m] "r" (val_m)              \
  );                                     \
}

#if (__mips == 64)
#define LD(psrc) ({                                 \
  const uint8_t *psrc_m = (const uint8_t *)(psrc);  \
  uint64_t val_m = 0;                               \
                                                    \
  __asm__ __volatile__ (                            \
      "uld  %[val_m],  %[psrc_m]  \n\t"             \
                                                    \
      : [val_m] "=r" (val_m)                        \
      : [psrc_m] "m" (*psrc_m)                      \
  );                                                \
                                                    \
  val_m;                                            \
})
#else  // !(__mips == 64)
#define LD(psrc) ({                                        \
  const uint8_t *psrc_m1 = (const uint8_t *)(psrc);        \
  uint32_t val0_m, val1_m;                                 \
  uint64_t val_m = 0;                                      \
                                                           \
  val0_m = LW(psrc_m1);                                    \
  val1_m = LW(psrc_m1 + 4);                                \
                                                           \
  val_m = (uint64_t)(val1_m);                              \
  val_m = (uint64_t)((val_m << 32) & 0xFFFFFFFF00000000);  \
  val_m = (uint64_t)(val_m | (uint64_t)val0_m);            \
                                                           \
  val_m;                                                   \
})
#endif  // (__mips == 64)

#define SD(val, pdst) {                                     \
  uint8_t *pdst_m1 = (uint8_t *)(pdst);                     \
  uint32_t val0_m, val1_m;                                  \
                                                            \
  val0_m = (uint32_t)((val) & 0x00000000FFFFFFFF);          \
  val1_m = (uint32_t)(((val) >> 32) & 0x00000000FFFFFFFF);  \
                                                            \
  SW(val0_m, pdst_m1);                                      \
  SW(val1_m, pdst_m1 + 4);                                  \
}
#endif  // (__mips_isa_rev >= 6)

/* Description : Load 4 words with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1, out2, out3
   Details     : Load word in 'out0' from (psrc)
                 Load word in 'out1' from (psrc + stride)
                 Load word in 'out2' from (psrc + 2 * stride)
                 Load word in 'out3' from (psrc + 3 * stride)
*/
#define LW4(psrc, stride, out0, out1, out2, out3) {  \
  out0 = LW((psrc));                                 \
  out1 = LW((psrc) + stride);                        \
  out2 = LW((psrc) + 2 * stride);                    \
  out3 = LW((psrc) + 3 * stride);                    \
}

/* Description : Load double words with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
   Details     : Load double word in 'out0' from (psrc)
                 Load double word in 'out1' from (psrc + stride)
*/
#define LD2(psrc, stride, out0, out1) {  \
  out0 = LD((psrc));                     \
  out1 = LD((psrc) + stride);            \
}
#define LD4(psrc, stride, out0, out1, out2, out3) {  \
  LD2((psrc), stride, out0, out1);                   \
  LD2((psrc) + 2 * stride, stride, out2, out3);      \
}

/* Description : Load vectors with 16 byte elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Load 16 byte elements in 'out0' from (psrc)
                 Load 16 byte elements in 'out1' from (psrc + stride)
*/
#define LD_B2(RTYPE, psrc, stride, out0, out1) {  \
  out0 = LD_B(RTYPE, (psrc));                     \
  out1 = LD_B(RTYPE, (psrc) + stride);            \
}
#define LD_UB2(...) LD_B2(v16u8, __VA_ARGS__)
#define LD_SB2(...) LD_B2(v16i8, __VA_ARGS__)

#define LD_B3(RTYPE, psrc, stride, out0, out1, out2) {  \
  LD_B2(RTYPE, (psrc), stride, out0, out1);             \
  out2 = LD_B(RTYPE, (psrc) + 2 * stride);              \
}
#define LD_UB3(...) LD_B3(v16u8, __VA_ARGS__)

#define LD_B4(RTYPE, psrc, stride, out0, out1, out2, out3) {  \
  LD_B2(RTYPE, (psrc), stride, out0, out1);                   \
  LD_B2(RTYPE, (psrc) + 2 * stride , stride, out2, out3);     \
}
#define LD_UB4(...) LD_B4(v16u8, __VA_ARGS__)
#define LD_SB4(...) LD_B4(v16i8, __VA_ARGS__)

#define LD_B5(RTYPE, psrc, stride, out0, out1, out2, out3, out4) {  \
  LD_B4(RTYPE, (psrc), stride, out0, out1, out2, out3);             \
  out4 = LD_B(RTYPE, (psrc) + 4 * stride);                          \
}
#define LD_UB5(...) LD_B5(v16u8, __VA_ARGS__)

#define LD_B8(RTYPE, psrc, stride,                                    \
              out0, out1, out2, out3, out4, out5, out6, out7) {       \
  LD_B4(RTYPE, (psrc), stride, out0, out1, out2, out3);               \
  LD_B4(RTYPE, (psrc) + 4 * stride, stride, out4, out5, out6, out7);  \
}
#define LD_UB8(...) LD_B8(v16u8, __VA_ARGS__)
#define LD_SB8(...) LD_B8(v16i8, __VA_ARGS__)

/* Description : Load vectors with 8 halfword elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
   Details     : Load 8 halfword elements in 'out0' from (psrc)
                 Load 8 halfword elements in 'out1' from (psrc + stride)
*/
#define LD_H2(RTYPE, psrc, stride, out0, out1) {  \
  out0 = LD_H(RTYPE, (psrc));                     \
  out1 = LD_H(RTYPE, (psrc) + (stride));          \
}

#define LD_H4(RTYPE, psrc, stride, out0, out1, out2, out3) {  \
  LD_H2(RTYPE, (psrc), stride, out0, out1);                   \
  LD_H2(RTYPE, (psrc) + 2 * stride, stride, out2, out3);      \
}
#define LD_SH4(...) LD_H4(v8i16, __VA_ARGS__)

/* Description : average with rounding (in0 + in1 + 1) / 2.
   Arguments   : Inputs  - in0, in1, in2, in3,
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Each unsigned byte element from 'in0' vector is added with
                 each unsigned byte element from 'in1' vector. Then the average
                 with rounding is calculated and written to 'out0'
*/
#define AVER_UB2(RTYPE, in0, in1, in2, in3, out0, out1) {  \
  out0 = (RTYPE)__msa_aver_u_b((v16u8)in0, (v16u8)in1);    \
  out1 = (RTYPE)__msa_aver_u_b((v16u8)in2, (v16u8)in3);    \
}
#define AVER_UB2_UB(...) AVER_UB2(v16u8, __VA_ARGS__)

#define AVER_UB4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                 out0, out1, out2, out3) {                       \
  AVER_UB2(RTYPE, in0, in1, in2, in3, out0, out1)                \
  AVER_UB2(RTYPE, in4, in5, in6, in7, out2, out3)                \
}
#define AVER_UB4_UB(...) AVER_UB4(v16u8, __VA_ARGS__)

/* Description : Immediate number of elements to slide
   Arguments   : Inputs  - in0_0, in0_1, in1_0, in1_1, slide_val
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Byte elements from 'in0_0' vector are slid into 'in1_0' by
                 value specified in the 'slide_val'
*/
#define SLDI_B2(RTYPE, in0_0, in0_1, in1_0, in1_1, out0, out1, slide_val) {  \
  out0 = (RTYPE)__msa_sldi_b((v16i8)in0_0, (v16i8)in1_0, slide_val);         \
  out1 = (RTYPE)__msa_sldi_b((v16i8)in0_1, (v16i8)in1_1, slide_val);         \
}
#define SLDI_B2_UB(...) SLDI_B2(v16u8, __VA_ARGS__)

/* Description : Dot product & addition of halfword vector elements
   Arguments   : Inputs  - mult0, mult1, cnst0, cnst1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Signed halfword elements from 'mult0' are multiplied with
                 signed halfword elements from 'cnst0' producing a result
                 twice the size of input i.e. signed word.
                 The multiplication result of adjacent odd-even elements
                 are added to the 'out0' vector
*/
#define DPADD_SH2(RTYPE, mult0, mult1, cnst0, cnst1, out0, out1) {         \
  out0 = (RTYPE)__msa_dpadd_s_w((v4i32)out0, (v8i16)mult0, (v8i16)cnst0);  \
  out1 = (RTYPE)__msa_dpadd_s_w((v4i32)out1, (v8i16)mult1, (v8i16)cnst1);  \
}
#define DPADD_SH2_SW(...) DPADD_SH2(v4i32, __VA_ARGS__)

/* Description : Dot product & addition of double word vector elements
   Arguments   : Inputs  - mult0, mult1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Each signed word element from 'mult0' is multiplied with itself
                 producing an intermediate result twice the size of it
                 i.e. signed double word
                 The multiplication result of adjacent odd-even elements
                 are added to the 'out0' vector
*/
#define DPADD_SD2(RTYPE, mult0, mult1, out0, out1) {                       \
  out0 = (RTYPE)__msa_dpadd_s_d((v2i64)out0, (v4i32)mult0, (v4i32)mult0);  \
  out1 = (RTYPE)__msa_dpadd_s_d((v2i64)out1, (v4i32)mult1, (v4i32)mult1);  \
}
#define DPADD_SD2_SD(...) DPADD_SD2(v2i64, __VA_ARGS__)

/* Description : Horizontal addition of 4 signed word elements of input vector
   Arguments   : Input  - in       (signed word vector)
                 Output - sum_m    (i32 sum)
                 Return Type - signed word (GP)
   Details     : 4 signed word elements of 'in' vector are added together and
                 the resulting integer sum is returned
*/
#define HADD_SW_S32(in) ({                        \
  v2i64 res0_m, res1_m;                           \
  int32_t sum_m;                                  \
                                                  \
  res0_m = __msa_hadd_s_d((v4i32)in, (v4i32)in);  \
  res1_m = __msa_splati_d(res0_m, 1);             \
  res0_m = res0_m + res1_m;                       \
  sum_m = __msa_copy_s_w((v4i32)res0_m, 0);       \
  sum_m;                                          \
})

/* Description : Horizontal addition of 8 unsigned halfword elements
   Arguments   : Inputs  - in       (unsigned halfword vector)
                 Outputs - sum_m    (u32 sum)
                 Return Type - unsigned word
   Details     : 8 unsigned halfword elements of input vector are added
                 together and the resulting integer sum is returned
*/
#define HADD_UH_U32(in) ({                           \
  v4u32 res_m;                                       \
  v2u64 res0_m, res1_m;                              \
  uint32_t sum_m;                                    \
                                                     \
  res_m = __msa_hadd_u_w((v8u16)in, (v8u16)in);      \
  res0_m = __msa_hadd_u_d(res_m, res_m);             \
  res1_m = (v2u64)__msa_splati_d((v2i64)res0_m, 1);  \
  res0_m = res0_m + res1_m;                          \
  sum_m = __msa_copy_u_w((v4i32)res0_m, 0);          \
  sum_m;                                             \
})

/* Description : Horizontal subtraction of unsigned byte vector elements
   Arguments   : Inputs  - in0, in1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Each unsigned odd byte element from 'in0' is subtracted from
                 even unsigned byte element from 'in0' (pairwise) and the
                 halfword result is written to 'out0'
*/
#define HSUB_UB2(RTYPE, in0, in1, out0, out1) {          \
  out0 = (RTYPE)__msa_hsub_u_h((v16u8)in0, (v16u8)in0);  \
  out1 = (RTYPE)__msa_hsub_u_h((v16u8)in1, (v16u8)in1);  \
}
#define HSUB_UB2_SH(...) HSUB_UB2(v8i16, __VA_ARGS__)

/* Description : SAD (Sum of Absolute Difference)
   Arguments   : Inputs  - in0, in1, ref0, ref1
                 Outputs - sad_m                 (halfword vector)
                 Return Type - unsigned halfword
   Details     : Absolute difference of all the byte elements from 'in0' with
                 'ref0' is calculated and preserved in 'diff0'. Then even-odd
                 pairs are added together to generate 8 halfword results.
*/
#define SAD_UB2_UH(in0, in1, ref0, ref1) ({                 \
  v16u8 diff0_m, diff1_m;                                   \
  v8u16 sad_m = { 0 };                                      \
                                                            \
  diff0_m = __msa_asub_u_b((v16u8)in0, (v16u8)ref0);        \
  diff1_m = __msa_asub_u_b((v16u8)in1, (v16u8)ref1);        \
                                                            \
  sad_m += __msa_hadd_u_h((v16u8)diff0_m, (v16u8)diff0_m);  \
  sad_m += __msa_hadd_u_h((v16u8)diff1_m, (v16u8)diff1_m);  \
                                                            \
  sad_m;                                                    \
})

/* Description : Set element n input vector to GPR value
   Arguments   : Inputs - in0, in1, in2, in3
                 Output - out
                 Return Type - as per RTYPE
   Details     : Set element 0 in vector 'out' to value specified in 'in0'
*/
#define INSERT_W4(RTYPE, in0, in1, in2, in3, out) {  \
  out = (RTYPE)__msa_insert_w((v4i32)out, 0, in0);   \
  out = (RTYPE)__msa_insert_w((v4i32)out, 1, in1);   \
  out = (RTYPE)__msa_insert_w((v4i32)out, 2, in2);   \
  out = (RTYPE)__msa_insert_w((v4i32)out, 3, in3);   \
}
#define INSERT_W4_UB(...) INSERT_W4(v16u8, __VA_ARGS__)
#define INSERT_W4_SB(...) INSERT_W4(v16i8, __VA_ARGS__)

#define INSERT_D2(RTYPE, in0, in1, out) {           \
  out = (RTYPE)__msa_insert_d((v2i64)out, 0, in0);  \
  out = (RTYPE)__msa_insert_d((v2i64)out, 1, in1);  \
}
#define INSERT_D2_UB(...) INSERT_D2(v16u8, __VA_ARGS__)
#define INSERT_D2_SB(...) INSERT_D2(v16i8, __VA_ARGS__)

/* Description : Interleave both left and right half of input vectors
   Arguments   : Inputs  - in0, in1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of byte elements from 'in0' and 'in1' are
                 interleaved and written to 'out0'
*/
#define ILVRL_B2(RTYPE, in0, in1, out0, out1) {        \
  out0 = (RTYPE)__msa_ilvr_b((v16i8)in0, (v16i8)in1);  \
  out1 = (RTYPE)__msa_ilvl_b((v16i8)in0, (v16i8)in1);  \
}
#define ILVRL_B2_UB(...) ILVRL_B2(v16u8, __VA_ARGS__)

#define ILVRL_H2(RTYPE, in0, in1, out0, out1) {        \
  out0 = (RTYPE)__msa_ilvr_h((v8i16)in0, (v8i16)in1);  \
  out1 = (RTYPE)__msa_ilvl_h((v8i16)in0, (v8i16)in1);  \
}
#define ILVRL_H2_SW(...) ILVRL_H2(v4i32, __VA_ARGS__)

/* Description : Pack even double word elements of vector pairs
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Even double elements of 'in0' are copied to the left half of
                 'out0' & even double elements of 'in1' are copied to the right
                 half of 'out0'.
*/
#define PCKEV_D2(RTYPE, in0, in1, in2, in3, out0, out1) {  \
  out0 = (RTYPE)__msa_pckev_d((v2i64)in0, (v2i64)in1);     \
  out1 = (RTYPE)__msa_pckev_d((v2i64)in2, (v2i64)in3);     \
}
#define PCKEV_D2_UB(...) PCKEV_D2(v16u8, __VA_ARGS__)

#define PCKEV_D4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                 out0, out1, out2, out3) {                       \
  PCKEV_D2(RTYPE, in0, in1, in2, in3, out0, out1);               \
  PCKEV_D2(RTYPE, in4, in5, in6, in7, out2, out3);               \
}
#define PCKEV_D4_UB(...) PCKEV_D4(v16u8, __VA_ARGS__)

/* Description : Sign extend halfword elements from input vector and return
                 the result in pair of vectors
   Arguments   : Input   - in            (halfword vector)
                 Outputs - out0, out1   (sign extended word vectors)
                 Return Type - signed word
   Details     : Sign bit of halfword elements from input vector 'in' is
                 extracted and interleaved right with same vector 'in0' to
                 generate 4 signed word elements in 'out0'
                 Then interleaved left with same vector 'in0' to
                 generate 4 signed word elements in 'out1'
*/
#define UNPCK_SH_SW(in, out0, out1) {    \
  v8i16 tmp_m;                           \
                                         \
  tmp_m = __msa_clti_s_h((v8i16)in, 0);  \
  ILVRL_H2_SW(tmp_m, in, out0, out1);    \
}

/* Description : Store 4 double words with stride
   Arguments   : Inputs - in0, in1, in2, in3, pdst, stride
   Details     : Store double word from 'in0' to (pdst)
                 Store double word from 'in1' to (pdst + stride)
                 Store double word from 'in2' to (pdst + 2 * stride)
                 Store double word from 'in3' to (pdst + 3 * stride)
*/
#define SD4(in0, in1, in2, in3, pdst, stride) {  \
  SD(in0, (pdst))                                \
  SD(in1, (pdst) + stride);                      \
  SD(in2, (pdst) + 2 * stride);                  \
  SD(in3, (pdst) + 3 * stride);                  \
}

/* Description : Store vectors of 8 halfword elements with stride
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Store 8 halfword elements from 'in0' to (pdst)
                 Store 8 halfword elements from 'in1' to (pdst + stride)
*/
#define ST_H2(RTYPE, in0, in1, pdst, stride) {  \
  ST_H(RTYPE, in0, (pdst));                     \
  ST_H(RTYPE, in1, (pdst) + stride);            \
}
#define ST_SH2(...) ST_H2(v8i16, __VA_ARGS__)

/* Description : Store 8x4 byte block to destination memory from input
                 vectors
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Index 0 double word element from 'in0' vector is copied to the
                 GP register and stored to (pdst)
                 Index 1 double word element from 'in0' vector is copied to the
                 GP register and stored to (pdst + stride)
                 Index 0 double word element from 'in1' vector is copied to the
                 GP register and stored to (pdst + 2 * stride)
                 Index 1 double word element from 'in1' vector is copied to the
                 GP register and stored to (pdst + 3 * stride)
*/
#define ST8x4_UB(in0, in1, pdst, stride) {                  \
  uint64_t out0_m, out1_m, out2_m, out3_m;                  \
  uint8_t *pblk_8x4_m = (uint8_t *)(pdst);                  \
                                                            \
  out0_m = __msa_copy_u_d((v2i64)in0, 0);                   \
  out1_m = __msa_copy_u_d((v2i64)in0, 1);                   \
  out2_m = __msa_copy_u_d((v2i64)in1, 0);                   \
  out3_m = __msa_copy_u_d((v2i64)in1, 1);                   \
                                                            \
  SD4(out0_m, out1_m, out2_m, out3_m, pblk_8x4_m, stride);  \
}
#endif  /* VPX_DSP_MIPS_MACROS_MSA_H_ */
