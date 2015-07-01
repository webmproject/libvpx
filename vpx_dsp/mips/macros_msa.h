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

#define LD_B4(RTYPE, psrc, stride, out0, out1, out2, out3) {  \
  LD_B2(RTYPE, (psrc), stride, out0, out1);                   \
  LD_B2(RTYPE, (psrc) + 2 * stride , stride, out2, out3);     \
}
#define LD_UB4(...) LD_B4(v16u8, __VA_ARGS__)

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
#endif  /* VPX_DSP_MIPS_MACROS_MSA_H_ */
