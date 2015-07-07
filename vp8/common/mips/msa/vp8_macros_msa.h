/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP8_COMMON_MIPS_MSA_VP8_MACROS_MSA_H_
#define VP8_COMMON_MIPS_MSA_VP8_MACROS_MSA_H_

#include <msa.h>

#include "./vpx_config.h"
#include "vpx/vpx_integer.h"

#define LD_B(RTYPE, psrc) *((const RTYPE *)(psrc))
#define LD_UB(...) LD_B(v16u8, __VA_ARGS__)
#define LD_SB(...) LD_B(v16i8, __VA_ARGS__)

#define LD_H(RTYPE, psrc) *((const RTYPE *)(psrc))
#define LD_UH(...) LD_H(v8u16, __VA_ARGS__)
#define LD_SH(...) LD_H(v8i16, __VA_ARGS__)

#define ST_B(RTYPE, in, pdst) *((RTYPE *)(pdst)) = (in)
#define ST_UB(...) ST_B(v16u8, __VA_ARGS__)
#define ST_SB(...) ST_B(v16i8, __VA_ARGS__)

#define ST_H(RTYPE, in, pdst) *((RTYPE *)(pdst)) = (in)
#define ST_UH(...) ST_H(v8u16, __VA_ARGS__)
#define ST_SH(...) ST_H(v8i16, __VA_ARGS__)

/* Description : Load vectors with 16 byte elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Load 16 byte elements in 'out0' from (psrc)
                 Load 16 byte elements in 'out1' from (psrc + stride)
*/
#define LD_B2(RTYPE, psrc, stride, out0, out1)  \
{                                               \
    out0 = LD_B(RTYPE, (psrc));                 \
    out1 = LD_B(RTYPE, (psrc) + stride);        \
}
#define LD_UB2(...) LD_B2(v16u8, __VA_ARGS__)
#define LD_SB2(...) LD_B2(v16i8, __VA_ARGS__)

#define LD_B4(RTYPE, psrc, stride, out0, out1, out2, out3)   \
{                                                            \
    LD_B2(RTYPE, (psrc), stride, out0, out1);                \
    LD_B2(RTYPE, (psrc) + 2 * stride , stride, out2, out3);  \
}
#define LD_UB4(...) LD_B4(v16u8, __VA_ARGS__)
#define LD_SB4(...) LD_B4(v16i8, __VA_ARGS__)

/* Description : Load vectors with 8 halfword elements with stride
   Arguments   : Inputs  - psrc, stride
                 Outputs - out0, out1
   Details     : Load 8 halfword elements in 'out0' from (psrc)
                 Load 8 halfword elements in 'out1' from (psrc + stride)
*/
#define LD_H2(RTYPE, psrc, stride, out0, out1)  \
{                                               \
    out0 = LD_H(RTYPE, (psrc));                 \
    out1 = LD_H(RTYPE, (psrc) + (stride));      \
}
#define LD_SH2(...) LD_H2(v8i16, __VA_ARGS__)

#define LD_H4(RTYPE, psrc, stride, out0, out1, out2, out3)  \
{                                                           \
    LD_H2(RTYPE, (psrc), stride, out0, out1);               \
    LD_H2(RTYPE, (psrc) + 2 * stride, stride, out2, out3);  \
}
#define LD_SH4(...) LD_H4(v8i16, __VA_ARGS__)

/* Description : Store vectors of 16 byte elements with stride
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Store 16 byte elements from 'in0' to (pdst)
                 Store 16 byte elements from 'in1' to (pdst + stride)
*/
#define ST_B2(RTYPE, in0, in1, pdst, stride)  \
{                                             \
    ST_B(RTYPE, in0, (pdst));                 \
    ST_B(RTYPE, in1, (pdst) + stride);        \
}

#define ST_B4(RTYPE, in0, in1, in2, in3, pdst, stride)    \
{                                                         \
    ST_B2(RTYPE, in0, in1, (pdst), stride);               \
    ST_B2(RTYPE, in2, in3, (pdst) + 2 * stride, stride);  \
}
#define ST_UB4(...) ST_B4(v16u8, __VA_ARGS__)
#define ST_SB4(...) ST_B4(v16i8, __VA_ARGS__)

/* Description : Store vectors of 8 halfword elements with stride
   Arguments   : Inputs - in0, in1, pdst, stride
   Details     : Store 8 halfword elements from 'in0' to (pdst)
                 Store 8 halfword elements from 'in1' to (pdst + stride)
*/
#define ST_H2(RTYPE, in0, in1, pdst, stride)  \
{                                             \
    ST_H(RTYPE, in0, (pdst));                 \
    ST_H(RTYPE, in1, (pdst) + stride);        \
}
#define ST_SH2(...) ST_H2(v8i16, __VA_ARGS__)

/* Description : Shuffle byte vector elements as per mask vector
   Arguments   : Inputs  - in0, in1, in2, in3, mask0, mask1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Byte elements from 'in0' & 'in1' are copied selectively to
                 'out0' as per control vector 'mask0'
*/
#define VSHF_B2(RTYPE, in0, in1, in2, in3, mask0, mask1, out0, out1)   \
{                                                                      \
    out0 = (RTYPE)__msa_vshf_b((v16i8)mask0, (v16i8)in1, (v16i8)in0);  \
    out1 = (RTYPE)__msa_vshf_b((v16i8)mask1, (v16i8)in3, (v16i8)in2);  \
}
#define VSHF_B2_SB(...) VSHF_B2(v16i8, __VA_ARGS__)

/* Description : Clips all signed halfword elements of input vector
                 between 0 & 255
   Arguments   : Input  - in
                 Output - out_m
                 Return Type - signed halfword
*/
#define CLIP_SH_0_255(in)                               \
({                                                      \
    v8i16 max_m = __msa_ldi_h(255);                     \
    v8i16 out_m;                                        \
                                                        \
    out_m = __msa_maxi_s_h((v8i16)in, 0);               \
    out_m = __msa_min_s_h((v8i16)max_m, (v8i16)out_m);  \
    out_m;                                              \
})
#define CLIP_SH2_0_255(in0, in1)  \
{                                 \
    in0 = CLIP_SH_0_255(in0);     \
    in1 = CLIP_SH_0_255(in1);     \
}
#define CLIP_SH4_0_255(in0, in1, in2, in3)  \
{                                           \
    CLIP_SH2_0_255(in0, in1);               \
    CLIP_SH2_0_255(in2, in3);               \
}

/* Description : Clips all signed word elements of input vector
                 between 0 & 255
   Arguments   : Input  - in
                 Output - out_m
                 Return Type - signed word
*/
#define CLIP_SW_0_255(in)                               \
({                                                      \
    v4i32 max_m = __msa_ldi_w(255);                     \
    v4i32 out_m;                                        \
                                                        \
    out_m = __msa_maxi_s_w((v4i32)in, 0);               \
    out_m = __msa_min_s_w((v4i32)max_m, (v4i32)out_m);  \
    out_m;                                              \
})

/* Description : Interleave left half of halfword elements from vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Left half of halfword elements of 'in0' and 'in1' are
                 interleaved and written to 'out0'.
*/
#define ILVL_H2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                        \
    out0 = (RTYPE)__msa_ilvl_h((v8i16)in0, (v8i16)in1);  \
    out1 = (RTYPE)__msa_ilvl_h((v8i16)in2, (v8i16)in3);  \
}
#define ILVL_H2_SH(...) ILVL_H2(v8i16, __VA_ARGS__)
#define ILVL_H2_SW(...) ILVL_H2(v4i32, __VA_ARGS__)

/* Description : Interleave left half of word elements from vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Left half of word elements of 'in0' and 'in1' are interleaved
                 and written to 'out0'.
*/
#define ILVL_W2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                        \
    out0 = (RTYPE)__msa_ilvl_w((v4i32)in0, (v4i32)in1);  \
    out1 = (RTYPE)__msa_ilvl_w((v4i32)in2, (v4i32)in3);  \
}
#define ILVL_W2_SH(...) ILVL_W2(v8i16, __VA_ARGS__)

/* Description : Interleave right half of byte elements from vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of byte elements of 'in0' and 'in1' are interleaved
                 and written to out0.
*/
#define ILVR_B2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                        \
    out0 = (RTYPE)__msa_ilvr_b((v16i8)in0, (v16i8)in1);  \
    out1 = (RTYPE)__msa_ilvr_b((v16i8)in2, (v16i8)in3);  \
}

#define ILVR_B4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                out0, out1, out2, out3)                         \
{                                                               \
    ILVR_B2(RTYPE, in0, in1, in2, in3, out0, out1);             \
    ILVR_B2(RTYPE, in4, in5, in6, in7, out2, out3);             \
}
#define ILVR_B4_SH(...) ILVR_B4(v8i16, __VA_ARGS__)
#define ILVR_B4_SW(...) ILVR_B4(v4i32, __VA_ARGS__)

/* Description : Interleave right half of halfword elements from vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of halfword elements of 'in0' and 'in1' are
                 interleaved and written to 'out0'.
*/
#define ILVR_H2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                        \
    out0 = (RTYPE)__msa_ilvr_h((v8i16)in0, (v8i16)in1);  \
    out1 = (RTYPE)__msa_ilvr_h((v8i16)in2, (v8i16)in3);  \
}
#define ILVR_H2_SH(...) ILVR_H2(v8i16, __VA_ARGS__)
#define ILVR_H2_SW(...) ILVR_H2(v4i32, __VA_ARGS__)

#define ILVR_H4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                out0, out1, out2, out3)                         \
{                                                               \
    ILVR_H2(RTYPE, in0, in1, in2, in3, out0, out1);             \
    ILVR_H2(RTYPE, in4, in5, in6, in7, out2, out3);             \
}
#define ILVR_H4_SH(...) ILVR_H4(v8i16, __VA_ARGS__)
#define ILVR_H4_SW(...) ILVR_H4(v4i32, __VA_ARGS__)

#define ILVR_W2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                        \
    out0 = (RTYPE)__msa_ilvr_w((v4i32)in0, (v4i32)in1);  \
    out1 = (RTYPE)__msa_ilvr_w((v4i32)in2, (v4i32)in3);  \
}
#define ILVR_W2_SH(...) ILVR_W2(v8i16, __VA_ARGS__)

/* Description : Interleave right half of double word elements from vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of double word elements of 'in0' and 'in1' are
                 interleaved and written to 'out0'.
*/
#define ILVR_D2(RTYPE, in0, in1, in2, in3, out0, out1)       \
{                                                            \
    out0 = (RTYPE)__msa_ilvr_d((v2i64)(in0), (v2i64)(in1));  \
    out1 = (RTYPE)__msa_ilvr_d((v2i64)(in2), (v2i64)(in3));  \
}
#define ILVR_D2_UB(...) ILVR_D2(v16u8, __VA_ARGS__)
#define ILVR_D2_SB(...) ILVR_D2(v16i8, __VA_ARGS__)
#define ILVR_D2_SH(...) ILVR_D2(v8i16, __VA_ARGS__)

/* Description : Interleave both left and right half of input vectors
   Arguments   : Inputs  - in0, in1
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Right half of byte elements from 'in0' and 'in1' are
                 interleaved and written to 'out0'
*/
#define ILVRL_H2(RTYPE, in0, in1, out0, out1)            \
{                                                        \
    out0 = (RTYPE)__msa_ilvr_h((v8i16)in0, (v8i16)in1);  \
    out1 = (RTYPE)__msa_ilvl_h((v8i16)in0, (v8i16)in1);  \
}
#define ILVRL_H2_SH(...) ILVRL_H2(v8i16, __VA_ARGS__)
#define ILVRL_H2_SW(...) ILVRL_H2(v4i32, __VA_ARGS__)

#define ILVRL_W2(RTYPE, in0, in1, out0, out1)            \
{                                                        \
    out0 = (RTYPE)__msa_ilvr_w((v4i32)in0, (v4i32)in1);  \
    out1 = (RTYPE)__msa_ilvl_w((v4i32)in0, (v4i32)in1);  \
}
#define ILVRL_W2_UB(...) ILVRL_W2(v16u8, __VA_ARGS__)
#define ILVRL_W2_SH(...) ILVRL_W2(v8i16, __VA_ARGS__)
#define ILVRL_W2_SW(...) ILVRL_W2(v4i32, __VA_ARGS__)

/* Description : Pack even byte elements of vector pairs
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Even byte elements of 'in0' are copied to the left half of
                 'out0' & even byte elements of 'in1' are copied to the right
                 half of 'out0'.
*/
#define PCKEV_B2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                         \
    out0 = (RTYPE)__msa_pckev_b((v16i8)in0, (v16i8)in1);  \
    out1 = (RTYPE)__msa_pckev_b((v16i8)in2, (v16i8)in3);  \
}
#define PCKEV_B2_SB(...) PCKEV_B2(v16i8, __VA_ARGS__)
#define PCKEV_B2_UB(...) PCKEV_B2(v16u8, __VA_ARGS__)

#define PCKEV_B4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                 out0, out1, out2, out3)                         \
{                                                                \
    PCKEV_B2(RTYPE, in0, in1, in2, in3, out0, out1);             \
    PCKEV_B2(RTYPE, in4, in5, in6, in7, out2, out3);             \
}
#define PCKEV_B4_SB(...) PCKEV_B4(v16i8, __VA_ARGS__)
#define PCKEV_B4_UB(...) PCKEV_B4(v16u8, __VA_ARGS__)
#define PCKEV_B4_SH(...) PCKEV_B4(v8i16, __VA_ARGS__)

/* Description : Pack even halfword elements of vector pairs
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Even halfword elements of 'in0' are copied to the left half of
                 'out0' & even halfword elements of 'in1' are copied to the
                 right half of 'out0'.
*/
#define PCKEV_H2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                         \
    out0 = (RTYPE)__msa_pckev_h((v8i16)in0, (v8i16)in1);  \
    out1 = (RTYPE)__msa_pckev_h((v8i16)in2, (v8i16)in3);  \
}
#define PCKEV_H2_SH(...) PCKEV_H2(v8i16, __VA_ARGS__)

#define PCKEV_H4(RTYPE, in0, in1, in2, in3, in4, in5, in6, in7,  \
                 out0, out1, out2, out3)                         \
{                                                                \
    PCKEV_H2(RTYPE, in0, in1, in2, in3, out0, out1);             \
    PCKEV_H2(RTYPE, in4, in5, in6, in7, out2, out3);             \
}
#define PCKEV_H4_SH(...) PCKEV_H4(v8i16, __VA_ARGS__)

/* Description : Pack even double word elements of vector pairs
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Even double elements of 'in0' are copied to the left half of
                 'out0' & even double elements of 'in1' are copied to the right
                 half of 'out0'.
*/
#define PCKEV_D2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                         \
    out0 = (RTYPE)__msa_pckev_d((v2i64)in0, (v2i64)in1);  \
    out1 = (RTYPE)__msa_pckev_d((v2i64)in2, (v2i64)in3);  \
}
#define PCKEV_D2_UB(...) PCKEV_D2(v16u8, __VA_ARGS__)
#define PCKEV_D2_SH(...) PCKEV_D2(v8i16, __VA_ARGS__)

/* Description : Pack odd double word elements of vector pairs
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
                 Return Type - as per RTYPE
   Details     : Odd double word elements of 'in0' are copied to the left half
                 of 'out0' & odd double word elements of 'in1' are copied to
                 the right half of 'out0'.
*/
#define PCKOD_D2(RTYPE, in0, in1, in2, in3, out0, out1)   \
{                                                         \
    out0 = (RTYPE)__msa_pckod_d((v2i64)in0, (v2i64)in1);  \
    out1 = (RTYPE)__msa_pckod_d((v2i64)in2, (v2i64)in3);  \
}
#define PCKOD_D2_UB(...) PCKOD_D2(v16u8, __VA_ARGS__)
#define PCKOD_D2_SH(...) PCKOD_D2(v8i16, __VA_ARGS__)

/* Description : Arithmetic shift right all elements of vector
                 (generic for all data types)
   Arguments   : Inputs  - in0, in1, in2, in3, shift
                 Outputs - in place operation
                 Return Type - as per input vector RTYPE
   Details     : Each element of vector 'in0' is right shifted by 'shift' and
                 the result is written in-place. 'shift' is a GP variable.
*/
#define SRA_4V(in0, in1, in2, in3, shift)  \
{                                          \
    in0 = in0 >> shift;                    \
    in1 = in1 >> shift;                    \
    in2 = in2 >> shift;                    \
    in3 = in3 >> shift;                    \
}

/* Description : Shift right arithmetic rounded (immediate)
   Arguments   : Inputs  - in0, in1, shift
                 Outputs - in place operation
                 Return Type - as per RTYPE
   Details     : Each element of vector 'in0' is shifted right arithmetically by
                 the value in 'shift'. The last discarded bit is added to the
                 shifted value for rounding and the result is written in-place.
                 'shift' is an immediate value.
*/
#define SRARI_H2(RTYPE, in0, in1, shift)            \
{                                                   \
    in0 = (RTYPE)__msa_srari_h((v8i16)in0, shift);  \
    in1 = (RTYPE)__msa_srari_h((v8i16)in1, shift);  \
}
#define SRARI_H2_UH(...) SRARI_H2(v8u16, __VA_ARGS__)
#define SRARI_H2_SH(...) SRARI_H2(v8i16, __VA_ARGS__)

#define SRARI_W2(RTYPE, in0, in1, shift)            \
{                                                   \
    in0 = (RTYPE)__msa_srari_w((v4i32)in0, shift);  \
    in1 = (RTYPE)__msa_srari_w((v4i32)in1, shift);  \
}

#define SRARI_W4(RTYPE, in0, in1, in2, in3, shift)  \
{                                                   \
    SRARI_W2(RTYPE, in0, in1, shift);               \
    SRARI_W2(RTYPE, in2, in3, shift);               \
}
#define SRARI_W4_SW(...) SRARI_W4(v4i32, __VA_ARGS__)

/* Description : Multiplication of pairs of vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
   Details     : Each element from 'in0' is multiplied with elements from 'in1'
                 and the result is written to 'out0'
*/
#define MUL2(in0, in1, in2, in3, out0, out1)  \
{                                             \
    out0 = in0 * in1;                         \
    out1 = in2 * in3;                         \
}
#define MUL4(in0, in1, in2, in3, in4, in5, in6, in7,  \
             out0, out1, out2, out3)                  \
{                                                     \
    MUL2(in0, in1, in2, in3, out0, out1);             \
    MUL2(in4, in5, in6, in7, out2, out3);             \
}

/* Description : Addition of 2 pairs of vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1
   Details     : Each element in 'in0' is added to 'in1' and result is written
                 to 'out0'.
*/
#define ADD2(in0, in1, in2, in3, out0, out1)  \
{                                             \
    out0 = in0 + in1;                         \
    out1 = in2 + in3;                         \
}
#define ADD4(in0, in1, in2, in3, in4, in5, in6, in7,  \
             out0, out1, out2, out3)                  \
{                                                     \
    ADD2(in0, in1, in2, in3, out0, out1);             \
    ADD2(in4, in5, in6, in7, out2, out3);             \
}

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
#define UNPCK_SH_SW(in, out0, out1)        \
{                                          \
    v8i16 tmp_m;                           \
                                           \
    tmp_m = __msa_clti_s_h((v8i16)in, 0);  \
    ILVRL_H2_SW(tmp_m, in, out0, out1);    \
}

/* Description : Butterfly of 4 input vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1, out2, out3
   Details     : Butterfly operation
*/
#define BUTTERFLY_4(in0, in1, in2, in3, out0, out1, out2, out3)  \
{                                                                \
    out0 = in0 + in3;                                            \
    out1 = in1 + in2;                                            \
                                                                 \
    out2 = in1 - in2;                                            \
    out3 = in0 - in3;                                            \
}

/* Description : Transpose 8x4 block with half word elements in vectors
   Arguments   : Inputs  - in0, in1, in2, in3, in4, in5, in6, in7
                 Outputs - out0, out1, out2, out3, out4, out5, out6, out7
                 Return Type - signed halfword
*/
#define TRANSPOSE8X4_SH_SH(in0, in1, in2, in3, out0, out1, out2, out3)  \
{                                                                       \
    v8i16 tmp0_m, tmp1_m, tmp2_m, tmp3_m;                               \
                                                                        \
    ILVR_H2_SH(in1, in0, in3, in2, tmp0_m, tmp1_m);                     \
    ILVL_H2_SH(in1, in0, in3, in2, tmp2_m, tmp3_m);                     \
    ILVR_W2_SH(tmp1_m, tmp0_m, tmp3_m, tmp2_m, out0, out2);             \
    ILVL_W2_SH(tmp1_m, tmp0_m, tmp3_m, tmp2_m, out1, out3);             \
}

/* Description : Transpose 4x4 block with word elements in vectors
   Arguments   : Inputs  - in0, in1, in2, in3
                 Outputs - out0, out1, out2, out3
                 Return Type - signed word
*/
#define TRANSPOSE4x4_SW_SW(in0, in1, in2, in3, out0, out1, out2, out3)  \
{                                                                       \
    v4i32 s0_m, s1_m, s2_m, s3_m;                                       \
                                                                        \
    ILVRL_W2_SW(in1, in0, s0_m, s1_m);                                  \
    ILVRL_W2_SW(in3, in2, s2_m, s3_m);                                  \
                                                                        \
    out0 = (v4i32)__msa_ilvr_d((v2i64)s2_m, (v2i64)s0_m);               \
    out1 = (v4i32)__msa_ilvl_d((v2i64)s2_m, (v2i64)s0_m);               \
    out2 = (v4i32)__msa_ilvr_d((v2i64)s3_m, (v2i64)s1_m);               \
    out3 = (v4i32)__msa_ilvl_d((v2i64)s3_m, (v2i64)s1_m);               \
}
#endif  /* VP8_COMMON_MIPS_MSA_VP8_MACROS_MSA_H_ */
