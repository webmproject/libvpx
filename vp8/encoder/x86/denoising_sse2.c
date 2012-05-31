/*
 *  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "vp8/encoder/denoising.h"

#include "vp8/common/reconinter.h"
#include "vpx/vpx_integer.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_rtcd.h"

#include <emmintrin.h>

union sum_union {
    __m128i v;
    short e[8];
};

int vp8_denoiser_filter_sse2(YV12_BUFFER_CONFIG *mc_running_avg,
                             YV12_BUFFER_CONFIG *running_avg,
                             MACROBLOCK *signal, unsigned int motion_magnitude,
                             int y_offset, int uv_offset)
{
    unsigned char filtered_buf[16*16];
    unsigned char *filtered = filtered_buf;
    unsigned char *sig = signal->thismb;
    int sig_stride = 16;
    unsigned char *mc_running_avg_y = mc_running_avg->y_buffer + y_offset;
    int mc_avg_y_stride = mc_running_avg->y_stride;
    unsigned char *running_avg_y = running_avg->y_buffer + y_offset;
    int avg_y_stride = running_avg->y_stride;
    const union coeff_pair *LUT = vp8_get_filter_coeff_LUT(motion_magnitude);
    int r, c;
    __m128i acc_diff = { 0 };

    for (r = 0; r < 16; ++r)
    {
        __m128i filter_coefficient_00, filter_coefficient_04;
        __m128i filter_coefficient_08, filter_coefficient_12;
        __m128i v_sig0, v_sig1;
        __m128i v_mc_running_avg_y0, v_mc_running_avg_y1;
        __m128i state0, state1, state2, state3;
        __m128i res0, res1, res2, res3;
        __m128i v_running_avg_y;
        __m128i diff0, diff1, diff0sq, diff1sq, diff_sq;
        const __m128i kNOISE_DIFF2_THRESHOLD =
                _mm_set1_epi8(NOISE_DIFF2_THRESHOLD);
        __m128i take_running, p0, p1, p2;
        const __m128i k_zero = _mm_set1_epi16(0);
        const __m128i k_128 = _mm_set1_epi32(128);

        // Calculate absolute differences
        DECLARE_ALIGNED_ARRAY(16,unsigned char,abs_diff,16);
        DECLARE_ALIGNED_ARRAY(16,uint32_t,filter_coefficient,16);
        __m128i v_sig = _mm_loadu_si128((__m128i *)(&sig[0]));
        __m128i v_mc_running_avg_y = _mm_loadu_si128(
                                         (__m128i *)(&mc_running_avg_y[0]));
        __m128i a_minus_b = _mm_subs_epu8(v_sig, v_mc_running_avg_y);
        __m128i b_minus_a = _mm_subs_epu8(v_mc_running_avg_y, v_sig);
        __m128i v_abs_diff = _mm_adds_epu8(a_minus_b, b_minus_a);
        _mm_store_si128((__m128i *)(&abs_diff[0]), v_abs_diff);

        // Use LUT to get filter coefficients (two 16b value; f and 256-f)
        for (c = 0; c < 16; ++c)
        {
            filter_coefficient[c] = LUT[abs_diff[c]].as_int;
        }

        // Filtering...
        // load filter coefficients (two 16b value; f and 256-f)
        filter_coefficient_00 = _mm_load_si128(
                (__m128i *)(&filter_coefficient[ 0]));
        filter_coefficient_04 = _mm_load_si128(
                (__m128i *)(&filter_coefficient[ 4]));
        filter_coefficient_08 = _mm_load_si128(
                (__m128i *)(&filter_coefficient[ 8]));
        filter_coefficient_12 = _mm_load_si128(
                (__m128i *)(&filter_coefficient[12]));

        // expand sig from 8b to 16b
        v_sig0 = _mm_unpacklo_epi8(v_sig, k_zero);
        v_sig1 = _mm_unpackhi_epi8(v_sig, k_zero);
        // expand mc_running_avg_y from 8b to 16b
        v_mc_running_avg_y0 = _mm_unpacklo_epi8(v_mc_running_avg_y, k_zero);
        v_mc_running_avg_y1 = _mm_unpackhi_epi8(v_mc_running_avg_y, k_zero);
        // interleave sig and mc_running_avg_y for upcoming multiply-add
        state0 = _mm_unpacklo_epi16(v_mc_running_avg_y0, v_sig0);
        state1 = _mm_unpackhi_epi16(v_mc_running_avg_y0, v_sig0);
        state2 = _mm_unpacklo_epi16(v_mc_running_avg_y1, v_sig1);
        state3 = _mm_unpackhi_epi16(v_mc_running_avg_y1, v_sig1);
        // blend values
        res0 = _mm_madd_epi16(filter_coefficient_00, state0);
        res1 = _mm_madd_epi16(filter_coefficient_04, state1);
        res2 = _mm_madd_epi16(filter_coefficient_08, state2);
        res3 = _mm_madd_epi16(filter_coefficient_12, state3);
        res0 = _mm_add_epi32(res0, k_128);
        res1 = _mm_add_epi32(res1, k_128);
        res2 = _mm_add_epi32(res2, k_128);
        res3 = _mm_add_epi32(res3, k_128);
        res0 = _mm_srai_epi32(res0, 8);
        res1 = _mm_srai_epi32(res1, 8);
        res2 = _mm_srai_epi32(res2, 8);
        res3 = _mm_srai_epi32(res3, 8);
        // combine the 32b results into a single 8b vector
        res0 = _mm_packs_epi32(res0, res1);
        res2 = _mm_packs_epi32(res2, res3);
        v_running_avg_y = _mm_packus_epi16(res0, res2);

        // Depending on the magnitude of the difference between the signal and
        // filtered version, either replace the signal by the filtered one or
        // update the filter state with the signal when the change in a pixel
        // isn't classified as noise.
        diff0 = _mm_sub_epi16(v_sig0, res0);
        diff1 = _mm_sub_epi16(v_sig1, res2);
        acc_diff = _mm_add_epi16(acc_diff, _mm_add_epi16(diff0, diff1));

        diff0sq = _mm_mullo_epi16(diff0, diff0);
        diff1sq = _mm_mullo_epi16(diff1, diff1);
        diff_sq = _mm_packus_epi16(diff0sq, diff1sq);
        take_running = _mm_cmplt_epi8(diff_sq, kNOISE_DIFF2_THRESHOLD);
        p0 = _mm_and_si128(take_running, v_running_avg_y);
        p1 = _mm_andnot_si128(take_running, v_sig);
        p2 = _mm_or_si128(p0, p1);
        _mm_storeu_si128((__m128i *)(&running_avg_y[0]), p2);
        _mm_storeu_si128((__m128i *)(&filtered[0]), p2);

        // Update pointers for next iteration.
        sig += sig_stride;
        filtered += 16;
        mc_running_avg_y += mc_avg_y_stride;
        running_avg_y += avg_y_stride;
    }
    {
        // Compute the sum of all pixel differences of this MB.
        union sum_union s;
        int sum_diff;
        s.v = acc_diff;
        sum_diff = s.e[0] + s.e[1] + s.e[2] + s.e[3] +
          s.e[4] + s.e[5] + s.e[6] + s.e[7];
        if (abs(sum_diff) > SUM_DIFF_THRESHOLD)
        {
            return COPY_BLOCK;
        }
    }
    vp8_copy_mem16x16(filtered_buf, 16, signal->thismb, sig_stride);
    return FILTER_BLOCK;
}
