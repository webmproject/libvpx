#include <tmmintrin.h>
#include <inttypes.h>
#include <emmintrin.h>
#include "vpx_ports/mem.h"
#include "vpx_ports/emmintrin_compat.h"
#pragma GCC push_options
#pragma GCC optimize("unroll-loops")

DECLARE_ALIGNED(16, const unsigned char, filt1_4_h8[16]) = {
  0, 1, 1, 2, 2, 3, 3, 4, 2, 3, 3, 4, 4, 5, 5, 6
};
DECLARE_ALIGNED(16, const unsigned char, filt2_4_h8[16]) = {
  4, 5, 5, 6, 6, 7, 7, 8, 6, 7, 7, 8, 8, 9, 9, 10
};
DECLARE_ALIGNED(16, const unsigned char, filt1_global[16]) = {
  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8
};
DECLARE_ALIGNED(16, const unsigned char, filt2_global[16]) = {
  2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10
};
DECLARE_ALIGNED(16, const unsigned char, filt3_global[16]) = {
  4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12
};
DECLARE_ALIGNED(16, const unsigned char, filt4_global[16]) = {
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14
};

void vp9_filter_block1d4_h8_intrin_ssse3(unsigned char *src_ptr,
                                         unsigned int src_pixels_per_line,
                                         unsigned char *output_ptr,
                                         unsigned int output_pitch,
                                         unsigned int output_height,
                                         short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters,
          thirdFilters, forthFilters, srcReg, srcRegFilt1, srcRegFilt2,
          srcRegFilt3, srcRegFilt4;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg = _mm_packs_epi16(filtersReg, filtersReg);
  firstFilters = _mm_shufflelo_epi16(filtersReg, 0);
  secondFilters = _mm_shufflelo_epi16(filtersReg, 0xAAu);
  firstFilters = _mm_shufflehi_epi16(firstFilters, 0x55u);
  secondFilters = _mm_shufflehi_epi16(secondFilters, 0xFFu);
  thirdFilters =_mm_load_si128((__m128i const *)filt1_4_h8);
  forthFilters = _mm_load_si128((__m128i const *)filt2_4_h8);

  for (i = 0; i < output_height; i++) {
    srcReg = _mm_loadu_si128((__m128i *)(src_ptr-3));
    srcRegFilt1= _mm_shuffle_epi8(srcReg, thirdFilters);
    srcRegFilt2= _mm_shuffle_epi8(srcReg, forthFilters);
    srcRegFilt1 = _mm_maddubs_epi16(srcRegFilt1, firstFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, secondFilters);
    srcRegFilt3 =  _mm_srli_si128(srcRegFilt1, 8);
    srcRegFilt4 =  _mm_srli_si128(srcRegFilt2, 8);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt3);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt4);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt2);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, addFilterReg127);
    srcRegFilt1 = _mm_srai_epi16(srcRegFilt1, 7);
    srcRegFilt1 = _mm_packus_epi16(srcRegFilt1, srcRegFilt1);
    src_ptr += src_pixels_per_line;
    *((int*)&output_ptr[0])= _mm_cvtsi128_si32(srcRegFilt1);
    output_ptr += output_pitch;
  }
}

void vp9_filter_block1d8_h8_intrin_ssse3(unsigned char *src_ptr,
                                         unsigned int src_pixels_per_line,
                                         unsigned char *output_ptr,
                                         unsigned int output_pitch,
                                         unsigned int output_height,
                                         short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters,
          thirdFilters, forthFilters, srcReg, srcRegFilt1, srcRegFilt2,
          srcRegFilt3, srcRegFilt4;
  __m128i filt1Reg, filt2Reg, filt3Reg, filt4Reg;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg =_mm_packs_epi16(filtersReg, filtersReg);

  firstFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x100u));
  secondFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x302u));
  thirdFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x504u));
  forthFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x706u));

  filt1Reg = _mm_load_si128((__m128i const *)filt1_global);
  filt2Reg = _mm_load_si128((__m128i const *)filt2_global);
  filt3Reg = _mm_load_si128((__m128i const *)filt3_global);
  filt4Reg = _mm_load_si128((__m128i const *)filt4_global);

  for (i = 0 ; i < output_height ; i++) {
    srcReg = _mm_loadu_si128((__m128i *)(src_ptr-3));
    srcRegFilt1= _mm_shuffle_epi8(srcReg, filt1Reg);
    srcRegFilt2= _mm_shuffle_epi8(srcReg, filt2Reg);
    srcRegFilt1 = _mm_maddubs_epi16(srcRegFilt1, firstFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, secondFilters);
    srcRegFilt3= _mm_shuffle_epi8(srcReg, filt3Reg);
    srcRegFilt4= _mm_shuffle_epi8(srcReg, filt4Reg);
    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, thirdFilters);
    srcRegFilt4 = _mm_maddubs_epi16(srcRegFilt4, forthFilters);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt2);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt4);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt3);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, addFilterReg127);
    srcRegFilt1 = _mm_srai_epi16(srcRegFilt1, 7);
    srcRegFilt1 = _mm_packus_epi16(srcRegFilt1, srcRegFilt1);
    src_ptr += src_pixels_per_line;
    _mm_storel_epi64((__m128i *)&output_ptr[0], srcRegFilt1);
    output_ptr += output_pitch;
  }
}

void vp9_filter_block1d16_h8_intrin_ssse3(unsigned char *src_ptr,
                                          unsigned int src_pixels_per_line,
                                          unsigned char *output_ptr,
                                          unsigned int output_pitch,
                                          unsigned int output_height,
                                          short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters,
          thirdFilters, forthFilters, srcReg1, srcReg2, srcRegFilt1_1,
          srcRegFilt2_1, srcRegFilt2, srcRegFilt3;
  __m128i filt1Reg, filt2Reg, filt3Reg, filt4Reg;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg =_mm_packs_epi16(filtersReg, filtersReg);

  firstFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x100u));
  secondFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x302u));
  thirdFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x504u));
  forthFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x706u));

  filt1Reg = _mm_load_si128((__m128i const *)filt1_global);
  filt2Reg = _mm_load_si128((__m128i const *)filt2_global);
  filt3Reg = _mm_load_si128((__m128i const *)filt3_global);
  filt4Reg = _mm_load_si128((__m128i const *)filt4_global);

  for (i = 0 ; i < output_height ; i++) {
    srcReg1 = _mm_loadu_si128((__m128i *)(src_ptr-3));
    srcRegFilt1_1 = _mm_shuffle_epi8(srcReg1, filt1Reg);
    srcRegFilt2 = _mm_shuffle_epi8(srcReg1, filt2Reg);

    srcRegFilt1_1 = _mm_maddubs_epi16(srcRegFilt1_1, firstFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, secondFilters);

    srcRegFilt1_1 = _mm_adds_epi16(srcRegFilt1_1, srcRegFilt2);

    srcRegFilt3 = _mm_shuffle_epi8(srcReg1, filt4Reg);
    srcRegFilt2 = _mm_shuffle_epi8(srcReg1, filt3Reg);

    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, forthFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, thirdFilters);

    srcRegFilt1_1 = _mm_adds_epi16(srcRegFilt1_1, srcRegFilt3);

    srcReg2 = _mm_loadu_si128((__m128i *)(src_ptr+5));

    srcRegFilt1_1 = _mm_adds_epi16(srcRegFilt1_1, srcRegFilt2);

    srcRegFilt2_1 = _mm_shuffle_epi8(srcReg2, filt1Reg);
    srcRegFilt2 = _mm_shuffle_epi8(srcReg2, filt2Reg);

    srcRegFilt2_1 = _mm_maddubs_epi16(srcRegFilt2_1, firstFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, secondFilters);

    srcRegFilt2_1 = _mm_adds_epi16(srcRegFilt2_1, srcRegFilt2);

    srcRegFilt3 = _mm_shuffle_epi8(srcReg2, filt4Reg);
    srcRegFilt2 = _mm_shuffle_epi8(srcReg2, filt3Reg);

    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, forthFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, thirdFilters);

    srcRegFilt2_1 = _mm_adds_epi16(srcRegFilt2_1, srcRegFilt3);
    srcRegFilt2_1 = _mm_adds_epi16(srcRegFilt2_1, srcRegFilt2);
    srcRegFilt1_1 = _mm_adds_epi16(srcRegFilt1_1, addFilterReg127);
    srcRegFilt2_1 = _mm_adds_epi16(srcRegFilt2_1, addFilterReg127);
    srcRegFilt1_1 = _mm_srai_epi16(srcRegFilt1_1, 7);
    srcRegFilt2_1 = _mm_srai_epi16(srcRegFilt2_1, 7);
    srcRegFilt1_1 = _mm_packus_epi16(srcRegFilt1_1, srcRegFilt2_1);

    src_ptr += src_pixels_per_line;
    _mm_store_si128((__m128i*)output_ptr, srcRegFilt1_1);
    output_ptr += output_pitch;
  }
}

void vp9_filter_block1d4_v8_intrin_ssse3(unsigned char *src_ptr,
                                         unsigned int src_pitch,
                                         unsigned char *output_ptr,
                                         unsigned int out_pitch,
                                         unsigned int output_height,
                                         short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters, srcRegFilt1,
          srcRegFilt2, srcRegFilt3, srcRegFilt4;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg =_mm_packs_epi16(filtersReg, filtersReg);
  firstFilters = _mm_shufflelo_epi16(filtersReg, 0);
  firstFilters = _mm_shufflehi_epi16(firstFilters, 0x55u);
  secondFilters = _mm_shufflelo_epi16(filtersReg, 0xAAu);
  secondFilters = _mm_shufflehi_epi16(secondFilters, 0xFFu);

  for (i = 0 ; i < output_height ; i++) {
    srcRegFilt1 = _mm_cvtsi32_si128(*((int*)&src_ptr[0]));
    srcRegFilt2 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch)[0]));

    srcRegFilt1 = _mm_unpacklo_epi8(srcRegFilt1, srcRegFilt2);

    srcRegFilt2 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*2)[0]));
    srcRegFilt3 =  _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*3)[0]));

    srcRegFilt2 = _mm_unpacklo_epi8(srcRegFilt2, srcRegFilt3);

    srcRegFilt3 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*4)[0]));
    srcRegFilt4 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*5)[0]));

    srcRegFilt3 = _mm_unpacklo_epi8(srcRegFilt3, srcRegFilt4);
    srcRegFilt1 = _mm_unpacklo_epi64(srcRegFilt1, srcRegFilt2);

    srcRegFilt4 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*6)[0]));
    srcRegFilt2 = _mm_cvtsi32_si128(*((int*)&(src_ptr+src_pitch*7)[0]));

    srcRegFilt4 = _mm_unpacklo_epi8(srcRegFilt4, srcRegFilt2);
    srcRegFilt3 = _mm_unpacklo_epi64(srcRegFilt3, srcRegFilt4);

    srcRegFilt1 = _mm_maddubs_epi16(srcRegFilt1, firstFilters);
    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, secondFilters);
    srcRegFilt2 = _mm_srli_si128(srcRegFilt1, 8);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, _mm_srli_si128(srcRegFilt3, 8));
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt2);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt3);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, addFilterReg127);
    srcRegFilt1 = _mm_srai_epi16(srcRegFilt1, 7);
    srcRegFilt1 = _mm_packus_epi16(srcRegFilt1, srcRegFilt1);

    src_ptr += src_pitch;
    *((int*)&output_ptr[0])= _mm_cvtsi128_si32(srcRegFilt1);
    output_ptr += out_pitch;
  }
}

void vp9_filter_block1d8_v8_intrin_ssse3(unsigned char *src_ptr,
                                         unsigned int src_pitch,
                                         unsigned char *output_ptr,
                                         unsigned int out_pitch,
                                         unsigned int output_height,
                                         short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters,
          thirdFilters, forthFilters, srcRegFilt1, srcRegFilt2, srcRegFilt3,
          srcRegFilt4, srcRegFilt5, srcRegFilt6;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg =_mm_packs_epi16(filtersReg, filtersReg);

  firstFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x100u));
  secondFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x302u));
  thirdFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x504u));
  forthFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x706u));

  for (i = 0 ; i < output_height ; i++) {
    srcRegFilt1 = _mm_loadl_epi64((__m128i*)&src_ptr[0]);
    srcRegFilt2 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch)[0]);

    srcRegFilt3 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*2)[0]);
    srcRegFilt4 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*3)[0]);

    srcRegFilt1 = _mm_unpacklo_epi8(srcRegFilt1, srcRegFilt2);

    srcRegFilt3 = _mm_unpacklo_epi8(srcRegFilt3, srcRegFilt4);

    srcRegFilt2 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*4)[0]);
    srcRegFilt4 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*5)[0]);

    srcRegFilt5 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*6)[0]);
    srcRegFilt6 = _mm_loadl_epi64((__m128i*)&(src_ptr+src_pitch*7)[0]);

    srcRegFilt2 = _mm_unpacklo_epi8(srcRegFilt2, srcRegFilt4);
    srcRegFilt5 = _mm_unpacklo_epi8(srcRegFilt5, srcRegFilt6);
    srcRegFilt1 = _mm_maddubs_epi16(srcRegFilt1, firstFilters);
    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, secondFilters);
    srcRegFilt2 = _mm_maddubs_epi16(srcRegFilt2, thirdFilters);
    srcRegFilt5 = _mm_maddubs_epi16(srcRegFilt5, forthFilters);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt5);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt3);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt2);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, addFilterReg127);
    srcRegFilt1 = _mm_srai_epi16(srcRegFilt1, 7);
    srcRegFilt1 = _mm_packus_epi16(srcRegFilt1, srcRegFilt1);

    src_ptr += src_pitch;
    _mm_storel_epi64((__m128i *)&output_ptr[0], srcRegFilt1);
    output_ptr += out_pitch;
  }
}

void vp9_filter_block1d16_v8_intrin_ssse3(unsigned char *src_ptr,
                                          unsigned int src_pitch,
                                          unsigned char *output_ptr,
                                          unsigned int out_pitch,
                                          unsigned int output_height,
                                          short *filter) {
  __m128i addFilterReg127, filtersReg, firstFilters, secondFilters,
          thirdFilters, forthFilters, srcRegFilt1, srcRegFilt2, srcRegFilt3,
          srcRegFilt4, srcRegFilt5, srcRegFilt6;
  unsigned int i;

  addFilterReg127 = _mm_shuffle_epi32(_mm_cvtsi32_si128((int)0x0400040u), 0);
  filtersReg = _mm_loadu_si128((__m128i *)filter);
  filtersReg =_mm_packs_epi16(filtersReg, filtersReg);

  firstFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x100u));
  secondFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x302u));
  thirdFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x504u));
  forthFilters = _mm_shuffle_epi8(filtersReg, _mm_set1_epi16(0x706u));

  for (i = 0 ; i < output_height ; i++) {
    srcRegFilt1 = _mm_loadu_si128((__m128i *)(src_ptr));
    srcRegFilt2 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch));

    srcRegFilt3 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*6));
    srcRegFilt4 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*7));

    srcRegFilt5 = _mm_unpacklo_epi8(srcRegFilt1, srcRegFilt2);
    srcRegFilt6 = _mm_unpacklo_epi8(srcRegFilt3, srcRegFilt4);
    srcRegFilt1 = _mm_unpackhi_epi8(srcRegFilt1, srcRegFilt2);
    srcRegFilt3 = _mm_unpackhi_epi8(srcRegFilt3, srcRegFilt4);

    srcRegFilt5 = _mm_maddubs_epi16(srcRegFilt5, firstFilters);
    srcRegFilt6 = _mm_maddubs_epi16(srcRegFilt6, forthFilters);
    srcRegFilt1 = _mm_maddubs_epi16(srcRegFilt1, firstFilters);
    srcRegFilt3 = _mm_maddubs_epi16(srcRegFilt3, forthFilters);

    srcRegFilt5 = _mm_adds_epi16(srcRegFilt5, srcRegFilt6);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt3);

    srcRegFilt2 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*2));
    srcRegFilt3 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*3));

    srcRegFilt4 = _mm_unpacklo_epi8(srcRegFilt2, srcRegFilt3);
    srcRegFilt6 = _mm_unpackhi_epi8(srcRegFilt2, srcRegFilt3);
    srcRegFilt4 = _mm_maddubs_epi16(srcRegFilt4, secondFilters);
    srcRegFilt6 = _mm_maddubs_epi16(srcRegFilt6, secondFilters);

    srcRegFilt5 = _mm_adds_epi16(srcRegFilt5, srcRegFilt4);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt6);

    srcRegFilt2 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*4));
    srcRegFilt3 = _mm_loadu_si128((__m128i *)(src_ptr+src_pitch*5));

    srcRegFilt4 = _mm_unpacklo_epi8(srcRegFilt2, srcRegFilt3);
    srcRegFilt6 = _mm_unpackhi_epi8(srcRegFilt2, srcRegFilt3);
    srcRegFilt4 = _mm_maddubs_epi16(srcRegFilt4, thirdFilters);
    srcRegFilt6 = _mm_maddubs_epi16(srcRegFilt6, thirdFilters);
    srcRegFilt5 = _mm_adds_epi16(srcRegFilt5, srcRegFilt4);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, srcRegFilt6);
    srcRegFilt5 = _mm_adds_epi16(srcRegFilt5, addFilterReg127);
    srcRegFilt1 = _mm_adds_epi16(srcRegFilt1, addFilterReg127);
    srcRegFilt5 = _mm_srai_epi16(srcRegFilt5, 7);
    srcRegFilt1 = _mm_srai_epi16(srcRegFilt1, 7);
    srcRegFilt1 = _mm_packus_epi16(srcRegFilt5, srcRegFilt1);

    src_ptr += src_pitch;
    _mm_store_si128((__m128i*)output_ptr, srcRegFilt1);
    output_ptr += out_pitch;
  }
}
