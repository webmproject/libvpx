#include "vp10/common/x86/vp10_txfm1d_sse4.h"

void vp10_fdct4_new_sse4_1(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 4;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[4];
  __m128i buf1[4];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[0], buf0[3]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[1], buf0[2]);

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[0], buf1[1], buf0[0],
                        buf0[1], bit);
    btf_32_sse4_1_type1(cospi[48], cospi[16], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[2];
    buf1[2] = buf0[1];
    buf1[3] = buf0[3];

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
  }
}

void vp10_fdct8_new_sse4_1(const __m128i* input, __m128i* output,
                           const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 8;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[8];
  __m128i buf1[8];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[0], buf0[7]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[1], buf0[6]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[2], buf0[5]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[4]);
    buf1[4] = _mm_sub_epi32(buf0[3], buf0[4]);

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[3]);
    buf0[3] = _mm_sub_epi32(buf1[0], buf1[3]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[2]);
    buf0[2] = _mm_sub_epi32(buf1[1], buf1[2]);
    buf0[4] = buf1[4];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[5], buf1[6], buf0[5],
                        buf0[6], bit);
    buf0[7] = buf1[7];

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf0[0], buf0[1], buf1[0],
                        buf1[1], bit);
    btf_32_sse4_1_type1(cospi[48], cospi[16], buf0[2], buf0[3], buf1[2],
                        buf1[3], bit);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[4], buf0[5]);
    buf1[6] = _mm_sub_epi32(buf0[7], buf0[6]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[6]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    btf_32_sse4_1_type1(cospi[56], cospi[8], buf1[4], buf1[7], buf0[4], buf0[7],
                        bit);
    btf_32_sse4_1_type1(cospi[24], cospi[40], buf1[5], buf1[6], buf0[5],
                        buf0[6], bit);

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[4];
    buf1[2] = buf0[2];
    buf1[3] = buf0[6];
    buf1[4] = buf0[1];
    buf1[5] = buf0[5];
    buf1[6] = buf0[3];
    buf1[7] = buf0[7];

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
  }
}

void vp10_fdct16_new_sse4_1(const __m128i* input, __m128i* output,
                            const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 16;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[16];
  __m128i buf1[16];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];
    buf0[8] = input[8 * col_num + col];
    buf0[9] = input[9 * col_num + col];
    buf0[10] = input[10 * col_num + col];
    buf0[11] = input[11 * col_num + col];
    buf0[12] = input[12 * col_num + col];
    buf0[13] = input[13 * col_num + col];
    buf0[14] = input[14 * col_num + col];
    buf0[15] = input[15 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[0], buf0[15]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[1], buf0[14]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[2], buf0[13]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[3], buf0[12]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[4], buf0[11]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[5], buf0[10]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[6], buf0[9]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[8]);
    buf1[8] = _mm_sub_epi32(buf0[7], buf0[8]);

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[7]);
    buf0[7] = _mm_sub_epi32(buf1[0], buf1[7]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[6]);
    buf0[6] = _mm_sub_epi32(buf1[1], buf1[6]);
    buf0[2] = _mm_add_epi32(buf1[2], buf1[5]);
    buf0[5] = _mm_sub_epi32(buf1[2], buf1[5]);
    buf0[3] = _mm_add_epi32(buf1[3], buf1[4]);
    buf0[4] = _mm_sub_epi32(buf1[3], buf1[4]);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[11], buf1[12], buf0[11],
                        buf0[12], bit);
    buf0[14] = buf1[14];
    buf0[15] = buf1[15];

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[0], buf0[3]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[1], buf0[2]);
    buf1[4] = buf0[4];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[5], buf0[6], buf1[5],
                        buf1[6], bit);
    buf1[7] = buf0[7];
    buf1[8] = _mm_add_epi32(buf0[8], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[8], buf0[11]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[9], buf0[10]);
    buf1[12] = _mm_sub_epi32(buf0[15], buf0[12]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[12]);
    buf1[13] = _mm_sub_epi32(buf0[14], buf0[13]);
    buf1[14] = _mm_add_epi32(buf0[14], buf0[13]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[0], buf1[1], buf0[0],
                        buf0[1], bit);
    btf_32_sse4_1_type1(cospi[48], cospi[16], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    buf0[4] = _mm_add_epi32(buf1[4], buf1[5]);
    buf0[5] = _mm_sub_epi32(buf1[4], buf1[5]);
    buf0[6] = _mm_sub_epi32(buf1[7], buf1[6]);
    buf0[7] = _mm_add_epi32(buf1[7], buf1[6]);
    buf0[8] = buf1[8];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[9], buf1[14], buf0[9],
                        buf0[14], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    buf0[11] = buf1[11];
    buf0[12] = buf1[12];
    buf0[15] = buf1[15];

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[1];
    buf1[2] = buf0[2];
    buf1[3] = buf0[3];
    btf_32_sse4_1_type1(cospi[56], cospi[8], buf0[4], buf0[7], buf1[4], buf1[7],
                        bit);
    btf_32_sse4_1_type1(cospi[24], cospi[40], buf0[5], buf0[6], buf1[5],
                        buf1[6], bit);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[8], buf0[9]);
    buf1[10] = _mm_sub_epi32(buf0[11], buf0[10]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[10]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[12], buf0[13]);
    buf1[14] = _mm_sub_epi32(buf0[15], buf0[14]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[14]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    btf_32_sse4_1_type1(cospi[60], cospi[4], buf1[8], buf1[15], buf0[8],
                        buf0[15], bit);
    btf_32_sse4_1_type1(cospi[28], cospi[36], buf1[9], buf1[14], buf0[9],
                        buf0[14], bit);
    btf_32_sse4_1_type1(cospi[44], cospi[20], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    btf_32_sse4_1_type1(cospi[12], cospi[52], buf1[11], buf1[12], buf0[11],
                        buf0[12], bit);

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[8];
    buf1[2] = buf0[4];
    buf1[3] = buf0[12];
    buf1[4] = buf0[2];
    buf1[5] = buf0[10];
    buf1[6] = buf0[6];
    buf1[7] = buf0[14];
    buf1[8] = buf0[1];
    buf1[9] = buf0[9];
    buf1[10] = buf0[5];
    buf1[11] = buf0[13];
    buf1[12] = buf0[3];
    buf1[13] = buf0[11];
    buf1[14] = buf0[7];
    buf1[15] = buf0[15];

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
    output[8 * col_num + col] = buf1[8];
    output[9 * col_num + col] = buf1[9];
    output[10 * col_num + col] = buf1[10];
    output[11 * col_num + col] = buf1[11];
    output[12 * col_num + col] = buf1[12];
    output[13 * col_num + col] = buf1[13];
    output[14 * col_num + col] = buf1[14];
    output[15 * col_num + col] = buf1[15];
  }
}

void vp10_fdct32_new_sse4_1(const __m128i* input, __m128i* output,
                            const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 32;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[32];
  __m128i buf1[32];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];
    buf0[8] = input[8 * col_num + col];
    buf0[9] = input[9 * col_num + col];
    buf0[10] = input[10 * col_num + col];
    buf0[11] = input[11 * col_num + col];
    buf0[12] = input[12 * col_num + col];
    buf0[13] = input[13 * col_num + col];
    buf0[14] = input[14 * col_num + col];
    buf0[15] = input[15 * col_num + col];
    buf0[16] = input[16 * col_num + col];
    buf0[17] = input[17 * col_num + col];
    buf0[18] = input[18 * col_num + col];
    buf0[19] = input[19 * col_num + col];
    buf0[20] = input[20 * col_num + col];
    buf0[21] = input[21 * col_num + col];
    buf0[22] = input[22 * col_num + col];
    buf0[23] = input[23 * col_num + col];
    buf0[24] = input[24 * col_num + col];
    buf0[25] = input[25 * col_num + col];
    buf0[26] = input[26 * col_num + col];
    buf0[27] = input[27 * col_num + col];
    buf0[28] = input[28 * col_num + col];
    buf0[29] = input[29 * col_num + col];
    buf0[30] = input[30 * col_num + col];
    buf0[31] = input[31 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[31]);
    buf1[31] = _mm_sub_epi32(buf0[0], buf0[31]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[30]);
    buf1[30] = _mm_sub_epi32(buf0[1], buf0[30]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[29]);
    buf1[29] = _mm_sub_epi32(buf0[2], buf0[29]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[28]);
    buf1[28] = _mm_sub_epi32(buf0[3], buf0[28]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[27]);
    buf1[27] = _mm_sub_epi32(buf0[4], buf0[27]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[26]);
    buf1[26] = _mm_sub_epi32(buf0[5], buf0[26]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[25]);
    buf1[25] = _mm_sub_epi32(buf0[6], buf0[25]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[24]);
    buf1[24] = _mm_sub_epi32(buf0[7], buf0[24]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[23]);
    buf1[23] = _mm_sub_epi32(buf0[8], buf0[23]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[22]);
    buf1[22] = _mm_sub_epi32(buf0[9], buf0[22]);
    buf1[10] = _mm_add_epi32(buf0[10], buf0[21]);
    buf1[21] = _mm_sub_epi32(buf0[10], buf0[21]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[20]);
    buf1[20] = _mm_sub_epi32(buf0[11], buf0[20]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[19]);
    buf1[19] = _mm_sub_epi32(buf0[12], buf0[19]);
    buf1[13] = _mm_add_epi32(buf0[13], buf0[18]);
    buf1[18] = _mm_sub_epi32(buf0[13], buf0[18]);
    buf1[14] = _mm_add_epi32(buf0[14], buf0[17]);
    buf1[17] = _mm_sub_epi32(buf0[14], buf0[17]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[16]);
    buf1[16] = _mm_sub_epi32(buf0[15], buf0[16]);

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[15]);
    buf0[15] = _mm_sub_epi32(buf1[0], buf1[15]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[14]);
    buf0[14] = _mm_sub_epi32(buf1[1], buf1[14]);
    buf0[2] = _mm_add_epi32(buf1[2], buf1[13]);
    buf0[13] = _mm_sub_epi32(buf1[2], buf1[13]);
    buf0[3] = _mm_add_epi32(buf1[3], buf1[12]);
    buf0[12] = _mm_sub_epi32(buf1[3], buf1[12]);
    buf0[4] = _mm_add_epi32(buf1[4], buf1[11]);
    buf0[11] = _mm_sub_epi32(buf1[4], buf1[11]);
    buf0[5] = _mm_add_epi32(buf1[5], buf1[10]);
    buf0[10] = _mm_sub_epi32(buf1[5], buf1[10]);
    buf0[6] = _mm_add_epi32(buf1[6], buf1[9]);
    buf0[9] = _mm_sub_epi32(buf1[6], buf1[9]);
    buf0[7] = _mm_add_epi32(buf1[7], buf1[8]);
    buf0[8] = _mm_sub_epi32(buf1[7], buf1[8]);
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    buf0[18] = buf1[18];
    buf0[19] = buf1[19];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[20], buf1[27], buf0[20],
                        buf0[27], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[21], buf1[26], buf0[21],
                        buf0[26], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[22], buf1[25], buf0[22],
                        buf0[25], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[23], buf1[24], buf0[23],
                        buf0[24], bit);
    buf0[28] = buf1[28];
    buf0[29] = buf1[29];
    buf0[30] = buf1[30];
    buf0[31] = buf1[31];

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[0], buf0[7]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[1], buf0[6]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[2], buf0[5]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[4]);
    buf1[4] = _mm_sub_epi32(buf0[3], buf0[4]);
    buf1[8] = buf0[8];
    buf1[9] = buf0[9];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[10], buf0[13], buf1[10],
                        buf1[13], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[11], buf0[12], buf1[11],
                        buf1[12], bit);
    buf1[14] = buf0[14];
    buf1[15] = buf0[15];
    buf1[16] = _mm_add_epi32(buf0[16], buf0[23]);
    buf1[23] = _mm_sub_epi32(buf0[16], buf0[23]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[22]);
    buf1[22] = _mm_sub_epi32(buf0[17], buf0[22]);
    buf1[18] = _mm_add_epi32(buf0[18], buf0[21]);
    buf1[21] = _mm_sub_epi32(buf0[18], buf0[21]);
    buf1[19] = _mm_add_epi32(buf0[19], buf0[20]);
    buf1[20] = _mm_sub_epi32(buf0[19], buf0[20]);
    buf1[24] = _mm_sub_epi32(buf0[31], buf0[24]);
    buf1[31] = _mm_add_epi32(buf0[31], buf0[24]);
    buf1[25] = _mm_sub_epi32(buf0[30], buf0[25]);
    buf1[30] = _mm_add_epi32(buf0[30], buf0[25]);
    buf1[26] = _mm_sub_epi32(buf0[29], buf0[26]);
    buf1[29] = _mm_add_epi32(buf0[29], buf0[26]);
    buf1[27] = _mm_sub_epi32(buf0[28], buf0[27]);
    buf1[28] = _mm_add_epi32(buf0[28], buf0[27]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[3]);
    buf0[3] = _mm_sub_epi32(buf1[0], buf1[3]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[2]);
    buf0[2] = _mm_sub_epi32(buf1[1], buf1[2]);
    buf0[4] = buf1[4];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[5], buf1[6], buf0[5],
                        buf0[6], bit);
    buf0[7] = buf1[7];
    buf0[8] = _mm_add_epi32(buf1[8], buf1[11]);
    buf0[11] = _mm_sub_epi32(buf1[8], buf1[11]);
    buf0[9] = _mm_add_epi32(buf1[9], buf1[10]);
    buf0[10] = _mm_sub_epi32(buf1[9], buf1[10]);
    buf0[12] = _mm_sub_epi32(buf1[15], buf1[12]);
    buf0[15] = _mm_add_epi32(buf1[15], buf1[12]);
    buf0[13] = _mm_sub_epi32(buf1[14], buf1[13]);
    buf0[14] = _mm_add_epi32(buf1[14], buf1[13]);
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[18], buf1[29], buf0[18],
                        buf0[29], bit);
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[19], buf1[28], buf0[19],
                        buf0[28], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[20], buf1[27], buf0[20],
                        buf0[27], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[21], buf1[26], buf0[21],
                        buf0[26], bit);
    buf0[22] = buf1[22];
    buf0[23] = buf1[23];
    buf0[24] = buf1[24];
    buf0[25] = buf1[25];
    buf0[30] = buf1[30];
    buf0[31] = buf1[31];

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf0[0], buf0[1], buf1[0],
                        buf1[1], bit);
    btf_32_sse4_1_type1(cospi[48], cospi[16], buf0[2], buf0[3], buf1[2],
                        buf1[3], bit);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[4], buf0[5]);
    buf1[6] = _mm_sub_epi32(buf0[7], buf0[6]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[6]);
    buf1[8] = buf0[8];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf0[9], buf0[14], buf1[9],
                        buf1[14], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf0[10], buf0[13], buf1[10],
                        buf1[13], bit);
    buf1[11] = buf0[11];
    buf1[12] = buf0[12];
    buf1[15] = buf0[15];
    buf1[16] = _mm_add_epi32(buf0[16], buf0[19]);
    buf1[19] = _mm_sub_epi32(buf0[16], buf0[19]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[18]);
    buf1[18] = _mm_sub_epi32(buf0[17], buf0[18]);
    buf1[20] = _mm_sub_epi32(buf0[23], buf0[20]);
    buf1[23] = _mm_add_epi32(buf0[23], buf0[20]);
    buf1[21] = _mm_sub_epi32(buf0[22], buf0[21]);
    buf1[22] = _mm_add_epi32(buf0[22], buf0[21]);
    buf1[24] = _mm_add_epi32(buf0[24], buf0[27]);
    buf1[27] = _mm_sub_epi32(buf0[24], buf0[27]);
    buf1[25] = _mm_add_epi32(buf0[25], buf0[26]);
    buf1[26] = _mm_sub_epi32(buf0[25], buf0[26]);
    buf1[28] = _mm_sub_epi32(buf0[31], buf0[28]);
    buf1[31] = _mm_add_epi32(buf0[31], buf0[28]);
    buf1[29] = _mm_sub_epi32(buf0[30], buf0[29]);
    buf1[30] = _mm_add_epi32(buf0[30], buf0[29]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    btf_32_sse4_1_type1(cospi[56], cospi[8], buf1[4], buf1[7], buf0[4], buf0[7],
                        bit);
    btf_32_sse4_1_type1(cospi[24], cospi[40], buf1[5], buf1[6], buf0[5],
                        buf0[6], bit);
    buf0[8] = _mm_add_epi32(buf1[8], buf1[9]);
    buf0[9] = _mm_sub_epi32(buf1[8], buf1[9]);
    buf0[10] = _mm_sub_epi32(buf1[11], buf1[10]);
    buf0[11] = _mm_add_epi32(buf1[11], buf1[10]);
    buf0[12] = _mm_add_epi32(buf1[12], buf1[13]);
    buf0[13] = _mm_sub_epi32(buf1[12], buf1[13]);
    buf0[14] = _mm_sub_epi32(buf1[15], buf1[14]);
    buf0[15] = _mm_add_epi32(buf1[15], buf1[14]);
    buf0[16] = buf1[16];
    btf_32_sse4_1_type0(-cospi[8], cospi[56], buf1[17], buf1[30], buf0[17],
                        buf0[30], bit);
    btf_32_sse4_1_type0(-cospi[56], -cospi[8], buf1[18], buf1[29], buf0[18],
                        buf0[29], bit);
    buf0[19] = buf1[19];
    buf0[20] = buf1[20];
    btf_32_sse4_1_type0(-cospi[40], cospi[24], buf1[21], buf1[26], buf0[21],
                        buf0[26], bit);
    btf_32_sse4_1_type0(-cospi[24], -cospi[40], buf1[22], buf1[25], buf0[22],
                        buf0[25], bit);
    buf0[23] = buf1[23];
    buf0[24] = buf1[24];
    buf0[27] = buf1[27];
    buf0[28] = buf1[28];
    buf0[31] = buf1[31];

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[1];
    buf1[2] = buf0[2];
    buf1[3] = buf0[3];
    buf1[4] = buf0[4];
    buf1[5] = buf0[5];
    buf1[6] = buf0[6];
    buf1[7] = buf0[7];
    btf_32_sse4_1_type1(cospi[60], cospi[4], buf0[8], buf0[15], buf1[8],
                        buf1[15], bit);
    btf_32_sse4_1_type1(cospi[28], cospi[36], buf0[9], buf0[14], buf1[9],
                        buf1[14], bit);
    btf_32_sse4_1_type1(cospi[44], cospi[20], buf0[10], buf0[13], buf1[10],
                        buf1[13], bit);
    btf_32_sse4_1_type1(cospi[12], cospi[52], buf0[11], buf0[12], buf1[11],
                        buf1[12], bit);
    buf1[16] = _mm_add_epi32(buf0[16], buf0[17]);
    buf1[17] = _mm_sub_epi32(buf0[16], buf0[17]);
    buf1[18] = _mm_sub_epi32(buf0[19], buf0[18]);
    buf1[19] = _mm_add_epi32(buf0[19], buf0[18]);
    buf1[20] = _mm_add_epi32(buf0[20], buf0[21]);
    buf1[21] = _mm_sub_epi32(buf0[20], buf0[21]);
    buf1[22] = _mm_sub_epi32(buf0[23], buf0[22]);
    buf1[23] = _mm_add_epi32(buf0[23], buf0[22]);
    buf1[24] = _mm_add_epi32(buf0[24], buf0[25]);
    buf1[25] = _mm_sub_epi32(buf0[24], buf0[25]);
    buf1[26] = _mm_sub_epi32(buf0[27], buf0[26]);
    buf1[27] = _mm_add_epi32(buf0[27], buf0[26]);
    buf1[28] = _mm_add_epi32(buf0[28], buf0[29]);
    buf1[29] = _mm_sub_epi32(buf0[28], buf0[29]);
    buf1[30] = _mm_sub_epi32(buf0[31], buf0[30]);
    buf1[31] = _mm_add_epi32(buf0[31], buf0[30]);

    // stage 8
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    buf0[10] = buf1[10];
    buf0[11] = buf1[11];
    buf0[12] = buf1[12];
    buf0[13] = buf1[13];
    buf0[14] = buf1[14];
    buf0[15] = buf1[15];
    btf_32_sse4_1_type1(cospi[62], cospi[2], buf1[16], buf1[31], buf0[16],
                        buf0[31], bit);
    btf_32_sse4_1_type1(cospi[30], cospi[34], buf1[17], buf1[30], buf0[17],
                        buf0[30], bit);
    btf_32_sse4_1_type1(cospi[46], cospi[18], buf1[18], buf1[29], buf0[18],
                        buf0[29], bit);
    btf_32_sse4_1_type1(cospi[14], cospi[50], buf1[19], buf1[28], buf0[19],
                        buf0[28], bit);
    btf_32_sse4_1_type1(cospi[54], cospi[10], buf1[20], buf1[27], buf0[20],
                        buf0[27], bit);
    btf_32_sse4_1_type1(cospi[22], cospi[42], buf1[21], buf1[26], buf0[21],
                        buf0[26], bit);
    btf_32_sse4_1_type1(cospi[38], cospi[26], buf1[22], buf1[25], buf0[22],
                        buf0[25], bit);
    btf_32_sse4_1_type1(cospi[6], cospi[58], buf1[23], buf1[24], buf0[23],
                        buf0[24], bit);

    // stage 9
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[16];
    buf1[2] = buf0[8];
    buf1[3] = buf0[24];
    buf1[4] = buf0[4];
    buf1[5] = buf0[20];
    buf1[6] = buf0[12];
    buf1[7] = buf0[28];
    buf1[8] = buf0[2];
    buf1[9] = buf0[18];
    buf1[10] = buf0[10];
    buf1[11] = buf0[26];
    buf1[12] = buf0[6];
    buf1[13] = buf0[22];
    buf1[14] = buf0[14];
    buf1[15] = buf0[30];
    buf1[16] = buf0[1];
    buf1[17] = buf0[17];
    buf1[18] = buf0[9];
    buf1[19] = buf0[25];
    buf1[20] = buf0[5];
    buf1[21] = buf0[21];
    buf1[22] = buf0[13];
    buf1[23] = buf0[29];
    buf1[24] = buf0[3];
    buf1[25] = buf0[19];
    buf1[26] = buf0[11];
    buf1[27] = buf0[27];
    buf1[28] = buf0[7];
    buf1[29] = buf0[23];
    buf1[30] = buf0[15];
    buf1[31] = buf0[31];

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
    output[8 * col_num + col] = buf1[8];
    output[9 * col_num + col] = buf1[9];
    output[10 * col_num + col] = buf1[10];
    output[11 * col_num + col] = buf1[11];
    output[12 * col_num + col] = buf1[12];
    output[13 * col_num + col] = buf1[13];
    output[14 * col_num + col] = buf1[14];
    output[15 * col_num + col] = buf1[15];
    output[16 * col_num + col] = buf1[16];
    output[17 * col_num + col] = buf1[17];
    output[18 * col_num + col] = buf1[18];
    output[19 * col_num + col] = buf1[19];
    output[20 * col_num + col] = buf1[20];
    output[21 * col_num + col] = buf1[21];
    output[22 * col_num + col] = buf1[22];
    output[23 * col_num + col] = buf1[23];
    output[24 * col_num + col] = buf1[24];
    output[25 * col_num + col] = buf1[25];
    output[26 * col_num + col] = buf1[26];
    output[27 * col_num + col] = buf1[27];
    output[28 * col_num + col] = buf1[28];
    output[29 * col_num + col] = buf1[29];
    output[30 * col_num + col] = buf1[30];
    output[31 * col_num + col] = buf1[31];
  }
}

void vp10_fadst4_new_sse4_1(const __m128i* input, __m128i* output,
                            const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 4;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[4];
  __m128i buf1[4];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[3];
    buf1[1] = buf0[0];
    buf1[2] = buf0[1];
    buf1[3] = buf0[2];

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[8], cospi[56], buf1[0], buf1[1], buf0[0], buf0[1],
                        bit);
    btf_32_sse4_1_type0(cospi[40], cospi[24], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[0], buf0[2]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[1], buf0[3]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[2]);
    buf1[2] = buf0[3];
    buf1[3] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[1]);

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
  }
}

void vp10_fadst8_new_sse4_1(const __m128i* input, __m128i* output,
                            const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 8;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[8];
  __m128i buf1[8];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[7];
    buf1[1] = buf0[0];
    buf1[2] = buf0[5];
    buf1[3] = buf0[2];
    buf1[4] = buf0[3];
    buf1[5] = buf0[4];
    buf1[6] = buf0[1];
    buf1[7] = buf0[6];

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[4], cospi[60], buf1[0], buf1[1], buf0[0], buf0[1],
                        bit);
    btf_32_sse4_1_type0(cospi[20], cospi[44], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    btf_32_sse4_1_type0(cospi[36], cospi[28], buf1[4], buf1[5], buf0[4],
                        buf0[5], bit);
    btf_32_sse4_1_type0(cospi[52], cospi[12], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[4]);
    buf1[4] = _mm_sub_epi32(buf0[0], buf0[4]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[1], buf0[5]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[2], buf0[6]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[3], buf0[7]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[4], buf1[5], buf0[4],
                        buf0[5], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[0], buf0[2]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[1], buf0[3]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[4], buf0[6]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[5], buf0[7]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[4]);
    buf1[2] = buf0[6];
    buf1[3] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[2]);
    buf1[4] = buf0[3];
    buf1[5] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[7]);
    buf1[6] = buf0[5];
    buf1[7] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[1]);

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
  }
}

void vp10_fadst16_new_sse4_1(const __m128i* input, __m128i* output,
                             const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 16;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[16];
  __m128i buf1[16];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];
    buf0[8] = input[8 * col_num + col];
    buf0[9] = input[9 * col_num + col];
    buf0[10] = input[10 * col_num + col];
    buf0[11] = input[11 * col_num + col];
    buf0[12] = input[12 * col_num + col];
    buf0[13] = input[13 * col_num + col];
    buf0[14] = input[14 * col_num + col];
    buf0[15] = input[15 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[15];
    buf1[1] = buf0[0];
    buf1[2] = buf0[13];
    buf1[3] = buf0[2];
    buf1[4] = buf0[11];
    buf1[5] = buf0[4];
    buf1[6] = buf0[9];
    buf1[7] = buf0[6];
    buf1[8] = buf0[7];
    buf1[9] = buf0[8];
    buf1[10] = buf0[5];
    buf1[11] = buf0[10];
    buf1[12] = buf0[3];
    buf1[13] = buf0[12];
    buf1[14] = buf0[1];
    buf1[15] = buf0[14];

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[2], cospi[62], buf1[0], buf1[1], buf0[0], buf0[1],
                        bit);
    btf_32_sse4_1_type0(cospi[10], cospi[54], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    btf_32_sse4_1_type0(cospi[18], cospi[46], buf1[4], buf1[5], buf0[4],
                        buf0[5], bit);
    btf_32_sse4_1_type0(cospi[26], cospi[38], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    btf_32_sse4_1_type0(cospi[34], cospi[30], buf1[8], buf1[9], buf0[8],
                        buf0[9], bit);
    btf_32_sse4_1_type0(cospi[42], cospi[22], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    btf_32_sse4_1_type0(cospi[50], cospi[14], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(cospi[58], cospi[6], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[8]);
    buf1[8] = _mm_sub_epi32(buf0[0], buf0[8]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[1], buf0[9]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[2], buf0[10]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[3], buf0[11]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[4], buf0[12]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[5], buf0[13]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[6], buf0[14]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[7], buf0[15]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    btf_32_sse4_1_type0(cospi[8], cospi[56], buf1[8], buf1[9], buf0[8], buf0[9],
                        bit);
    btf_32_sse4_1_type0(cospi[40], cospi[24], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    btf_32_sse4_1_type0(-cospi[56], cospi[8], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[24], cospi[40], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[4]);
    buf1[4] = _mm_sub_epi32(buf0[0], buf0[4]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[1], buf0[5]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[2], buf0[6]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[3], buf0[7]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[8], buf0[12]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[9], buf0[13]);
    buf1[10] = _mm_add_epi32(buf0[10], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[10], buf0[14]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[11], buf0[15]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[4], buf1[5], buf0[4],
                        buf0[5], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    buf0[10] = buf1[10];
    buf0[11] = buf1[11];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[0], buf0[2]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[1], buf0[3]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[4], buf0[6]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[5], buf0[7]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[8], buf0[10]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[9], buf0[11]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[12], buf0[14]);
    buf1[13] = _mm_add_epi32(buf0[13], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[13], buf0[15]);

    // stage 8
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    buf0[12] = buf1[12];
    buf0[13] = buf1[13];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);

    // stage 9
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[8]);
    buf1[2] = buf0[12];
    buf1[3] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[4]);
    buf1[4] = buf0[6];
    buf1[5] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[14]);
    buf1[6] = buf0[10];
    buf1[7] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[2]);
    buf1[8] = buf0[3];
    buf1[9] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[11]);
    buf1[10] = buf0[15];
    buf1[11] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[7]);
    buf1[12] = buf0[5];
    buf1[13] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[13]);
    buf1[14] = buf0[9];
    buf1[15] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[1]);

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
    output[8 * col_num + col] = buf1[8];
    output[9 * col_num + col] = buf1[9];
    output[10 * col_num + col] = buf1[10];
    output[11 * col_num + col] = buf1[11];
    output[12 * col_num + col] = buf1[12];
    output[13 * col_num + col] = buf1[13];
    output[14 * col_num + col] = buf1[14];
    output[15 * col_num + col] = buf1[15];
  }
}

void vp10_fadst32_new_sse4_1(const __m128i* input, __m128i* output,
                             const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 32;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[32];
  __m128i buf1[32];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];
    buf0[8] = input[8 * col_num + col];
    buf0[9] = input[9 * col_num + col];
    buf0[10] = input[10 * col_num + col];
    buf0[11] = input[11 * col_num + col];
    buf0[12] = input[12 * col_num + col];
    buf0[13] = input[13 * col_num + col];
    buf0[14] = input[14 * col_num + col];
    buf0[15] = input[15 * col_num + col];
    buf0[16] = input[16 * col_num + col];
    buf0[17] = input[17 * col_num + col];
    buf0[18] = input[18 * col_num + col];
    buf0[19] = input[19 * col_num + col];
    buf0[20] = input[20 * col_num + col];
    buf0[21] = input[21 * col_num + col];
    buf0[22] = input[22 * col_num + col];
    buf0[23] = input[23 * col_num + col];
    buf0[24] = input[24 * col_num + col];
    buf0[25] = input[25 * col_num + col];
    buf0[26] = input[26 * col_num + col];
    buf0[27] = input[27 * col_num + col];
    buf0[28] = input[28 * col_num + col];
    buf0[29] = input[29 * col_num + col];
    buf0[30] = input[30 * col_num + col];
    buf0[31] = input[31 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[31];
    buf1[1] = buf0[0];
    buf1[2] = buf0[29];
    buf1[3] = buf0[2];
    buf1[4] = buf0[27];
    buf1[5] = buf0[4];
    buf1[6] = buf0[25];
    buf1[7] = buf0[6];
    buf1[8] = buf0[23];
    buf1[9] = buf0[8];
    buf1[10] = buf0[21];
    buf1[11] = buf0[10];
    buf1[12] = buf0[19];
    buf1[13] = buf0[12];
    buf1[14] = buf0[17];
    buf1[15] = buf0[14];
    buf1[16] = buf0[15];
    buf1[17] = buf0[16];
    buf1[18] = buf0[13];
    buf1[19] = buf0[18];
    buf1[20] = buf0[11];
    buf1[21] = buf0[20];
    buf1[22] = buf0[9];
    buf1[23] = buf0[22];
    buf1[24] = buf0[7];
    buf1[25] = buf0[24];
    buf1[26] = buf0[5];
    buf1[27] = buf0[26];
    buf1[28] = buf0[3];
    buf1[29] = buf0[28];
    buf1[30] = buf0[1];
    buf1[31] = buf0[30];

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[1], cospi[63], buf1[0], buf1[1], buf0[0], buf0[1],
                        bit);
    btf_32_sse4_1_type0(cospi[5], cospi[59], buf1[2], buf1[3], buf0[2], buf0[3],
                        bit);
    btf_32_sse4_1_type0(cospi[9], cospi[55], buf1[4], buf1[5], buf0[4], buf0[5],
                        bit);
    btf_32_sse4_1_type0(cospi[13], cospi[51], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    btf_32_sse4_1_type0(cospi[17], cospi[47], buf1[8], buf1[9], buf0[8],
                        buf0[9], bit);
    btf_32_sse4_1_type0(cospi[21], cospi[43], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    btf_32_sse4_1_type0(cospi[25], cospi[39], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(cospi[29], cospi[35], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);
    btf_32_sse4_1_type0(cospi[33], cospi[31], buf1[16], buf1[17], buf0[16],
                        buf0[17], bit);
    btf_32_sse4_1_type0(cospi[37], cospi[27], buf1[18], buf1[19], buf0[18],
                        buf0[19], bit);
    btf_32_sse4_1_type0(cospi[41], cospi[23], buf1[20], buf1[21], buf0[20],
                        buf0[21], bit);
    btf_32_sse4_1_type0(cospi[45], cospi[19], buf1[22], buf1[23], buf0[22],
                        buf0[23], bit);
    btf_32_sse4_1_type0(cospi[49], cospi[15], buf1[24], buf1[25], buf0[24],
                        buf0[25], bit);
    btf_32_sse4_1_type0(cospi[53], cospi[11], buf1[26], buf1[27], buf0[26],
                        buf0[27], bit);
    btf_32_sse4_1_type0(cospi[57], cospi[7], buf1[28], buf1[29], buf0[28],
                        buf0[29], bit);
    btf_32_sse4_1_type0(cospi[61], cospi[3], buf1[30], buf1[31], buf0[30],
                        buf0[31], bit);

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[16]);
    buf1[16] = _mm_sub_epi32(buf0[0], buf0[16]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[17]);
    buf1[17] = _mm_sub_epi32(buf0[1], buf0[17]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[18]);
    buf1[18] = _mm_sub_epi32(buf0[2], buf0[18]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[19]);
    buf1[19] = _mm_sub_epi32(buf0[3], buf0[19]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[20]);
    buf1[20] = _mm_sub_epi32(buf0[4], buf0[20]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[21]);
    buf1[21] = _mm_sub_epi32(buf0[5], buf0[21]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[22]);
    buf1[22] = _mm_sub_epi32(buf0[6], buf0[22]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[23]);
    buf1[23] = _mm_sub_epi32(buf0[7], buf0[23]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[24]);
    buf1[24] = _mm_sub_epi32(buf0[8], buf0[24]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[25]);
    buf1[25] = _mm_sub_epi32(buf0[9], buf0[25]);
    buf1[10] = _mm_add_epi32(buf0[10], buf0[26]);
    buf1[26] = _mm_sub_epi32(buf0[10], buf0[26]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[27]);
    buf1[27] = _mm_sub_epi32(buf0[11], buf0[27]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[28]);
    buf1[28] = _mm_sub_epi32(buf0[12], buf0[28]);
    buf1[13] = _mm_add_epi32(buf0[13], buf0[29]);
    buf1[29] = _mm_sub_epi32(buf0[13], buf0[29]);
    buf1[14] = _mm_add_epi32(buf0[14], buf0[30]);
    buf1[30] = _mm_sub_epi32(buf0[14], buf0[30]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[31]);
    buf1[31] = _mm_sub_epi32(buf0[15], buf0[31]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    buf0[10] = buf1[10];
    buf0[11] = buf1[11];
    buf0[12] = buf1[12];
    buf0[13] = buf1[13];
    buf0[14] = buf1[14];
    buf0[15] = buf1[15];
    btf_32_sse4_1_type0(cospi[4], cospi[60], buf1[16], buf1[17], buf0[16],
                        buf0[17], bit);
    btf_32_sse4_1_type0(cospi[20], cospi[44], buf1[18], buf1[19], buf0[18],
                        buf0[19], bit);
    btf_32_sse4_1_type0(cospi[36], cospi[28], buf1[20], buf1[21], buf0[20],
                        buf0[21], bit);
    btf_32_sse4_1_type0(cospi[52], cospi[12], buf1[22], buf1[23], buf0[22],
                        buf0[23], bit);
    btf_32_sse4_1_type0(-cospi[60], cospi[4], buf1[24], buf1[25], buf0[24],
                        buf0[25], bit);
    btf_32_sse4_1_type0(-cospi[44], cospi[20], buf1[26], buf1[27], buf0[26],
                        buf0[27], bit);
    btf_32_sse4_1_type0(-cospi[28], cospi[36], buf1[28], buf1[29], buf0[28],
                        buf0[29], bit);
    btf_32_sse4_1_type0(-cospi[12], cospi[52], buf1[30], buf1[31], buf0[30],
                        buf0[31], bit);

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[8]);
    buf1[8] = _mm_sub_epi32(buf0[0], buf0[8]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[1], buf0[9]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[2], buf0[10]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[3], buf0[11]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[4], buf0[12]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[5], buf0[13]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[6], buf0[14]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[7], buf0[15]);
    buf1[16] = _mm_add_epi32(buf0[16], buf0[24]);
    buf1[24] = _mm_sub_epi32(buf0[16], buf0[24]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[25]);
    buf1[25] = _mm_sub_epi32(buf0[17], buf0[25]);
    buf1[18] = _mm_add_epi32(buf0[18], buf0[26]);
    buf1[26] = _mm_sub_epi32(buf0[18], buf0[26]);
    buf1[19] = _mm_add_epi32(buf0[19], buf0[27]);
    buf1[27] = _mm_sub_epi32(buf0[19], buf0[27]);
    buf1[20] = _mm_add_epi32(buf0[20], buf0[28]);
    buf1[28] = _mm_sub_epi32(buf0[20], buf0[28]);
    buf1[21] = _mm_add_epi32(buf0[21], buf0[29]);
    buf1[29] = _mm_sub_epi32(buf0[21], buf0[29]);
    buf1[22] = _mm_add_epi32(buf0[22], buf0[30]);
    buf1[30] = _mm_sub_epi32(buf0[22], buf0[30]);
    buf1[23] = _mm_add_epi32(buf0[23], buf0[31]);
    buf1[31] = _mm_sub_epi32(buf0[23], buf0[31]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    btf_32_sse4_1_type0(cospi[8], cospi[56], buf1[8], buf1[9], buf0[8], buf0[9],
                        bit);
    btf_32_sse4_1_type0(cospi[40], cospi[24], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    btf_32_sse4_1_type0(-cospi[56], cospi[8], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[24], cospi[40], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    buf0[18] = buf1[18];
    buf0[19] = buf1[19];
    buf0[20] = buf1[20];
    buf0[21] = buf1[21];
    buf0[22] = buf1[22];
    buf0[23] = buf1[23];
    btf_32_sse4_1_type0(cospi[8], cospi[56], buf1[24], buf1[25], buf0[24],
                        buf0[25], bit);
    btf_32_sse4_1_type0(cospi[40], cospi[24], buf1[26], buf1[27], buf0[26],
                        buf0[27], bit);
    btf_32_sse4_1_type0(-cospi[56], cospi[8], buf1[28], buf1[29], buf0[28],
                        buf0[29], bit);
    btf_32_sse4_1_type0(-cospi[24], cospi[40], buf1[30], buf1[31], buf0[30],
                        buf0[31], bit);

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[4]);
    buf1[4] = _mm_sub_epi32(buf0[0], buf0[4]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[5]);
    buf1[5] = _mm_sub_epi32(buf0[1], buf0[5]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[2], buf0[6]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[3], buf0[7]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[8], buf0[12]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[9], buf0[13]);
    buf1[10] = _mm_add_epi32(buf0[10], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[10], buf0[14]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[11], buf0[15]);
    buf1[16] = _mm_add_epi32(buf0[16], buf0[20]);
    buf1[20] = _mm_sub_epi32(buf0[16], buf0[20]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[21]);
    buf1[21] = _mm_sub_epi32(buf0[17], buf0[21]);
    buf1[18] = _mm_add_epi32(buf0[18], buf0[22]);
    buf1[22] = _mm_sub_epi32(buf0[18], buf0[22]);
    buf1[19] = _mm_add_epi32(buf0[19], buf0[23]);
    buf1[23] = _mm_sub_epi32(buf0[19], buf0[23]);
    buf1[24] = _mm_add_epi32(buf0[24], buf0[28]);
    buf1[28] = _mm_sub_epi32(buf0[24], buf0[28]);
    buf1[25] = _mm_add_epi32(buf0[25], buf0[29]);
    buf1[29] = _mm_sub_epi32(buf0[25], buf0[29]);
    buf1[26] = _mm_add_epi32(buf0[26], buf0[30]);
    buf1[30] = _mm_sub_epi32(buf0[26], buf0[30]);
    buf1[27] = _mm_add_epi32(buf0[27], buf0[31]);
    buf1[31] = _mm_sub_epi32(buf0[27], buf0[31]);

    // stage 8
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[4], buf1[5], buf0[4],
                        buf0[5], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    buf0[10] = buf1[10];
    buf0[11] = buf1[11];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[12], buf1[13], buf0[12],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    buf0[18] = buf1[18];
    buf0[19] = buf1[19];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[20], buf1[21], buf0[20],
                        buf0[21], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[22], buf1[23], buf0[22],
                        buf0[23], bit);
    buf0[24] = buf1[24];
    buf0[25] = buf1[25];
    buf0[26] = buf1[26];
    buf0[27] = buf1[27];
    btf_32_sse4_1_type0(cospi[16], cospi[48], buf1[28], buf1[29], buf0[28],
                        buf0[29], bit);
    btf_32_sse4_1_type0(-cospi[48], cospi[16], buf1[30], buf1[31], buf0[30],
                        buf0[31], bit);

    // stage 9
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[0], buf0[2]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[1], buf0[3]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[6]);
    buf1[6] = _mm_sub_epi32(buf0[4], buf0[6]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[7]);
    buf1[7] = _mm_sub_epi32(buf0[5], buf0[7]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[8], buf0[10]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[9], buf0[11]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[12], buf0[14]);
    buf1[13] = _mm_add_epi32(buf0[13], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[13], buf0[15]);
    buf1[16] = _mm_add_epi32(buf0[16], buf0[18]);
    buf1[18] = _mm_sub_epi32(buf0[16], buf0[18]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[19]);
    buf1[19] = _mm_sub_epi32(buf0[17], buf0[19]);
    buf1[20] = _mm_add_epi32(buf0[20], buf0[22]);
    buf1[22] = _mm_sub_epi32(buf0[20], buf0[22]);
    buf1[21] = _mm_add_epi32(buf0[21], buf0[23]);
    buf1[23] = _mm_sub_epi32(buf0[21], buf0[23]);
    buf1[24] = _mm_add_epi32(buf0[24], buf0[26]);
    buf1[26] = _mm_sub_epi32(buf0[24], buf0[26]);
    buf1[25] = _mm_add_epi32(buf0[25], buf0[27]);
    buf1[27] = _mm_sub_epi32(buf0[25], buf0[27]);
    buf1[28] = _mm_add_epi32(buf0[28], buf0[30]);
    buf1[30] = _mm_sub_epi32(buf0[28], buf0[30]);
    buf1[29] = _mm_add_epi32(buf0[29], buf0[31]);
    buf1[31] = _mm_sub_epi32(buf0[29], buf0[31]);

    // stage 10
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[6], buf1[7], buf0[6],
                        buf0[7], bit);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[10], buf1[11], buf0[10],
                        buf0[11], bit);
    buf0[12] = buf1[12];
    buf0[13] = buf1[13];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[14], buf1[15], buf0[14],
                        buf0[15], bit);
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[18], buf1[19], buf0[18],
                        buf0[19], bit);
    buf0[20] = buf1[20];
    buf0[21] = buf1[21];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[22], buf1[23], buf0[22],
                        buf0[23], bit);
    buf0[24] = buf1[24];
    buf0[25] = buf1[25];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[26], buf1[27], buf0[26],
                        buf0[27], bit);
    buf0[28] = buf1[28];
    buf0[29] = buf1[29];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[30], buf1[31], buf0[30],
                        buf0[31], bit);

    // stage 11
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[16]);
    buf1[2] = buf0[24];
    buf1[3] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[8]);
    buf1[4] = buf0[12];
    buf1[5] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[28]);
    buf1[6] = buf0[20];
    buf1[7] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[4]);
    buf1[8] = buf0[6];
    buf1[9] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[22]);
    buf1[10] = buf0[30];
    buf1[11] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[14]);
    buf1[12] = buf0[10];
    buf1[13] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[26]);
    buf1[14] = buf0[18];
    buf1[15] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[2]);
    buf1[16] = buf0[3];
    buf1[17] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[19]);
    buf1[18] = buf0[27];
    buf1[19] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[11]);
    buf1[20] = buf0[15];
    buf1[21] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[31]);
    buf1[22] = buf0[23];
    buf1[23] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[7]);
    buf1[24] = buf0[5];
    buf1[25] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[21]);
    buf1[26] = buf0[29];
    buf1[27] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[13]);
    buf1[28] = buf0[9];
    buf1[29] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[25]);
    buf1[30] = buf0[17];
    buf1[31] = _mm_sub_epi32(_mm_set1_epi32(0), buf0[1]);

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
    output[8 * col_num + col] = buf1[8];
    output[9 * col_num + col] = buf1[9];
    output[10 * col_num + col] = buf1[10];
    output[11 * col_num + col] = buf1[11];
    output[12 * col_num + col] = buf1[12];
    output[13 * col_num + col] = buf1[13];
    output[14 * col_num + col] = buf1[14];
    output[15 * col_num + col] = buf1[15];
    output[16 * col_num + col] = buf1[16];
    output[17 * col_num + col] = buf1[17];
    output[18 * col_num + col] = buf1[18];
    output[19 * col_num + col] = buf1[19];
    output[20 * col_num + col] = buf1[20];
    output[21 * col_num + col] = buf1[21];
    output[22 * col_num + col] = buf1[22];
    output[23 * col_num + col] = buf1[23];
    output[24 * col_num + col] = buf1[24];
    output[25 * col_num + col] = buf1[25];
    output[26 * col_num + col] = buf1[26];
    output[27 * col_num + col] = buf1[27];
    output[28 * col_num + col] = buf1[28];
    output[29 * col_num + col] = buf1[29];
    output[30 * col_num + col] = buf1[30];
    output[31 * col_num + col] = buf1[31];
  }
}

void vp10_fdct64_new_sse4_1(const __m128i* input, __m128i* output,
                            const int8_t* cos_bit, const int8_t* stage_range) {
  const int txfm_size = 64;
  const int num_per_128 = 4;
  const int32_t* cospi;
  __m128i buf0[64];
  __m128i buf1[64];
  int col_num = txfm_size / num_per_128;
  int bit;
  int col;
  (void)stage_range;
  for (col = 0; col < col_num; col++) {
    // stage 0;
    int32_t stage_idx = 0;
    buf0[0] = input[0 * col_num + col];
    buf0[1] = input[1 * col_num + col];
    buf0[2] = input[2 * col_num + col];
    buf0[3] = input[3 * col_num + col];
    buf0[4] = input[4 * col_num + col];
    buf0[5] = input[5 * col_num + col];
    buf0[6] = input[6 * col_num + col];
    buf0[7] = input[7 * col_num + col];
    buf0[8] = input[8 * col_num + col];
    buf0[9] = input[9 * col_num + col];
    buf0[10] = input[10 * col_num + col];
    buf0[11] = input[11 * col_num + col];
    buf0[12] = input[12 * col_num + col];
    buf0[13] = input[13 * col_num + col];
    buf0[14] = input[14 * col_num + col];
    buf0[15] = input[15 * col_num + col];
    buf0[16] = input[16 * col_num + col];
    buf0[17] = input[17 * col_num + col];
    buf0[18] = input[18 * col_num + col];
    buf0[19] = input[19 * col_num + col];
    buf0[20] = input[20 * col_num + col];
    buf0[21] = input[21 * col_num + col];
    buf0[22] = input[22 * col_num + col];
    buf0[23] = input[23 * col_num + col];
    buf0[24] = input[24 * col_num + col];
    buf0[25] = input[25 * col_num + col];
    buf0[26] = input[26 * col_num + col];
    buf0[27] = input[27 * col_num + col];
    buf0[28] = input[28 * col_num + col];
    buf0[29] = input[29 * col_num + col];
    buf0[30] = input[30 * col_num + col];
    buf0[31] = input[31 * col_num + col];
    buf0[32] = input[32 * col_num + col];
    buf0[33] = input[33 * col_num + col];
    buf0[34] = input[34 * col_num + col];
    buf0[35] = input[35 * col_num + col];
    buf0[36] = input[36 * col_num + col];
    buf0[37] = input[37 * col_num + col];
    buf0[38] = input[38 * col_num + col];
    buf0[39] = input[39 * col_num + col];
    buf0[40] = input[40 * col_num + col];
    buf0[41] = input[41 * col_num + col];
    buf0[42] = input[42 * col_num + col];
    buf0[43] = input[43 * col_num + col];
    buf0[44] = input[44 * col_num + col];
    buf0[45] = input[45 * col_num + col];
    buf0[46] = input[46 * col_num + col];
    buf0[47] = input[47 * col_num + col];
    buf0[48] = input[48 * col_num + col];
    buf0[49] = input[49 * col_num + col];
    buf0[50] = input[50 * col_num + col];
    buf0[51] = input[51 * col_num + col];
    buf0[52] = input[52 * col_num + col];
    buf0[53] = input[53 * col_num + col];
    buf0[54] = input[54 * col_num + col];
    buf0[55] = input[55 * col_num + col];
    buf0[56] = input[56 * col_num + col];
    buf0[57] = input[57 * col_num + col];
    buf0[58] = input[58 * col_num + col];
    buf0[59] = input[59 * col_num + col];
    buf0[60] = input[60 * col_num + col];
    buf0[61] = input[61 * col_num + col];
    buf0[62] = input[62 * col_num + col];
    buf0[63] = input[63 * col_num + col];

    // stage 1
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[63]);
    buf1[63] = _mm_sub_epi32(buf0[0], buf0[63]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[62]);
    buf1[62] = _mm_sub_epi32(buf0[1], buf0[62]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[61]);
    buf1[61] = _mm_sub_epi32(buf0[2], buf0[61]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[60]);
    buf1[60] = _mm_sub_epi32(buf0[3], buf0[60]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[59]);
    buf1[59] = _mm_sub_epi32(buf0[4], buf0[59]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[58]);
    buf1[58] = _mm_sub_epi32(buf0[5], buf0[58]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[57]);
    buf1[57] = _mm_sub_epi32(buf0[6], buf0[57]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[56]);
    buf1[56] = _mm_sub_epi32(buf0[7], buf0[56]);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[55]);
    buf1[55] = _mm_sub_epi32(buf0[8], buf0[55]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[54]);
    buf1[54] = _mm_sub_epi32(buf0[9], buf0[54]);
    buf1[10] = _mm_add_epi32(buf0[10], buf0[53]);
    buf1[53] = _mm_sub_epi32(buf0[10], buf0[53]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[52]);
    buf1[52] = _mm_sub_epi32(buf0[11], buf0[52]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[51]);
    buf1[51] = _mm_sub_epi32(buf0[12], buf0[51]);
    buf1[13] = _mm_add_epi32(buf0[13], buf0[50]);
    buf1[50] = _mm_sub_epi32(buf0[13], buf0[50]);
    buf1[14] = _mm_add_epi32(buf0[14], buf0[49]);
    buf1[49] = _mm_sub_epi32(buf0[14], buf0[49]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[48]);
    buf1[48] = _mm_sub_epi32(buf0[15], buf0[48]);
    buf1[16] = _mm_add_epi32(buf0[16], buf0[47]);
    buf1[47] = _mm_sub_epi32(buf0[16], buf0[47]);
    buf1[17] = _mm_add_epi32(buf0[17], buf0[46]);
    buf1[46] = _mm_sub_epi32(buf0[17], buf0[46]);
    buf1[18] = _mm_add_epi32(buf0[18], buf0[45]);
    buf1[45] = _mm_sub_epi32(buf0[18], buf0[45]);
    buf1[19] = _mm_add_epi32(buf0[19], buf0[44]);
    buf1[44] = _mm_sub_epi32(buf0[19], buf0[44]);
    buf1[20] = _mm_add_epi32(buf0[20], buf0[43]);
    buf1[43] = _mm_sub_epi32(buf0[20], buf0[43]);
    buf1[21] = _mm_add_epi32(buf0[21], buf0[42]);
    buf1[42] = _mm_sub_epi32(buf0[21], buf0[42]);
    buf1[22] = _mm_add_epi32(buf0[22], buf0[41]);
    buf1[41] = _mm_sub_epi32(buf0[22], buf0[41]);
    buf1[23] = _mm_add_epi32(buf0[23], buf0[40]);
    buf1[40] = _mm_sub_epi32(buf0[23], buf0[40]);
    buf1[24] = _mm_add_epi32(buf0[24], buf0[39]);
    buf1[39] = _mm_sub_epi32(buf0[24], buf0[39]);
    buf1[25] = _mm_add_epi32(buf0[25], buf0[38]);
    buf1[38] = _mm_sub_epi32(buf0[25], buf0[38]);
    buf1[26] = _mm_add_epi32(buf0[26], buf0[37]);
    buf1[37] = _mm_sub_epi32(buf0[26], buf0[37]);
    buf1[27] = _mm_add_epi32(buf0[27], buf0[36]);
    buf1[36] = _mm_sub_epi32(buf0[27], buf0[36]);
    buf1[28] = _mm_add_epi32(buf0[28], buf0[35]);
    buf1[35] = _mm_sub_epi32(buf0[28], buf0[35]);
    buf1[29] = _mm_add_epi32(buf0[29], buf0[34]);
    buf1[34] = _mm_sub_epi32(buf0[29], buf0[34]);
    buf1[30] = _mm_add_epi32(buf0[30], buf0[33]);
    buf1[33] = _mm_sub_epi32(buf0[30], buf0[33]);
    buf1[31] = _mm_add_epi32(buf0[31], buf0[32]);
    buf1[32] = _mm_sub_epi32(buf0[31], buf0[32]);

    // stage 2
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[31]);
    buf0[31] = _mm_sub_epi32(buf1[0], buf1[31]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[30]);
    buf0[30] = _mm_sub_epi32(buf1[1], buf1[30]);
    buf0[2] = _mm_add_epi32(buf1[2], buf1[29]);
    buf0[29] = _mm_sub_epi32(buf1[2], buf1[29]);
    buf0[3] = _mm_add_epi32(buf1[3], buf1[28]);
    buf0[28] = _mm_sub_epi32(buf1[3], buf1[28]);
    buf0[4] = _mm_add_epi32(buf1[4], buf1[27]);
    buf0[27] = _mm_sub_epi32(buf1[4], buf1[27]);
    buf0[5] = _mm_add_epi32(buf1[5], buf1[26]);
    buf0[26] = _mm_sub_epi32(buf1[5], buf1[26]);
    buf0[6] = _mm_add_epi32(buf1[6], buf1[25]);
    buf0[25] = _mm_sub_epi32(buf1[6], buf1[25]);
    buf0[7] = _mm_add_epi32(buf1[7], buf1[24]);
    buf0[24] = _mm_sub_epi32(buf1[7], buf1[24]);
    buf0[8] = _mm_add_epi32(buf1[8], buf1[23]);
    buf0[23] = _mm_sub_epi32(buf1[8], buf1[23]);
    buf0[9] = _mm_add_epi32(buf1[9], buf1[22]);
    buf0[22] = _mm_sub_epi32(buf1[9], buf1[22]);
    buf0[10] = _mm_add_epi32(buf1[10], buf1[21]);
    buf0[21] = _mm_sub_epi32(buf1[10], buf1[21]);
    buf0[11] = _mm_add_epi32(buf1[11], buf1[20]);
    buf0[20] = _mm_sub_epi32(buf1[11], buf1[20]);
    buf0[12] = _mm_add_epi32(buf1[12], buf1[19]);
    buf0[19] = _mm_sub_epi32(buf1[12], buf1[19]);
    buf0[13] = _mm_add_epi32(buf1[13], buf1[18]);
    buf0[18] = _mm_sub_epi32(buf1[13], buf1[18]);
    buf0[14] = _mm_add_epi32(buf1[14], buf1[17]);
    buf0[17] = _mm_sub_epi32(buf1[14], buf1[17]);
    buf0[15] = _mm_add_epi32(buf1[15], buf1[16]);
    buf0[16] = _mm_sub_epi32(buf1[15], buf1[16]);
    buf0[32] = buf1[32];
    buf0[33] = buf1[33];
    buf0[34] = buf1[34];
    buf0[35] = buf1[35];
    buf0[36] = buf1[36];
    buf0[37] = buf1[37];
    buf0[38] = buf1[38];
    buf0[39] = buf1[39];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[40], buf1[55], buf0[40],
                        buf0[55], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[41], buf1[54], buf0[41],
                        buf0[54], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[42], buf1[53], buf0[42],
                        buf0[53], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[43], buf1[52], buf0[43],
                        buf0[52], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[44], buf1[51], buf0[44],
                        buf0[51], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[45], buf1[50], buf0[45],
                        buf0[50], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[46], buf1[49], buf0[46],
                        buf0[49], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[47], buf1[48], buf0[47],
                        buf0[48], bit);
    buf0[56] = buf1[56];
    buf0[57] = buf1[57];
    buf0[58] = buf1[58];
    buf0[59] = buf1[59];
    buf0[60] = buf1[60];
    buf0[61] = buf1[61];
    buf0[62] = buf1[62];
    buf0[63] = buf1[63];

    // stage 3
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[15]);
    buf1[15] = _mm_sub_epi32(buf0[0], buf0[15]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[14]);
    buf1[14] = _mm_sub_epi32(buf0[1], buf0[14]);
    buf1[2] = _mm_add_epi32(buf0[2], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[2], buf0[13]);
    buf1[3] = _mm_add_epi32(buf0[3], buf0[12]);
    buf1[12] = _mm_sub_epi32(buf0[3], buf0[12]);
    buf1[4] = _mm_add_epi32(buf0[4], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[4], buf0[11]);
    buf1[5] = _mm_add_epi32(buf0[5], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[5], buf0[10]);
    buf1[6] = _mm_add_epi32(buf0[6], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[6], buf0[9]);
    buf1[7] = _mm_add_epi32(buf0[7], buf0[8]);
    buf1[8] = _mm_sub_epi32(buf0[7], buf0[8]);
    buf1[16] = buf0[16];
    buf1[17] = buf0[17];
    buf1[18] = buf0[18];
    buf1[19] = buf0[19];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[20], buf0[27], buf1[20],
                        buf1[27], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[21], buf0[26], buf1[21],
                        buf1[26], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[22], buf0[25], buf1[22],
                        buf1[25], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[23], buf0[24], buf1[23],
                        buf1[24], bit);
    buf1[28] = buf0[28];
    buf1[29] = buf0[29];
    buf1[30] = buf0[30];
    buf1[31] = buf0[31];
    buf1[32] = _mm_add_epi32(buf0[32], buf0[47]);
    buf1[47] = _mm_sub_epi32(buf0[32], buf0[47]);
    buf1[33] = _mm_add_epi32(buf0[33], buf0[46]);
    buf1[46] = _mm_sub_epi32(buf0[33], buf0[46]);
    buf1[34] = _mm_add_epi32(buf0[34], buf0[45]);
    buf1[45] = _mm_sub_epi32(buf0[34], buf0[45]);
    buf1[35] = _mm_add_epi32(buf0[35], buf0[44]);
    buf1[44] = _mm_sub_epi32(buf0[35], buf0[44]);
    buf1[36] = _mm_add_epi32(buf0[36], buf0[43]);
    buf1[43] = _mm_sub_epi32(buf0[36], buf0[43]);
    buf1[37] = _mm_add_epi32(buf0[37], buf0[42]);
    buf1[42] = _mm_sub_epi32(buf0[37], buf0[42]);
    buf1[38] = _mm_add_epi32(buf0[38], buf0[41]);
    buf1[41] = _mm_sub_epi32(buf0[38], buf0[41]);
    buf1[39] = _mm_add_epi32(buf0[39], buf0[40]);
    buf1[40] = _mm_sub_epi32(buf0[39], buf0[40]);
    buf1[48] = _mm_sub_epi32(buf0[63], buf0[48]);
    buf1[63] = _mm_add_epi32(buf0[63], buf0[48]);
    buf1[49] = _mm_sub_epi32(buf0[62], buf0[49]);
    buf1[62] = _mm_add_epi32(buf0[62], buf0[49]);
    buf1[50] = _mm_sub_epi32(buf0[61], buf0[50]);
    buf1[61] = _mm_add_epi32(buf0[61], buf0[50]);
    buf1[51] = _mm_sub_epi32(buf0[60], buf0[51]);
    buf1[60] = _mm_add_epi32(buf0[60], buf0[51]);
    buf1[52] = _mm_sub_epi32(buf0[59], buf0[52]);
    buf1[59] = _mm_add_epi32(buf0[59], buf0[52]);
    buf1[53] = _mm_sub_epi32(buf0[58], buf0[53]);
    buf1[58] = _mm_add_epi32(buf0[58], buf0[53]);
    buf1[54] = _mm_sub_epi32(buf0[57], buf0[54]);
    buf1[57] = _mm_add_epi32(buf0[57], buf0[54]);
    buf1[55] = _mm_sub_epi32(buf0[56], buf0[55]);
    buf1[56] = _mm_add_epi32(buf0[56], buf0[55]);

    // stage 4
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = _mm_add_epi32(buf1[0], buf1[7]);
    buf0[7] = _mm_sub_epi32(buf1[0], buf1[7]);
    buf0[1] = _mm_add_epi32(buf1[1], buf1[6]);
    buf0[6] = _mm_sub_epi32(buf1[1], buf1[6]);
    buf0[2] = _mm_add_epi32(buf1[2], buf1[5]);
    buf0[5] = _mm_sub_epi32(buf1[2], buf1[5]);
    buf0[3] = _mm_add_epi32(buf1[3], buf1[4]);
    buf0[4] = _mm_sub_epi32(buf1[3], buf1[4]);
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf1[11], buf1[12], buf0[11],
                        buf0[12], bit);
    buf0[14] = buf1[14];
    buf0[15] = buf1[15];
    buf0[16] = _mm_add_epi32(buf1[16], buf1[23]);
    buf0[23] = _mm_sub_epi32(buf1[16], buf1[23]);
    buf0[17] = _mm_add_epi32(buf1[17], buf1[22]);
    buf0[22] = _mm_sub_epi32(buf1[17], buf1[22]);
    buf0[18] = _mm_add_epi32(buf1[18], buf1[21]);
    buf0[21] = _mm_sub_epi32(buf1[18], buf1[21]);
    buf0[19] = _mm_add_epi32(buf1[19], buf1[20]);
    buf0[20] = _mm_sub_epi32(buf1[19], buf1[20]);
    buf0[24] = _mm_sub_epi32(buf1[31], buf1[24]);
    buf0[31] = _mm_add_epi32(buf1[31], buf1[24]);
    buf0[25] = _mm_sub_epi32(buf1[30], buf1[25]);
    buf0[30] = _mm_add_epi32(buf1[30], buf1[25]);
    buf0[26] = _mm_sub_epi32(buf1[29], buf1[26]);
    buf0[29] = _mm_add_epi32(buf1[29], buf1[26]);
    buf0[27] = _mm_sub_epi32(buf1[28], buf1[27]);
    buf0[28] = _mm_add_epi32(buf1[28], buf1[27]);
    buf0[32] = buf1[32];
    buf0[33] = buf1[33];
    buf0[34] = buf1[34];
    buf0[35] = buf1[35];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[36], buf1[59], buf0[36],
                        buf0[59], bit);
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[37], buf1[58], buf0[37],
                        buf0[58], bit);
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[38], buf1[57], buf0[38],
                        buf0[57], bit);
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[39], buf1[56], buf0[39],
                        buf0[56], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[40], buf1[55], buf0[40],
                        buf0[55], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[41], buf1[54], buf0[41],
                        buf0[54], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[42], buf1[53], buf0[42],
                        buf0[53], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[43], buf1[52], buf0[43],
                        buf0[52], bit);
    buf0[44] = buf1[44];
    buf0[45] = buf1[45];
    buf0[46] = buf1[46];
    buf0[47] = buf1[47];
    buf0[48] = buf1[48];
    buf0[49] = buf1[49];
    buf0[50] = buf1[50];
    buf0[51] = buf1[51];
    buf0[60] = buf1[60];
    buf0[61] = buf1[61];
    buf0[62] = buf1[62];
    buf0[63] = buf1[63];

    // stage 5
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = _mm_add_epi32(buf0[0], buf0[3]);
    buf1[3] = _mm_sub_epi32(buf0[0], buf0[3]);
    buf1[1] = _mm_add_epi32(buf0[1], buf0[2]);
    buf1[2] = _mm_sub_epi32(buf0[1], buf0[2]);
    buf1[4] = buf0[4];
    btf_32_sse4_1_type0(-cospi[32], cospi[32], buf0[5], buf0[6], buf1[5],
                        buf1[6], bit);
    buf1[7] = buf0[7];
    buf1[8] = _mm_add_epi32(buf0[8], buf0[11]);
    buf1[11] = _mm_sub_epi32(buf0[8], buf0[11]);
    buf1[9] = _mm_add_epi32(buf0[9], buf0[10]);
    buf1[10] = _mm_sub_epi32(buf0[9], buf0[10]);
    buf1[12] = _mm_sub_epi32(buf0[15], buf0[12]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[12]);
    buf1[13] = _mm_sub_epi32(buf0[14], buf0[13]);
    buf1[14] = _mm_add_epi32(buf0[14], buf0[13]);
    buf1[16] = buf0[16];
    buf1[17] = buf0[17];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf0[18], buf0[29], buf1[18],
                        buf1[29], bit);
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf0[19], buf0[28], buf1[19],
                        buf1[28], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf0[20], buf0[27], buf1[20],
                        buf1[27], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf0[21], buf0[26], buf1[21],
                        buf1[26], bit);
    buf1[22] = buf0[22];
    buf1[23] = buf0[23];
    buf1[24] = buf0[24];
    buf1[25] = buf0[25];
    buf1[30] = buf0[30];
    buf1[31] = buf0[31];
    buf1[32] = _mm_add_epi32(buf0[32], buf0[39]);
    buf1[39] = _mm_sub_epi32(buf0[32], buf0[39]);
    buf1[33] = _mm_add_epi32(buf0[33], buf0[38]);
    buf1[38] = _mm_sub_epi32(buf0[33], buf0[38]);
    buf1[34] = _mm_add_epi32(buf0[34], buf0[37]);
    buf1[37] = _mm_sub_epi32(buf0[34], buf0[37]);
    buf1[35] = _mm_add_epi32(buf0[35], buf0[36]);
    buf1[36] = _mm_sub_epi32(buf0[35], buf0[36]);
    buf1[40] = _mm_sub_epi32(buf0[47], buf0[40]);
    buf1[47] = _mm_add_epi32(buf0[47], buf0[40]);
    buf1[41] = _mm_sub_epi32(buf0[46], buf0[41]);
    buf1[46] = _mm_add_epi32(buf0[46], buf0[41]);
    buf1[42] = _mm_sub_epi32(buf0[45], buf0[42]);
    buf1[45] = _mm_add_epi32(buf0[45], buf0[42]);
    buf1[43] = _mm_sub_epi32(buf0[44], buf0[43]);
    buf1[44] = _mm_add_epi32(buf0[44], buf0[43]);
    buf1[48] = _mm_add_epi32(buf0[48], buf0[55]);
    buf1[55] = _mm_sub_epi32(buf0[48], buf0[55]);
    buf1[49] = _mm_add_epi32(buf0[49], buf0[54]);
    buf1[54] = _mm_sub_epi32(buf0[49], buf0[54]);
    buf1[50] = _mm_add_epi32(buf0[50], buf0[53]);
    buf1[53] = _mm_sub_epi32(buf0[50], buf0[53]);
    buf1[51] = _mm_add_epi32(buf0[51], buf0[52]);
    buf1[52] = _mm_sub_epi32(buf0[51], buf0[52]);
    buf1[56] = _mm_sub_epi32(buf0[63], buf0[56]);
    buf1[63] = _mm_add_epi32(buf0[63], buf0[56]);
    buf1[57] = _mm_sub_epi32(buf0[62], buf0[57]);
    buf1[62] = _mm_add_epi32(buf0[62], buf0[57]);
    buf1[58] = _mm_sub_epi32(buf0[61], buf0[58]);
    buf1[61] = _mm_add_epi32(buf0[61], buf0[58]);
    buf1[59] = _mm_sub_epi32(buf0[60], buf0[59]);
    buf1[60] = _mm_add_epi32(buf0[60], buf0[59]);

    // stage 6
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    btf_32_sse4_1_type0(cospi[32], cospi[32], buf1[0], buf1[1], buf0[0],
                        buf0[1], bit);
    btf_32_sse4_1_type1(cospi[48], cospi[16], buf1[2], buf1[3], buf0[2],
                        buf0[3], bit);
    buf0[4] = _mm_add_epi32(buf1[4], buf1[5]);
    buf0[5] = _mm_sub_epi32(buf1[4], buf1[5]);
    buf0[6] = _mm_sub_epi32(buf1[7], buf1[6]);
    buf0[7] = _mm_add_epi32(buf1[7], buf1[6]);
    buf0[8] = buf1[8];
    btf_32_sse4_1_type0(-cospi[16], cospi[48], buf1[9], buf1[14], buf0[9],
                        buf0[14], bit);
    btf_32_sse4_1_type0(-cospi[48], -cospi[16], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    buf0[11] = buf1[11];
    buf0[12] = buf1[12];
    buf0[15] = buf1[15];
    buf0[16] = _mm_add_epi32(buf1[16], buf1[19]);
    buf0[19] = _mm_sub_epi32(buf1[16], buf1[19]);
    buf0[17] = _mm_add_epi32(buf1[17], buf1[18]);
    buf0[18] = _mm_sub_epi32(buf1[17], buf1[18]);
    buf0[20] = _mm_sub_epi32(buf1[23], buf1[20]);
    buf0[23] = _mm_add_epi32(buf1[23], buf1[20]);
    buf0[21] = _mm_sub_epi32(buf1[22], buf1[21]);
    buf0[22] = _mm_add_epi32(buf1[22], buf1[21]);
    buf0[24] = _mm_add_epi32(buf1[24], buf1[27]);
    buf0[27] = _mm_sub_epi32(buf1[24], buf1[27]);
    buf0[25] = _mm_add_epi32(buf1[25], buf1[26]);
    buf0[26] = _mm_sub_epi32(buf1[25], buf1[26]);
    buf0[28] = _mm_sub_epi32(buf1[31], buf1[28]);
    buf0[31] = _mm_add_epi32(buf1[31], buf1[28]);
    buf0[29] = _mm_sub_epi32(buf1[30], buf1[29]);
    buf0[30] = _mm_add_epi32(buf1[30], buf1[29]);
    buf0[32] = buf1[32];
    buf0[33] = buf1[33];
    btf_32_sse4_1_type0(-cospi[8], cospi[56], buf1[34], buf1[61], buf0[34],
                        buf0[61], bit);
    btf_32_sse4_1_type0(-cospi[8], cospi[56], buf1[35], buf1[60], buf0[35],
                        buf0[60], bit);
    btf_32_sse4_1_type0(-cospi[56], -cospi[8], buf1[36], buf1[59], buf0[36],
                        buf0[59], bit);
    btf_32_sse4_1_type0(-cospi[56], -cospi[8], buf1[37], buf1[58], buf0[37],
                        buf0[58], bit);
    buf0[38] = buf1[38];
    buf0[39] = buf1[39];
    buf0[40] = buf1[40];
    buf0[41] = buf1[41];
    btf_32_sse4_1_type0(-cospi[40], cospi[24], buf1[42], buf1[53], buf0[42],
                        buf0[53], bit);
    btf_32_sse4_1_type0(-cospi[40], cospi[24], buf1[43], buf1[52], buf0[43],
                        buf0[52], bit);
    btf_32_sse4_1_type0(-cospi[24], -cospi[40], buf1[44], buf1[51], buf0[44],
                        buf0[51], bit);
    btf_32_sse4_1_type0(-cospi[24], -cospi[40], buf1[45], buf1[50], buf0[45],
                        buf0[50], bit);
    buf0[46] = buf1[46];
    buf0[47] = buf1[47];
    buf0[48] = buf1[48];
    buf0[49] = buf1[49];
    buf0[54] = buf1[54];
    buf0[55] = buf1[55];
    buf0[56] = buf1[56];
    buf0[57] = buf1[57];
    buf0[62] = buf1[62];
    buf0[63] = buf1[63];

    // stage 7
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[1];
    buf1[2] = buf0[2];
    buf1[3] = buf0[3];
    btf_32_sse4_1_type1(cospi[56], cospi[8], buf0[4], buf0[7], buf1[4], buf1[7],
                        bit);
    btf_32_sse4_1_type1(cospi[24], cospi[40], buf0[5], buf0[6], buf1[5],
                        buf1[6], bit);
    buf1[8] = _mm_add_epi32(buf0[8], buf0[9]);
    buf1[9] = _mm_sub_epi32(buf0[8], buf0[9]);
    buf1[10] = _mm_sub_epi32(buf0[11], buf0[10]);
    buf1[11] = _mm_add_epi32(buf0[11], buf0[10]);
    buf1[12] = _mm_add_epi32(buf0[12], buf0[13]);
    buf1[13] = _mm_sub_epi32(buf0[12], buf0[13]);
    buf1[14] = _mm_sub_epi32(buf0[15], buf0[14]);
    buf1[15] = _mm_add_epi32(buf0[15], buf0[14]);
    buf1[16] = buf0[16];
    btf_32_sse4_1_type0(-cospi[8], cospi[56], buf0[17], buf0[30], buf1[17],
                        buf1[30], bit);
    btf_32_sse4_1_type0(-cospi[56], -cospi[8], buf0[18], buf0[29], buf1[18],
                        buf1[29], bit);
    buf1[19] = buf0[19];
    buf1[20] = buf0[20];
    btf_32_sse4_1_type0(-cospi[40], cospi[24], buf0[21], buf0[26], buf1[21],
                        buf1[26], bit);
    btf_32_sse4_1_type0(-cospi[24], -cospi[40], buf0[22], buf0[25], buf1[22],
                        buf1[25], bit);
    buf1[23] = buf0[23];
    buf1[24] = buf0[24];
    buf1[27] = buf0[27];
    buf1[28] = buf0[28];
    buf1[31] = buf0[31];
    buf1[32] = _mm_add_epi32(buf0[32], buf0[35]);
    buf1[35] = _mm_sub_epi32(buf0[32], buf0[35]);
    buf1[33] = _mm_add_epi32(buf0[33], buf0[34]);
    buf1[34] = _mm_sub_epi32(buf0[33], buf0[34]);
    buf1[36] = _mm_sub_epi32(buf0[39], buf0[36]);
    buf1[39] = _mm_add_epi32(buf0[39], buf0[36]);
    buf1[37] = _mm_sub_epi32(buf0[38], buf0[37]);
    buf1[38] = _mm_add_epi32(buf0[38], buf0[37]);
    buf1[40] = _mm_add_epi32(buf0[40], buf0[43]);
    buf1[43] = _mm_sub_epi32(buf0[40], buf0[43]);
    buf1[41] = _mm_add_epi32(buf0[41], buf0[42]);
    buf1[42] = _mm_sub_epi32(buf0[41], buf0[42]);
    buf1[44] = _mm_sub_epi32(buf0[47], buf0[44]);
    buf1[47] = _mm_add_epi32(buf0[47], buf0[44]);
    buf1[45] = _mm_sub_epi32(buf0[46], buf0[45]);
    buf1[46] = _mm_add_epi32(buf0[46], buf0[45]);
    buf1[48] = _mm_add_epi32(buf0[48], buf0[51]);
    buf1[51] = _mm_sub_epi32(buf0[48], buf0[51]);
    buf1[49] = _mm_add_epi32(buf0[49], buf0[50]);
    buf1[50] = _mm_sub_epi32(buf0[49], buf0[50]);
    buf1[52] = _mm_sub_epi32(buf0[55], buf0[52]);
    buf1[55] = _mm_add_epi32(buf0[55], buf0[52]);
    buf1[53] = _mm_sub_epi32(buf0[54], buf0[53]);
    buf1[54] = _mm_add_epi32(buf0[54], buf0[53]);
    buf1[56] = _mm_add_epi32(buf0[56], buf0[59]);
    buf1[59] = _mm_sub_epi32(buf0[56], buf0[59]);
    buf1[57] = _mm_add_epi32(buf0[57], buf0[58]);
    buf1[58] = _mm_sub_epi32(buf0[57], buf0[58]);
    buf1[60] = _mm_sub_epi32(buf0[63], buf0[60]);
    buf1[63] = _mm_add_epi32(buf0[63], buf0[60]);
    buf1[61] = _mm_sub_epi32(buf0[62], buf0[61]);
    buf1[62] = _mm_add_epi32(buf0[62], buf0[61]);

    // stage 8
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    btf_32_sse4_1_type1(cospi[60], cospi[4], buf1[8], buf1[15], buf0[8],
                        buf0[15], bit);
    btf_32_sse4_1_type1(cospi[28], cospi[36], buf1[9], buf1[14], buf0[9],
                        buf0[14], bit);
    btf_32_sse4_1_type1(cospi[44], cospi[20], buf1[10], buf1[13], buf0[10],
                        buf0[13], bit);
    btf_32_sse4_1_type1(cospi[12], cospi[52], buf1[11], buf1[12], buf0[11],
                        buf0[12], bit);
    buf0[16] = _mm_add_epi32(buf1[16], buf1[17]);
    buf0[17] = _mm_sub_epi32(buf1[16], buf1[17]);
    buf0[18] = _mm_sub_epi32(buf1[19], buf1[18]);
    buf0[19] = _mm_add_epi32(buf1[19], buf1[18]);
    buf0[20] = _mm_add_epi32(buf1[20], buf1[21]);
    buf0[21] = _mm_sub_epi32(buf1[20], buf1[21]);
    buf0[22] = _mm_sub_epi32(buf1[23], buf1[22]);
    buf0[23] = _mm_add_epi32(buf1[23], buf1[22]);
    buf0[24] = _mm_add_epi32(buf1[24], buf1[25]);
    buf0[25] = _mm_sub_epi32(buf1[24], buf1[25]);
    buf0[26] = _mm_sub_epi32(buf1[27], buf1[26]);
    buf0[27] = _mm_add_epi32(buf1[27], buf1[26]);
    buf0[28] = _mm_add_epi32(buf1[28], buf1[29]);
    buf0[29] = _mm_sub_epi32(buf1[28], buf1[29]);
    buf0[30] = _mm_sub_epi32(buf1[31], buf1[30]);
    buf0[31] = _mm_add_epi32(buf1[31], buf1[30]);
    buf0[32] = buf1[32];
    btf_32_sse4_1_type0(-cospi[4], cospi[60], buf1[33], buf1[62], buf0[33],
                        buf0[62], bit);
    btf_32_sse4_1_type0(-cospi[60], -cospi[4], buf1[34], buf1[61], buf0[34],
                        buf0[61], bit);
    buf0[35] = buf1[35];
    buf0[36] = buf1[36];
    btf_32_sse4_1_type0(-cospi[36], cospi[28], buf1[37], buf1[58], buf0[37],
                        buf0[58], bit);
    btf_32_sse4_1_type0(-cospi[28], -cospi[36], buf1[38], buf1[57], buf0[38],
                        buf0[57], bit);
    buf0[39] = buf1[39];
    buf0[40] = buf1[40];
    btf_32_sse4_1_type0(-cospi[20], cospi[44], buf1[41], buf1[54], buf0[41],
                        buf0[54], bit);
    btf_32_sse4_1_type0(-cospi[44], -cospi[20], buf1[42], buf1[53], buf0[42],
                        buf0[53], bit);
    buf0[43] = buf1[43];
    buf0[44] = buf1[44];
    btf_32_sse4_1_type0(-cospi[52], cospi[12], buf1[45], buf1[50], buf0[45],
                        buf0[50], bit);
    btf_32_sse4_1_type0(-cospi[12], -cospi[52], buf1[46], buf1[49], buf0[46],
                        buf0[49], bit);
    buf0[47] = buf1[47];
    buf0[48] = buf1[48];
    buf0[51] = buf1[51];
    buf0[52] = buf1[52];
    buf0[55] = buf1[55];
    buf0[56] = buf1[56];
    buf0[59] = buf1[59];
    buf0[60] = buf1[60];
    buf0[63] = buf1[63];

    // stage 9
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[1];
    buf1[2] = buf0[2];
    buf1[3] = buf0[3];
    buf1[4] = buf0[4];
    buf1[5] = buf0[5];
    buf1[6] = buf0[6];
    buf1[7] = buf0[7];
    buf1[8] = buf0[8];
    buf1[9] = buf0[9];
    buf1[10] = buf0[10];
    buf1[11] = buf0[11];
    buf1[12] = buf0[12];
    buf1[13] = buf0[13];
    buf1[14] = buf0[14];
    buf1[15] = buf0[15];
    btf_32_sse4_1_type1(cospi[62], cospi[2], buf0[16], buf0[31], buf1[16],
                        buf1[31], bit);
    btf_32_sse4_1_type1(cospi[30], cospi[34], buf0[17], buf0[30], buf1[17],
                        buf1[30], bit);
    btf_32_sse4_1_type1(cospi[46], cospi[18], buf0[18], buf0[29], buf1[18],
                        buf1[29], bit);
    btf_32_sse4_1_type1(cospi[14], cospi[50], buf0[19], buf0[28], buf1[19],
                        buf1[28], bit);
    btf_32_sse4_1_type1(cospi[54], cospi[10], buf0[20], buf0[27], buf1[20],
                        buf1[27], bit);
    btf_32_sse4_1_type1(cospi[22], cospi[42], buf0[21], buf0[26], buf1[21],
                        buf1[26], bit);
    btf_32_sse4_1_type1(cospi[38], cospi[26], buf0[22], buf0[25], buf1[22],
                        buf1[25], bit);
    btf_32_sse4_1_type1(cospi[6], cospi[58], buf0[23], buf0[24], buf1[23],
                        buf1[24], bit);
    buf1[32] = _mm_add_epi32(buf0[32], buf0[33]);
    buf1[33] = _mm_sub_epi32(buf0[32], buf0[33]);
    buf1[34] = _mm_sub_epi32(buf0[35], buf0[34]);
    buf1[35] = _mm_add_epi32(buf0[35], buf0[34]);
    buf1[36] = _mm_add_epi32(buf0[36], buf0[37]);
    buf1[37] = _mm_sub_epi32(buf0[36], buf0[37]);
    buf1[38] = _mm_sub_epi32(buf0[39], buf0[38]);
    buf1[39] = _mm_add_epi32(buf0[39], buf0[38]);
    buf1[40] = _mm_add_epi32(buf0[40], buf0[41]);
    buf1[41] = _mm_sub_epi32(buf0[40], buf0[41]);
    buf1[42] = _mm_sub_epi32(buf0[43], buf0[42]);
    buf1[43] = _mm_add_epi32(buf0[43], buf0[42]);
    buf1[44] = _mm_add_epi32(buf0[44], buf0[45]);
    buf1[45] = _mm_sub_epi32(buf0[44], buf0[45]);
    buf1[46] = _mm_sub_epi32(buf0[47], buf0[46]);
    buf1[47] = _mm_add_epi32(buf0[47], buf0[46]);
    buf1[48] = _mm_add_epi32(buf0[48], buf0[49]);
    buf1[49] = _mm_sub_epi32(buf0[48], buf0[49]);
    buf1[50] = _mm_sub_epi32(buf0[51], buf0[50]);
    buf1[51] = _mm_add_epi32(buf0[51], buf0[50]);
    buf1[52] = _mm_add_epi32(buf0[52], buf0[53]);
    buf1[53] = _mm_sub_epi32(buf0[52], buf0[53]);
    buf1[54] = _mm_sub_epi32(buf0[55], buf0[54]);
    buf1[55] = _mm_add_epi32(buf0[55], buf0[54]);
    buf1[56] = _mm_add_epi32(buf0[56], buf0[57]);
    buf1[57] = _mm_sub_epi32(buf0[56], buf0[57]);
    buf1[58] = _mm_sub_epi32(buf0[59], buf0[58]);
    buf1[59] = _mm_add_epi32(buf0[59], buf0[58]);
    buf1[60] = _mm_add_epi32(buf0[60], buf0[61]);
    buf1[61] = _mm_sub_epi32(buf0[60], buf0[61]);
    buf1[62] = _mm_sub_epi32(buf0[63], buf0[62]);
    buf1[63] = _mm_add_epi32(buf0[63], buf0[62]);

    // stage 10
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf0[0] = buf1[0];
    buf0[1] = buf1[1];
    buf0[2] = buf1[2];
    buf0[3] = buf1[3];
    buf0[4] = buf1[4];
    buf0[5] = buf1[5];
    buf0[6] = buf1[6];
    buf0[7] = buf1[7];
    buf0[8] = buf1[8];
    buf0[9] = buf1[9];
    buf0[10] = buf1[10];
    buf0[11] = buf1[11];
    buf0[12] = buf1[12];
    buf0[13] = buf1[13];
    buf0[14] = buf1[14];
    buf0[15] = buf1[15];
    buf0[16] = buf1[16];
    buf0[17] = buf1[17];
    buf0[18] = buf1[18];
    buf0[19] = buf1[19];
    buf0[20] = buf1[20];
    buf0[21] = buf1[21];
    buf0[22] = buf1[22];
    buf0[23] = buf1[23];
    buf0[24] = buf1[24];
    buf0[25] = buf1[25];
    buf0[26] = buf1[26];
    buf0[27] = buf1[27];
    buf0[28] = buf1[28];
    buf0[29] = buf1[29];
    buf0[30] = buf1[30];
    buf0[31] = buf1[31];
    btf_32_sse4_1_type1(cospi[63], cospi[1], buf1[32], buf1[63], buf0[32],
                        buf0[63], bit);
    btf_32_sse4_1_type1(cospi[31], cospi[33], buf1[33], buf1[62], buf0[33],
                        buf0[62], bit);
    btf_32_sse4_1_type1(cospi[47], cospi[17], buf1[34], buf1[61], buf0[34],
                        buf0[61], bit);
    btf_32_sse4_1_type1(cospi[15], cospi[49], buf1[35], buf1[60], buf0[35],
                        buf0[60], bit);
    btf_32_sse4_1_type1(cospi[55], cospi[9], buf1[36], buf1[59], buf0[36],
                        buf0[59], bit);
    btf_32_sse4_1_type1(cospi[23], cospi[41], buf1[37], buf1[58], buf0[37],
                        buf0[58], bit);
    btf_32_sse4_1_type1(cospi[39], cospi[25], buf1[38], buf1[57], buf0[38],
                        buf0[57], bit);
    btf_32_sse4_1_type1(cospi[7], cospi[57], buf1[39], buf1[56], buf0[39],
                        buf0[56], bit);
    btf_32_sse4_1_type1(cospi[59], cospi[5], buf1[40], buf1[55], buf0[40],
                        buf0[55], bit);
    btf_32_sse4_1_type1(cospi[27], cospi[37], buf1[41], buf1[54], buf0[41],
                        buf0[54], bit);
    btf_32_sse4_1_type1(cospi[43], cospi[21], buf1[42], buf1[53], buf0[42],
                        buf0[53], bit);
    btf_32_sse4_1_type1(cospi[11], cospi[53], buf1[43], buf1[52], buf0[43],
                        buf0[52], bit);
    btf_32_sse4_1_type1(cospi[51], cospi[13], buf1[44], buf1[51], buf0[44],
                        buf0[51], bit);
    btf_32_sse4_1_type1(cospi[19], cospi[45], buf1[45], buf1[50], buf0[45],
                        buf0[50], bit);
    btf_32_sse4_1_type1(cospi[35], cospi[29], buf1[46], buf1[49], buf0[46],
                        buf0[49], bit);
    btf_32_sse4_1_type1(cospi[3], cospi[61], buf1[47], buf1[48], buf0[47],
                        buf0[48], bit);

    // stage 11
    stage_idx++;
    bit = cos_bit[stage_idx];
    cospi = cospi_arr[bit - cos_bit_min];
    buf1[0] = buf0[0];
    buf1[1] = buf0[32];
    buf1[2] = buf0[16];
    buf1[3] = buf0[48];
    buf1[4] = buf0[8];
    buf1[5] = buf0[40];
    buf1[6] = buf0[24];
    buf1[7] = buf0[56];
    buf1[8] = buf0[4];
    buf1[9] = buf0[36];
    buf1[10] = buf0[20];
    buf1[11] = buf0[52];
    buf1[12] = buf0[12];
    buf1[13] = buf0[44];
    buf1[14] = buf0[28];
    buf1[15] = buf0[60];
    buf1[16] = buf0[2];
    buf1[17] = buf0[34];
    buf1[18] = buf0[18];
    buf1[19] = buf0[50];
    buf1[20] = buf0[10];
    buf1[21] = buf0[42];
    buf1[22] = buf0[26];
    buf1[23] = buf0[58];
    buf1[24] = buf0[6];
    buf1[25] = buf0[38];
    buf1[26] = buf0[22];
    buf1[27] = buf0[54];
    buf1[28] = buf0[14];
    buf1[29] = buf0[46];
    buf1[30] = buf0[30];
    buf1[31] = buf0[62];
    buf1[32] = buf0[1];
    buf1[33] = buf0[33];
    buf1[34] = buf0[17];
    buf1[35] = buf0[49];
    buf1[36] = buf0[9];
    buf1[37] = buf0[41];
    buf1[38] = buf0[25];
    buf1[39] = buf0[57];
    buf1[40] = buf0[5];
    buf1[41] = buf0[37];
    buf1[42] = buf0[21];
    buf1[43] = buf0[53];
    buf1[44] = buf0[13];
    buf1[45] = buf0[45];
    buf1[46] = buf0[29];
    buf1[47] = buf0[61];
    buf1[48] = buf0[3];
    buf1[49] = buf0[35];
    buf1[50] = buf0[19];
    buf1[51] = buf0[51];
    buf1[52] = buf0[11];
    buf1[53] = buf0[43];
    buf1[54] = buf0[27];
    buf1[55] = buf0[59];
    buf1[56] = buf0[7];
    buf1[57] = buf0[39];
    buf1[58] = buf0[23];
    buf1[59] = buf0[55];
    buf1[60] = buf0[15];
    buf1[61] = buf0[47];
    buf1[62] = buf0[31];
    buf1[63] = buf0[63];

    output[0 * col_num + col] = buf1[0];
    output[1 * col_num + col] = buf1[1];
    output[2 * col_num + col] = buf1[2];
    output[3 * col_num + col] = buf1[3];
    output[4 * col_num + col] = buf1[4];
    output[5 * col_num + col] = buf1[5];
    output[6 * col_num + col] = buf1[6];
    output[7 * col_num + col] = buf1[7];
    output[8 * col_num + col] = buf1[8];
    output[9 * col_num + col] = buf1[9];
    output[10 * col_num + col] = buf1[10];
    output[11 * col_num + col] = buf1[11];
    output[12 * col_num + col] = buf1[12];
    output[13 * col_num + col] = buf1[13];
    output[14 * col_num + col] = buf1[14];
    output[15 * col_num + col] = buf1[15];
    output[16 * col_num + col] = buf1[16];
    output[17 * col_num + col] = buf1[17];
    output[18 * col_num + col] = buf1[18];
    output[19 * col_num + col] = buf1[19];
    output[20 * col_num + col] = buf1[20];
    output[21 * col_num + col] = buf1[21];
    output[22 * col_num + col] = buf1[22];
    output[23 * col_num + col] = buf1[23];
    output[24 * col_num + col] = buf1[24];
    output[25 * col_num + col] = buf1[25];
    output[26 * col_num + col] = buf1[26];
    output[27 * col_num + col] = buf1[27];
    output[28 * col_num + col] = buf1[28];
    output[29 * col_num + col] = buf1[29];
    output[30 * col_num + col] = buf1[30];
    output[31 * col_num + col] = buf1[31];
    output[32 * col_num + col] = buf1[32];
    output[33 * col_num + col] = buf1[33];
    output[34 * col_num + col] = buf1[34];
    output[35 * col_num + col] = buf1[35];
    output[36 * col_num + col] = buf1[36];
    output[37 * col_num + col] = buf1[37];
    output[38 * col_num + col] = buf1[38];
    output[39 * col_num + col] = buf1[39];
    output[40 * col_num + col] = buf1[40];
    output[41 * col_num + col] = buf1[41];
    output[42 * col_num + col] = buf1[42];
    output[43 * col_num + col] = buf1[43];
    output[44 * col_num + col] = buf1[44];
    output[45 * col_num + col] = buf1[45];
    output[46 * col_num + col] = buf1[46];
    output[47 * col_num + col] = buf1[47];
    output[48 * col_num + col] = buf1[48];
    output[49 * col_num + col] = buf1[49];
    output[50 * col_num + col] = buf1[50];
    output[51 * col_num + col] = buf1[51];
    output[52 * col_num + col] = buf1[52];
    output[53 * col_num + col] = buf1[53];
    output[54 * col_num + col] = buf1[54];
    output[55 * col_num + col] = buf1[55];
    output[56 * col_num + col] = buf1[56];
    output[57 * col_num + col] = buf1[57];
    output[58 * col_num + col] = buf1[58];
    output[59 * col_num + col] = buf1[59];
    output[60 * col_num + col] = buf1[60];
    output[61 * col_num + col] = buf1[61];
    output[62 * col_num + col] = buf1[62];
    output[63 * col_num + col] = buf1[63];
  }
}
