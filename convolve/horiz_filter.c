#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <tmmintrin.h>

static inline unsigned int readtsc(void) {
  unsigned int tsc;
  __asm__ __volatile__("rdtsc\n\t":"=a"(tsc):);
  return tsc;
}

#define FILTER_BITS         (7)
#define TEST_NUM            (32)
#define HD_WIDTH            (1920)
#define SUPER_BLOCK_HEIGHT  (128)
#define RAND_SEED           (0xabc)
unsigned int seed = RAND_SEED;

int round_power_of_two(int x, int n) {
  int ret = (x + (1 << (n - 1))) >> n;
  return ret;
}

uint8_t inline clip_pixel(int x) {
  uint8_t ret = x;
  if (x < 0) {
    ret = 0;
  }
  if (x > 255) {
    ret = 255;
  }
  return ret;
}

static int filtering(const uint8_t *src, const int16_t *filter, int flen) {
  int k;
  int sum = 0;
  int prod;
  for (k = 0; k < flen; ++k) {
    prod = src[k] * filter[k];
    sum += prod;
  }
  return sum;
}

void convolve(const uint8_t *src, int w, const int16_t *filter, int flen,
                uint8_t *buffer) {
  int i;
  int sum;

  for (i = 0; i < w; ++i) {
    sum = filtering(src, filter, flen);
    buffer[i] = clip_pixel(round_power_of_two(sum, FILTER_BITS));
    src += 1;
  }
}

void init_state(uint8_t *buf, uint8_t *pixel,
                int width, int height, int stride,
                const unsigned int random) {

  int row, col;
  int block = HD_WIDTH * SUPER_BLOCK_HEIGHT;

  memset(buf, 0, sizeof(buf[0]) * block);
  memset(pixel, 0, sizeof(pixel[0]) * block);

  seed = random;
  for (row = 0; row < height; ++row) {
    for (col = 0; col < width; ++col) {
      pixel[col] = clip_pixel(rand_r(&seed) % 255);
    }
    pixel += stride;
  }
}

void check_buffer(const uint8_t *buf1, const uint8_t *buf2,
                  int width, int height, int stride) {
  int row, col;

  for (row = 0; row < height; ++row) {
    for (col = 0; col < width; ++col) {
      if (buf1[col] != buf2[col]) {
        printf("Not bit-exact at col: %d row: %d\n", col, row);
        printf("Expected: 0x%x, Actual: 0x%x\n", buf1[col], buf2[col]);
        return;
      }
    }
    buf1 += stride;
    buf2 += stride;
  }
}

static const int16_t filter12[12] __attribute__ ((aligned(16))) = {
  -1,   3,  -4,   8, -18, 120,  28, -12,   7,  -4,   2, -1};

static const int16_t filter10[10] __attribute__ ((aligned(16))) = {
  1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1};

// SSSE3

struct Filter {
  const int8_t (*coeffs)[16];
  int tapsNum;
};

// for Horiz4 method (subpixel)
// 12-tap filter
static const int8_t filter12_subpixel_ns[6][16] __attribute__ ((aligned(16))) = {
  {-1, 3, -1, 3, -1, 3, -1, 3, -1, 3, -1, 3, -1, 3, -1, 3},
  {-4, 8, -4, 8, -4, 8, -4, 8, -4, 8, -4, 8, -4, 8, -4, 8},
  {-18, 120, -18, 120, -18, 120, -18, 120, -18, 120, -18, 120, -18, 120, -18, 120},
  {28, -12, 28, -12, 28, -12, 28, -12, 28, -12, 28, -12, 28, -12, 28, -12},
  {7, -4, 7, -4, 7, -4, 7, -4, 7, -4, 7, -4, 7, -4, 7, -4},
  {2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1},
};

// 10-tap filter
static const int8_t filter10_subpixel_ns[6][16] __attribute__ ((aligned(16))) = {
  { 0, 1,  0, 1,  0, 1,  0, 1,  0, 1,  0, 1,  0, 1,  0, 1},
  {-3, 7, -3, 7, -3, 7, -3, 7, -3, 7, -3, 7, -3, 7, -3, 7},
  {-17, 119, -17, 119, -17, 119, -17, 119, -17, 119, -17, 119, -17, 119, -17, 119},
  {28, -11, 28, -11, 28, -11, 28, -11, 28, -11, 28, -11, 28, -11, 28, -11},
  {5, -2, 5, -2, 5, -2, 5, -2, 5, -2, 5, -2, 5, -2, 5, -2},
  {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},
};

static const struct Filter pfilter_12tap_subpixel = {
  filter12_subpixel_ns, 12
};

static const struct Filter pfilter_10tap_subpixel = {
  filter10_subpixel_ns, 10
};

// for HorizP method
// 12-tap filter
const int8_t pfilter12[2][16] __attribute__ ((aligned(16))) = {
  {-1,  3, -4,  8, -18, 120,  28, -12,   7,  -4,   2,  -1,  0,  0,  0,  0},
  { 0,  0, -1,  3,  -4,   8, -18, 120,  28, -12,   7,  -4,  2, -1,  0,  0},
};

// 10-tap filter
const int8_t pfilter10[2][16] __attribute__ ((aligned(16))) = {
  {0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0, 0},
  {0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0},
};

const struct Filter pfilter_12tap = {
  pfilter12, 12
};

const struct Filter pfilter_10tap = {
  pfilter10, 10
};

void inline transpose_4x8(const __m128i *in, __m128i *out) {
  __m128i t0, t1;

  t0 = _mm_unpacklo_epi16(in[0], in[1]);
  t1 = _mm_unpacklo_epi16(in[2], in[3]);

  out[0] = _mm_unpacklo_epi32(t0, t1);
  out[1] = _mm_srli_si128(out[0], 8);
  out[2] = _mm_unpackhi_epi32(t0, t1);
  out[3] = _mm_srli_si128(out[2], 8);

  t0 = _mm_unpackhi_epi16(in[0], in[1]);
  t1 = _mm_unpackhi_epi16(in[2], in[3]);

  out[4] = _mm_unpacklo_epi32(t0, t1);
  out[5] = _mm_srli_si128(out[4], 8);
  // Note: We ignore out[6] and out[7] because
  // they're zero vectors.
}

void horiz_w4_ssse3(const uint8_t *src, const __m128i *f,
                    int tapsNum, uint8_t *buffer) {
  __m128i sumPairRow[4];
  __m128i sumPairCol[8];
  __m128i pixel;
  const __m128i k_256 = _mm_set1_epi16(1 << 8);

  if (10 == tapsNum) {
    src -= 1;
  }

  pixel = _mm_loadu_si128((__m128i const *)src);
  sumPairRow[0] = _mm_maddubs_epi16(pixel, f[0]);
  sumPairRow[2] = _mm_maddubs_epi16(pixel, f[1]);
  sumPairRow[2] = _mm_srli_si128(sumPairRow[2], 2);

  pixel = _mm_loadu_si128((__m128i const *)(src + 1));
  sumPairRow[1] = _mm_maddubs_epi16(pixel, f[0]);
  sumPairRow[3] = _mm_maddubs_epi16(pixel, f[1]);
  sumPairRow[3] = _mm_srli_si128(sumPairRow[3], 2);

  transpose_4x8(sumPairRow, sumPairCol);

  sumPairRow[0] = _mm_adds_epi16(sumPairCol[0], sumPairCol[1]);
  sumPairRow[1] = _mm_adds_epi16(sumPairCol[4], sumPairCol[5]);

  sumPairRow[2] = _mm_min_epi16(sumPairCol[2], sumPairCol[3]);
  sumPairRow[3] = _mm_max_epi16(sumPairCol[2], sumPairCol[3]);

  sumPairRow[0] = _mm_adds_epi16(sumPairRow[0], sumPairRow[1]);
  sumPairRow[0] = _mm_adds_epi16(sumPairRow[0], sumPairRow[2]);
  sumPairRow[0] = _mm_adds_epi16(sumPairRow[0], sumPairRow[3]);

  sumPairRow[1] = _mm_mulhrs_epi16(sumPairRow[0], k_256);
  sumPairRow[2] = _mm_packus_epi16(sumPairRow[1], sumPairRow[1]);

  *(int *)buffer = _mm_cvtsi128_si32(sumPairRow[2]);
}

void horiz_w8_ssse3(const uint8_t *src, const __m128i *f, int tapsNum,
                    uint8_t *buf) {
  horiz_w4_ssse3(src, f, tapsNum, buf);
  src += 4;
  buf += 4;
  horiz_w4_ssse3(src, f, tapsNum, buf);
}

void horiz_w16_ssse3(const uint8_t *src, const __m128i *f, int tapsNum,
                     uint8_t *buf) {
  horiz_w8_ssse3(src, f, tapsNum, buf);
  src += 8;
  buf += 8;
  horiz_w8_ssse3(src, f, tapsNum, buf);
}

void horiz_w32_ssse3(const uint8_t *src, const __m128i *f, int tapsNum,
                     uint8_t *buf) {
  horiz_w16_ssse3(src, f, tapsNum, buf);
  src += 16;
  buf += 16;
  horiz_w16_ssse3(src, f, tapsNum, buf);
}

void horiz_w64_ssse3(const uint8_t *src, const __m128i *f, int tapsNum,
                     uint8_t *buf) {
  horiz_w32_ssse3(src, f, tapsNum, buf);
  src += 32;
  buf += 32;
  horiz_w32_ssse3(src, f, tapsNum, buf);
}

void horiz_w128_ssse3(const uint8_t *src, const __m128i *f, int tapsNum,
                      uint8_t *buf) {
  horiz_w64_ssse3(src, f, tapsNum, buf);
  src += 64;
  buf += 64;
  horiz_w64_ssse3(src, f, tapsNum, buf);
}

void (*horizTab[6])(const uint8_t *, const __m128i *, int, uint8_t *) = {
   horiz_w4_ssse3,
   horiz_w8_ssse3,
   horiz_w16_ssse3,
   horiz_w32_ssse3,
   horiz_w64_ssse3,
   horiz_w128_ssse3,
};

void horiz_filter_ssse3(const uint8_t *src, const struct Filter fData,
                        int width, uint8_t *buffer) {
  const int16_t *filter = (const int16_t *) fData.coeffs;
  __m128i f[2];

  f[0] = *((__m128i *)(fData.coeffs));
  f[1] = *((__m128i *)(fData.coeffs + 1));

  switch (width) {
    case 4:
      horizTab[0](src, f, fData.tapsNum, buffer);
      break;
    case 8:
      horizTab[1](src, f, fData.tapsNum, buffer);
      break;
    case 16:
      horizTab[2](src, f, fData.tapsNum, buffer);
      break;
    case 32:
      horizTab[3](src, f, fData.tapsNum, buffer);
      break;
    case 64:
      horizTab[4](src, f, fData.tapsNum, buffer);
      break;
    case 128:
      horizTab[5](src, f, fData.tapsNum, buffer);
      break;
    default:
      assert(0);
  }
}

// sub-pixel 4x4 method
static void transpose4x4_to_dst(const uint8_t *src, ptrdiff_t src_stride,
                                uint8_t *dst, ptrdiff_t dst_stride) {
  __m128i A = _mm_cvtsi32_si128(*(const int *)src);
  __m128i B = _mm_cvtsi32_si128(*(const int *)(src + src_stride));
  __m128i C = _mm_cvtsi32_si128(*(const int *)(src + src_stride * 2));
  __m128i D = _mm_cvtsi32_si128(*(const int *)(src + src_stride * 3));
  // 00 10 01 11 02 12 03 13
  const __m128i tr0_0 = _mm_unpacklo_epi8(A, B);
  // 20 30 21 31 22 32 23 33
  const __m128i tr0_1 = _mm_unpacklo_epi8(C, D);
  // 00 10 20 30 01 11 21 31 02 12 22 32 03 13 23 33
  A = _mm_unpacklo_epi16(tr0_0, tr0_1);
  B = _mm_srli_si128(A, 4);
  C = _mm_srli_si128(A, 8);
  D = _mm_srli_si128(A, 12);

  *(int *)(dst) =  _mm_cvtsi128_si32(A);
  *(int *)(dst + dst_stride) =  _mm_cvtsi128_si32(B);
  *(int *)(dst + dst_stride * 2) =  _mm_cvtsi128_si32(C);
  *(int *)(dst + dst_stride * 3) =  _mm_cvtsi128_si32(D);
}

static void filter_horiz_w4_ssse3(const uint8_t *src_ptr, ptrdiff_t src_pitch,
                                  __m128i *f, int tapsNum, uint8_t *dst) {
  const __m128i k_256 = _mm_set1_epi16(1 << 8);
#if 0
  // pack and duplicate the filter values
  const __m128i f1f0 = *((__m128i *)(filter.coeffs + 0));
  const __m128i f3f2 = *((__m128i *)(filter.coeffs + 1));
  const __m128i f5f4 = *((__m128i *)(filter.coeffs + 2));
  const __m128i f7f6 = *((__m128i *)(filter.coeffs + 3));
  const __m128i f9f8 = *((__m128i *)(filter.coeffs + 4));
  const __m128i fbfa = *((__m128i *)(filter.coeffs + 5));
#endif
  if (tapsNum == 10) {
    src_ptr -= 1;
  }
  const __m128i A = _mm_loadu_si128((const __m128i *)src_ptr);
  const __m128i B = _mm_loadu_si128((const __m128i *)(src_ptr + src_pitch));
  const __m128i C = _mm_loadu_si128((const __m128i *)(src_ptr + src_pitch * 2));
  const __m128i D = _mm_loadu_si128((const __m128i *)(src_ptr + src_pitch * 3));
  // TRANSPOSE...
  // 00 01 02 03 04 05 06 07
  // 10 11 12 13 14 15 16 17
  // 20 21 22 23 24 25 26 27
  // 30 31 32 33 34 35 36 37
  //
  // TO
  //
  // 00 10 20 30
  // 01 11 21 31
  // 02 12 22 32
  // 03 13 23 33
  // 04 14 24 34
  // 05 15 25 35
  // 06 16 26 36
  // 07 17 27 37
  //
  // 00 01 10 11 02 03 12 13 04 05 14 15 06 07 16 17
  __m128i tr0_0 = _mm_unpacklo_epi16(A, B);
  // 20 21 30 31 22 23 32 33 24 25 34 35 26 27 36 37
  __m128i tr0_1 = _mm_unpacklo_epi16(C, D);
  // 00 01 10 11 20 21 30 31 02 03 12 13 22 23 32 33
  const __m128i s1s0  = _mm_unpacklo_epi32(tr0_0, tr0_1);
  // 04 05 14 15 24 25 34 35 06 07 16 17 26 27 36 37
  const __m128i s5s4 = _mm_unpackhi_epi32(tr0_0, tr0_1);
  // 02 03 12 13 22 23 32 33
  const __m128i s3s2 = _mm_srli_si128(s1s0, 8);
  // 06 07 16 17 26 27 36 37
  const __m128i s7s6 = _mm_srli_si128(s5s4, 8);

  tr0_0 = _mm_unpackhi_epi16(A, B);
  tr0_1 = _mm_unpackhi_epi16(C, D);
  const __m128i s9s8  = _mm_unpacklo_epi32(tr0_0, tr0_1);
  const __m128i sbsa = _mm_srli_si128(s9s8, 8);

  // multiply 2 adjacent elements with the filter and add the result
  const __m128i x0 = _mm_maddubs_epi16(s1s0, f[0]);
  const __m128i x1 = _mm_maddubs_epi16(s3s2, f[1]);
  const __m128i x2 = _mm_maddubs_epi16(s5s4, f[2]);
  const __m128i x3 = _mm_maddubs_epi16(s7s6, f[3]);
  const __m128i x4 = _mm_maddubs_epi16(s9s8, f[4]);
  const __m128i x5 = _mm_maddubs_epi16(sbsa, f[5]);
  // add and saturate the results together
  const __m128i min_x2x3 = _mm_min_epi16(x2, x3);
  const __m128i max_x2x3 = _mm_max_epi16(x2, x3);
  __m128i temp = _mm_adds_epi16(x0, x1);
  temp = _mm_adds_epi16(temp, x5);
  temp = _mm_adds_epi16(temp, x4);

  temp = _mm_adds_epi16(temp, min_x2x3);
  temp = _mm_adds_epi16(temp, max_x2x3);
  // round and shift by 7 bit each 16 bit
  temp = _mm_mulhrs_epi16(temp, k_256);
  // shrink to 8 bit each 16 bits
  temp = _mm_packus_epi16(temp, temp);
  // save only 4 bytes
  *(int *)dst = _mm_cvtsi128_si32(temp);
}

// Testing wrapper functions

void run_prototype_filter(uint8_t *src, int width, int height, int stride,
                          const int16_t *filter, int flen, uint8_t *dst) {
  uint32_t start, end;
  int count = 0;

  start = readtsc();
  do {
    convolve(src, width, filter, flen, dst);
    src += stride;
    dst += stride;
    count++;
  } while (count < height);
  end = readtsc();

  printf("C version cycles:\t%d\n", end - start);
}

void run_target_filter(uint8_t *src, int width, int height, int stride,
                       struct Filter filter, uint8_t *dst) {
  uint32_t start, end;
  int count = 0;

  start = readtsc();
  do {
    horiz_filter_ssse3(src, filter, width, dst);
    src += stride;
    dst += stride;
    count++;
  } while (count < height);
  end = readtsc();

  printf("SIMD HorizP cycles:\t%d\n\n", end - start);
}

// for Horiz4 method (subpixel)
void run_subpixel_filter(uint8_t *src, int width, int height, int stride,
                         const struct Filter filter, uint8_t *dst) {
  uint8_t temp[4 * 4] __attribute__ ((aligned(16)));
  __m128i f[6];
  int tapsNum;
  uint8_t *src_ptr;
  uint32_t start, end;
  int count;
  int block_height;
  int col, i;

  start = readtsc();

  tapsNum = filter.tapsNum;
  count = 0;
  block_height = height >> 2;
  src_ptr = src;
  f[0] = *((__m128i *)(filter.coeffs));
  f[1] = *((__m128i *)(filter.coeffs + 1));
  f[2] = *((__m128i *)(filter.coeffs + 2));
  f[3] = *((__m128i *)(filter.coeffs + 3));
  f[4] = *((__m128i *)(filter.coeffs + 4));
  f[5] = *((__m128i *)(filter.coeffs + 5));

  do {
    for (col = 0; col < width; col += 4) {
      for (i = 0; i < 4; ++i) {
        filter_horiz_w4_ssse3(src_ptr, stride, f, tapsNum, temp + (i * 4));
        src_ptr += 1;
      }
      transpose4x4_to_dst(temp, 4, dst + col, stride);
    }
    count++;
    src_ptr = src + count * stride * 4;
    dst += stride * 4;
  } while (count < block_height);
  end = readtsc();

  printf("SIMD Horiz4 cycles:\t%d\n\n", end - start);
}

// Test driver main()
int main(int argc, char **argv)
{
  // We simulate HD width (1920) and max super block height (128)
  const size_t block_size = HD_WIDTH * SUPER_BLOCK_HEIGHT;

  if (argc != 4) {
    printf("Usage: horiz_filter <seed> <width> <height>\n");
    printf("width/height = 4, 8, 16, 32, 64, 128, seed = random seed number.\n");
    return -1;
  }

  const int width = atoi(argv[2]);
  const unsigned int random_seed = atoi(argv[1]);
  const int height = atoi(argv[3]);
  const int stride = HD_WIDTH;

  assert(0 == width % 4);
  assert(0 == height % 4);

  uint8_t *buffer = (uint8_t *) malloc(2 * sizeof(buffer[0]) * block_size);
  uint8_t *pixel = (uint8_t *) malloc(2 * sizeof(pixel[0]) * block_size);
  uint8_t *ppixel = pixel + block_size;
  uint8_t *pbuffer = buffer + block_size;

  init_state(buffer, pixel, 8 + width, height, stride, random_seed);
  init_state(pbuffer, ppixel, 8 + width, height, stride, random_seed);
#if 1
  run_prototype_filter(pixel, width, height, stride, filter12, 12, buffer);
  run_target_filter(ppixel, width, height, stride, pfilter_12tap, pbuffer);
  check_buffer(buffer, pbuffer, width, height, stride);

  run_subpixel_filter(ppixel, width, height, stride,
                      pfilter_12tap_subpixel, pbuffer);
  check_buffer(buffer, pbuffer, width, height, stride);
#endif
  run_prototype_filter(pixel, width, height, stride, filter10, 10, buffer);
  run_target_filter(ppixel, width, height, stride, pfilter_10tap, pbuffer);
  check_buffer(buffer, pbuffer, width, height, stride);

  run_subpixel_filter(ppixel, width, height, stride,
                      pfilter_10tap_subpixel, pbuffer);
  check_buffer(buffer, pbuffer, width, height, stride);

  free(buffer);
  free(pixel);
  return 0;
}
