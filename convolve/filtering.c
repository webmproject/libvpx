#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <tmmintrin.h>

#define RAND_SEED (0xbeef)
unsigned int seed = RAND_SEED;

static inline unsigned int readtsc(void) {
  unsigned int tsc;
  __asm__ __volatile__("rdtsc\n\t":"=a"(tsc):);
  return tsc;
}

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
              int *buffer) {
  int i;
  int sum;

  for (i = 0; i < w; ++i) {
    sum = filtering(src, filter, flen);
    buffer[i] = sum;
    src += 1;
  }
}

void init_state(int *buf, uint8_t *pixel, int w, int block) {
  int i;

  memset(buf, 0, sizeof(buf[0]) * block);
  memset(pixel, 0, sizeof(pixel[0]) * block);
  seed = RAND_SEED;
  for (i = 0; i < w; ++i) {
    pixel[i] = clip_pixel(rand_r(&seed) % RAND_SEED);
  }
}

void check_buffer(const int *buf1, const int *buf2, int width) {
  int i;
  for (i = 0; i < width; ++i) {
    if (buf1[i] != buf2[i]) {
      printf("Not bit-exact on index %d\n", i);
      printf("Expected: %d, Actual: %d\n", buf1[i], buf2[i]);
      return;
    }
  }
}

static const int16_t filter12[12] = {
  -1,   3,  -4,   8, -18, 120,  28, -12,   7,  -4,   2, -1};

static const int16_t filter10[10] = {
  1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1};

// SSE4.1

const int16_t pfilter12[5][16] __attribute__ ((aligned(16))) = {
  {-1,  3, -4,  8, -18, 120,  28, -12,   7,  -4,   2,  -1,  0,  0,  0,  0},
  { 0, -1,  3, -4,   8, -18, 120,  28, -12,   7,  -4,   2, -1,  0,  0,  0},
  { 0,  0, -1,  3,  -4,   8, -18, 120,  28, -12,   7,  -4,  2, -1,  0,  0},
  { 0,  0,  0, -1,   3,  -4,   8, -18, 120,  28, -12,   7, -4,  2, -1,  0},
  { 0,  0,  0,  0,  -1,   3,  -4,   8, -18, 120,  28, -12,  7, -4,  2, -1},
};

const int16_t pfilter12_2[3][16] __attribute__ ((aligned(16))) = {
  {-1,    3,  -1,   3,  -1,   3,  -1,   3, -4,   8, -4,   8, -4,   8, -4,   8},
  {-18, 120, -18, 120, -18, 120, -18, 120, 28, -12, 28, -12, 28, -12, 28, -12},
  {  7,  -4,   7,  -4,   7,  -4,   7,  -4,  2,  -1,  2,  -1,  2,  -1,  2,  -1},
};

const int16_t pfilter10[7][16] __attribute__ ((aligned(16))) = {
  {1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0, 0, 0},
  {0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0, 0},
  {0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0},
  {0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0},
  {0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0},
  {0, 0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0},
  {0, 0, 0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1},
};

// Note:
//  We can do this:
//  typedef const int16_t filter_array[16];
//
//  struct Filter {
//    filter_array *coeffs;
//  };

struct Filter {
  const int16_t (*coeffs)[16];
  int tapsNum;
  int signalSpan;
};

const struct Filter pfilter_12tap = {
  pfilter12, 12, 5
};

const struct Filter pfilter_12tap_2 = {
  pfilter12_2, 12, 5
};

const struct Filter pfilter_10tap = {
  pfilter10, 10, 7
};

static inline __m128i multiply_add(const __m128i ps0, const __m128i ps1,
                                   const __m128i fl0, const __m128i fl1) {
  __m128i sum;
  __m128i sum0 = _mm_madd_epi16(ps0, fl0);
  __m128i sum1 = _mm_madd_epi16(ps1, fl1);

  sum = _mm_hadd_epi32(sum0, sum1);
  sum = _mm_hadd_epi32(sum, sum);
  return _mm_hadd_epi32(sum, sum);
}

// Note:
//  Instead of, const int16_t (*filter)[16],
//  we can also, const int16_t filter[][16]

// Note:
// This function memory access is 4 (12-tap filter) or 6 (10-tap filter) bytes
// more than C version
void inline convolve_w4_sse4_1(const uint8_t *src, const __m128i *f,
                               int *buffer) {
  __m128i u0, u1, x, s[4];
  __m128i pixel, pixello, pixelhi;
  const __m128i zero = _mm_setzero_si128();

  pixel = _mm_loadu_si128((__m128i const *)src);
  pixello = _mm_unpacklo_epi8(pixel, zero);
  pixelhi = _mm_unpackhi_epi8(pixel, zero);

  s[0] = multiply_add(pixello, pixelhi, f[0], f[1]);
  s[1] = multiply_add(pixello, pixelhi, f[2], f[3]);
  s[2] = multiply_add(pixello, pixelhi, f[4], f[5]);
  s[3] = multiply_add(pixello, pixelhi, f[6], f[7]);

  u0 = _mm_unpacklo_epi32(s[0], s[1]);
  u1 = _mm_unpacklo_epi32(s[2], s[3]);
  x = _mm_unpacklo_epi64(u0, u1);
  _mm_storeu_si128((__m128i *)buffer, x);
}

void convolve_w8_sse4_1(const uint8_t *src, const __m128i *f, int *buffer) {
  convolve_w4_sse4_1(src, f, buffer);
  src += 4;
  buffer += 4;
  convolve_w4_sse4_1(src, f, buffer);
}

void convolve_w16_sse4_1(const uint8_t *src, const __m128i *f, int *buffer) {
  convolve_w8_sse4_1(src, f, buffer);
  src += 8;
  buffer += 8;
  convolve_w8_sse4_1(src, f, buffer);
}

void convolve_w32_sse4_1(const uint8_t *src, const __m128i *f, int *buffer) {
  convolve_w16_sse4_1(src, f, buffer);
  src += 16;
  buffer += 16;
  convolve_w16_sse4_1(src, f, buffer);
}

void convolve_w64_sse4_1(const uint8_t *src, const __m128i *f, int *buffer) {
  convolve_w32_sse4_1(src, f, buffer);
  src += 32;
  buffer += 32;
  convolve_w32_sse4_1(src, f, buffer);
}

void (*convolveTab[5])(const uint8_t *, const __m128i *, int *) = {
   convolve_w4_sse4_1,
   convolve_w8_sse4_1,
   convolve_w16_sse4_1,
   convolve_w32_sse4_1,
   convolve_w64_sse4_1,
};


void convolve_sse4_1(const uint8_t *src, const struct Filter fData,
                     int width, int *buffer) {
  const int16_t *filter = (const int16_t *) fData.coeffs;
  __m128i f[10];

  f[0] = *((__m128i *)(filter + 0 * 8));
  f[1] = *((__m128i *)(filter + 1 * 8));
  f[2] = *((__m128i *)(filter + 2 * 8));
  f[3] = *((__m128i *)(filter + 3 * 8));
  f[4] = *((__m128i *)(filter + 4 * 8));
  f[5] = *((__m128i *)(filter + 5 * 8));
  f[6] = *((__m128i *)(filter + 6 * 8));
  f[7] = *((__m128i *)(filter + 7 * 8));
  f[8] = *((__m128i *)(filter + 8 * 8));
  f[9] = *((__m128i *)(filter + 9 * 8));

  switch (width) {
    case 4:
      convolveTab[0](src, f, buffer);
      break;
    case 8:
      convolveTab[1](src, f, buffer);
      break;
    case 16:
      convolveTab[2](src, f, buffer);
      break;
    case 32:
      convolveTab[3](src, f, buffer);
      break;
    case 64:
      convolveTab[4](src, f, buffer);
      break;
    default:
      assert(0);
  }
}

// Solution 2
void inline convolve2_w4_sse4_1(const uint8_t *src, const __m128i *f,
                                const __m128i *cm, int *buffer) {
  __m128i pixel;
  __m128i s0, s1, s2, s3, s4, s5;
  const __m128i zero = _mm_setzero_si128();

  pixel = _mm_loadu_si128((__m128i const *)src);

  __m128i y0 = _mm_shuffle_epi8(pixel, cm[0]);
  __m128i y1 = _mm_shuffle_epi8(pixel, cm[1]);
  __m128i y2 = _mm_shuffle_epi8(pixel, cm[2]);
  __m128i y3 = _mm_shuffle_epi8(pixel, cm[3]);
  __m128i y4 = _mm_shuffle_epi8(pixel, cm[4]);
  __m128i y5 = _mm_shuffle_epi8(pixel, cm[5]);

  y0 = _mm_unpacklo_epi8(y0, zero);
  y1 = _mm_unpacklo_epi8(y1, zero);
  y2 = _mm_unpacklo_epi8(y2, zero);
  y3 = _mm_unpacklo_epi8(y3, zero);
  y4 = _mm_unpacklo_epi8(y4, zero);
  y5 = _mm_unpacklo_epi8(y5, zero);

  s0 = _mm_madd_epi16(y0, f[0]);
  s1 = _mm_madd_epi16(y1, f[1]);
  s2 = _mm_madd_epi16(y2, f[2]);
  s3 = _mm_madd_epi16(y3, f[3]);
  s4 = _mm_madd_epi16(y4, f[4]);
  s5 = _mm_madd_epi16(y5, f[5]);

  s0 = _mm_add_epi32(s0, s1);
  s0 = _mm_add_epi32(s0, s2);
  s0 = _mm_add_epi32(s0, s3);
  s0 = _mm_add_epi32(s0, s4);
  s0 = _mm_add_epi32(s0, s5);

  _mm_storeu_si128((__m128i *)buffer, s0);
}

void convolve2_w8_sse4_1(const uint8_t *src, const __m128i *f,
                         const __m128i *cm, int *buffer) {
  convolve2_w4_sse4_1(src, f, cm, buffer);
  src += 4;
  buffer += 4;
  convolve2_w4_sse4_1(src, f, cm, buffer);
}

void convolve2_w16_sse4_1(const uint8_t *src, const __m128i *f,
                          const __m128i *cm, int *buffer) {
  convolve2_w8_sse4_1(src, f, cm, buffer);
  src += 8;
  buffer += 8;
  convolve2_w8_sse4_1(src, f, cm, buffer);
}

void convolve2_w32_sse4_1(const uint8_t *src, const __m128i *f,
                          const __m128i *cm, int *buffer) {
  convolve2_w16_sse4_1(src, f, cm, buffer);
  src += 16;
  buffer += 16;
  convolve2_w16_sse4_1(src, f, cm, buffer);
}

void convolve2_w64_sse4_1(const uint8_t *src, const __m128i *f,
                          const __m128i *cm, int *buffer) {
  convolve2_w32_sse4_1(src, f, cm, buffer);
  src += 32;
  buffer += 32;
  convolve2_w32_sse4_1(src, f, cm, buffer);
}

void (*convolveTab2[5])(const uint8_t *, const __m128i *, const __m128i *,
                        int *) = {
   convolve2_w4_sse4_1,
   convolve2_w8_sse4_1,
   convolve2_w16_sse4_1,
   convolve2_w32_sse4_1,
   convolve2_w64_sse4_1,
};

void convolve2_sse4_1(const uint8_t *src, const struct Filter fData,
                     int width, int *buffer) {
  __m128i f[6];
  __m128i cm[6];
  const int16_t *filter = (const int16_t *) fData.coeffs;

  f[0] = *((__m128i *)(filter + 0 * 8));
  f[1] = *((__m128i *)(filter + 1 * 8));
  f[2] = *((__m128i *)(filter + 2 * 8));
  f[3] = *((__m128i *)(filter + 3 * 8));
  f[4] = *((__m128i *)(filter + 4 * 8));
  f[5] = *((__m128i *)(filter + 5 * 8));

  cm[0] = _mm_setr_epi16(0x100, 0x201, 0x302, 0x403, 0, 0, 0, 0);
  cm[1] = _mm_setr_epi16(0x302, 0x403, 0x504, 0x605, 0, 0, 0, 0);
  cm[2] = _mm_setr_epi16(0x504, 0x605, 0x706, 0x807, 0, 0, 0, 0);
  cm[3] = _mm_setr_epi16(0x706, 0x807, 0x908, 0xa09, 0, 0, 0, 0);
  cm[4] = _mm_setr_epi16(0x908, 0xa09, 0xb0a, 0xc0b, 0, 0, 0, 0);
  cm[5] = _mm_setr_epi16(0xb0a, 0xc0b, 0xd0c, 0xe0d, 0, 0, 0, 0);

  switch (width) {
    case 4:
      convolveTab2[0](src, f, cm, buffer);
      break;
    case 8:
      convolveTab2[1](src, f, cm, buffer);
      break;
    case 16:
      convolveTab2[2](src, f, cm, buffer);
      break;
    case 32:
      convolveTab2[3](src, f, cm, buffer);
      break;
    case 64:
      convolveTab2[4](src, f, cm, buffer);
      break;
    default:
      assert(0);
  }
}

#define TEST_NUM (32)

int main(int argc, char **argv)
{
  const size_t block_size = 256;

  if (argc != 2) {
    printf("Usage: filtering <width>, where width = 4, 8, 16, 32, 64\n");
    return -1;
  }

  const int width = atoi(argv[1]);

  int *buffer = (int *) malloc(3 * sizeof(buffer[0]) * block_size);
  uint8_t *pixel = (uint8_t *) malloc(3 * sizeof(pixel[0]) * block_size);
  uint8_t *ppixel = pixel + block_size;
  int *pbuffer = buffer + block_size;
  uint8_t *ppixel2 = pixel + 2 * block_size;
  int *pbuffer2 = buffer + 2 * block_size;
  uint32_t start, end;
  int count;

  init_state(buffer, pixel, width, block_size);
  init_state(pbuffer, ppixel, width, block_size);
  init_state(pbuffer2, ppixel2, width, block_size);

  count = 0;
  start = readtsc();
  do {
    convolve(pixel, width, filter12, 12, buffer);
    count++;
  } while (count < TEST_NUM);
  end = readtsc();

  printf("C version cycles: %d\n", end - start);

  // Solution 1
  count = 0;
  start = readtsc();
  do {
    convolve_sse4_1(ppixel, pfilter_12tap, width, pbuffer);
    count++;
  } while (count < TEST_NUM);
  end = readtsc();

  printf("SIMD version cycles: %d\n", end - start);

  check_buffer(buffer, pbuffer, width);

  // Solution 2
  count = 0;
  start = readtsc();
  do {
    convolve2_sse4_1(ppixel2, pfilter_12tap_2, width, pbuffer2);
    count++;
  } while (count < TEST_NUM);
  end = readtsc();

  printf("SIMD version 2 cycles: %d\n", end - start);

  check_buffer(buffer, pbuffer2, width);

  free(buffer);
  free(pixel);
  return 0;
}
