#include <assert.h>
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

#define FILTER_BITS 7
#define RAND_SEED (0xabcd)
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

void init_state(uint8_t *buf, uint8_t *pixel, int w, int block) {
  int i;

  memset(buf, 0, sizeof(buf[0]) * block);
  memset(pixel, 0, sizeof(pixel[0]) * block);

  seed = RAND_SEED;
  for (i = 0; i < w; ++i) {
    pixel[i] = clip_pixel(rand_r(&seed) % 255);
  }
}

void check_buffer(const uint8_t *buf1, const uint8_t *buf2, int width) {
  int i;
  for (i = 0; i < width; ++i) {
    if (buf1[i] != buf2[i]) {
      printf("Not bit-exact on index %d\n", i);
      printf("Expected: 0x%x, Actual: 0x%x\n", buf1[i], buf2[i]);
      return;
    }
  }
}

static const int16_t filter12[12] = {
  -1,   3,  -4,   8, -18, 120,  28, -12,   7,  -4,   2, -1};

static const int16_t filter10[10] = {
  1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1};

// SSSE3

const int8_t pfilter12[3][16] __attribute__ ((aligned(16))) = {
  {-1,  3, -4,  8, -18, 120,  28, -12,   7,  -4,   2,  -1,  0,  0,  0,  0},
  { 0,  0, -1,  3,  -4,   8, -18, 120,  28, -12,   7,  -4,  2, -1,  0,  0},
  { 0,  0,  0,  0,  -1,   3,  -4,   8, -18, 120,  28, -12,  7, -4,  2, -1},
};

const int8_t pfilter10[3][16] __attribute__ ((aligned(16))) = {
  //{1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0, 0, 0},
  {0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0, 0},
  //{0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0, 0},
  {0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0, 0},
  //{0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0, 0},
  {0, 0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1, 0},
  //{0, 0, 0, 0, 0, 0, 1,  -3,   7, -17, 119,  28, -11,   5,  -2, 1},
};

struct Filter {
  const int8_t (*coeffs)[16];
  int tapsNum;
  int signalSpan;
};

const struct Filter pfilter_12tap = {
  pfilter12, 12, 5
};

const struct Filter pfilter_10tap = {
  pfilter10, 10, 7
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
  out[6] = _mm_unpackhi_epi32(t0, t1);  // 0?
  out[7] = _mm_srli_si128(out[6], 8);   // 0?
}

void horiz_w4_ssse3(const uint8_t *src, const __m128i *f,
                    uint8_t *buffer) {
  __m128i sr[4];
  __m128i sc[8];
  __m128i pixel;
  const __m128i k_256 = _mm_set1_epi16(1 << 8);

  pixel = _mm_loadu_si128((__m128i const *)src);
  sr[0] = _mm_maddubs_epi16(pixel, f[0]);
  sr[2] = _mm_maddubs_epi16(pixel, f[1]);

  pixel = _mm_loadu_si128((__m128i const *)(src + 1));
  sr[1] = _mm_maddubs_epi16(pixel, f[0]);
  sr[3] = _mm_maddubs_epi16(pixel, f[1]);

  transpose_4x8(sr, sc);

  sr[0] = _mm_adds_epi16(sc[0], sc[1]);
  sr[0] = _mm_adds_epi16(sr[0], sc[2]);

  sr[1] = _mm_adds_epi16(sc[7], sc[6]);
  sr[1] = _mm_adds_epi16(sr[1], sc[5]);

  sr[2] = _mm_min_epi16(sc[3], sc[4]);
  sr[3] = _mm_max_epi16(sc[3], sc[4]);

  sr[0] = _mm_adds_epi16(sr[0], sr[1]);
  sr[0] = _mm_adds_epi16(sr[0], sr[2]);
  sr[0] = _mm_adds_epi16(sr[0], sr[3]);

  sr[1] = _mm_mulhrs_epi16(sr[0], k_256);
  sr[2] = _mm_packus_epi16(sr[1], sr[1]);

  *(int *)buffer = _mm_cvtsi128_si32(sr[2]);
}

void horiz_w8_ssse3(const uint8_t *src, const __m128i *f, uint8_t *buf) {
horiz_w4_ssse3(src, f, buf);
src += 4;
buf += 4;
horiz_w4_ssse3(src, f, buf);
}

void horiz_w16_ssse3(const uint8_t *src, const __m128i *f, uint8_t *buf) {
horiz_w8_ssse3(src, f, buf);
src += 8;
buf += 8;
horiz_w8_ssse3(src, f, buf);
}

void horiz_w32_ssse3(const uint8_t *src, const __m128i *f, uint8_t *buf) {
horiz_w16_ssse3(src, f, buf);
src += 16;
buf += 16;
horiz_w16_ssse3(src, f, buf);
}

void horiz_w64_ssse3(const uint8_t *src, const __m128i *f, uint8_t *buf) {
horiz_w32_ssse3(src, f, buf);
src += 32;
buf += 32;
horiz_w32_ssse3(src, f, buf);
}

void (*horizTab[5])(const uint8_t *, const __m128i *, uint8_t *) = {
   horiz_w4_ssse3,
   horiz_w8_ssse3,
   horiz_w16_ssse3,
   horiz_w32_ssse3,
   horiz_w64_ssse3,
};

void horiz_filter_ssse3(const uint8_t *src, const struct Filter fData,
                        int width, uint8_t *buffer) {
  const int16_t *filter = (const int16_t *) fData.coeffs;
  __m128i f[7];

  if (fData.tapsNum == 12) {
    f[0] = *((__m128i *)(fData.coeffs));
    f[1] = *((__m128i *)(fData.coeffs + 1));
    f[2] = *((__m128i *)(fData.coeffs + 2));
    f[3] = *((__m128i *)(fData.coeffs + 3));
  } else {
    f[0] = *((__m128i *)(fData.coeffs));
    f[1] = *((__m128i *)(fData.coeffs + 1));
    f[2] = *((__m128i *)(fData.coeffs + 2));
    f[3] = *((__m128i *)(fData.coeffs + 3));
    f[4] = *((__m128i *)(fData.coeffs + 4));
    f[5] = *((__m128i *)(fData.coeffs + 5));
    f[6] = *((__m128i *)(fData.coeffs + 6));
  }

  switch (width) {
    case 4:
      horizTab[0](src, f, buffer);
      break;
    case 8:
      horizTab[1](src, f, buffer);
      break;
    case 16:
      horizTab[2](src, f, buffer);
      break;
    case 32:
      horizTab[3](src, f, buffer);
      break;
    case 64:
      horizTab[4](src, f, buffer);
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

  uint8_t *buffer = (uint8_t *) malloc(2 * sizeof(buffer[0]) * block_size);
  uint8_t *pixel = (uint8_t *) malloc(2 * sizeof(pixel[0]) * block_size);
  uint8_t *ppixel = pixel + block_size;
  uint8_t *pbuffer = buffer + block_size;

  uint32_t start, end;
  int count;

  init_state(buffer, pixel, width, block_size);
  init_state(pbuffer, ppixel, width, block_size);

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
    horiz_filter_ssse3(ppixel, pfilter_12tap, width, pbuffer);
    count++;
  } while (count < TEST_NUM);
  end = readtsc();

  printf("SIMD version cycles: %d\n", end - start);

  check_buffer(buffer, pbuffer, width);

  free(buffer);
  free(pixel);
  return 0;
}
