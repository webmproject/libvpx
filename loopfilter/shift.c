#include <stdio.h>
#include <stdint.h>

//----------------------------------------------------------------------------
// Note:
// Establish 64x64 block, contructed by 256 (16x16) 4x4 sub-block.
// Every 4 rows would be represented by one uint64_t mask. Hence,
// there are 4 uint64_t bitmask[4] to represent the whole 64x64.
//
// Given a location by (idx, idy), This function returns the index
// 0, 1, 2, 3 to select which bitmask[] to use.
// Then the pointer y_shift contains the shift value in the bit mask.
// Function returns y_shift; y_index contains the index.
//
int get_y_index_shift(int idx, int idy, int *y_index) {
  *y_index = idy >> 4;
  const int y_idy = (idy >> 2) % 4;
  return (y_idy << 4) + (idx >> 2);
}

// Note:
// For 4:2:0 format sampling, establish 32x32 block, constructed by
// 64 (8x8), 4x4 sub-block. We need one uint64_t bitmask to present
// all edge information
// Function returns uv_shift.
//
int get_uv_index_shift(int idx, int idy) {
  return ((idy >> 3) << 3) + (idx >> 3);
}

//----------------------------------------------------------------------------

// AV1 has 4x4 coding block
// I use 4 uint64_t integer to describe block edge information by a bit mask
//
void get_y_shift_value(int size) {
  int x, y;

  for (y = 0; y < size; y += 4) {
    for (x = 0; x < size; x += 4) {
      printf("[%02d,%02d] ", x, y);
    }
    printf("\n");
  }

  printf("\n");
  int v_index;
  int y_shift;
  for (y = 0; y < size; y += 4) {
    for (x = 0; x < size; x += 4) {
      // cb8x8
      //int shift = ((y >> 3) << 3) + (x >> 3);
      //printf("%02d ", shift);

      // cb4x4
      y_shift = get_y_index_shift(x, y, &v_index);
      printf("%02d ", y_shift);
    }
    printf("Index %d\n", v_index);
  }
}

void get_uv_shift_value(int size) {
  int x, y;

  int uv_shift = 0;
  const int step = 4;
  for (y = 0; y < size; y += step) {
    for (x = 0; x < size; x += step) {
      // cb8x8
      // int uv_shift = ((y >> 3) << 2) + (x >> 3);

      // cb4x4
      uv_shift = get_uv_index_shift(x, y);
      printf("%02d ", uv_shift);
    }
    printf("\n");
  }
}

//---------------------------------------------------------------------------
int get_uv_shift(int idx, int idy) {
  return (((idy - 2) >> 2) << 3) + (idx >> 2);
}
//---------------------------------------------------------------------------

// AV1: AV1=1
// VP9: AV1=0
#define AV1 1

#if AV1
#define MAX_MIB_SIZE_LOG2 (4)
const int num = 16;
typedef struct {
  uint64_t bits[4];
} FilterMaskY;
#else
#define MAX_MIB_SIZE_LOG2 (3)
const int num = 8;
#endif


int main() {
  get_y_shift_value(64);
  printf("\n");
  get_uv_shift_value(64);
  printf("\n");

  int y_index = 0;
  const int x = 0;
  int y;
  int i;

  // Remaining rows are 1, 2, ..., num - 1
  // VP9 : 1-7
  // AV1 : 1-15
  for (i = 1; i < num; ++i) {
#if AV1
    y = i << 2;
    int y_shift = get_y_index_shift(x, y, &y_index);
    int uv_shift = get_uv_shift(x >> 1, y >> 1);

    printf("[%02d,%02d] index=%d y_shift=%02d uv_shift=%02d mask_y ",
           x, y, y_index, y_shift, uv_shift);

    FilterMaskY mask = {0, 0, 0, 0};
    int j;
    for (j = 0; j < y_index; ++j) {
      mask.bits[j] = 0xffffffffffffffffULL;
    }
    mask.bits[y_index] = ((uint64_t)1 << y_shift) - 1;
    for (j = 0; j < 4; ++j) {
      printf("0x%016llx ", (unsigned long long int)mask.bits[j]);
    }

    uint64_t mask_uv = (((uint64_t)1 << (uv_shift + 8)) - 1);
    if (uv_shift + 8 == 64) mask_uv = 0xffffffffffffffffULL;

    printf("mask_uv 0x%016llx", (unsigned long long int)mask_uv);
    printf("\n");
#else
    const uint64_t mask_y = (((uint64_t)1 << (i << MAX_MIB_SIZE_LOG2)) - 1);
    const uint16_t mask_uv =
        (((uint16_t)1 << (((i + 1) >> 1) << (MAX_MIB_SIZE_LOG2 - 1))) - 1);
    printf("mask_y=%016llx, mask_uv=%04x\n", (long long unsigned int)mask_y, mask_uv);
#endif
  }

  return 0;
}
