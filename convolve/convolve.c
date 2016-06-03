#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int i;
  int height = 27;
  int x_q4 = 8;
  int x_step_q4 = 16;
  const int subpel_bits = 4;
  const int subpel_mask = (1 << subpel_bits) - 1;

  if (3 != argc) {
    printf("Usage: convolve <q4> <step_q4>\n");
    return -1;
  }

  x_q4 = atoi(argv[1]);
  x_step_q4 = atoi(argv[2]);

  for (i = 0; i < height; ++i) {
    printf("x_q4: 0x%X, sigidx: 0x%X, filidx: 0x%X\n",
           x_q4, x_q4 >> subpel_bits, x_q4 & subpel_mask);
    x_q4 += x_step_q4;
  }

  return 0;
}
