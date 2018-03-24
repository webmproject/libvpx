#include "vpx_dsp/ppc/types_vsx.h"

void vpx_round_store4x4_vsx(int16x8_t *in, int16x8_t *out, uint8_t *dest,
                            int stride);
void vpx_idct4_vsx(int16x8_t *in, int16x8_t *out);
void vp9_iadst4_vsx(int16x8_t *in, int16x8_t *out);
