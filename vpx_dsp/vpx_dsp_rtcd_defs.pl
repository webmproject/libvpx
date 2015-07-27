sub vpx_dsp_forward_decls() {
print <<EOF
/*
 * DSP
 */

#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"

EOF
}
forward_decls qw/vpx_dsp_forward_decls/;

# x86inc.asm had specific constraints. break it out so it's easy to disable.
# zero all the variables to avoid tricky else conditions.
$mmx_x86inc = $sse_x86inc = $sse2_x86inc = $ssse3_x86inc = $avx_x86inc =
  $avx2_x86inc = '';
$mmx_x86_64_x86inc = $sse_x86_64_x86inc = $sse2_x86_64_x86inc =
  $ssse3_x86_64_x86inc = $avx_x86_64_x86inc = $avx2_x86_64_x86inc = '';
if (vpx_config("CONFIG_USE_X86INC") eq "yes") {
  $mmx_x86inc = 'mmx';
  $sse_x86inc = 'sse';
  $sse2_x86inc = 'sse2';
  $ssse3_x86inc = 'ssse3';
  $avx_x86inc = 'avx';
  $avx2_x86inc = 'avx2';
  if ($opts{arch} eq "x86_64") {
    $mmx_x86_64_x86inc = 'mmx';
    $sse_x86_64_x86inc = 'sse';
    $sse2_x86_64_x86inc = 'sse2';
    $ssse3_x86_64_x86inc = 'ssse3';
    $avx_x86_64_x86inc = 'avx';
    $avx2_x86_64_x86inc = 'avx2';
  }
}

# functions that are 64 bit only.
$mmx_x86_64 = $sse2_x86_64 = $ssse3_x86_64 = $avx_x86_64 = $avx2_x86_64 = '';
if ($opts{arch} eq "x86_64") {
  $mmx_x86_64 = 'mmx';
  $sse2_x86_64 = 'sse2';
  $ssse3_x86_64 = 'ssse3';
  $avx_x86_64 = 'avx';
  $avx2_x86_64 = 'avx2';
}

#
# Loopfilter
#
add_proto qw/void vpx_lpf_vertical_16/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_16 sse2 neon_asm msa/;
$vpx_lpf_vertical_16_neon_asm=vpx_lpf_vertical_16_neon;

add_proto qw/void vpx_lpf_vertical_16_dual/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_16_dual sse2 neon_asm msa/;
$vpx_lpf_vertical_16_dual_neon_asm=vpx_lpf_vertical_16_dual_neon;

add_proto qw/void vpx_lpf_vertical_8/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count";
specialize qw/vpx_lpf_vertical_8 sse2 neon msa/;

add_proto qw/void vpx_lpf_vertical_8_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_vertical_8_dual sse2 neon_asm msa/;
$vpx_lpf_vertical_8_dual_neon_asm=vpx_lpf_vertical_8_dual_neon;

add_proto qw/void vpx_lpf_vertical_4/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count";
specialize qw/vpx_lpf_vertical_4 mmx neon msa/;

add_proto qw/void vpx_lpf_vertical_4_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_vertical_4_dual sse2 neon msa/;

add_proto qw/void vpx_lpf_horizontal_16/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count";
specialize qw/vpx_lpf_horizontal_16 sse2 avx2 neon_asm msa/;
$vpx_lpf_horizontal_16_neon_asm=vpx_lpf_horizontal_16_neon;

add_proto qw/void vpx_lpf_horizontal_8/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count";
specialize qw/vpx_lpf_horizontal_8 sse2 neon msa/;

add_proto qw/void vpx_lpf_horizontal_8_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_horizontal_8_dual sse2 neon_asm msa/;
$vpx_lpf_horizontal_8_dual_neon_asm=vpx_lpf_horizontal_8_dual_neon;

add_proto qw/void vpx_lpf_horizontal_4/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count";
specialize qw/vpx_lpf_horizontal_4 mmx neon msa/;

add_proto qw/void vpx_lpf_horizontal_4_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_horizontal_4_dual sse2 neon msa/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_lpf_vertical_16/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_16 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_16_dual/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_16_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_8/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count, int bd";
  specialize qw/vpx_highbd_lpf_vertical_8 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_8_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_vertical_8_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_4/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count, int bd";
  specialize qw/vpx_highbd_lpf_vertical_4 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_4_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_vertical_4_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_16/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_16 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_8/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_8 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_8_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_8_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_4/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_4 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_4_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_4_dual sse2/;
}  # CONFIG_VP9_HIGHBITDEPTH

#
# Encoder functions.
#

#
# Forward transform
#
if (vpx_config("CONFIG_VP9_ENCODER") eq "yes") {
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vp9_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct4x4 sse2/;

  add_proto qw/void vp9_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct8x8 sse2/;

  add_proto qw/void vp9_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct16x16 sse2/;

  add_proto qw/void vp9_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct32x32 sse2/;

  add_proto qw/void vp9_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct32x32_rd sse2/;

  add_proto qw/void vp9_highbd_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_highbd_fdct4x4 sse2/;

  add_proto qw/void vp9_highbd_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_highbd_fdct8x8 sse2/;

  add_proto qw/void vp9_highbd_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_highbd_fdct16x16 sse2/;

  add_proto qw/void vp9_highbd_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_highbd_fdct32x32 sse2/;

  add_proto qw/void vp9_highbd_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_highbd_fdct32x32_rd sse2/;
} else {
  add_proto qw/void vp9_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct4x4 sse2 msa/;

  add_proto qw/void vp9_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct8x8 sse2 neon msa/, "$ssse3_x86_64_x86inc";

  add_proto qw/void vp9_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct16x16 sse2 msa/;

  add_proto qw/void vp9_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct32x32 sse2 avx2 msa/;

  add_proto qw/void vp9_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vp9_fdct32x32_rd sse2 avx2 msa/;
}  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_VP9_ENCODER

#
# Quantization
#
if (vpx_config("CONFIG_VP9_ENCODER") eq "yes") {
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vp9_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_quantize_b/;

  add_proto qw/void vp9_quantize_b_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_quantize_b_32x32/;

  add_proto qw/void vp9_highbd_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_highbd_quantize_b sse2/;

  add_proto qw/void vp9_highbd_quantize_b_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_highbd_quantize_b_32x32 sse2/;
} else {
  add_proto qw/void vp9_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_quantize_b sse2/, "$ssse3_x86_64_x86inc";

  add_proto qw/void vp9_quantize_b_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vp9_quantize_b_32x32/, "$ssse3_x86_64_x86inc";
}  # CONFIG_VP9_ENCODER
}  # CONFIG_VP9_HIGHBITDEPTH

if (vpx_config("CONFIG_ENCODERS") eq "yes") {
#
# Block subtraction
#
add_proto qw/void vpx_subtract_block/, "int rows, int cols, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr, ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride";
specialize qw/vpx_subtract_block neon msa/, "$sse2_x86inc";

#
# Single block SAD
#
add_proto qw/unsigned int vpx_sad64x64/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad64x64 avx2 neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad64x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad64x32 avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x64/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad32x64 avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad32x32 avx2 neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad32x16 avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad16x32 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad16x16 mmx media neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad16x8 mmx neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad8x16 mmx neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad8x8 mmx neon msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x4/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad8x4 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad4x8 msa/, "$sse_x86inc";

add_proto qw/unsigned int vpx_sad4x4/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
specialize qw/vpx_sad4x4 mmx neon msa/, "$sse_x86inc";

#
# Avg
#
add_proto qw/unsigned int vpx_sad64x64_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad64x64_avg avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad64x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad64x32_avg avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x64_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad32x64_avg avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad32x32_avg avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad32x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad32x16_avg avx2 msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad16x32_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad16x16_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad16x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad16x8_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad8x16_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad8x8_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad8x4_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad8x4_avg msa/, "$sse2_x86inc";

add_proto qw/unsigned int vpx_sad4x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad4x8_avg msa/, "$sse_x86inc";

add_proto qw/unsigned int vpx_sad4x4_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
specialize qw/vpx_sad4x4_avg msa/, "$sse_x86inc";

#
# Multi-block SAD, comparing a reference to N blocks 1 pixel apart horizontally
#
# Blocks of 3
add_proto qw/void vpx_sad64x64x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad64x64x3 msa/;

add_proto qw/void vpx_sad32x32x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad32x32x3 msa/;

add_proto qw/void vpx_sad16x16x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x16x3 sse3 ssse3 msa/;

add_proto qw/void vpx_sad16x8x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x8x3 sse3 ssse3 msa/;

add_proto qw/void vpx_sad8x16x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x16x3 sse3 msa/;

add_proto qw/void vpx_sad8x8x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x8x3 sse3 msa/;

add_proto qw/void vpx_sad4x4x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x4x3 sse3 msa/;

# Blocks of 8
add_proto qw/void vpx_sad64x64x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad64x64x8 msa/;

add_proto qw/void vpx_sad32x32x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad32x32x8 msa/;

add_proto qw/void vpx_sad16x16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x16x8 sse4_1 msa/;

add_proto qw/void vpx_sad16x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x8x8 sse4_1 msa/;

add_proto qw/void vpx_sad8x16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x16x8 sse4_1 msa/;

add_proto qw/void vpx_sad8x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x8x8 sse4_1 msa/;

add_proto qw/void vpx_sad8x4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x4x8 msa/;

add_proto qw/void vpx_sad4x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x8x8 msa/;

add_proto qw/void vpx_sad4x4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x4x8 sse4_1 msa/;

#
# Multi-block SAD, comparing a reference to N independent blocks
#
add_proto qw/void vpx_sad64x64x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad64x64x4d avx2 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad64x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad64x32x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad32x64x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad32x64x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad32x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad32x32x4d avx2 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad32x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad32x16x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad16x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x32x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad16x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x16x4d neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad16x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x8x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad8x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x16x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad8x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x8x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad8x4x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x4x4d msa/, "$sse2_x86inc";

add_proto qw/void vpx_sad4x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x8x4d msa/, "$sse_x86inc";

add_proto qw/void vpx_sad4x4x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x4x4d msa/, "$sse_x86inc";

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  #
  # Block subtraction
  #
  add_proto qw/void vpx_highbd_subtract_block/, "int rows, int cols, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr, ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride, int bd";
  specialize qw/vpx_highbd_subtract_block/;

  #
  # Single block SAD
  #
  add_proto qw/unsigned int vpx_highbd_sad64x64/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad64x64/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad64x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad64x32/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x64/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad32x64/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad32x32/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad32x16/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x32/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad16x32/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad16x16/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad16x8/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x16/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad8x16/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad8x8/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x4/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad8x4/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad4x8/;

  add_proto qw/unsigned int vpx_highbd_sad4x4/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  specialize qw/vpx_highbd_sad4x4/;

  #
  # Avg
  #
  add_proto qw/unsigned int vpx_highbd_sad64x64_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad64x64_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad64x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad64x32_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x64_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad32x64_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad32x32_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad32x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad32x16_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x32_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad16x32_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad16x16_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad16x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad16x8_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x16_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad8x16_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad8x8_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad8x4_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad8x4_avg/, "$sse2_x86inc";

  add_proto qw/unsigned int vpx_highbd_sad4x8_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad4x8_avg/;

  add_proto qw/unsigned int vpx_highbd_sad4x4_avg/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
  specialize qw/vpx_highbd_sad4x4_avg/;

  #
  # Multi-block SAD, comparing a reference to N blocks 1 pixel apart horizontally
  #
  # Blocks of 3
  add_proto qw/void vpx_highbd_sad64x64x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad64x64x3/;

  add_proto qw/void vpx_highbd_sad32x32x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad32x32x3/;

  add_proto qw/void vpx_highbd_sad16x16x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x16x3/;

  add_proto qw/void vpx_highbd_sad16x8x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x8x3/;

  add_proto qw/void vpx_highbd_sad8x16x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x16x3/;

  add_proto qw/void vpx_highbd_sad8x8x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x8x3/;

  add_proto qw/void vpx_highbd_sad4x4x3/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad4x4x3/;

  # Blocks of 8
  add_proto qw/void vpx_highbd_sad64x64x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad64x64x8/;

  add_proto qw/void vpx_highbd_sad32x32x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad32x32x8/;

  add_proto qw/void vpx_highbd_sad16x16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x16x8/;

  add_proto qw/void vpx_highbd_sad16x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x8x8/;

  add_proto qw/void vpx_highbd_sad8x16x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x16x8/;

  add_proto qw/void vpx_highbd_sad8x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x8x8/;

  add_proto qw/void vpx_highbd_sad8x4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x4x8/;

  add_proto qw/void vpx_highbd_sad4x8x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad4x8x8/;

  add_proto qw/void vpx_highbd_sad4x4x8/, "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad4x4x8/;

  #
  # Multi-block SAD, comparing a reference to N independent blocks
  #
  add_proto qw/void vpx_highbd_sad64x64x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad64x64x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad64x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad64x32x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad32x64x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad32x64x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad32x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad32x32x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad32x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad32x16x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad16x32x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x32x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad16x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x16x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad16x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad16x8x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad8x16x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x16x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad8x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x8x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad8x4x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad8x4x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad4x8x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad4x8x4d/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_sad4x4x4d/, "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, uint32_t *sad_array";
  specialize qw/vpx_highbd_sad4x4x4d/, "$sse2_x86inc";

}  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_ENCODERS

if (vpx_config("CONFIG_ENCODERS") eq "yes" || vpx_config("CONFIG_POSTPROC") eq "yes" || vpx_config("CONFIG_VP9_POSTPROC") eq "yes") {

#
# Variance
#
add_proto qw/unsigned int vpx_variance64x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance64x64 sse2 avx2 neon msa/;

add_proto qw/unsigned int vpx_variance64x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance64x32 sse2 avx2 neon msa/;

add_proto qw/unsigned int vpx_variance32x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance32x64 sse2 neon msa/;

add_proto qw/unsigned int vpx_variance32x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance32x32 sse2 avx2 neon msa/;

add_proto qw/unsigned int vpx_variance32x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance32x16 sse2 avx2 msa/;

add_proto qw/unsigned int vpx_variance16x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance16x32 sse2 msa/;

add_proto qw/unsigned int vpx_variance16x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance16x16 mmx sse2 avx2 media neon msa/;

add_proto qw/unsigned int vpx_variance16x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance16x8 mmx sse2 neon msa/;

add_proto qw/unsigned int vpx_variance8x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance8x16 mmx sse2 neon msa/;

add_proto qw/unsigned int vpx_variance8x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance8x8 mmx sse2 media neon msa/;

add_proto qw/unsigned int vpx_variance8x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance8x4 sse2 msa/;

add_proto qw/unsigned int vpx_variance4x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance4x8 sse2 msa/;

add_proto qw/unsigned int vpx_variance4x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_variance4x4 mmx sse2 msa/;

#
# Specialty Variance
#
add_proto qw/void vpx_get16x16var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
  specialize qw/vpx_get16x16var sse2 avx2 neon msa/;

add_proto qw/void vpx_get8x8var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
  specialize qw/vpx_get8x8var mmx sse2 neon msa/;

add_proto qw/unsigned int vpx_mse16x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_mse16x16 mmx sse2 avx2 media neon msa/;

add_proto qw/unsigned int vpx_mse16x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_mse16x8 sse2 msa/;

add_proto qw/unsigned int vpx_mse8x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_mse8x16 sse2 msa/;

add_proto qw/unsigned int vpx_mse8x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_mse8x8 sse2 msa/;

add_proto qw/unsigned int vpx_get_mb_ss/, "const int16_t *";
  specialize qw/vpx_get_mb_ss mmx sse2 msa/;

add_proto qw/unsigned int vpx_get4x4sse_cs/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride";
  specialize qw/vpx_get4x4sse_cs neon msa/;

add_proto qw/void vpx_comp_avg_pred/, "uint8_t *comp_pred, const uint8_t *pred, int width, int height, const uint8_t *ref, int ref_stride";

#
# Subpixel Variance
#
add_proto qw/uint32_t vpx_sub_pixel_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance64x64 avx2 neon msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance64x32 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance32x64 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance32x32 avx2 neon msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance32x16 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance16x32 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance16x16 mmx media neon msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance16x8 mmx msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance8x16 mmx msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance8x8 mmx media neon msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance8x4 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance4x8 msa/, "$sse_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_sub_pixel_variance4x4 mmx msa/, "$sse_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance64x64 avx2 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance64x32 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance32x64 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance32x32 avx2 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance32x16 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance16x32 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance16x16 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance16x8 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance8x16 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance8x8 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance8x4 msa/, "$sse2_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance4x8 msa/, "$sse_x86inc", "$ssse3_x86inc";

add_proto qw/uint32_t vpx_sub_pixel_avg_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_sub_pixel_avg_variance4x4 msa/, "$sse_x86inc", "$ssse3_x86inc";

#
# Specialty Subpixel
#
add_proto qw/uint32_t vpx_variance_halfpixvar16x16_h/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, uint32_t *sse";
  specialize qw/vpx_variance_halfpixvar16x16_h mmx media/;

add_proto qw/uint32_t vpx_variance_halfpixvar16x16_v/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, uint32_t *sse";
  specialize qw/vpx_variance_halfpixvar16x16_v mmx media/;

add_proto qw/uint32_t vpx_variance_halfpixvar16x16_hv/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, uint32_t *sse";
  specialize qw/vpx_variance_halfpixvar16x16_hv mmx media/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/unsigned int vpx_highbd_12_variance64x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance64x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance64x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance64x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance32x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance32x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance32x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance32x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance32x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance32x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance16x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance16x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance16x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance16x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance16x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance8x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance8x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance8x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_variance8x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_variance8x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_12_variance4x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_12_variance4x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";

  add_proto qw/unsigned int vpx_highbd_10_variance64x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance64x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance64x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance64x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance32x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance32x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance32x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance32x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance32x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance32x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance16x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance16x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance16x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance16x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance16x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance8x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance8x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance8x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_variance8x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_variance8x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_10_variance4x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_10_variance4x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";

  add_proto qw/unsigned int vpx_highbd_8_variance64x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance64x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance64x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance64x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance32x64/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance32x64 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance32x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance32x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance32x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance32x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance16x32/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance16x32 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance16x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance16x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance16x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance8x16/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance8x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance8x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_variance8x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_variance8x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_8_variance4x8/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_8_variance4x4/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";

  add_proto qw/void vpx_highbd_8_get16x16var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
  add_proto qw/void vpx_highbd_8_get8x8var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";

  add_proto qw/void vpx_highbd_10_get16x16var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
  add_proto qw/void vpx_highbd_10_get8x8var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";

  add_proto qw/void vpx_highbd_12_get16x16var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
  add_proto qw/void vpx_highbd_12_get8x8var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";

  add_proto qw/unsigned int vpx_highbd_8_mse16x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_mse16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_8_mse16x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_8_mse8x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_8_mse8x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_8_mse8x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_mse16x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_mse16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_10_mse16x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_10_mse8x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_10_mse8x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_10_mse8x8 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_mse16x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_mse16x16 sse2/;

  add_proto qw/unsigned int vpx_highbd_12_mse16x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_12_mse8x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  add_proto qw/unsigned int vpx_highbd_12_mse8x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
  specialize qw/vpx_highbd_12_mse8x8 sse2/;

  add_proto qw/void vpx_highbd_comp_avg_pred/, "uint16_t *comp_pred, const uint8_t *pred8, int width, int height, const uint8_t *ref8, int ref_stride";

  #
  # Subpixel Variance
  #
  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_12_sub_pixel_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_10_sub_pixel_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  specialize qw/vpx_highbd_8_sub_pixel_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_12_sub_pixel_avg_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  add_proto qw/uint32_t vpx_highbd_12_sub_pixel_avg_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_10_sub_pixel_avg_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  add_proto qw/uint32_t vpx_highbd_10_sub_pixel_avg_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance64x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance64x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance64x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance64x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance32x64/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance32x64/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance32x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance32x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance32x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance32x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance16x32/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance16x32/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance16x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance16x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance16x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance16x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance8x16/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance8x16/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance8x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance8x8/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance8x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  specialize qw/vpx_highbd_8_sub_pixel_avg_variance8x4/, "$sse2_x86inc";

  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance4x8/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
  add_proto qw/uint32_t vpx_highbd_8_sub_pixel_avg_variance4x4/, "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";

}  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_ENCODERS || CONFIG_POSTPROC || CONFIG_VP9_POSTPROC

1;
