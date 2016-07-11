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

# optimizations which depend on multiple features
$avx2_ssse3 = '';
if ((vpx_config("HAVE_AVX2") eq "yes") && (vpx_config("HAVE_SSSE3") eq "yes")) {
  $avx2_ssse3 = 'avx2';
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

if (vpx_config("CONFIG_EXT_PARTITION") eq "yes") {
  @block_widths = (4, 8, 16, 32, 64, 128)
} else {
  @block_widths = (4, 8, 16, 32, 64)
}

@block_sizes = ();
foreach $w (@block_widths) {
  foreach $h (@block_widths) {
    push @block_sizes, [$w, $h] if ($w <= 2*$h && $h <= 2*$w) ;
  }
}

#
# Intra prediction
#

add_proto qw/void vpx_d207_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207_predictor_4x4/, "$ssse3_x86inc";

add_proto qw/void vpx_d207e_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207e_predictor_4x4/;

add_proto qw/void vpx_d45_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45_predictor_4x4 neon/, "$ssse3_x86inc";

add_proto qw/void vpx_d45e_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45e_predictor_4x4/;

add_proto qw/void vpx_d63_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63_predictor_4x4/, "$ssse3_x86inc";

add_proto qw/void vpx_d63e_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63e_predictor_4x4/;

add_proto qw/void vpx_d63f_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63f_predictor_4x4/;

add_proto qw/void vpx_h_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_h_predictor_4x4 neon dspr2 msa/, "$sse2_x86inc";

add_proto qw/void vpx_he_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_he_predictor_4x4/;

add_proto qw/void vpx_d117_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d117_predictor_4x4/;

add_proto qw/void vpx_d135_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d135_predictor_4x4 neon/;

add_proto qw/void vpx_d153_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d153_predictor_4x4/, "$ssse3_x86inc";

add_proto qw/void vpx_v_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_v_predictor_4x4 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_ve_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_ve_predictor_4x4/;

add_proto qw/void vpx_tm_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_tm_predictor_4x4 neon dspr2 msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_predictor_4x4 dspr2 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_top_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_top_predictor_4x4 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_left_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_left_predictor_4x4 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_128_predictor_4x4/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_128_predictor_4x4 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_d207_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207_predictor_8x8/, "$ssse3_x86inc";

add_proto qw/void vpx_d207e_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207e_predictor_8x8/;

add_proto qw/void vpx_d45_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45_predictor_8x8 neon/, "$ssse3_x86inc";

add_proto qw/void vpx_d45e_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45e_predictor_8x8/;

add_proto qw/void vpx_d63_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63_predictor_8x8/, "$ssse3_x86inc";

add_proto qw/void vpx_d63e_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63e_predictor_8x8/;

add_proto qw/void vpx_h_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_h_predictor_8x8 neon dspr2 msa/, "$sse2_x86inc";

add_proto qw/void vpx_d117_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d117_predictor_8x8/;

add_proto qw/void vpx_d135_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d135_predictor_8x8/;

add_proto qw/void vpx_d153_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d153_predictor_8x8/, "$ssse3_x86inc";

add_proto qw/void vpx_v_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_v_predictor_8x8 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_tm_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_tm_predictor_8x8 neon dspr2 msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_predictor_8x8 dspr2 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_top_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_top_predictor_8x8 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_left_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_left_predictor_8x8 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_128_predictor_8x8/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_128_predictor_8x8 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_d207_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207_predictor_16x16/, "$ssse3_x86inc";

add_proto qw/void vpx_d207e_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207e_predictor_16x16/;

add_proto qw/void vpx_d45_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45_predictor_16x16 neon/, "$ssse3_x86inc";

add_proto qw/void vpx_d45e_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45e_predictor_16x16/;

add_proto qw/void vpx_d63_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63_predictor_16x16/, "$ssse3_x86inc";

add_proto qw/void vpx_d63e_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63e_predictor_16x16/;

add_proto qw/void vpx_h_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_h_predictor_16x16 neon dspr2 msa/, "$sse2_x86inc";

add_proto qw/void vpx_d117_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d117_predictor_16x16/;

add_proto qw/void vpx_d135_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d135_predictor_16x16/;

add_proto qw/void vpx_d153_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d153_predictor_16x16/, "$ssse3_x86inc";

add_proto qw/void vpx_v_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_v_predictor_16x16 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_tm_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_tm_predictor_16x16 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_predictor_16x16 dspr2 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_top_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_top_predictor_16x16 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_left_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_left_predictor_16x16 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_128_predictor_16x16/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_128_predictor_16x16 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_d207_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207_predictor_32x32/, "$ssse3_x86inc";

add_proto qw/void vpx_d207e_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d207e_predictor_32x32/;

add_proto qw/void vpx_d45_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45_predictor_32x32/, "$ssse3_x86inc";

add_proto qw/void vpx_d45e_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d45e_predictor_32x32/;

add_proto qw/void vpx_d63_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63_predictor_32x32/, "$ssse3_x86inc";

add_proto qw/void vpx_d63e_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d63e_predictor_32x32/;

add_proto qw/void vpx_h_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_h_predictor_32x32 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_d117_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d117_predictor_32x32/;

add_proto qw/void vpx_d135_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d135_predictor_32x32/;

add_proto qw/void vpx_d153_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_d153_predictor_32x32/, "$ssse3_x86inc";

add_proto qw/void vpx_v_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_v_predictor_32x32 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_tm_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_tm_predictor_32x32 neon msa/, "$sse2_x86inc";

add_proto qw/void vpx_dc_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_predictor_32x32 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_top_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_top_predictor_32x32 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_left_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_left_predictor_32x32 msa neon/, "$sse2_x86inc";

add_proto qw/void vpx_dc_128_predictor_32x32/, "uint8_t *dst, ptrdiff_t y_stride, const uint8_t *above, const uint8_t *left";
specialize qw/vpx_dc_128_predictor_32x32 msa neon/, "$sse2_x86inc";

# High bitdepth functions
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_d207_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207_predictor_4x4/;

  add_proto qw/void vpx_highbd_d207e_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207e_predictor_4x4/;

  add_proto qw/void vpx_highbd_d45_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45_predictor_4x4/;

  add_proto qw/void vpx_highbd_d45e_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45e_predictor_4x4/;

  add_proto qw/void vpx_highbd_d63_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63_predictor_4x4/;

  add_proto qw/void vpx_highbd_d63e_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63e_predictor_4x4/;

  add_proto qw/void vpx_highbd_h_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_h_predictor_4x4/;

  add_proto qw/void vpx_highbd_d117_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d117_predictor_4x4/;

  add_proto qw/void vpx_highbd_d135_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d135_predictor_4x4/;

  add_proto qw/void vpx_highbd_d153_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d153_predictor_4x4/;

  add_proto qw/void vpx_highbd_v_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_v_predictor_4x4/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_tm_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_tm_predictor_4x4/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_predictor_4x4/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_top_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_top_predictor_4x4/;

  add_proto qw/void vpx_highbd_dc_left_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_left_predictor_4x4/;

  add_proto qw/void vpx_highbd_dc_128_predictor_4x4/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_128_predictor_4x4/;

  add_proto qw/void vpx_highbd_d207_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207_predictor_8x8/;

  add_proto qw/void vpx_highbd_d207e_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207e_predictor_8x8/;

  add_proto qw/void vpx_highbd_d45_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45_predictor_8x8/;

  add_proto qw/void vpx_highbd_d45e_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45e_predictor_8x8/;

  add_proto qw/void vpx_highbd_d63_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63_predictor_8x8/;

  add_proto qw/void vpx_highbd_d63e_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63e_predictor_8x8/;

  add_proto qw/void vpx_highbd_h_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_h_predictor_8x8/;

  add_proto qw/void vpx_highbd_d117_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d117_predictor_8x8/;

  add_proto qw/void vpx_highbd_d135_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d135_predictor_8x8/;

  add_proto qw/void vpx_highbd_d153_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d153_predictor_8x8/;

  add_proto qw/void vpx_highbd_v_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_v_predictor_8x8/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_tm_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_tm_predictor_8x8/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_predictor_8x8/, "$sse2_x86inc";;

  add_proto qw/void vpx_highbd_dc_top_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_top_predictor_8x8/;

  add_proto qw/void vpx_highbd_dc_left_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_left_predictor_8x8/;

  add_proto qw/void vpx_highbd_dc_128_predictor_8x8/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_128_predictor_8x8/;

  add_proto qw/void vpx_highbd_d207_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207_predictor_16x16/;

  add_proto qw/void vpx_highbd_d207e_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207e_predictor_16x16/;

  add_proto qw/void vpx_highbd_d45_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45_predictor_16x16/;

  add_proto qw/void vpx_highbd_d45e_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45e_predictor_16x16/;

  add_proto qw/void vpx_highbd_d63_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63_predictor_16x16/;

  add_proto qw/void vpx_highbd_d63e_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63e_predictor_16x16/;

  add_proto qw/void vpx_highbd_h_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_h_predictor_16x16/;

  add_proto qw/void vpx_highbd_d117_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d117_predictor_16x16/;

  add_proto qw/void vpx_highbd_d135_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d135_predictor_16x16/;

  add_proto qw/void vpx_highbd_d153_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d153_predictor_16x16/;

  add_proto qw/void vpx_highbd_v_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_v_predictor_16x16/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_tm_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_tm_predictor_16x16/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_predictor_16x16/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_top_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_top_predictor_16x16/;

  add_proto qw/void vpx_highbd_dc_left_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_left_predictor_16x16/;

  add_proto qw/void vpx_highbd_dc_128_predictor_16x16/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_128_predictor_16x16/;

  add_proto qw/void vpx_highbd_d207_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207_predictor_32x32/;

  add_proto qw/void vpx_highbd_d207e_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d207e_predictor_32x32/;

  add_proto qw/void vpx_highbd_d45_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45_predictor_32x32/;

  add_proto qw/void vpx_highbd_d45e_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d45e_predictor_32x32/;

  add_proto qw/void vpx_highbd_d63_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63_predictor_32x32/;

  add_proto qw/void vpx_highbd_d63e_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d63e_predictor_32x32/;

  add_proto qw/void vpx_highbd_h_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_h_predictor_32x32/;

  add_proto qw/void vpx_highbd_d117_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d117_predictor_32x32/;

  add_proto qw/void vpx_highbd_d135_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d135_predictor_32x32/;

  add_proto qw/void vpx_highbd_d153_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_d153_predictor_32x32/;

  add_proto qw/void vpx_highbd_v_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_v_predictor_32x32/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_tm_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_tm_predictor_32x32/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_predictor_32x32/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_dc_top_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_top_predictor_32x32/;

  add_proto qw/void vpx_highbd_dc_left_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_left_predictor_32x32/;

  add_proto qw/void vpx_highbd_dc_128_predictor_32x32/, "uint16_t *dst, ptrdiff_t y_stride, const uint16_t *above, const uint16_t *left, int bd";
  specialize qw/vpx_highbd_dc_128_predictor_32x32/;
}  # CONFIG_VP9_HIGHBITDEPTH

#
# Sub Pixel Filters
#
add_proto qw/void vpx_convolve_copy/,       "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve_avg/,        "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8/,           "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8_horiz/,     "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8_vert/,      "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8_avg/,       "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8_avg_horiz/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_convolve8_avg_vert/,  "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_2d/,           "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_horiz/,        "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_vert/,         "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_avg_2d/,       "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_avg_horiz/,    "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";
add_proto qw/void vpx_scaled_avg_vert/,     "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h";

specialize qw/vpx_convolve_copy                 /, "$sse2_x86inc";
specialize qw/vpx_convolve_avg                  /, "$sse2_x86inc";
specialize qw/vpx_convolve8           sse2 ssse3/, "$avx2_ssse3";
specialize qw/vpx_convolve8_horiz     sse2 ssse3/, "$avx2_ssse3";
specialize qw/vpx_convolve8_vert      sse2 ssse3/, "$avx2_ssse3";
specialize qw/vpx_convolve8_avg       sse2 ssse3/;
specialize qw/vpx_convolve8_avg_horiz sse2 ssse3/;
specialize qw/vpx_convolve8_avg_vert  sse2 ssse3/;
specialize qw/vpx_scaled_2d                ssse3/;

# TODO(any): These need to be extended to up to 128x128 block sizes
if (!(vpx_config("CONFIG_VP10") eq "yes" && vpx_config("CONFIG_EXT_PARTITION") eq "yes")) {
  specialize qw/vpx_convolve_copy       neon dspr2 msa/;
  specialize qw/vpx_convolve_avg        neon dspr2 msa/;
  specialize qw/vpx_convolve8           neon dspr2 msa/;
  specialize qw/vpx_convolve8_horiz     neon dspr2 msa/;
  specialize qw/vpx_convolve8_vert      neon dspr2 msa/;
  specialize qw/vpx_convolve8_avg       neon dspr2 msa/;
  specialize qw/vpx_convolve8_avg_horiz neon dspr2 msa/;
  specialize qw/vpx_convolve8_avg_vert  neon dspr2 msa/;
}

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_convolve_copy/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve_copy/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_convolve_avg/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve_avg/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_convolve8/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8/, "$sse2_x86_64";

  add_proto qw/void vpx_highbd_convolve8_horiz/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8_horiz/, "$sse2_x86_64";

  add_proto qw/void vpx_highbd_convolve8_vert/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8_vert/, "$sse2_x86_64";

  add_proto qw/void vpx_highbd_convolve8_avg/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8_avg/, "$sse2_x86_64";

  add_proto qw/void vpx_highbd_convolve8_avg_horiz/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8_avg_horiz/, "$sse2_x86_64";

  add_proto qw/void vpx_highbd_convolve8_avg_vert/, "const uint8_t *src, ptrdiff_t src_stride, uint8_t *dst, ptrdiff_t dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h, int bps";
  specialize qw/vpx_highbd_convolve8_avg_vert/, "$sse2_x86_64";
}  # CONFIG_VP9_HIGHBITDEPTH

#
# Loopfilter
#
add_proto qw/void vpx_lpf_vertical_16/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_16 sse2 neon_asm dspr2 msa/;
$vpx_lpf_vertical_16_neon_asm=vpx_lpf_vertical_16_neon;

add_proto qw/void vpx_lpf_vertical_16_dual/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_16_dual sse2 neon_asm dspr2 msa/;
$vpx_lpf_vertical_16_dual_neon_asm=vpx_lpf_vertical_16_dual_neon;

add_proto qw/void vpx_lpf_vertical_8/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_8 sse2 neon dspr2 msa/;

add_proto qw/void vpx_lpf_vertical_8_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_vertical_8_dual sse2 neon_asm dspr2 msa/;
$vpx_lpf_vertical_8_dual_neon_asm=vpx_lpf_vertical_8_dual_neon;

add_proto qw/void vpx_lpf_vertical_4/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_vertical_4 neon dspr2 msa/, "$mmx_x86inc";

add_proto qw/void vpx_lpf_vertical_4_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_vertical_4_dual sse2 neon dspr2 msa/;

add_proto qw/void vpx_lpf_horizontal_edge_8/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_horizontal_edge_8 sse2 avx2 neon_asm dspr2 msa/;
$vpx_lpf_horizontal_edge_8_neon_asm=vpx_lpf_horizontal_edge_8_neon;

add_proto qw/void vpx_lpf_horizontal_edge_16/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_horizontal_edge_16 sse2 avx2 neon_asm dspr2 msa/;
$vpx_lpf_horizontal_edge_16_neon_asm=vpx_lpf_horizontal_edge_16_neon;

add_proto qw/void vpx_lpf_horizontal_8/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_horizontal_8 sse2 neon dspr2 msa/;

add_proto qw/void vpx_lpf_horizontal_8_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_horizontal_8_dual sse2 neon_asm dspr2 msa/;
$vpx_lpf_horizontal_8_dual_neon_asm=vpx_lpf_horizontal_8_dual_neon;

add_proto qw/void vpx_lpf_horizontal_4/, "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh";
specialize qw/vpx_lpf_horizontal_4 neon dspr2 msa/, "$mmx_x86inc";

add_proto qw/void vpx_lpf_horizontal_4_dual/, "uint8_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1";
specialize qw/vpx_lpf_horizontal_4_dual sse2 neon dspr2 msa/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_lpf_vertical_16/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_16 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_16_dual/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_16_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_8/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_8 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_8_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_vertical_8_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_4/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_vertical_4 sse2/;

  add_proto qw/void vpx_highbd_lpf_vertical_4_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_vertical_4_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_edge_8/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_edge_8 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_edge_16/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_edge_16 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_8/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_8 sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_8_dual/, "uint16_t *s, int pitch, const uint8_t *blimit0, const uint8_t *limit0, const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1, const uint8_t *thresh1, int bd";
  specialize qw/vpx_highbd_lpf_horizontal_8_dual sse2/;

  add_proto qw/void vpx_highbd_lpf_horizontal_4/, "uint16_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int bd";
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
if ((vpx_config("CONFIG_VP9_ENCODER") eq "yes") || (vpx_config("CONFIG_VP10_ENCODER") eq "yes")) {
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct4x4 sse2/;

  add_proto qw/void vpx_fdct4x4_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct4x4_1 sse2/;

  add_proto qw/void vpx_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct8x8 sse2/;

  add_proto qw/void vpx_fdct8x8_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct8x8_1 sse2/;

  add_proto qw/void vpx_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct16x16 sse2/;

  add_proto qw/void vpx_fdct16x16_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct16x16_1 sse2/;

  add_proto qw/void vpx_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32 sse2/;

  add_proto qw/void vpx_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32_rd sse2/;

  add_proto qw/void vpx_fdct32x32_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32_1 sse2/;

  add_proto qw/void vpx_highbd_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct4x4 sse2/;

  add_proto qw/void vpx_highbd_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct8x8 sse2/;

  add_proto qw/void vpx_highbd_fdct8x8_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct8x8_1/;

  add_proto qw/void vpx_highbd_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct16x16 sse2/;

  add_proto qw/void vpx_highbd_fdct16x16_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct16x16_1/;

  add_proto qw/void vpx_highbd_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct32x32 sse2/;

  add_proto qw/void vpx_highbd_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct32x32_rd sse2/;

  add_proto qw/void vpx_highbd_fdct32x32_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_highbd_fdct32x32_1/;
} else {
  add_proto qw/void vpx_fdct4x4/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct4x4 sse2 msa/;

  add_proto qw/void vpx_fdct4x4_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct4x4_1 sse2/;

  add_proto qw/void vpx_fdct8x8/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct8x8 sse2 neon msa/, "$ssse3_x86_64_x86inc";

  add_proto qw/void vpx_fdct8x8_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct8x8_1 sse2 neon msa/;

  add_proto qw/void vpx_fdct16x16/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct16x16 sse2 msa/;

  add_proto qw/void vpx_fdct16x16_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct16x16_1 sse2 msa/;

  add_proto qw/void vpx_fdct32x32/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32 sse2 avx2 msa/;

  add_proto qw/void vpx_fdct32x32_rd/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32_rd sse2 avx2 msa/;

  add_proto qw/void vpx_fdct32x32_1/, "const int16_t *input, tran_low_t *output, int stride";
  specialize qw/vpx_fdct32x32_1 sse2 msa/;
}  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_VP9_ENCODER || CONFIG_VP10_ENCODER

#
# Inverse transform
if ((vpx_config("CONFIG_VP9") eq "yes") || (vpx_config("CONFIG_VP10") eq "yes")) {
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  # Note as optimized versions of these functions are added we need to add a check to ensure
  # that when CONFIG_EMULATE_HARDWARE is on, it defaults to the C versions only.
  add_proto qw/void vpx_iwht4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
  specialize qw/vpx_iwht4x4_1_add/;

  add_proto qw/void vpx_iwht4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
  specialize qw/vpx_iwht4x4_16_add/, "$sse2_x86inc";

  add_proto qw/void vpx_highbd_idct4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct4x4_1_add/;

  add_proto qw/void vpx_highbd_idct8x8_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct8x8_1_add/;

  add_proto qw/void vpx_highbd_idct16x16_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct16x16_1_add/;

  add_proto qw/void vpx_highbd_idct32x32_1024_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct32x32_1024_add/;

  add_proto qw/void vpx_highbd_idct32x32_34_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct32x32_34_add/;

  add_proto qw/void vpx_highbd_idct32x32_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_idct32x32_1_add/;

  add_proto qw/void vpx_highbd_iwht4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_iwht4x4_1_add/;

  add_proto qw/void vpx_highbd_iwht4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
  specialize qw/vpx_highbd_iwht4x4_16_add/;

  # Force C versions if CONFIG_EMULATE_HARDWARE is 1
  if (vpx_config("CONFIG_EMULATE_HARDWARE") eq "yes") {
    add_proto qw/void vpx_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_16_add/;

    add_proto qw/void vpx_idct4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_1_add/;

    add_proto qw/void vpx_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_64_add/;

    add_proto qw/void vpx_idct8x8_12_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_12_add/;

    add_proto qw/void vpx_idct8x8_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_1_add/;

    add_proto qw/void vpx_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_256_add/;

    add_proto qw/void vpx_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_10_add/;

    add_proto qw/void vpx_idct16x16_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_1_add/;

    add_proto qw/void vpx_idct32x32_1024_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1024_add/;

    add_proto qw/void vpx_idct32x32_135_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_135_add/;

    add_proto qw/void vpx_idct32x32_34_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_34_add/;

    add_proto qw/void vpx_idct32x32_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1_add/;

    add_proto qw/void vpx_highbd_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct4x4_16_add/;

    add_proto qw/void vpx_highbd_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct8x8_64_add/;

    add_proto qw/void vpx_highbd_idct8x8_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct8x8_10_add/;

    add_proto qw/void vpx_highbd_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct16x16_256_add/;

    add_proto qw/void vpx_highbd_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct16x16_10_add/;
  } else {
    add_proto qw/void vpx_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_16_add sse2/;

    add_proto qw/void vpx_idct4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_1_add sse2/;

    add_proto qw/void vpx_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_64_add sse2/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct8x8_12_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_12_add sse2/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct8x8_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_1_add sse2/;

    add_proto qw/void vpx_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_256_add sse2/;

    add_proto qw/void vpx_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_10_add sse2/;

    add_proto qw/void vpx_idct16x16_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_1_add sse2/;

    add_proto qw/void vpx_idct32x32_1024_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1024_add sse2/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct32x32_135_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_135_add sse2/, "$ssse3_x86_64_x86inc";
    # Need to add 135 eob idct32x32 implementations.
    $vpx_idct32x32_135_add_sse2=vpx_idct32x32_1024_add_sse2;

    add_proto qw/void vpx_idct32x32_34_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_34_add sse2/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct32x32_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1_add sse2/;

    add_proto qw/void vpx_highbd_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct4x4_16_add sse2/;

    add_proto qw/void vpx_highbd_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct8x8_64_add sse2/;

    add_proto qw/void vpx_highbd_idct8x8_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct8x8_10_add sse2/;

    add_proto qw/void vpx_highbd_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct16x16_256_add sse2/;

    add_proto qw/void vpx_highbd_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride, int bd";
    specialize qw/vpx_highbd_idct16x16_10_add sse2/;
  }  # CONFIG_EMULATE_HARDWARE
} else {
  # Force C versions if CONFIG_EMULATE_HARDWARE is 1
  if (vpx_config("CONFIG_EMULATE_HARDWARE") eq "yes") {
    add_proto qw/void vpx_idct4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_1_add/;

    add_proto qw/void vpx_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_16_add/;

    add_proto qw/void vpx_idct8x8_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_1_add/;

    add_proto qw/void vpx_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_64_add/;

    add_proto qw/void vpx_idct8x8_12_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_12_add/;

    add_proto qw/void vpx_idct16x16_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_1_add/;

    add_proto qw/void vpx_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_256_add/;

    add_proto qw/void vpx_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_10_add/;

    add_proto qw/void vpx_idct32x32_1024_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1024_add/;

    add_proto qw/void vpx_idct32x32_135_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_135_add/;

    add_proto qw/void vpx_idct32x32_34_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_34_add/;

    add_proto qw/void vpx_idct32x32_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1_add/;

    add_proto qw/void vpx_iwht4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_iwht4x4_1_add/;

    add_proto qw/void vpx_iwht4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_iwht4x4_16_add/;
  } else {
    add_proto qw/void vpx_idct4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_1_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct4x4_16_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct8x8_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_1_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct8x8_64_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_64_add sse2 neon dspr2 msa/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct8x8_12_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct8x8_12_add sse2 neon dspr2 msa/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct16x16_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_1_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct16x16_256_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_256_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct16x16_10_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct16x16_10_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_idct32x32_1024_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1024_add sse2 neon dspr2 msa/, "$ssse3_x86_64_x86inc";

    add_proto qw/void vpx_idct32x32_135_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_135_add sse2 neon dspr2 msa/, "$ssse3_x86_64_x86inc";
    # Need to add 135 eob idct32x32 implementations.
    $vpx_idct32x32_135_add_sse2=vpx_idct32x32_1024_add_sse2;
    $vpx_idct32x32_135_add_neon=vpx_idct32x32_1024_add_neon;
    $vpx_idct32x32_135_add_dspr2=vpx_idct32x32_1024_add_dspr2;
    $vpx_idct32x32_135_add_msa=vpx_idct32x32_1024_add_msa;

    add_proto qw/void vpx_idct32x32_34_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_34_add sse2 neon dspr2 msa/, "$ssse3_x86_64_x86inc";
    # Need to add 34 eob idct32x32 neon implementation.
    $vpx_idct32x32_34_add_neon=vpx_idct32x32_1024_add_neon;

    add_proto qw/void vpx_idct32x32_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_idct32x32_1_add sse2 neon dspr2 msa/;

    add_proto qw/void vpx_iwht4x4_1_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_iwht4x4_1_add msa/;

    add_proto qw/void vpx_iwht4x4_16_add/, "const tran_low_t *input, uint8_t *dest, int dest_stride";
    specialize qw/vpx_iwht4x4_16_add msa/, "$sse2_x86inc";
  }  # CONFIG_EMULATE_HARDWARE
}  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_VP9 || CONFIG_VP10

#
# Quantization
#
if ((vpx_config("CONFIG_VP9_ENCODER") eq "yes") || (vpx_config("CONFIG_VP10_ENCODER") eq "yes")) {
  add_proto qw/void vpx_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vpx_quantize_b sse2/, "$ssse3_x86_64_x86inc", "$avx_x86_64_x86inc";

  add_proto qw/void vpx_quantize_b_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
  specialize qw/vpx_quantize_b_32x32/, "$ssse3_x86_64_x86inc", "$avx_x86_64_x86inc";

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    add_proto qw/void vpx_highbd_quantize_b/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    specialize qw/vpx_highbd_quantize_b sse2/;

    add_proto qw/void vpx_highbd_quantize_b_32x32/, "const tran_low_t *coeff_ptr, intptr_t n_coeffs, int skip_block, const int16_t *zbin_ptr, const int16_t *round_ptr, const int16_t *quant_ptr, const int16_t *quant_shift_ptr, tran_low_t *qcoeff_ptr, tran_low_t *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan";
    specialize qw/vpx_highbd_quantize_b_32x32 sse2/;
  }  # CONFIG_VP9_HIGHBITDEPTH
}  # CONFIG_VP9_ENCODER || CONFIG_VP10_ENCODER

if (vpx_config("CONFIG_VP10") eq "yes") {
  #
  # Alpha blending with mask
  #
  add_proto qw/void vpx_blend_a64_mask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int h, int w, int suby, int subx";
  add_proto qw/void vpx_blend_a64_hmask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int h, int w";
  add_proto qw/void vpx_blend_a64_vmask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int h, int w";
  specialize "vpx_blend_a64_mask", qw/sse4_1/;
  specialize "vpx_blend_a64_hmask", qw/sse4_1/;
  specialize "vpx_blend_a64_vmask", qw/sse4_1/;

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    add_proto qw/void vpx_highbd_blend_a64_mask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride, int h, int w, int suby, int subx, int bd";
    add_proto qw/void vpx_highbd_blend_a64_hmask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int h, int w, int bd";
    add_proto qw/void vpx_highbd_blend_a64_vmask/, "uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int h, int w, int bd";
    specialize "vpx_highbd_blend_a64_mask", qw/sse4_1/;
    specialize "vpx_highbd_blend_a64_hmask", qw/sse4_1/;
    specialize "vpx_highbd_blend_a64_vmask", qw/sse4_1/;
  }
}  # CONFIG_VP10

if (vpx_config("CONFIG_ENCODERS") eq "yes") {
#
# Block subtraction
#
add_proto qw/void vpx_subtract_block/, "int rows, int cols, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr, ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride";
specialize qw/vpx_subtract_block neon msa/, "$sse2_x86inc";

if (vpx_config("CONFIG_VP10_ENCODER") eq "yes") {
  #
  # Sum of Squares
  #
  add_proto qw/uint64_t vpx_sum_squares_2d_i16/, "const int16_t *src, int stride, int size";
  specialize qw/vpx_sum_squares_2d_i16 sse2/;

  add_proto qw/uint64_t vpx_sum_squares_i16/, "const int16_t *src, uint32_t N";
  specialize qw/vpx_sum_squares_i16 sse2/;
}

if ((vpx_config("CONFIG_VP9_ENCODER") eq "yes") || (vpx_config("CONFIG_VP10_ENCODER") eq "yes")) {
  #
  # Avg
  #
  add_proto qw/unsigned int vpx_avg_8x8/, "const uint8_t *, int p";
  specialize qw/vpx_avg_8x8 sse2 neon msa/;
  add_proto qw/unsigned int vpx_avg_4x4/, "const uint8_t *, int p";
  specialize qw/vpx_avg_4x4 sse2 neon msa/;
  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    add_proto qw/unsigned int vpx_highbd_avg_8x8/, "const uint8_t *, int p";
    specialize qw/vpx_highbd_avg_8x8/;
    add_proto qw/unsigned int vpx_highbd_avg_4x4/, "const uint8_t *, int p";
    specialize qw/vpx_highbd_avg_4x4/;
    add_proto qw/void vpx_highbd_subtract_block/, "int rows, int cols, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr, ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride, int bd";
    specialize qw/vpx_highbd_subtract_block sse2/;
  }

  #
  # Minmax
  #
  add_proto qw/void vpx_minmax_8x8/, "const uint8_t *s, int p, const uint8_t *d, int dp, int *min, int *max";
  specialize qw/vpx_minmax_8x8 sse2/;
  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    add_proto qw/void vpx_highbd_minmax_8x8/, "const uint8_t *s, int p, const uint8_t *d, int dp, int *min, int *max";
    specialize qw/vpx_highbd_minmax_8x8/;
  }

  add_proto qw/void vpx_hadamard_8x8/, "const int16_t *src_diff, int src_stride, int16_t *coeff";
  specialize qw/vpx_hadamard_8x8 sse2/, "$ssse3_x86_64_x86inc";

  add_proto qw/void vpx_hadamard_16x16/, "const int16_t *src_diff, int src_stride, int16_t *coeff";
  specialize qw/vpx_hadamard_16x16 sse2/;

  add_proto qw/int vpx_satd/, "const int16_t *coeff, int length";
  specialize qw/vpx_satd sse2 neon/;

  add_proto qw/void vpx_int_pro_row/, "int16_t *hbuf, const uint8_t *ref, const int ref_stride, const int height";
  specialize qw/vpx_int_pro_row sse2 neon/;

  add_proto qw/int16_t vpx_int_pro_col/, "const uint8_t *ref, const int width";
  specialize qw/vpx_int_pro_col sse2 neon/;

  add_proto qw/int vpx_vector_var/, "const int16_t *ref, const int16_t *src, const int bwl";
  specialize qw/vpx_vector_var neon sse2/;
}  # CONFIG_VP9_ENCODER || CONFIG_VP10_ENCODER

#
# Single block SAD / Single block Avg SAD
#
foreach (@block_sizes) {
  ($w, $h) = @$_;
  add_proto qw/unsigned int/, "vpx_sad${w}x${h}", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
  add_proto qw/unsigned int/, "vpx_sad${w}x${h}_avg", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
}

specialize qw/vpx_sad128x128                       /, "$sse2_x86inc";
specialize qw/vpx_sad128x64                        /, "$sse2_x86inc";
specialize qw/vpx_sad64x128                        /, "$sse2_x86inc";
specialize qw/vpx_sad64x64      avx2       neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad64x32      avx2            msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x64      avx2            msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x32      avx2       neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x16      avx2            msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x32                      msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x16   mmx     media neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x8    mmx           neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x16    mmx           neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x8     mmx           neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x4                        msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x8                        msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x4     mmx           neon msa/, "$sse2_x86inc";

specialize qw/vpx_sad128x128_avg         /, "$sse2_x86inc";
specialize qw/vpx_sad128x64_avg          /, "$sse2_x86inc";
specialize qw/vpx_sad64x128_avg          /, "$sse2_x86inc";
specialize qw/vpx_sad64x64_avg   avx2 msa/, "$sse2_x86inc";
specialize qw/vpx_sad64x32_avg   avx2 msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x64_avg   avx2 msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x32_avg   avx2 msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x16_avg   avx2 msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x32_avg        msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x16_avg        msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x8_avg         msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x16_avg         msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x8_avg          msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x4_avg          msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x8_avg          msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x4_avg          msa/, "$sse2_x86inc";

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/unsigned int/, "vpx_highbd_sad${w}x${h}", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride";
    add_proto qw/unsigned int/, "vpx_highbd_sad${w}x${h}_avg", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *second_pred";
    if ($w != 128 && $h != 128 && $w != 4) {
      specialize "vpx_highbd_sad${w}x${h}", "$sse2_x86inc";
      specialize "vpx_highbd_sad${w}x${h}_avg", "$sse2_x86inc";
    }
  }
}

#
# Masked SAD
#
if (vpx_config("CONFIG_EXT_INTER") eq "yes") {
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/unsigned int/, "vpx_masked_sad${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *mask, int mask_stride";
    specialize "vpx_masked_sad${w}x${h}", qw/ssse3/;
  }

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    foreach (@block_sizes) {
      ($w, $h) = @$_;
      add_proto qw/unsigned int/, "vpx_highbd_masked_sad${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *mask, int mask_stride";
      specialize "vpx_highbd_masked_sad${w}x${h}", qw/ssse3/;
    }
  }
}

#
# OBMC SAD
#
if (vpx_config("CONFIG_OBMC") eq "yes") {
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/unsigned int/, "vpx_obmc_sad${w}x${h}", "const uint8_t *ref_ptr, int ref_stride, const int32_t *wsrc_ptr, const int32_t *mask";
    specialize "vpx_obmc_sad${w}x${h}", qw/sse4_1/;
  }

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    foreach (@block_sizes) {
      ($w, $h) = @$_;
      add_proto qw/unsigned int/, "vpx_highbd_obmc_sad${w}x${h}", "const uint8_t *ref_ptr, int ref_stride, const int32_t *wsrc_ptr, const int32_t *mask";
      specialize "vpx_highbd_obmc_sad${w}x${h}", qw/sse4_1/;
    }
  }
}

#
# Multi-block SAD, comparing a reference to N blocks 1 pixel apart horizontally
#
# Blocks of 3
foreach $s (@block_widths) {
  add_proto qw/void/, "vpx_sad${s}x${s}x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
}
specialize qw/vpx_sad64x64x3            msa/;
specialize qw/vpx_sad32x32x3            msa/;
specialize qw/vpx_sad16x16x3 sse3 ssse3 msa/;
specialize qw/vpx_sad8x8x3   sse3       msa/;
specialize qw/vpx_sad4x4x3   sse3       msa/;

add_proto qw/void/, "vpx_sad16x8x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x8x3 sse3 ssse3 msa/;
add_proto qw/void/, "vpx_sad8x16x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x16x3 sse3 msa/;

# Blocks of 8
foreach $s (@block_widths) {
  add_proto qw/void/, "vpx_sad${s}x${s}x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
}
specialize qw/vpx_sad64x64x8        msa/;
specialize qw/vpx_sad32x32x8        msa/;
specialize qw/vpx_sad16x16x8 sse4_1 msa/;
specialize qw/vpx_sad8x8x8   sse4_1 msa/;
specialize qw/vpx_sad4x4x8   sse4_1 msa/;

add_proto qw/void/, "vpx_sad16x8x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad16x8x8 sse4_1 msa/;
add_proto qw/void/, "vpx_sad8x16x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x16x8 sse4_1 msa/;
add_proto qw/void/, "vpx_sad8x4x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad8x4x8 msa/;
add_proto qw/void/, "vpx_sad4x8x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
specialize qw/vpx_sad4x8x8 msa/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  foreach $s (@block_widths) {
    # Blocks of 3
    add_proto qw/void/, "vpx_highbd_sad${s}x${s}x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
    # Blocks of 8
    add_proto qw/void/, "vpx_highbd_sad${s}x${s}x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  }
  # Blocks of 3
  add_proto qw/void/, "vpx_highbd_sad16x8x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  add_proto qw/void/, "vpx_highbd_sad8x16x3", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  # Blocks of 8
  add_proto qw/void/, "vpx_highbd_sad16x8x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  add_proto qw/void/, "vpx_highbd_sad8x16x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  add_proto qw/void/, "vpx_highbd_sad8x4x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
  add_proto qw/void/, "vpx_highbd_sad4x8x8", "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array";
}

#
# Multi-block SAD, comparing a reference to N independent blocks
#
foreach (@block_sizes) {
  ($w, $h) = @$_;
  add_proto qw/void/, "vpx_sad${w}x${h}x4d", "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
}

specialize qw/vpx_sad128x128x4d              /, "$sse2_x86inc";
specialize qw/vpx_sad128x64x4d               /, "$sse2_x86inc";
specialize qw/vpx_sad64x128x4d               /, "$sse2_x86inc";
specialize qw/vpx_sad64x64x4d   avx2 neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad64x32x4d             msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x64x4d             msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x32x4d   avx2 neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad32x16x4d             msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x32x4d             msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x16x4d        neon msa/, "$sse2_x86inc";
specialize qw/vpx_sad16x8x4d              msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x16x4d              msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x8x4d               msa/, "$sse2_x86inc";
specialize qw/vpx_sad8x4x4d               msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x8x4d               msa/, "$sse2_x86inc";
specialize qw/vpx_sad4x4x4d               msa/, "$sse2_x86inc";

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  #
  # Multi-block SAD, comparing a reference to N independent blocks
  #
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/void/, "vpx_highbd_sad${w}x${h}x4d", "const uint8_t *src_ptr, int src_stride, const uint8_t * const ref_ptr[], int ref_stride, uint32_t *sad_array";
    if ($w != 128 && $h != 128) {
      specialize "vpx_highbd_sad${w}x${h}x4d", "$sse2_x86inc";
    }
  }
}

#
# Structured Similarity (SSIM)
#
if (vpx_config("CONFIG_INTERNAL_STATS") eq "yes") {
  add_proto qw/void vpx_ssim_parms_8x8/, "const uint8_t *s, int sp, const uint8_t *r, int rp, uint32_t *sum_s, uint32_t *sum_r, uint32_t *sum_sq_s, uint32_t *sum_sq_r, uint32_t *sum_sxr";
  specialize qw/vpx_ssim_parms_8x8/, "$sse2_x86_64";

  add_proto qw/void vpx_ssim_parms_16x16/, "const uint8_t *s, int sp, const uint8_t *r, int rp, uint32_t *sum_s, uint32_t *sum_r, uint32_t *sum_sq_s, uint32_t *sum_sq_r, uint32_t *sum_sxr";
  specialize qw/vpx_ssim_parms_16x16/, "$sse2_x86_64";

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    add_proto qw/void vpx_highbd_ssim_parms_8x8/, "const uint16_t *s, int sp, const uint16_t *r, int rp, uint32_t *sum_s, uint32_t *sum_r, uint32_t *sum_sq_s, uint32_t *sum_sq_r, uint32_t *sum_sxr";
  }
}
}  # CONFIG_ENCODERS

if (vpx_config("CONFIG_ENCODERS") eq "yes" || vpx_config("CONFIG_POSTPROC") eq "yes" || vpx_config("CONFIG_VP9_POSTPROC") eq "yes") {

#
# Specialty Variance
#
add_proto qw/void vpx_get16x16var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
add_proto qw/void vpx_get8x8var/, "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";

specialize qw/vpx_get16x16var     avx2 sse2 neon msa/;
specialize qw/vpx_get8x8var   mmx      sse2 neon msa/;

add_proto qw/unsigned int vpx_mse16x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
add_proto qw/unsigned int vpx_mse16x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
add_proto qw/unsigned int vpx_mse8x16/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
add_proto qw/unsigned int vpx_mse8x8/, "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";

specialize qw/vpx_mse16x16 mmx avx2 sse2 media neon msa/;
specialize qw/vpx_mse16x8           sse2            msa/;
specialize qw/vpx_mse8x16           sse2            msa/;
specialize qw/vpx_mse8x8            sse2            msa/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  foreach $bd (8, 10, 12) {
    add_proto qw/void/, "vpx_highbd_${bd}_get16x16var", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";
    add_proto qw/void/, "vpx_highbd_${bd}_get8x8var", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum";

    add_proto qw/unsigned int/, "vpx_highbd_${bd}_mse16x16", "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
    add_proto qw/unsigned int/, "vpx_highbd_${bd}_mse16x8", "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
    add_proto qw/unsigned int/, "vpx_highbd_${bd}_mse8x16", "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";
    add_proto qw/unsigned int/, "vpx_highbd_${bd}_mse8x8", "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse";

    specialize "vpx_highbd_${bd}_mse16x16", qw/sse2/;
    specialize "vpx_highbd_${bd}_mse8x8", qw/sse2/;
  }
}

#
# ...
#
add_proto qw/void vpx_upsampled_pred/, "uint8_t *comp_pred, int width, int height, const uint8_t *ref, int ref_stride";
specialize qw/vpx_upsampled_pred sse2/;
add_proto qw/void vpx_comp_avg_upsampled_pred/, "uint8_t *comp_pred, const uint8_t *pred, int width, int height, const uint8_t *ref, int ref_stride";
specialize qw/vpx_comp_avg_upsampled_pred sse2/;

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_upsampled_pred/, "uint16_t *comp_pred, int width, int height, const uint8_t *ref8, int ref_stride";
  specialize qw/vpx_highbd_upsampled_pred sse2/;
  add_proto qw/void vpx_highbd_comp_avg_upsampled_pred/, "uint16_t *comp_pred, const uint8_t *pred8, int width, int height, const uint8_t *ref8, int ref_stride";
  specialize qw/vpx_highbd_comp_avg_upsampled_pred sse2/;
}

#
# ...
#
add_proto qw/unsigned int vpx_get_mb_ss/, "const int16_t *";
add_proto qw/unsigned int vpx_get4x4sse_cs/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride";

specialize qw/vpx_get_mb_ss mmx sse2 msa/;
specialize qw/vpx_get4x4sse_cs neon msa/;

#
# Variance / Subpixel Variance / Subpixel Avg Variance
#
foreach (@block_sizes) {
  ($w, $h) = @$_;
  add_proto qw/unsigned int/, "vpx_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse";
  add_proto qw/uint32_t/, "vpx_sub_pixel_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
  add_proto qw/uint32_t/, "vpx_sub_pixel_avg_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
}

specialize qw/vpx_variance64x64     sse2 avx2       neon msa/;
specialize qw/vpx_variance64x32     sse2 avx2       neon msa/;
specialize qw/vpx_variance32x64     sse2            neon msa/;
specialize qw/vpx_variance32x32     sse2 avx2       neon msa/;
specialize qw/vpx_variance32x16     sse2 avx2            msa/;
specialize qw/vpx_variance16x32     sse2                 msa/;
specialize qw/vpx_variance16x16 mmx sse2 avx2 media neon msa/;
specialize qw/vpx_variance16x8  mmx sse2            neon msa/;
specialize qw/vpx_variance8x16  mmx sse2            neon msa/;
specialize qw/vpx_variance8x8   mmx sse2      media neon msa/;
specialize qw/vpx_variance8x4       sse2                 msa/;
specialize qw/vpx_variance4x8       sse2                 msa/;
specialize qw/vpx_variance4x4   mmx sse2                 msa/;

specialize qw/vpx_sub_pixel_variance64x64     avx2       neon msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance64x32                     msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance32x64                     msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance32x32     avx2       neon msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance32x16                     msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance16x32                     msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance16x16 mmx      media neon msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance16x8  mmx                 msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance8x16  mmx                 msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance8x8   mmx      media neon msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance8x4                       msa/,                 "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance4x8                       msa/, "$sse_x86inc",                  "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_variance4x4   mmx                 msa/, "$sse_x86inc",                  "$ssse3_x86inc";

specialize qw/vpx_sub_pixel_avg_variance64x64 avx2 msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance64x32      msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance32x64      msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance32x32 avx2 msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance32x16      msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance16x32      msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance16x16      msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance16x8       msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance8x16       msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance8x8        msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance8x4        msa/,                "$sse2_x86inc", "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance4x8        msa/, "$sse_x86inc",                 "$ssse3_x86inc";
specialize qw/vpx_sub_pixel_avg_variance4x4        msa/, "$sse_x86inc",                 "$ssse3_x86inc";

if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  foreach $bd (8, 10, 12) {
    foreach (@block_sizes) {
      ($w, $h) = @$_;
      add_proto qw/unsigned int/, "vpx_highbd_${bd}_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
      add_proto qw/uint32_t/, "vpx_highbd_${bd}_sub_pixel_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse";
      add_proto qw/uint32_t/, "vpx_highbd_${bd}_sub_pixel_avg_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, uint32_t *sse, const uint8_t *second_pred";
      if ($w != 128 && $h != 128 && $w != 4 && $h != 4) {
        specialize "vpx_highbd_${bd}_variance${w}x${h}", "sse2";
      }
      if ($w == 4 && $h == 4) {
        specialize "vpx_highbd_${bd}_variance${w}x${h}", "sse4_1";
      }
      if ($w != 128 && $h != 128 && $w != 4) {
        specialize "vpx_highbd_${bd}_sub_pixel_variance${w}x${h}", $sse2_x86inc;
        specialize "vpx_highbd_${bd}_sub_pixel_avg_variance${w}x${h}", $sse2_x86inc;
      }
      if ($w == 4 && $h == 4) {
        specialize "vpx_highbd_${bd}_sub_pixel_variance${w}x${h}", "sse4_1";
        specialize "vpx_highbd_${bd}_sub_pixel_avg_variance${w}x${h}", "sse4_1";
      }
    }
  }
}  # CONFIG_VP9_HIGHBITDEPTH

if (vpx_config("CONFIG_EXT_INTER") eq "yes") {
#
# Masked Variance / Masked Subpixel Variance
#
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/unsigned int/, "vpx_masked_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *mask, int mask_stride, unsigned int *sse";
    add_proto qw/unsigned int/, "vpx_masked_sub_pixel_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, const uint8_t *mask, int mask_stride, unsigned int *sse";
    specialize "vpx_masked_variance${w}x${h}", qw/ssse3/;
    specialize "vpx_masked_sub_pixel_variance${w}x${h}", qw/ssse3/;
  }

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    foreach $bd ("_", "_10_", "_12_") {
      foreach (@block_sizes) {
        ($w, $h) = @$_;
        add_proto qw/unsigned int/, "vpx_highbd${bd}masked_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, const uint8_t *mask, int mask_stride, unsigned int *sse";
        add_proto qw/unsigned int/, "vpx_highbd${bd}masked_sub_pixel_variance${w}x${h}", "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, const uint8_t *m, int m_stride, unsigned int *sse";
        specialize "vpx_highbd${bd}masked_variance${w}x${h}", qw/ssse3/;
        specialize "vpx_highbd${bd}masked_sub_pixel_variance${w}x${h}", qw/ssse3/;
      }
    }
  }
}

#
# OBMC Variance / OBMC Subpixel Variance
#
if (vpx_config("CONFIG_OBMC") eq "yes") {
  foreach (@block_sizes) {
    ($w, $h) = @$_;
    add_proto qw/unsigned int/, "vpx_obmc_variance${w}x${h}", "const uint8_t *pre_ptr, int pre_stride, const int32_t *wsrc_ptr, const int32_t *mask, unsigned int *sse";
    add_proto qw/unsigned int/, "vpx_obmc_sub_pixel_variance${w}x${h}", "const uint8_t *pre_ptr, int pre_stride, int xoffset, int  yoffset, const int32_t *wsrc_ptr, const int32_t *mask, unsigned int *sse";
    specialize "vpx_obmc_variance${w}x${h}";
    specialize "vpx_obmc_sub_pixel_variance${w}x${h}";
  }

  if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
    foreach $bd ("_", "_10_", "_12_") {
      foreach (@block_sizes) {
        ($w, $h) = @$_;
        add_proto qw/unsigned int/, "vpx_highbd${bd}obmc_variance${w}x${h}", "const uint8_t *pre_ptr, int pre_stride, const int32_t *wsrc_ptr, const int32_t *mask, unsigned int *sse";
        add_proto qw/unsigned int/, "vpx_highbd${bd}obmc_sub_pixel_variance${w}x${h}", "const uint8_t *pre_ptr, int pre_stride, int xoffset, int  yoffset, const int32_t *wsrc_ptr, const int32_t *mask, unsigned int *sse";
        specialize "vpx_highbd${bd}obmc_variance${w}x${h}";
        specialize "vpx_highbd${bd}obmc_sub_pixel_variance${w}x${h}";
      }
    }
  }
}

#
# Specialty Subpixel
#
add_proto qw/uint32_t vpx_variance_halfpixvar16x16_h/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, uint32_t *sse";
specialize qw/vpx_variance_halfpixvar16x16_h mmx sse2 media/;

add_proto qw/uint32_t vpx_variance_halfpixvar16x16_v/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, uint32_t *sse";
specialize qw/vpx_variance_halfpixvar16x16_v mmx sse2 media/;

add_proto qw/uint32_t vpx_variance_halfpixvar16x16_hv/, "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, uint32_t *sse";
specialize qw/vpx_variance_halfpixvar16x16_hv mmx sse2 media/;

#
# Comp Avg
#
add_proto qw/void vpx_comp_avg_pred/, "uint8_t *comp_pred, const uint8_t *pred, int width, int height, const uint8_t *ref, int ref_stride";
if (vpx_config("CONFIG_VP9_HIGHBITDEPTH") eq "yes") {
  add_proto qw/void vpx_highbd_comp_avg_pred/, "uint16_t *comp_pred, const uint8_t *pred8, int width, int height, const uint8_t *ref8, int ref_stride";
}

}  # CONFIG_ENCODERS || CONFIG_POSTPROC || CONFIG_VP9_POSTPROC

1;
