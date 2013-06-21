vp9_common_forward_decls() {
cat <<EOF
/*
 * VP9
 */

#include "vpx/vpx_integer.h"
#include "vp9/common/vp9_enums.h"

struct loop_filter_info;
struct macroblockd;
struct loop_filter_info;

/* Encoder forward decls */
struct macroblock;
struct vp9_variance_vtable;

#define DEC_MVCOSTS int *mvjcost, int *mvcost[2]
union int_mv;
struct yv12_buffer_config;
EOF
}
forward_decls vp9_common_forward_decls

#
# Dequant
#

prototype void vp9_idct_add_16x16 "int16_t *input, uint8_t *dest, int stride, int eob"
specialize vp9_idct_add_16x16

prototype void vp9_idct_add_8x8 "int16_t *input, uint8_t *dest, int stride, int eob"
specialize vp9_idct_add_8x8

prototype void vp9_idct_add "int16_t *input, uint8_t *dest, int stride, int eob"
specialize vp9_idct_add



prototype void vp9_idct_add_32x32 "int16_t *q, uint8_t *dst, int stride, int eob"
specialize vp9_idct_add_32x32

#
# RECON
#
prototype void vp9_copy_mem16x16 "const uint8_t *src, int src_pitch, uint8_t *dst, int dst_pitch"
specialize vp9_copy_mem16x16 mmx sse2 dspr2
vp9_copy_mem16x16_dspr2=vp9_copy_mem16x16_dspr2

prototype void vp9_copy_mem8x8 "const uint8_t *src, int src_pitch, uint8_t *dst, int dst_pitch"
specialize vp9_copy_mem8x8 mmx dspr2
vp9_copy_mem8x8_dspr2=vp9_copy_mem8x8_dspr2

prototype void vp9_copy_mem8x4 "const uint8_t *src, int src_pitch, uint8_t *dst, int dst_pitch"
specialize vp9_copy_mem8x4 mmx

prototype void vp9_build_intra_predictors "uint8_t *src, int src_stride, uint8_t *pred, int y_stride, int mode, int bw, int bh, int up_available, int left_available, int right_available"
specialize void vp9_build_intra_predictors

prototype void vp9_build_intra_predictors_sby_s "struct macroblockd *x, enum BLOCK_SIZE_TYPE bsize"
specialize vp9_build_intra_predictors_sby_s

prototype void vp9_build_intra_predictors_sbuv_s "struct macroblockd *x, enum BLOCK_SIZE_TYPE bsize"
specialize vp9_build_intra_predictors_sbuv_s

prototype void vp9_intra4x4_predict "struct macroblockd *xd, int block, enum BLOCK_SIZE_TYPE bsize, int b_mode, uint8_t *predictor, int pre_stride"
specialize vp9_intra4x4_predict;

if [ "$CONFIG_VP9_DECODER" = "yes" ]; then
prototype void vp9_add_constant_residual_8x8 "const int16_t diff, uint8_t *dest, int stride"
specialize vp9_add_constant_residual_8x8 sse2

prototype void vp9_add_constant_residual_16x16 "const int16_t diff, uint8_t *dest, int stride"
specialize vp9_add_constant_residual_16x16 sse2

prototype void vp9_add_constant_residual_32x32 "const int16_t diff, uint8_t *dest, int stride"
specialize vp9_add_constant_residual_32x32 sse2
fi

#
# Loopfilter
#
prototype void vp9_mb_lpf_vertical_edge_w "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh"
specialize vp9_mb_lpf_vertical_edge_w sse2

prototype void vp9_mbloop_filter_vertical_edge "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count"
specialize vp9_mbloop_filter_vertical_edge sse2

prototype void vp9_loop_filter_vertical_edge "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count"
specialize vp9_loop_filter_vertical_edge mmx

prototype void vp9_mb_lpf_horizontal_edge_w "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh"
specialize vp9_mb_lpf_horizontal_edge_w sse2

prototype void vp9_mbloop_filter_horizontal_edge "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count"
specialize vp9_mbloop_filter_horizontal_edge sse2

prototype void vp9_loop_filter_horizontal_edge "uint8_t *s, int pitch, const uint8_t *blimit, const uint8_t *limit, const uint8_t *thresh, int count"
specialize vp9_loop_filter_horizontal_edge mmx

#
# post proc
#
if [ "$CONFIG_POSTPROC" = "yes" ]; then
prototype void vp9_mbpost_proc_down "uint8_t *dst, int pitch, int rows, int cols, int flimit"
specialize vp9_mbpost_proc_down mmx sse2
vp9_mbpost_proc_down_sse2=vp9_mbpost_proc_down_xmm

prototype void vp9_mbpost_proc_across_ip "uint8_t *src, int pitch, int rows, int cols, int flimit"
specialize vp9_mbpost_proc_across_ip sse2
vp9_mbpost_proc_across_ip_sse2=vp9_mbpost_proc_across_ip_xmm

prototype void vp9_post_proc_down_and_across "const uint8_t *src_ptr, uint8_t *dst_ptr, int src_pixels_per_line, int dst_pixels_per_line, int rows, int cols, int flimit"
specialize vp9_post_proc_down_and_across mmx sse2
vp9_post_proc_down_and_across_sse2=vp9_post_proc_down_and_across_xmm

prototype void vp9_plane_add_noise "uint8_t *Start, char *noise, char blackclamp[16], char whiteclamp[16], char bothclamp[16], unsigned int Width, unsigned int Height, int Pitch"
specialize vp9_plane_add_noise mmx sse2
vp9_plane_add_noise_sse2=vp9_plane_add_noise_wmt
fi

prototype void vp9_blend_mb_inner "uint8_t *y, uint8_t *u, uint8_t *v, int y1, int u1, int v1, int alpha, int stride"
specialize vp9_blend_mb_inner

prototype void vp9_blend_mb_outer "uint8_t *y, uint8_t *u, uint8_t *v, int y1, int u1, int v1, int alpha, int stride"
specialize vp9_blend_mb_outer

prototype void vp9_blend_b "uint8_t *y, uint8_t *u, uint8_t *v, int y1, int u1, int v1, int alpha, int stride"
specialize vp9_blend_b

#
# Sub Pixel Filters
#
prototype void vp9_convolve8 "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8 ssse3

prototype void vp9_convolve8_horiz "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8_horiz ssse3

prototype void vp9_convolve8_vert "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8_vert ssse3

prototype void vp9_convolve8_avg "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8_avg ssse3

prototype void vp9_convolve8_avg_horiz "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8_avg_horiz ssse3

prototype void vp9_convolve8_avg_vert "const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, const int16_t *filter_x, int x_step_q4, const int16_t *filter_y, int y_step_q4, int w, int h"
specialize vp9_convolve8_avg_vert ssse3

#
# dct
#
prototype void vp9_short_idct4x4_1_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct4x4_1_add

prototype void vp9_short_idct4x4_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct4x4_add sse2

prototype void vp9_short_idct8x8_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct8x8_add sse2

prototype void vp9_short_idct10_8x8_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct10_8x8_add sse2

prototype void vp9_short_idct1_8x8 "int16_t *input, int16_t *output"
specialize vp9_short_idct1_8x8

prototype void vp9_short_idct16x16_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct16x16_add sse2

prototype void vp9_short_idct10_16x16_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct10_16x16_add sse2

prototype void vp9_short_idct1_16x16 "int16_t *input, int16_t *output"
specialize vp9_short_idct1_16x16

prototype void vp9_short_idct32x32_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct32x32_add sse2

prototype void vp9_short_idct1_32x32 "int16_t *input, int16_t *output"
specialize vp9_short_idct1_32x32

prototype void vp9_short_idct10_32x32_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_idct10_32x32_add

prototype void vp9_short_iht4x4_add "int16_t *input, uint8_t *dest, int dest_stride, int tx_type"
specialize vp9_short_iht4x4_add

prototype void vp9_short_iht8x8_add "int16_t *input, uint8_t *dest, int dest_stride, int tx_type"
specialize vp9_short_iht8x8_add

prototype void vp9_short_iht16x16_add "int16_t *input, uint8_t *output, int pitch, int tx_type"
specialize vp9_short_iht16x16_add

prototype void vp9_idct4_1d "int16_t *input, int16_t *output"
specialize vp9_idct4_1d sse2
# dct and add

prototype void vp9_dc_only_idct_add "int input_dc, uint8_t *pred_ptr, uint8_t *dst_ptr, int pitch, int stride"
specialize vp9_dc_only_idct_add sse2

prototype void vp9_short_iwalsh4x4_1_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_iwalsh4x4_1_add

prototype void vp9_short_iwalsh4x4_add "int16_t *input, uint8_t *dest, int dest_stride"
specialize vp9_short_iwalsh4x4_add

prototype unsigned int vp9_sad32x3 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int ref_stride, int max_sad"
specialize vp9_sad32x3

prototype unsigned int vp9_sad3x32 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int ref_stride, int max_sad"
specialize vp9_sad3x32

#
# Encoder functions below this point.
#
if [ "$CONFIG_VP9_ENCODER" = "yes" ]; then


# variance
[ $arch = "x86_64" ] && mmx_x86_64=mmx && sse2_x86_64=sse2

prototype unsigned int vp9_variance32x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance32x16 sse2

prototype unsigned int vp9_variance16x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance16x32 sse2

prototype unsigned int vp9_variance64x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance64x32 sse2

prototype unsigned int vp9_variance32x64 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance32x64 sse2

prototype unsigned int vp9_variance32x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance32x32 sse2

prototype unsigned int vp9_variance64x64 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance64x64 sse2

prototype unsigned int vp9_variance16x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance16x16 mmx sse2

prototype unsigned int vp9_variance16x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance16x8 mmx sse2

prototype unsigned int vp9_variance8x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance8x16 mmx sse2

prototype unsigned int vp9_variance8x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance8x8 mmx sse2

prototype void vp9_get_sse_sum_8x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, int *sum"
specialize vp9_get_sse_sum_8x8 sse2
vp9_get_sse_sum_8x8_sse2=vp9_get8x8var_sse2

prototype unsigned int vp9_variance8x4 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance8x4 sse2

prototype unsigned int vp9_variance4x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance4x8 sse2

prototype unsigned int vp9_variance4x4 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance4x4 mmx sse2

prototype unsigned int vp9_sub_pixel_variance64x64 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance64x64 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance64x64 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance64x64 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance32x64 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance32x64 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance32x64 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance32x64 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance64x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance64x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance64x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance64x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance32x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance32x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance32x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance32x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance16x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance16x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance16x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance16x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance32x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance32x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance32x32 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance32x32 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance16x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance16x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance16x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance16x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance8x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance8x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance8x16 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance8x16 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance16x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance16x8 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance16x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance16x8 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance8x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance8x8 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance8x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance8x8 sse2 ssse3

# TODO(jingning): need to convert 8x4/4x8 functions into mmx/sse form
prototype unsigned int vp9_sub_pixel_variance8x4 "const uint8_t *src_ptr, int source_stride, int xoffset, int yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance8x4 sse2 ssse3

prototype unsigned int vp9_sub_pixel_avg_variance8x4 "const uint8_t *src_ptr, int source_stride, int xoffset, int yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance8x4 sse2 ssse3

prototype unsigned int vp9_sub_pixel_variance4x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance4x8 sse ssse3

prototype unsigned int vp9_sub_pixel_avg_variance4x8 "const uint8_t *src_ptr, int source_stride, int xoffset, int yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance4x8 sse ssse3

prototype unsigned int vp9_sub_pixel_variance4x4 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_variance4x4 sse ssse3
#vp9_sub_pixel_variance4x4_sse2=vp9_sub_pixel_variance4x4_wmt

prototype unsigned int vp9_sub_pixel_avg_variance4x4 "const uint8_t *src_ptr, int source_stride, int xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse, const uint8_t *second_pred"
specialize vp9_sub_pixel_avg_variance4x4 sse ssse3

prototype unsigned int vp9_sad64x64 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad64x64 sse2

prototype unsigned int vp9_sad32x64 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad32x64 sse2

prototype unsigned int vp9_sad64x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad64x32 sse2

prototype unsigned int vp9_sad32x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad32x16 sse2

prototype unsigned int vp9_sad16x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad16x32 sse2

prototype unsigned int vp9_sad32x32 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad32x32 sse2

prototype unsigned int vp9_sad16x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad16x16 mmx sse2

prototype unsigned int vp9_sad16x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad16x8 mmx sse2

prototype unsigned int vp9_sad8x16 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad8x16 mmx sse2

prototype unsigned int vp9_sad8x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad8x8 mmx sse2

# TODO(jingning): need to covert these functions into mmx/sse2 form
prototype unsigned int vp9_sad8x4 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad8x4 sse2

prototype unsigned int vp9_sad4x8 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int max_sad"
specialize vp9_sad4x8 sse

prototype unsigned int vp9_sad4x4 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad4x4 mmx sse

prototype unsigned int vp9_variance_halfpixvar16x16_h "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_h sse2
vp9_variance_halfpixvar16x16_h_sse2=vp9_variance_halfpixvar16x16_h_wmt

prototype unsigned int vp9_variance_halfpixvar16x16_v "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_v sse2
vp9_variance_halfpixvar16x16_v_sse2=vp9_variance_halfpixvar16x16_v_wmt

prototype unsigned int vp9_variance_halfpixvar16x16_hv "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_hv sse2
vp9_variance_halfpixvar16x16_hv_sse2=vp9_variance_halfpixvar16x16_hv_wmt

prototype unsigned int vp9_variance_halfpixvar64x64_h "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar64x64_h

prototype unsigned int vp9_variance_halfpixvar64x64_v "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar64x64_v

prototype unsigned int vp9_variance_halfpixvar64x64_hv "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar64x64_hv

prototype unsigned int vp9_variance_halfpixvar32x32_h "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_h

prototype unsigned int vp9_variance_halfpixvar32x32_v "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_v

prototype unsigned int vp9_variance_halfpixvar32x32_hv "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_hv

prototype void vp9_sad64x64x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad64x64x3

prototype void vp9_sad32x32x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x32x3

prototype void vp9_sad16x16x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x16x3 sse3 ssse3

prototype void vp9_sad16x8x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x8x3 sse3 ssse3

prototype void vp9_sad8x16x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x16x3 sse3

prototype void vp9_sad8x8x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x8x3 sse3

prototype void vp9_sad4x4x3 "const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad4x4x3 sse3

prototype void vp9_sad64x64x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad64x64x8

prototype void vp9_sad32x32x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad32x32x8

prototype void vp9_sad16x16x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad16x16x8 sse4

prototype void vp9_sad16x8x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad16x8x8 sse4

prototype void vp9_sad8x16x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad8x16x8 sse4

prototype void vp9_sad8x8x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad8x8x8 sse4

prototype void vp9_sad8x4x8 "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array"
specialize vp9_sad8x4x8

prototype void vp9_sad4x8x8 "const uint8_t *src_ptr, int src_stride, const uint8_t *ref_ptr, int ref_stride, uint32_t *sad_array"
specialize vp9_sad4x8x8

prototype void vp9_sad4x4x8 "const uint8_t *src_ptr, int  src_stride, const uint8_t *ref_ptr, int  ref_stride, uint32_t *sad_array"
specialize vp9_sad4x4x8 sse4

prototype void vp9_sad64x64x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad64x64x4d sse2

prototype void vp9_sad32x64x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x64x4d sse2

prototype void vp9_sad64x32x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad64x32x4d sse2

prototype void vp9_sad32x16x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x16x4d sse2

prototype void vp9_sad16x32x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x32x4d sse2

prototype void vp9_sad32x32x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x32x4d sse2

prototype void vp9_sad16x16x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x16x4d sse2

prototype void vp9_sad16x8x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x8x4d sse2

prototype void vp9_sad8x16x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x16x4d sse2

prototype void vp9_sad8x8x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x8x4d sse2

# TODO(jingning): need to convert these 4x8/8x4 functions into sse2 form
prototype void vp9_sad8x4x4d "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, unsigned int *sad_array"
specialize vp9_sad8x4x4d sse2

prototype void vp9_sad4x8x4d "const uint8_t *src_ptr, int src_stride, const uint8_t* const ref_ptr[], int ref_stride, unsigned int *sad_array"
specialize vp9_sad4x8x4d sse

prototype void vp9_sad4x4x4d "const uint8_t *src_ptr, int  src_stride, const uint8_t* const ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad4x4x4d sse

#prototype unsigned int vp9_sub_pixel_mse16x16 "const uint8_t *src_ptr, int  src_pixels_per_line, int  xoffset, int  yoffset, const uint8_t *dst_ptr, int dst_pixels_per_line, unsigned int *sse"
#specialize vp9_sub_pixel_mse16x16 sse2 mmx

prototype unsigned int vp9_mse16x16 "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp9_mse16x16 mmx sse2
vp9_mse16x16_sse2=vp9_mse16x16_wmt

prototype unsigned int vp9_mse8x16 "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp9_mse8x16

prototype unsigned int vp9_mse16x8 "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp9_mse16x8

prototype unsigned int vp9_mse8x8 "const uint8_t *src_ptr, int  source_stride, const uint8_t *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp9_mse8x8

prototype unsigned int vp9_sub_pixel_mse64x64 "const uint8_t *src_ptr, int  source_stride, int  xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_mse64x64

prototype unsigned int vp9_sub_pixel_mse32x32 "const uint8_t *src_ptr, int  source_stride, int  xoffset, int  yoffset, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_sub_pixel_mse32x32

prototype unsigned int vp9_get_mb_ss "const int16_t *"
specialize vp9_get_mb_ss mmx sse2
# ENCODEMB INVOKE

prototype int64_t vp9_block_error "int16_t *coeff, int16_t *dqcoeff, intptr_t block_size"
specialize vp9_block_error sse2

prototype void vp9_subtract_block "int rows, int cols, int16_t *diff_ptr, ptrdiff_t diff_stride, const uint8_t *src_ptr, ptrdiff_t src_stride, const uint8_t *pred_ptr, ptrdiff_t pred_stride"
specialize vp9_subtract_block sse2

#
# Structured Similarity (SSIM)
#
if [ "$CONFIG_INTERNAL_STATS" = "yes" ]; then
    [ $arch = "x86_64" ] && sse2_on_x86_64=sse2

    prototype void vp9_ssim_parms_8x8 "uint8_t *s, int sp, uint8_t *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp9_ssim_parms_8x8 $sse2_on_x86_64

    prototype void vp9_ssim_parms_16x16 "uint8_t *s, int sp, uint8_t *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp9_ssim_parms_16x16 $sse2_on_x86_64
fi

# fdct functions
prototype void vp9_short_fht4x4 "int16_t *InputData, int16_t *OutputData, int pitch, int tx_type"
specialize vp9_short_fht4x4

prototype void vp9_short_fht8x8 "int16_t *InputData, int16_t *OutputData, int pitch, int tx_type"
specialize vp9_short_fht8x8

prototype void vp9_short_fht16x16 "int16_t *InputData, int16_t *OutputData, int pitch, int tx_type"
specialize vp9_short_fht16x16

prototype void vp9_short_fdct8x8 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct8x8 sse2

prototype void vp9_short_fdct4x4 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct4x4 sse2

prototype void vp9_short_fdct8x4 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct8x4 sse2

prototype void vp9_short_fdct32x32 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct32x32

prototype void vp9_short_fdct32x32_rd "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct32x32_rd

prototype void vp9_short_fdct16x16 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_fdct16x16 sse2

prototype void vp9_short_walsh4x4 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_walsh4x4

prototype void vp9_short_walsh8x4 "int16_t *InputData, int16_t *OutputData, int pitch"
specialize vp9_short_walsh8x4

#
# Motion search
#
prototype int vp9_full_search_sad "struct macroblock *x, union int_mv *ref_mv, int sad_per_bit, int distance, struct vp9_variance_vtable *fn_ptr, DEC_MVCOSTS, union int_mv *center_mv, int n"
specialize vp9_full_search_sad sse3 sse4_1
vp9_full_search_sad_sse3=vp9_full_search_sadx3
vp9_full_search_sad_sse4_1=vp9_full_search_sadx8

prototype int vp9_refining_search_sad "struct macroblock *x, union int_mv *ref_mv, int sad_per_bit, int distance, struct vp9_variance_vtable *fn_ptr, DEC_MVCOSTS, union int_mv *center_mv"
specialize vp9_refining_search_sad sse3
vp9_refining_search_sad_sse3=vp9_refining_search_sadx4

prototype int vp9_diamond_search_sad "struct macroblock *x, union int_mv *ref_mv, union int_mv *best_mv, int search_param, int sad_per_bit, int *num00, struct vp9_variance_vtable *fn_ptr, DEC_MVCOSTS, union int_mv *center_mv"
specialize vp9_diamond_search_sad sse3
vp9_diamond_search_sad_sse3=vp9_diamond_search_sadx4

prototype void vp9_temporal_filter_apply "uint8_t *frame1, unsigned int stride, uint8_t *frame2, unsigned int block_size, int strength, int filter_weight, unsigned int *accumulator, uint16_t *count"
specialize vp9_temporal_filter_apply sse2

prototype void vp9_yv12_copy_partial_frame "struct yv12_buffer_config *src_ybc, struct yv12_buffer_config *dst_ybc, int fraction"
specialize vp9_yv12_copy_partial_frame


fi
# end encoder functions
