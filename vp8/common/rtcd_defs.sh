common_forward_decls() {
cat <<EOF
#include "vp8/common/blockd.h"

struct loop_filter_info;

/* Encoder forward decls */
struct variance_vtable;
union int_mv;
struct yv12_buffer_config;
EOF
}
forward_decls common_forward_decls

prototype void vp8_filter_block2d_4x4_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp8_filter_block2d_8x4_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp8_filter_block2d_8x8_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp8_filter_block2d_16x16_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"

# At the very least, MSVC 2008 has compiler bug exhibited by this code; code
# compiles warning free but a dissassembly of generated code show bugs. To be
# on the safe side, only enabled when compiled with 'gcc'.
if [ "$CONFIG_GCC" = "yes" ]; then
    specialize vp8_filter_block2d_4x4_8 sse4_1 sse2
    specialize vp8_filter_block2d_8x4_8 sse4_1 sse2
    specialize vp8_filter_block2d_8x8_8 sse4_1 sse2
    specialize vp8_filter_block2d_16x16_8 sse4_1 sse2
fi


#
# RECON
#
prototype void vp8_copy_mem16x16 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem16x16 mmx sse2 media neon dspr2
vp8_copy_mem16x16_media=vp8_copy_mem16x16_v6
vp8_copy_mem16x16_dspr2=vp8_copy_mem16x16_dspr2

prototype void vp8_copy_mem8x8 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem8x8 mmx media neon dspr2
vp8_copy_mem8x8_media=vp8_copy_mem8x8_v6
vp8_copy_mem8x8_dspr2=vp8_copy_mem8x8_dspr2

prototype void vp8_copy_mem8x4 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem8x4 mmx

prototype void vp8_intra4x4_predict "unsigned char *Above, unsigned char *yleft, int left_stride, B_PREDICTION_MODE b_mode, unsigned char *dst, int dst_stride, unsigned char top_left"
specialize vp8_intra4x4_predict

prototype void vp8_avg_mem16x16 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_avg_mem16x16

prototype void vp8_avg_mem8x8 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_avg_mem8x8

prototype void vp8_copy_mem8x4 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem8x4 mmx media neon dspr2
vp8_copy_mem8x4_media=vp8_copy_mem8x4_v6
vp8_copy_mem8x4_dspr2=vp8_copy_mem8x4_dspr2

prototype void vp8_recon_b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp8_recon_b

prototype void vp8_recon_uv_b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp8_recon_uv_b

prototype void vp8_recon2b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp8_recon2b sse2

prototype void vp8_recon4b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp8_recon4b sse2

prototype void vp8_recon_mb "MACROBLOCKD *x"
specialize vp8_recon_mb

prototype void vp8_recon_mby "MACROBLOCKD *x"
specialize vp8_recon_mby

prototype void vp8_build_intra_predictors_mby_s "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_mby_s

prototype void vp8_build_intra_predictors_sby_s "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_sby_s;

prototype void vp8_build_intra_predictors_sbuv_s "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_sbuv_s;

prototype void vp8_build_intra_predictors_mby "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_mby;

prototype void vp8_build_comp_intra_predictors_mby "MACROBLOCKD *x"
specialize vp8_build_comp_intra_predictors_mby;

prototype void vp8_build_intra_predictors_mby_s "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_mby_s;

prototype void vp8_build_intra_predictors_mbuv "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_mbuv;

prototype void vp8_build_intra_predictors_mbuv_s "MACROBLOCKD *x"
specialize vp8_build_intra_predictors_mbuv_s;

prototype void vp8_build_comp_intra_predictors_mbuv "MACROBLOCKD *x"
specialize vp8_build_comp_intra_predictors_mbuv;

prototype void vp8_intra4x4_predict "BLOCKD *x, int b_mode, unsigned char *predictor"
specialize vp8_intra4x4_predict;

prototype void vp8_comp_intra4x4_predict "BLOCKD *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp8_comp_intra4x4_predict;

prototype void vp8_intra8x8_predict "BLOCKD *x, int b_mode, unsigned char *predictor"
specialize vp8_intra8x8_predict;

prototype void vp8_comp_intra8x8_predict "BLOCKD *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp8_comp_intra8x8_predict;

prototype void vp8_intra_uv4x4_predict "BLOCKD *x, int b_mode, unsigned char *predictor"
specialize vp8_intra_uv4x4_predict;

prototype void vp8_comp_intra_uv4x4_predict "BLOCKD *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp8_comp_intra_uv4x4_predict;

#
# Loopfilter
#
prototype void vp8_loop_filter_mbv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_mbv sse2

prototype void vp8_loop_filter_bv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bv sse2

prototype void vp8_loop_filter_bv8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bv8x8 sse2

prototype void vp8_loop_filter_mbh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_mbh sse2

prototype void vp8_loop_filter_bh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bh sse2

prototype void vp8_loop_filter_bh8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bh8x8 sse2

prototype void vp8_loop_filter_simple_mbv "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp8_loop_filter_simple_mbv mmx sse2 media neon
vp8_loop_filter_simple_mbv_c=vp8_loop_filter_simple_vertical_edge_c
vp8_loop_filter_simple_mbv_mmx=vp8_loop_filter_simple_vertical_edge_mmx
vp8_loop_filter_simple_mbv_sse2=vp8_loop_filter_simple_vertical_edge_sse2
vp8_loop_filter_simple_mbv_media=vp8_loop_filter_simple_vertical_edge_armv6
vp8_loop_filter_simple_mbv_neon=vp8_loop_filter_mbvs_neon

prototype void vp8_loop_filter_simple_mbh "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp8_loop_filter_simple_mbh mmx sse2 media neon
vp8_loop_filter_simple_mbh_c=vp8_loop_filter_simple_horizontal_edge_c
vp8_loop_filter_simple_mbh_mmx=vp8_loop_filter_simple_horizontal_edge_mmx
vp8_loop_filter_simple_mbh_sse2=vp8_loop_filter_simple_horizontal_edge_sse2
vp8_loop_filter_simple_mbh_media=vp8_loop_filter_simple_horizontal_edge_armv6
vp8_loop_filter_simple_mbh_neon=vp8_loop_filter_mbhs_neon

prototype void vp8_loop_filter_simple_bv "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp8_loop_filter_simple_bv mmx sse2 media neon
vp8_loop_filter_simple_bv_c=vp8_loop_filter_bvs_c
vp8_loop_filter_simple_bv_mmx=vp8_loop_filter_bvs_mmx
vp8_loop_filter_simple_bv_sse2=vp8_loop_filter_bvs_sse2
vp8_loop_filter_simple_bv_media=vp8_loop_filter_bvs_armv6
vp8_loop_filter_simple_bv_neon=vp8_loop_filter_bvs_neon

prototype void vp8_loop_filter_simple_bh "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp8_loop_filter_simple_bh mmx sse2 media neon
vp8_loop_filter_simple_bh_c=vp8_loop_filter_bhs_c
vp8_loop_filter_simple_bh_mmx=vp8_loop_filter_bhs_mmx
vp8_loop_filter_simple_bh_sse2=vp8_loop_filter_bhs_sse2
vp8_loop_filter_simple_bh_media=vp8_loop_filter_bhs_armv6
vp8_loop_filter_simple_bh_neon=vp8_loop_filter_bhs_neon

#
# Encoder functions below this point.
#
if [ "$CONFIG_VP8_ENCODER" = "yes" ]; then


# variance
[ $arch = "x86_64" ] && mmx_x86_64=mmx && sse2_x86_64=sse2

prototype unsigned int vp8_variance32x32 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance32x32

prototype unsigned int vp8_variance16x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance16x16 mmx sse2
vp8_variance16x16_sse2=vp8_variance16x16_wmt
vp8_variance16x16_mmx=vp8_variance16x16_mmx

prototype unsigned int vp8_variance16x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance16x8 mmx sse2
vp8_variance16x8_sse2=vp8_variance16x8_wmt
vp8_variance16x8_mmx=vp8_variance16x8_mmx

prototype unsigned int vp8_variance8x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance8x16 mmx sse2
vp8_variance8x16_sse2=vp8_variance8x16_wmt
vp8_variance8x16_mmx=vp8_variance8x16_mmx

prototype unsigned int vp8_variance8x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance8x8 mmx sse2
vp8_variance8x8_sse2=vp8_variance8x8_wmt
vp8_variance8x8_mmx=vp8_variance8x8_mmx

prototype unsigned int vp8_variance4x4 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance4x4 mmx sse2
vp8_variance4x4_sse2=vp8_variance4x4_wmt
vp8_variance4x4_mmx=vp8_variance4x4_mmx

prototype unsigned int vp8_sub_pixel_variance32x32 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance32x32

prototype unsigned int vp8_sub_pixel_variance16x16 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance16x16 sse2 mmx ssse3
vp8_sub_pixel_variance16x16_sse2=vp8_sub_pixel_variance16x16_wmt

prototype unsigned int vp8_sub_pixel_variance8x16 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance8x16 sse2 mmx
vp8_sub_pixel_variance8x16_sse2=vp8_sub_pixel_variance8x16_wmt

prototype unsigned int vp8_sub_pixel_variance16x8 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance16x8 sse2 mmx ssse3
vp8_sub_pixel_variance16x8_sse2=vp8_sub_pixel_variance16x8_ssse3;
vp8_sub_pixel_variance16x8_sse2=vp8_sub_pixel_variance16x8_wmt

prototype unsigned int vp8_sub_pixel_variance8x8 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance8x8 sse2 mmx
vp8_sub_pixel_variance8x8_sse2=vp8_sub_pixel_variance8x8_wmt

prototype unsigned int vp8_sub_pixel_variance4x4 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_variance4x4 sse2 mmx
vp8_sub_pixel_variance4x4_sse2=vp8_sub_pixel_variance4x4_wmt

prototype unsigned int vp8_sad32x32 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad32x32

prototype unsigned int vp8_sad16x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad16x16 mmx sse2 sse3
vp8_sad16x16_sse2=vp8_sad16x16_wmt

prototype unsigned int vp8_sad16x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad16x8 mmx sse2
vp8_sad16x8_sse2=vp8_sad16x8_wmt

prototype unsigned int vp8_sad8x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad8x16 mmx sse2
vp8_sad8x16_sse2=vp8_sad8x16_wmt

prototype unsigned int vp8_sad8x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad8x8 mmx sse2
vp8_sad8x8_sse2=vp8_sad8x8_wmt

prototype unsigned int vp8_sad4x4 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp8_sad4x4 mmx sse2
vp8_sad4x4_sse2=vp8_sad4x4_wmt

prototype unsigned int vp8_variance_halfpixvar16x16_h "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar16x16_h mmx sse2
vp8_variance_halfpixvar16x16_h_sse2=vp8_variance_halfpixvar16x16_h_wmt

prototype unsigned int vp8_variance_halfpixvar16x16_v "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar16x16_v mmx sse2
vp8_variance_halfpixvar16x16_v_sse2=vp8_variance_halfpixvar16x16_v_wmt

prototype unsigned int vp8_variance_halfpixvar16x16_hv "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar16x16_hv mmx sse2
vp8_variance_halfpixvar16x16_hv_sse2=vp8_variance_halfpixvar16x16_hv_wmt

prototype unsigned int vp8_variance_halfpixvar32x32_h "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar32x32_h

prototype unsigned int vp8_variance_halfpixvar32x32_v "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar32x32_v

prototype unsigned int vp8_variance_halfpixvar32x32_hv "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp8_variance_halfpixvar32x32_hv

prototype void vp8_sad32x32x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad32x32x3

prototype void vp8_sad16x16x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad16x16x3 sse3 ssse3

prototype void vp8_sad16x8x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad16x8x3 sse3 ssse3

prototype void vp8_sad8x16x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad8x16x3 sse3

prototype void vp8_sad8x8x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad8x8x3 sse3

prototype void vp8_sad4x4x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp8_sad4x4x3 sse3

prototype void vp8_sad32x32x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad32x32x8

prototype void vp8_sad16x16x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad16x16x8 sse4

prototype void vp8_sad16x8x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad16x8x8 sse4

prototype void vp8_sad8x16x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad8x16x8 sse4

prototype void vp8_sad8x8x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad8x8x8 sse4

prototype void vp8_sad4x4x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp8_sad4x4x8 sse4

prototype void vp8_sad32x32x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad32x32x4d

prototype void vp8_sad16x16x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad16x16x4d sse3

prototype void vp8_sad16x8x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad16x8x4d sse3

prototype void vp8_sad8x16x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad8x16x4d sse3

prototype void vp8_sad8x8x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad8x8x4d sse3

prototype void vp8_sad4x4x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp8_sad4x4x4d sse3

#
# Block copy
#
case $arch in
    x86*)
    prototype void vp8_copy32xn "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, int n"
    specialize vp8_copy32xn sse2 sse3
    ;;
esac

prototype unsigned int vp8_sub_pixel_mse16x16 "const unsigned char  *src_ptr, int  src_pixels_per_line, int  xoffset, int  yoffset, const unsigned char *dst_ptr, int dst_pixels_per_line, unsigned int *sse"
specialize vp8_sub_pixel_mse16x16 sse2 mmx
vp8_sub_pixel_mse16x16_sse2=vp8_sub_pixel_mse16x16_wmt

prototype unsigned int vp8_mse16x16 "const unsigned char *src_ptr, int  source_stride, const unsigned char *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp8_mse16x16 mmx sse2
vp8_mse16x16_sse2=vp8_mse16x16_wmt

prototype unsigned int vp8_sub_pixel_mse32x32 "const unsigned char  *src_ptr, int  source_stride, int  xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp8_sub_pixel_mse32x32

prototype unsigned int vp8_get_mb_ss "const short *"
specialize vp8_get_mb_ss mmx sse2

#
# Structured Similarity (SSIM)
#
if [ "$CONFIG_INTERNAL_STATS" = "yes" ]; then
    [ $arch = "x86_64" ] && sse2_on_x86_64=sse2

    prototype void vp8_ssim_parms_8x8 "unsigned char *s, int sp, unsigned char *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp8_ssim_parms_8x8 $sse2_on_x86_64

    prototype void vp8_ssim_parms_16x16 "unsigned char *s, int sp, unsigned char *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp8_ssim_parms_16x16 $sse2_on_x86_64
fi

fi
# end encoder functions
