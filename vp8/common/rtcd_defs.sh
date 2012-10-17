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
specialize vp8_loop_filter_mbv;

prototype void vp8_loop_filter_bv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bv;

prototype void vp8_loop_filter_bv8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bv8x8;

prototype void vp8_loop_filter_mbh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_mbh;

prototype void vp8_loop_filter_bh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bh;

prototype void vp8_loop_filter_bh8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bh8x8;

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

