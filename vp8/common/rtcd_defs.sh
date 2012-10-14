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
prototype void vp8_intra4x4_predict "unsigned char *Above, unsigned char *yleft, int left_stride, B_PREDICTION_MODE b_mode, unsigned char *dst, int dst_stride, unsigned char top_left"
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
specialize vp8_recon2b
prototype void vp8_recon4b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp8_recon4b
prototype void vp8_recon_mb "MACROBLOCKD *x"
specialize vp8_recon_mb
prototype void vp8_recon_mby "MACROBLOCKD *x"
specialize vp8_recon_mby
prototype void vp8_build_intra_predictors_mby_s "MACROBLOCKD *x"
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