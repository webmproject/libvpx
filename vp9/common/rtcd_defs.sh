common_forward_decls() {
cat <<EOF

struct loop_filter_info;
struct blockd;
struct macroblockd;
struct loop_filter_info;

/* Encoder forward decls */
struct block;
struct macroblock;
struct variance_vtable;

/* Encoder forward decls */
struct variance_vtable;
union int_mv;
struct yv12_buffer_config;
EOF
}
forward_decls common_forward_decls

prototype void vp9_filter_block2d_4x4_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp9_filter_block2d_8x4_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp9_filter_block2d_8x8_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"
prototype void vp9_filter_block2d_16x16_8 "const unsigned char *src_ptr, const unsigned int src_stride, const short *HFilter_aligned16, const short *VFilter_aligned16, unsigned char *dst_ptr, unsigned int dst_stride"

# At the very least, MSVC 2008 has compiler bug exhibited by this code; code
# compiles warning free but a dissassembly of generated code show bugs. To be
# on the safe side, only enabled when compiled with 'gcc'.
if [ "$CONFIG_GCC" = "yes" ]; then
    specialize vp9_filter_block2d_4x4_8 sse4_1 sse2
fi
    specialize vp9_filter_block2d_8x4_8 ssse3 #sse4_1 sse2
    specialize vp9_filter_block2d_8x8_8 ssse3 #sse4_1 sse2
    specialize vp9_filter_block2d_16x16_8 ssse3 #sse4_1 sse2

#
# Dequant
#
prototype void vp9_dequantize_b "struct blockd *x"
specialize vp9_dequantize_b mmx

prototype void vp9_dequantize_b_2x2 "struct blockd *x"
specialize vp9_dequantize_b_2x2

prototype void vp9_dequant_dc_idct_add_y_block_8x8 "short *q, short *dq, unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs, short *dc, struct macroblockd *xd"
specialize vp9_dequant_dc_idct_add_y_block_8x8

prototype void vp9_dequant_idct_add_y_block_8x8 "short *q, short *dq, unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs, struct macroblockd *xd"
specialize vp9_dequant_idct_add_y_block_8x8

prototype void vp9_dequant_idct_add_uv_block_8x8 "short *q, short *dq, unsigned char *pre, unsigned char *dstu, unsigned char *dstv, int stride, unsigned short *eobs, struct macroblockd *xd"
specialize vp9_dequant_idct_add_uv_block_8x8

prototype void vp9_dequant_idct_add_16x16 "short *input, short *dq, unsigned char *pred, unsigned char *dest, int pitch, int stride"
specialize vp9_dequant_idct_add_16x16

prototype void vp9_dequant_idct_add "short *input, short *dq, unsigned char *pred, unsigned char *dest, int pitch, int stride"
specialize vp9_dequant_idct_add

prototype void vp9_dequant_dc_idct_add "short *input, short *dq, unsigned char *pred, unsigned char *dest, int pitch, int stride, int Dc"
specialize vp9_dequant_dc_idct_add

prototype void vp9_dequant_dc_idct_add_y_block "short *q, short *dq, unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs, short *dc"
specialize vp9_dequant_dc_idct_add_y_block mmx

prototype void vp9_dequant_idct_add_y_block "short *q, short *dq, unsigned char *pre, unsigned char *dst, int stride, unsigned short *eobs"
specialize vp9_dequant_idct_add_y_block mmx

prototype void vp9_dequant_idct_add_uv_block "short *q, short *dq, unsigned char *pre, unsigned char *dstu, unsigned char *dstv, int stride, unsigned short *eobs"
specialize vp9_dequant_idct_add_uv_block mmx

#
# RECON
#
prototype void vp9_copy_mem16x16 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_copy_mem16x16 mmx sse2 media neon dspr2
vp9_copy_mem16x16_media=vp9_copy_mem16x16_v6
vp9_copy_mem16x16_dspr2=vp9_copy_mem16x16_dspr2

prototype void vp9_copy_mem8x8 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_copy_mem8x8 mmx media neon dspr2
vp9_copy_mem8x8_media=vp9_copy_mem8x8_v6
vp9_copy_mem8x8_dspr2=vp9_copy_mem8x8_dspr2

prototype void vp9_copy_mem8x4 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_copy_mem8x4 mmx

prototype void vp9_intra4x4_predict "unsigned char *Above, unsigned char *yleft, int left_stride, B_PREDICTION_MODE b_mode, unsigned char *dst, int dst_stride, unsigned char top_left"
specialize vp9_intra4x4_predict

prototype void vp9_avg_mem16x16 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_avg_mem16x16

prototype void vp9_avg_mem8x8 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_avg_mem8x8

prototype void vp9_copy_mem8x4 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp9_copy_mem8x4 mmx media neon dspr2
vp9_copy_mem8x4_media=vp9_copy_mem8x4_v6
vp9_copy_mem8x4_dspr2=vp9_copy_mem8x4_dspr2

prototype void vp9_recon_b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp9_recon_b

prototype void vp9_recon_uv_b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp9_recon_uv_b

prototype void vp9_recon2b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp9_recon2b sse2

prototype void vp9_recon4b "unsigned char *pred_ptr, short *diff_ptr, unsigned char *dst_ptr, int stride"
specialize vp9_recon4b sse2

prototype void vp9_recon_mb "struct macroblockd *x"
specialize vp9_recon_mb

prototype void vp9_recon_mby "struct macroblockd *x"
specialize vp9_recon_mby

prototype void vp9_build_intra_predictors_mby_s "struct macroblockd *x"
specialize vp9_build_intra_predictors_mby_s

prototype void vp9_build_intra_predictors_sby_s "struct macroblockd *x"
specialize vp9_build_intra_predictors_sby_s;

prototype void vp9_build_intra_predictors_sbuv_s "struct macroblockd *x"
specialize vp9_build_intra_predictors_sbuv_s;

prototype void vp9_build_intra_predictors_mby "struct macroblockd *x"
specialize vp9_build_intra_predictors_mby;

prototype void vp9_build_comp_intra_predictors_mby "struct macroblockd *x"
specialize vp9_build_comp_intra_predictors_mby;

prototype void vp9_build_intra_predictors_mby_s "struct macroblockd *x"
specialize vp9_build_intra_predictors_mby_s;

prototype void vp9_build_intra_predictors_mbuv "struct macroblockd *x"
specialize vp9_build_intra_predictors_mbuv;

prototype void vp9_build_intra_predictors_mbuv_s "struct macroblockd *x"
specialize vp9_build_intra_predictors_mbuv_s;

prototype void vp9_build_comp_intra_predictors_mbuv "struct macroblockd *x"
specialize vp9_build_comp_intra_predictors_mbuv;

prototype void vp9_intra4x4_predict "struct blockd *x, int b_mode, unsigned char *predictor"
specialize vp9_intra4x4_predict;

prototype void vp9_comp_intra4x4_predict "struct blockd *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp9_comp_intra4x4_predict;

prototype void vp9_intra8x8_predict "struct blockd *x, int b_mode, unsigned char *predictor"
specialize vp9_intra8x8_predict;

prototype void vp9_comp_intra8x8_predict "struct blockd *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp9_comp_intra8x8_predict;

prototype void vp9_intra_uv4x4_predict "struct blockd *x, int b_mode, unsigned char *predictor"
specialize vp9_intra_uv4x4_predict;

prototype void vp9_comp_intra_uv4x4_predict "struct blockd *x, int b_mode, int second_mode, unsigned char *predictor"
specialize vp9_comp_intra_uv4x4_predict;

#
# Loopfilter
#
prototype void vp9_loop_filter_mbv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_mbv sse2

prototype void vp9_loop_filter_bv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_bv sse2

prototype void vp9_loop_filter_bv8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_bv8x8 sse2

prototype void vp9_loop_filter_mbh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_mbh sse2

prototype void vp9_loop_filter_bh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_bh sse2

prototype void vp9_loop_filter_bh8x8 "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp9_loop_filter_bh8x8 sse2

prototype void vp9_loop_filter_simple_mbv "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp9_loop_filter_simple_mbv mmx sse2 media neon
vp9_loop_filter_simple_mbv_c=vp9_loop_filter_simple_vertical_edge_c
vp9_loop_filter_simple_mbv_mmx=vp9_loop_filter_simple_vertical_edge_mmx
vp9_loop_filter_simple_mbv_sse2=vp9_loop_filter_simple_vertical_edge_sse2
vp9_loop_filter_simple_mbv_media=vp9_loop_filter_simple_vertical_edge_armv6
vp9_loop_filter_simple_mbv_neon=vp9_loop_filter_mbvs_neon

prototype void vp9_loop_filter_simple_mbh "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp9_loop_filter_simple_mbh mmx sse2 media neon
vp9_loop_filter_simple_mbh_c=vp9_loop_filter_simple_horizontal_edge_c
vp9_loop_filter_simple_mbh_mmx=vp9_loop_filter_simple_horizontal_edge_mmx
vp9_loop_filter_simple_mbh_sse2=vp9_loop_filter_simple_horizontal_edge_sse2
vp9_loop_filter_simple_mbh_media=vp9_loop_filter_simple_horizontal_edge_armv6
vp9_loop_filter_simple_mbh_neon=vp9_loop_filter_mbhs_neon

prototype void vp9_loop_filter_simple_bv "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp9_loop_filter_simple_bv mmx sse2 media neon
vp9_loop_filter_simple_bv_c=vp9_loop_filter_bvs_c
vp9_loop_filter_simple_bv_mmx=vp9_loop_filter_bvs_mmx
vp9_loop_filter_simple_bv_sse2=vp9_loop_filter_bvs_sse2
vp9_loop_filter_simple_bv_media=vp9_loop_filter_bvs_armv6
vp9_loop_filter_simple_bv_neon=vp9_loop_filter_bvs_neon

prototype void vp9_loop_filter_simple_bh "unsigned char *y, int ystride, const unsigned char *blimit"
specialize vp9_loop_filter_simple_bh mmx sse2 media neon
vp9_loop_filter_simple_bh_c=vp9_loop_filter_bhs_c
vp9_loop_filter_simple_bh_mmx=vp9_loop_filter_bhs_mmx
vp9_loop_filter_simple_bh_sse2=vp9_loop_filter_bhs_sse2
vp9_loop_filter_simple_bh_media=vp9_loop_filter_bhs_armv6
vp9_loop_filter_simple_bh_neon=vp9_loop_filter_bhs_neon

#
# sad 16x3, 3x16
#
if [ "$CONFIG_NEWBESTREFMV" = "yes" ]; then
prototype unsigned int vp9_sad16x3 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int ref_stride, int max_sad"
specialize vp9_sad16x3 sse2

prototype unsigned int vp9_sad3x16 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int ref_stride, int max_sad"
specialize vp9_sad3x16 sse2
fi

#
# Encoder functions below this point.
#
if [ "$CONFIG_VP9_ENCODER" = "yes" ]; then


# variance
[ $arch = "x86_64" ] && mmx_x86_64=mmx && sse2_x86_64=sse2

prototype unsigned int vp9_variance32x32 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance32x32

prototype unsigned int vp9_variance16x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance16x16 mmx sse2
vp9_variance16x16_sse2=vp9_variance16x16_wmt
vp9_variance16x16_mmx=vp9_variance16x16_mmx

prototype unsigned int vp9_variance16x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance16x8 mmx sse2
vp9_variance16x8_sse2=vp9_variance16x8_wmt
vp9_variance16x8_mmx=vp9_variance16x8_mmx

prototype unsigned int vp9_variance8x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance8x16 mmx sse2
vp9_variance8x16_sse2=vp9_variance8x16_wmt
vp9_variance8x16_mmx=vp9_variance8x16_mmx

prototype unsigned int vp9_variance8x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance8x8 mmx sse2
vp9_variance8x8_sse2=vp9_variance8x8_wmt
vp9_variance8x8_mmx=vp9_variance8x8_mmx

prototype unsigned int vp9_variance4x4 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance4x4 mmx sse2
vp9_variance4x4_sse2=vp9_variance4x4_wmt
vp9_variance4x4_mmx=vp9_variance4x4_mmx

prototype unsigned int vp9_sub_pixel_variance32x32 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance32x32

prototype unsigned int vp9_sub_pixel_variance16x16 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance16x16 sse2 mmx ssse3
vp9_sub_pixel_variance16x16_sse2=vp9_sub_pixel_variance16x16_wmt

prototype unsigned int vp9_sub_pixel_variance8x16 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance8x16 sse2 mmx
vp9_sub_pixel_variance8x16_sse2=vp9_sub_pixel_variance8x16_wmt

prototype unsigned int vp9_sub_pixel_variance16x8 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance16x8 sse2 mmx ssse3
vp9_sub_pixel_variance16x8_sse2=vp9_sub_pixel_variance16x8_ssse3;
vp9_sub_pixel_variance16x8_sse2=vp9_sub_pixel_variance16x8_wmt

prototype unsigned int vp9_sub_pixel_variance8x8 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance8x8 sse2 mmx
vp9_sub_pixel_variance8x8_sse2=vp9_sub_pixel_variance8x8_wmt

prototype unsigned int vp9_sub_pixel_variance4x4 "const unsigned char *src_ptr, int source_stride, int xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_variance4x4 sse2 mmx
vp9_sub_pixel_variance4x4_sse2=vp9_sub_pixel_variance4x4_wmt

prototype unsigned int vp9_sad32x32 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad32x32

prototype unsigned int vp9_sad16x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad16x16 mmx sse2 sse3
vp9_sad16x16_sse2=vp9_sad16x16_wmt

prototype unsigned int vp9_sad16x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad16x8 mmx sse2
vp9_sad16x8_sse2=vp9_sad16x8_wmt

prototype unsigned int vp9_sad8x16 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad8x16 mmx sse2
vp9_sad8x16_sse2=vp9_sad8x16_wmt

prototype unsigned int vp9_sad8x8 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad8x8 mmx sse2
vp9_sad8x8_sse2=vp9_sad8x8_wmt

prototype unsigned int vp9_sad4x4 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int max_sad"
specialize vp9_sad4x4 mmx sse2
vp9_sad4x4_sse2=vp9_sad4x4_wmt

prototype unsigned int vp9_variance_halfpixvar16x16_h "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_h mmx sse2
vp9_variance_halfpixvar16x16_h_sse2=vp9_variance_halfpixvar16x16_h_wmt

prototype unsigned int vp9_variance_halfpixvar16x16_v "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_v mmx sse2
vp9_variance_halfpixvar16x16_v_sse2=vp9_variance_halfpixvar16x16_v_wmt

prototype unsigned int vp9_variance_halfpixvar16x16_hv "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar16x16_hv mmx sse2
vp9_variance_halfpixvar16x16_hv_sse2=vp9_variance_halfpixvar16x16_hv_wmt

prototype unsigned int vp9_variance_halfpixvar32x32_h "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_h

prototype unsigned int vp9_variance_halfpixvar32x32_v "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_v

prototype unsigned int vp9_variance_halfpixvar32x32_hv "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, unsigned int *sse"
specialize vp9_variance_halfpixvar32x32_hv

prototype void vp9_sad32x32x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x32x3

prototype void vp9_sad16x16x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x16x3 sse3 ssse3

prototype void vp9_sad16x8x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x8x3 sse3 ssse3

prototype void vp9_sad8x16x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x16x3 sse3

prototype void vp9_sad8x8x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x8x3 sse3

prototype void vp9_sad4x4x3 "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned int *sad_array"
specialize vp9_sad4x4x3 sse3

prototype void vp9_sad32x32x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad32x32x8

prototype void vp9_sad16x16x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad16x16x8 sse4

prototype void vp9_sad16x8x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad16x8x8 sse4

prototype void vp9_sad8x16x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad8x16x8 sse4

prototype void vp9_sad8x8x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad8x8x8 sse4

prototype void vp9_sad4x4x8 "const unsigned char *src_ptr, int  src_stride, const unsigned char *ref_ptr, int  ref_stride, unsigned short *sad_array"
specialize vp9_sad4x4x8 sse4

prototype void vp9_sad32x32x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad32x32x4d

prototype void vp9_sad16x16x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x16x4d sse3

prototype void vp9_sad16x8x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad16x8x4d sse3

prototype void vp9_sad8x16x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x16x4d sse3

prototype void vp9_sad8x8x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad8x8x4d sse3

prototype void vp9_sad4x4x4d "const unsigned char *src_ptr, int  src_stride, unsigned char *ref_ptr[], int  ref_stride, unsigned int *sad_array"
specialize vp9_sad4x4x4d sse3

#
# Block copy
#
case $arch in
    x86*)
    prototype void vp9_copy32xn "const unsigned char *src_ptr, int source_stride, const unsigned char *ref_ptr, int ref_stride, int n"
    specialize vp9_copy32xn sse2 sse3
    ;;
esac

prototype unsigned int vp9_sub_pixel_mse16x16 "const unsigned char  *src_ptr, int  src_pixels_per_line, int  xoffset, int  yoffset, const unsigned char *dst_ptr, int dst_pixels_per_line, unsigned int *sse"
specialize vp9_sub_pixel_mse16x16 sse2 mmx
vp9_sub_pixel_mse16x16_sse2=vp9_sub_pixel_mse16x16_wmt

prototype unsigned int vp9_mse16x16 "const unsigned char *src_ptr, int  source_stride, const unsigned char *ref_ptr, int  recon_stride, unsigned int *sse"
specialize vp9_mse16x16 mmx sse2
vp9_mse16x16_sse2=vp9_mse16x16_wmt

prototype unsigned int vp9_sub_pixel_mse32x32 "const unsigned char  *src_ptr, int  source_stride, int  xoffset, int  yoffset, const unsigned char *ref_ptr, int Refstride, unsigned int *sse"
specialize vp9_sub_pixel_mse32x32

prototype unsigned int vp9_get_mb_ss "const short *"
specialize vp9_get_mb_ss mmx sse2
# ENCODEMB INVOKE
prototype int vp9_mbblock_error "struct macroblock *mb, int dc"
specialize vp9_mbblock_error mmx sse2
vp9_mbblock_error_sse2=vp9_mbblock_error_xmm

prototype int vp9_block_error "short *coeff, short *dqcoeff, int block_size"
specialize vp9_block_error mmx sse2
vp9_block_error_sse2=vp9_block_error_xmm

prototype void vp9_subtract_b "struct block *be, struct blockd *bd, int pitch"
specialize vp9_subtract_b mmx sse2

prototype int vp9_mbuverror "struct macroblock *mb"
specialize vp9_mbuverror mmx sse2
vp9_mbuverror_sse2=vp9_mbuverror_xmm

prototype void vp9_subtract_b "struct block *be, struct blockd *bd, int pitch"
specialize vp9_subtract_b mmx sse2

prototype void vp9_subtract_mby "short *diff, unsigned char *src, unsigned char *pred, int stride"
specialize vp9_subtract_mby mmx sse2

prototype void vp9_subtract_mbuv "short *diff, unsigned char *usrc, unsigned char *vsrc, unsigned char *pred, int stride"
specialize vp9_subtract_mbuv mmx sse2

#
# Structured Similarity (SSIM)
#
if [ "$CONFIG_INTERNAL_STATS" = "yes" ]; then
    [ $arch = "x86_64" ] && sse2_on_x86_64=sse2

    prototype void vp9_ssim_parms_8x8 "unsigned char *s, int sp, unsigned char *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp9_ssim_parms_8x8 $sse2_on_x86_64

    prototype void vp9_ssim_parms_16x16 "unsigned char *s, int sp, unsigned char *r, int rp, unsigned long *sum_s, unsigned long *sum_r, unsigned long *sum_sq_s, unsigned long *sum_sq_r, unsigned long *sum_sxr"
    specialize vp9_ssim_parms_16x16 $sse2_on_x86_64
fi

# fdct functions
prototype void vp9_fht "const short *input, int pitch, short *output, int tx_type, int tx_dim"
specialize vp9_fht

prototype void vp9_short_fdct8x8 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_fdct8x8

prototype void vp9_short_fhaar2x2 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_fhaar2x2

prototype void vp9_short_fdct4x4 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_fdct4x4

prototype void vp9_short_fdct8x4 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_fdct8x4

prototype void vp9_short_walsh4x4 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_walsh4x4

prototype void vp9_short_fdct16x16 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_fdct16x16

prototype void vp9_short_walsh4x4_lossless "short *InputData, short *OutputData, int pitch"
specialize vp9_short_walsh4x4_lossless

prototype void vp9_short_walsh4x4_x8 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_walsh4x4_x8

prototype void vp9_short_walsh8x4_x8 "short *InputData, short *OutputData, int pitch"
specialize vp9_short_walsh8x4_x8

fi
# end encoder functions

# Scaler functions
if [ "CONFIG_SPATIAL_RESAMPLING" != "yes" ]; then
    prototype void vp8_horizontal_line_4_5_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_4_5_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_last_vertical_band_4_5_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_2_3_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_2_3_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_last_vertical_band_2_3_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_3_5_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_3_5_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_last_vertical_band_3_5_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_3_4_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_3_4_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_last_vertical_band_3_4_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_1_2_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_1_2_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_last_vertical_band_1_2_scale "unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_5_4_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_5_4_scale "unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_5_3_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_5_3_scale "unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_horizontal_line_2_1_scale "const unsigned char *source, unsigned int source_width, unsigned char *dest, unsigned int dest_width"
    prototype void vp8_vertical_band_2_1_scale "unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
    prototype void vp8_vertical_band_2_1_scale_i "unsigned char *source, unsigned int src_pitch, unsigned char *dest, unsigned int dest_pitch, unsigned int dest_width"
fi

prototype void vp8_yv12_extend_frame_borders "struct yv12_buffer_config *ybf"
specialize vp8_yv12_extend_frame_borders neon

prototype void vp8_yv12_copy_frame "struct yv12_buffer_config *src_ybc, struct yv12_buffer_config *dst_ybc"
specialize vp8_yv12_copy_frame neon

prototype void vp8_yv12_copy_y "struct yv12_buffer_config *src_ybc, struct yv12_buffer_config *dst_ybc"
specialize vp8_yv12_copy_y neon

