common_forward_decls() {
cat <<EOF
struct blockd;
struct macroblockd;
struct loop_filter_info;
EOF
}
forward_decls common_forward_decls

#
# Dequant
#
prototype void vp8_dequantize_b "struct blockd*, short *dqc"
specialize vp8_dequantize_b mmx media neon
vp8_dequantize_b_media=vp8_dequantize_b_v6

prototype void vp8_dequant_idct_add "short *input, short *dq, unsigned char *output, int stride"
specialize vp8_dequant_idct_add mmx media neon
vp8_dequant_idct_add_media=vp8_dequant_idct_add_v6

prototype void vp8_dequant_idct_add_y_block "short *q, short *dq, unsigned char *dst, int stride, char *eobs"
specialize vp8_dequant_idct_add_y_block mmx sse2 media neon
vp8_dequant_idct_add_y_block_media=vp8_dequant_idct_add_y_block_v6

prototype void vp8_dequant_idct_add_uv_block "short *q, short *dq, unsigned char *dst_u, unsigned char *dst_v, int stride, char *eobs"
specialize vp8_dequant_idct_add_uv_block mmx sse2 media neon
vp8_dequant_idct_add_uv_block_media=vp8_dequant_idct_add_uv_block_v6

#
# Loopfilter
#
prototype void vp8_loop_filter_mbv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_mbv mmx sse2 media neon
vp8_loop_filter_mbv_media=vp8_loop_filter_mbv_armv6

prototype void vp8_loop_filter_bv "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bv mmx sse2 media neon
vp8_loop_filter_bv_media=vp8_loop_filter_bv_armv6

prototype void vp8_loop_filter_mbh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_mbh mmx sse2 media neon
vp8_loop_filter_mbh_media=vp8_loop_filter_mbh_armv6

prototype void vp8_loop_filter_bh "unsigned char *y, unsigned char *u, unsigned char *v, int ystride, int uv_stride, struct loop_filter_info *lfi"
specialize vp8_loop_filter_bh mmx sse2 media neon
vp8_loop_filter_bh_media=vp8_loop_filter_bh_armv6


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
# IDCT
#
#idct16
prototype void vp8_short_idct4x4llm "short *input, unsigned char *pred, int pitch, unsigned char *dst, int dst_stride"
specialize vp8_short_idct4x4llm mmx media neon
vp8_short_idct4x4llm_media=vp8_short_idct4x4llm_v6_dual

#iwalsh1
prototype void vp8_short_inv_walsh4x4_1 "short *input, short *output"
# no asm yet

#iwalsh16
prototype void vp8_short_inv_walsh4x4 "short *input, short *output"
specialize vp8_short_inv_walsh4x4 mmx sse2 media neon
vp8_short_inv_walsh4x4_media=vp8_short_inv_walsh4x4_v6

#idct1_scalar_add
prototype void vp8_dc_only_idct_add "short input, unsigned char *pred, int pred_stride, unsigned char *dst, int dst_stride"
specialize vp8_dc_only_idct_add	mmx media neon
vp8_dc_only_idct_add_media=vp8_dc_only_idct_add_v6

#
# RECON
#
prototype void vp8_copy_mem16x16 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem16x16 mmx sse2 media neon
vp8_copy_mem16x16_media=vp8_copy_mem16x16_v6

prototype void vp8_copy_mem8x8 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem8x8 mmx media neon
vp8_copy_mem8x8_media=vp8_copy_mem8x8_v6

prototype void vp8_copy_mem8x4 "unsigned char *src, int src_pitch, unsigned char *dst, int dst_pitch"
specialize vp8_copy_mem8x4 mmx media neon
vp8_copy_mem8x4_media=vp8_copy_mem8x4_v6

prototype void vp8_build_intra_predictors_mby "struct macroblockd *x"
specialize vp8_build_intra_predictors_mby sse2 ssse3 neon

prototype void vp8_build_intra_predictors_mby_s "struct macroblockd *x"
specialize vp8_build_intra_predictors_mby_s sse2 ssse3 neon

prototype void vp8_build_intra_predictors_mbuv "struct macroblockd *x"
specialize vp8_build_intra_predictors_mbuv sse2 ssse3

prototype void vp8_build_intra_predictors_mbuv_s "struct macroblockd *x"
specialize vp8_build_intra_predictors_mbuv_s sse2 ssse3

prototype void vp8_intra4x4_predict "unsigned char *src, int src_stride, int b_mode, unsigned char *dst, int dst_stride"
specialize vp8_intra4x4_predict media
vp8_intra4x4_predict_media=vp8_intra4x4_predict_armv6

#
# Postproc
#
if [ "$CONFIG_POSTPROC" = "yes" ]; then
    prototype void vp8_mbpost_proc_down "unsigned char *dst, int pitch, int rows, int cols,int flimit"
    specialize vp8_mbpost_proc_down mmx sse2
    vp8_mbpost_proc_down_sse2=vp8_mbpost_proc_down_xmm

    prototype void vp8_mbpost_proc_across_ip "unsigned char *dst, int pitch, int rows, int cols,int flimit"
    specialize vp8_mbpost_proc_across_ip sse2
    vp8_mbpost_proc_across_ip_sse2=vp8_mbpost_proc_across_ip_xmm

    prototype void vp8_post_proc_down_and_across "unsigned char *src, unsigned char *dst, int src_pitch, int dst_pitch, int rows, int cols, int flimit"
    specialize vp8_post_proc_down_and_across mmx sse2
    vp8_post_proc_down_and_across_sse2=vp8_post_proc_down_and_across_xmm

    prototype void vp8_plane_add_noise "unsigned char *s, char *noise, char blackclamp[16], char whiteclamp[16], char bothclamp[16], unsigned int w, unsigned int h, int pitch"
    specialize vp8_plane_add_noise mmx sse2
    vp8_plane_add_noise_sse2=vp8_plane_add_noise_wmt

    prototype void vp8_blend_mb_inner "unsigned char *y, unsigned char *u, unsigned char *v, int y1, int u1, int v1, int alpha, int stride"
    # no asm yet

    prototype void vp8_blend_mb_outer "unsigned char *y, unsigned char *u, unsigned char *v, int y1, int u1, int v1, int alpha, int stride"
    # no asm yet

    prototype void vp8_blend_b "unsigned char *y, unsigned char *u, unsigned char *v, int y1, int u1, int v1, int alpha, int stride"
    # no asm yet
fi
