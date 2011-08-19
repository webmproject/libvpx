common_forward_decls() {
cat <<EOF
struct blockd;
EOF
}
forward_decls common_forward_decls

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
