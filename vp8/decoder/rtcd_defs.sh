decoder_forward_decls() {
cat <<EOF
struct blockd;
EOF
}
forward_decls decoder_forward_decls

prototype void vp8_dequantize_b "struct blockd*"
specialize vp8_dequantize_b mmx v6 neon

prototype void vp8_dequant_idct_add "short *input, short *dq, unsigned char *output, int stride"
specialize vp8_dequant_idct_add mmx v6 neon

prototype void vp8_dequant_idct_add_y_block "short *q, short *dq, unsigned char *dst, int stride, char *eobs"
specialize vp8_dequant_idct_add_y_block mmx sse2 v6 neon

prototype void vp8_dequant_idct_add_uv_block "short *q, short *dq, unsigned char *dst_u, unsigned char *dst_v, int stride, char *eobs"
specialize vp8_dequant_idct_add_uv_block mmx sse2 v6 neon
