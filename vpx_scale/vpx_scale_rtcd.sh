vpx_scale_forward_decls() {
cat <<EOF
struct yv12_buffer_config;
EOF
}
forward_decls vpx_scale_forward_decls

# Scaler functions
if [ "CONFIG_SPATIAL_RESAMPLING" != "yes" ]; then
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

prototype void vp8_yv12_copy_frame "const struct yv12_buffer_config *src_ybc, struct yv12_buffer_config *dst_ybc"
specialize vp8_yv12_copy_frame neon

prototype void vpx_yv12_copy_y "const struct yv12_buffer_config *src_ybc, struct yv12_buffer_config *dst_ybc"
specialize vpx_yv12_copy_y neon

if [ "$CONFIG_VP9" = "yes" ]; then
    prototype void vp9_extend_frame_borders "struct yv12_buffer_config *ybf, int subsampling_x, int subsampling_y"
    specialize vp9_extend_frame_borders dspr2

    prototype void vp9_extend_frame_inner_borders "struct yv12_buffer_config *ybf, int subsampling_x, int subsampling_y"
    specialize vp9_extend_frame_inner_borders dspr2
fi
