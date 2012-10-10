common_forward_decls() {
cat <<EOF
struct blockd;
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
