##
## Copyright (c) 2015 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

DSP_SRCS-yes += vpx_dsp.mk
DSP_SRCS-yes += vpx_dsp_common.h

DSP_SRCS-$(HAVE_MSA)    += mips/macros_msa.h

# bit reader
DSP_SRCS-yes += vp9_prob.h
DSP_SRCS-yes += vp9_prob.c

ifeq ($(CONFIG_DECODERS),yes)
DSP_SRCS-yes += vp9_reader.h
DSP_SRCS-yes += vp9_reader.c
DSP_SRCS-yes += vp9_read_bit_buffer.c
DSP_SRCS-yes += vp9_read_bit_buffer.h
endif

# loop filters
DSP_SRCS-yes += loopfilter.c

DSP_SRCS-$(ARCH_X86)$(ARCH_X86_64)   += x86/loopfilter_sse2.c
DSP_SRCS-$(HAVE_AVX2)                += x86/loopfilter_avx2.c
DSP_SRCS-$(HAVE_MMX)                 += x86/loopfilter_mmx.asm

DSP_SRCS-$(HAVE_NEON)   += arm/loopfilter_neon.c
ifeq ($(HAVE_NEON_ASM),yes)
DSP_SRCS-yes  += arm/loopfilter_mb_neon$(ASM)
DSP_SRCS-yes  += arm/loopfilter_16_neon$(ASM)
DSP_SRCS-yes  += arm/loopfilter_8_neon$(ASM)
DSP_SRCS-yes  += arm/loopfilter_4_neon$(ASM)
else
ifeq ($(HAVE_NEON),yes)
DSP_SRCS-yes   += arm/loopfilter_16_neon.c
DSP_SRCS-yes   += arm/loopfilter_8_neon.c
DSP_SRCS-yes   += arm/loopfilter_4_neon.c
endif  # HAVE_NEON
endif  # HAVE_NEON_ASM

DSP_SRCS-$(HAVE_MSA)    += mips/loopfilter_msa.h
DSP_SRCS-$(HAVE_MSA)    += mips/loopfilter_16_msa.c
DSP_SRCS-$(HAVE_MSA)    += mips/loopfilter_8_msa.c
DSP_SRCS-$(HAVE_MSA)    += mips/loopfilter_4_msa.c

ifeq ($(CONFIG_VP9_HIGHBITDEPTH),yes)
DSP_SRCS-$(HAVE_SSE2)   += x86/highbd_loopfilter_sse2.c
endif  # CONFIG_VP9_HIGHBITDEPTH

ifeq ($(CONFIG_ENCODERS),yes)
DSP_SRCS-yes            += sad.c
DSP_SRCS-yes            += subtract.c

DSP_SRCS-$(HAVE_MEDIA)  += arm/sad_media$(ASM)
DSP_SRCS-$(HAVE_NEON)   += arm/sad4d_neon.c
DSP_SRCS-$(HAVE_NEON)   += arm/sad_neon.c
DSP_SRCS-$(HAVE_NEON)   += arm/subtract_neon.c

DSP_SRCS-$(HAVE_MSA)    += mips/sad_msa.c
DSP_SRCS-$(HAVE_MSA)    += mips/subtract_msa.c

DSP_SRCS-$(HAVE_MMX)    += x86/sad_mmx.asm
DSP_SRCS-$(HAVE_SSE3)   += x86/sad_sse3.asm
DSP_SRCS-$(HAVE_SSSE3)  += x86/sad_ssse3.asm
DSP_SRCS-$(HAVE_SSE4_1) += x86/sad_sse4.asm
DSP_SRCS-$(HAVE_AVX2)   += x86/sad4d_avx2.c
DSP_SRCS-$(HAVE_AVX2)   += x86/sad_avx2.c

ifeq ($(CONFIG_USE_X86INC),yes)
DSP_SRCS-$(HAVE_SSE2)   += x86/sad4d_sse2.asm
DSP_SRCS-$(HAVE_SSE2)   += x86/sad_sse2.asm
DSP_SRCS-$(HAVE_SSE2)   += x86/subtract_sse2.asm

ifeq ($(CONFIG_VP9_HIGHBITDEPTH),yes)
DSP_SRCS-$(HAVE_SSE2) += x86/highbd_sad4d_sse2.asm
DSP_SRCS-$(HAVE_SSE2) += x86/highbd_sad_sse2.asm
endif  # CONFIG_VP9_HIGHBITDEPTH
endif  # CONFIG_USE_X86INC

endif  # CONFIG_ENCODERS

ifneq ($(filter yes,$(CONFIG_ENCODERS) $(CONFIG_POSTPROC) $(CONFIG_VP9_POSTPROC)),)
DSP_SRCS-yes            += variance.c
DSP_SRCS-yes            += variance.h

DSP_SRCS-$(HAVE_MEDIA)  += arm/bilinear_filter_media$(ASM)
DSP_SRCS-$(HAVE_MEDIA)  += arm/subpel_variance_media.c
DSP_SRCS-$(HAVE_MEDIA)  += arm/variance_halfpixvar16x16_h_media$(ASM)
DSP_SRCS-$(HAVE_MEDIA)  += arm/variance_halfpixvar16x16_hv_media$(ASM)
DSP_SRCS-$(HAVE_MEDIA)  += arm/variance_halfpixvar16x16_v_media$(ASM)
DSP_SRCS-$(HAVE_MEDIA)  += arm/variance_media$(ASM)
DSP_SRCS-$(HAVE_NEON)   += arm/subpel_variance_neon.c
DSP_SRCS-$(HAVE_NEON)   += arm/variance_neon.c

DSP_SRCS-$(HAVE_MSA)    += mips/variance_msa.c
DSP_SRCS-$(HAVE_MSA)    += mips/sub_pixel_variance_msa.c

DSP_SRCS-$(HAVE_MMX)    += x86/variance_mmx.c
DSP_SRCS-$(HAVE_MMX)    += x86/variance_impl_mmx.asm
DSP_SRCS-$(HAVE_SSE2)   += x86/variance_sse2.c  # Contains SSE2 and SSSE3
DSP_SRCS-$(HAVE_AVX2)   += x86/variance_avx2.c
DSP_SRCS-$(HAVE_AVX2)   += x86/variance_impl_avx2.c

ifeq ($(CONFIG_USE_X86INC),yes)
DSP_SRCS-$(HAVE_SSE2)   += x86/subpel_variance_sse2.asm  # Contains SSE2 and SSSE3
endif  # CONFIG_USE_X86INC

ifeq ($(CONFIG_VP9_HIGHBITDEPTH),yes)
DSP_SRCS-$(HAVE_SSE2)   += x86/highbd_variance_sse2.c
DSP_SRCS-$(HAVE_SSE2)   += x86/highbd_variance_impl_sse2.asm
ifeq ($(CONFIG_USE_X86INC),yes)
DSP_SRCS-$(HAVE_SSE2)   += x86/highbd_subpel_variance_impl_sse2.asm
endif  # CONFIG_USE_X86INC
endif  # CONFIG_VP9_HIGHBITDEPTH
endif  # CONFIG_ENCODERS || CONFIG_POSTPROC || CONFIG_VP9_POSTPROC

DSP_SRCS-no += $(DSP_SRCS_REMOVE-yes)

DSP_SRCS-yes += vpx_dsp_rtcd.c
DSP_SRCS-yes += vpx_dsp_rtcd_defs.pl

$(eval $(call rtcd_h_template,vpx_dsp_rtcd,vpx_dsp/vpx_dsp_rtcd_defs.pl))
