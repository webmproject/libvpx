##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

VP8_COMMON_SRCS-yes += vp8_common.mk
VP8_COMMON_SRCS-yes += common/pragmas.h
VP8_COMMON_SRCS-yes += common/ppflags.h
VP8_COMMON_SRCS-yes += common/onyx.h
VP8_COMMON_SRCS-yes += common/onyxd.h
VP8_COMMON_SRCS-yes += common/alloccommon.c
VP8_COMMON_SRCS-yes += common/asm_com_offsets.c
VP8_COMMON_SRCS-yes += common/blockd.c
VP8_COMMON_SRCS-yes += common/coefupdateprobs.h
VP8_COMMON_SRCS-yes += common/debugmodes.c
VP8_COMMON_SRCS-yes += common/default_coef_probs.h
VP8_COMMON_SRCS-yes += common/dequantize.c
VP8_COMMON_SRCS-yes += common/dequantize.h
VP8_COMMON_SRCS-yes += common/entropy.c
VP8_COMMON_SRCS-yes += common/entropymode.c
VP8_COMMON_SRCS-yes += common/entropymv.c
VP8_COMMON_SRCS-yes += common/extend.c
VP8_COMMON_SRCS-yes += common/filter.c
VP8_COMMON_SRCS-yes += common/filter.h
VP8_COMMON_SRCS-yes += common/findnearmv.c
VP8_COMMON_SRCS-yes += common/generic/systemdependent.c
VP8_COMMON_SRCS-yes += common/idct_blk.c
VP8_COMMON_SRCS-yes += common/idctllm.c
VP8_COMMON_SRCS-yes += common/alloccommon.h
VP8_COMMON_SRCS-yes += common/blockd.h
VP8_COMMON_SRCS-yes += common/common.h
VP8_COMMON_SRCS-yes += common/entropy.h
VP8_COMMON_SRCS-yes += common/entropymode.h
VP8_COMMON_SRCS-yes += common/entropymv.h
VP8_COMMON_SRCS-yes += common/extend.h
VP8_COMMON_SRCS-yes += common/findnearmv.h
VP8_COMMON_SRCS-yes += common/header.h
VP8_COMMON_SRCS-yes += common/idct.h
VP8_COMMON_SRCS-yes += common/invtrans.h
VP8_COMMON_SRCS-yes += common/loopfilter.h
VP8_COMMON_SRCS-yes += common/modecont.h
VP8_COMMON_SRCS-yes += common/mv.h
VP8_COMMON_SRCS-yes += common/onyxc_int.h
VP8_COMMON_SRCS-yes += common/quant_common.h
VP8_COMMON_SRCS-yes += common/recon.h
VP8_COMMON_SRCS-yes += common/reconinter.h
VP8_COMMON_SRCS-yes += common/reconintra.h
VP8_COMMON_SRCS-yes += common/reconintra4x4.h
VP8_COMMON_SRCS-yes += common/setupintrarecon.h
VP8_COMMON_SRCS-yes += common/subpixel.h
VP8_COMMON_SRCS-yes += common/swapyv12buffer.h
VP8_COMMON_SRCS-yes += common/systemdependent.h
VP8_COMMON_SRCS-yes += common/threading.h
VP8_COMMON_SRCS-yes += common/treecoder.h
VP8_COMMON_SRCS-yes += common/loopfilter.c
VP8_COMMON_SRCS-yes += common/loopfilter_filters.c
VP8_COMMON_SRCS-yes += common/mbpitch.c
VP8_COMMON_SRCS-yes += common/modecont.c
VP8_COMMON_SRCS-yes += common/modecontext.c
VP8_COMMON_SRCS-yes += common/quant_common.c
VP8_COMMON_SRCS-yes += common/reconinter.c
VP8_COMMON_SRCS-yes += common/reconintra.c
VP8_COMMON_SRCS-yes += common/reconintra4x4.c
VP8_COMMON_SRCS-yes += common/setupintrarecon.c
VP8_COMMON_SRCS-yes += common/swapyv12buffer.c



VP8_COMMON_SRCS-$(CONFIG_POSTPROC_VISUALIZER) += common/textblit.c
VP8_COMMON_SRCS-yes += common/treecoder.c

VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/dequantize_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/filter_x86.c
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/filter_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/idct_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/subpixel_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/recon_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/loopfilter_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/postproc_x86.h
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/x86_systemdependent.c
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/vp8_asm_stubs.c
VP8_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/loopfilter_x86.c
VP8_COMMON_SRCS-$(CONFIG_POSTPROC) += common/postproc.h
VP8_COMMON_SRCS-$(CONFIG_POSTPROC) += common/postproc.c
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/dequantize_mmx.asm
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/idct_blk_mmx.c
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/idctllm_mmx.asm
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/iwalsh_mmx.asm
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/recon_mmx.asm
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/subpixel_mmx.asm
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/loopfilter_mmx.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/idct_blk_sse2.c
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/idctllm_sse2.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/recon_sse2.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/recon_wrapper_sse2.c
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/subpixel_sse2.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/loopfilter_sse2.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/iwalsh_sse2.asm
VP8_COMMON_SRCS-$(HAVE_SSSE3) += common/x86/subpixel_ssse3.asm
ifeq ($(CONFIG_POSTPROC),yes)
VP8_COMMON_SRCS-$(HAVE_MMX) += common/x86/postproc_mmx.asm
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/postproc_sse2.asm
endif
ifeq ($(ARCH_X86_64),yes)
VP8_COMMON_SRCS-$(HAVE_SSE2) += common/x86/loopfilter_block_sse2.asm
endif

# common (c)
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/arm_systemdependent.c
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/filter_arm.c
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/idct_arm.h
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/loopfilter_arm.c
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/loopfilter_arm.h
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/recon_arm.h
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/reconintra_arm.c
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/subpixel_arm.h
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/dequantize_arm.c
VP8_COMMON_SRCS-$(ARCH_ARM)  += common/arm/dequantize_arm.h

# common (armv6)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/bilinearfilter_arm.c
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/bilinearfilter_arm.h
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/bilinearfilter_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem8x4_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem8x8_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem16x16_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/dc_only_idct_add_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/iwalsh_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/filter_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/idct_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/loopfilter_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/simpleloopfilter_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/sixtappredict8x4_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/intra4x4_predict_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/dequant_idct_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/dequantize_v6$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/idct_blk_v6.c

# common (neon)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict4x4_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict8x4_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict8x8_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict16x16_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem8x4_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem8x8_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem16x16_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/dc_only_idct_add_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/iwalsh_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfilter_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfiltersimplehorizontaledge_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfiltersimpleverticaledge_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/mbloopfilter_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/shortidct4x4llm_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict4x4_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict8x4_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict8x8_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict16x16_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/buildintrapredictorsmby_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/save_neon_reg$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/dequant_idct_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/idct_dequant_full_2x_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/idct_dequant_0_2x_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/dequantizeb_neon$(ASM)
VP8_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/idct_blk_neon.c
