##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

VP9_COMMON_SRCS-yes += vp9_common.mk
VP9_COMMON_SRCS-yes += common/type_aliases.h
VP9_COMMON_SRCS-yes += common/pragmas.h
VP9_COMMON_SRCS-yes += common/ppflags.h
VP9_COMMON_SRCS-yes += common/onyx.h
VP9_COMMON_SRCS-yes += common/onyxd.h
VP9_COMMON_SRCS-yes += common/alloccommon.c
VP9_COMMON_SRCS-yes += common/asm_com_offsets.c
VP9_COMMON_SRCS-yes += common/blockd.c
VP9_COMMON_SRCS-yes += common/coefupdateprobs.h
VP9_COMMON_SRCS-yes += common/debugmodes.c
VP9_COMMON_SRCS-yes += common/entropy.c
VP9_COMMON_SRCS-yes += common/entropymode.c
VP9_COMMON_SRCS-yes += common/entropymv.c
VP9_COMMON_SRCS-yes += common/extend.c
VP9_COMMON_SRCS-yes += common/filter.c
VP9_COMMON_SRCS-yes += common/filter.h
VP9_COMMON_SRCS-yes += common/findnearmv.c
VP9_COMMON_SRCS-yes += common/generic/systemdependent.c
VP9_COMMON_SRCS-yes += common/idctllm.c
VP9_COMMON_SRCS-yes += common/alloccommon.h
VP9_COMMON_SRCS-yes += common/blockd.h
VP9_COMMON_SRCS-yes += common/common.h
VP9_COMMON_SRCS-yes += common/common_types.h
VP9_COMMON_SRCS-yes += common/entropy.h
VP9_COMMON_SRCS-yes += common/entropymode.h
VP9_COMMON_SRCS-yes += common/entropymv.h
VP9_COMMON_SRCS-yes += common/extend.h
VP9_COMMON_SRCS-yes += common/findnearmv.h
VP9_COMMON_SRCS-yes += common/header.h
VP9_COMMON_SRCS-yes += common/invtrans.h
VP9_COMMON_SRCS-yes += common/loopfilter.h
VP9_COMMON_SRCS-yes += common/modecont.h
VP9_COMMON_SRCS-yes += common/mv.h
VP9_COMMON_SRCS-yes += common/onyxc_int.h
VP9_COMMON_SRCS-yes += common/pred_common.h
VP9_COMMON_SRCS-yes += common/pred_common.c
VP9_COMMON_SRCS-yes += common/quant_common.h
VP9_COMMON_SRCS-yes += common/reconinter.h
VP9_COMMON_SRCS-yes += common/reconintra.h
VP9_COMMON_SRCS-yes += common/reconintra4x4.h
VP9_COMMON_SRCS-yes += common/rtcd.c
VP9_COMMON_SRCS-yes += common/rtcd_defs.sh
VP9_COMMON_SRCS-yes += common/sadmxn.h
VP9_COMMON_SRCS-yes += common/subpelvar.h
VP9_COMMON_SRCS-yes += common/seg_common.h
VP9_COMMON_SRCS-yes += common/seg_common.c
VP9_COMMON_SRCS-yes += common/setupintrarecon.h
VP9_COMMON_SRCS-yes += common/subpixel.h
VP9_COMMON_SRCS-yes += common/swapyv12buffer.h
VP9_COMMON_SRCS-yes += common/systemdependent.h
VP9_COMMON_SRCS-yes += common/treecoder.h
VP9_COMMON_SRCS-yes += common/invtrans.c
VP9_COMMON_SRCS-yes += common/loopfilter.c
VP9_COMMON_SRCS-yes += common/loopfilter_filters.c
VP9_COMMON_SRCS-yes += common/mbpitch.c
VP9_COMMON_SRCS-yes += common/modecont.c
VP9_COMMON_SRCS-yes += common/modecontext.c
VP9_COMMON_SRCS-yes += common/mvref_common.c
VP9_COMMON_SRCS-yes += common/mvref_common.h
VP9_COMMON_SRCS-yes += common/quant_common.c
VP9_COMMON_SRCS-yes += common/recon.c
VP9_COMMON_SRCS-yes += common/reconinter.c
VP9_COMMON_SRCS-yes += common/reconintra.c
VP9_COMMON_SRCS-yes += common/reconintra4x4.c
VP9_COMMON_SRCS-yes += common/setupintrarecon.c
VP9_COMMON_SRCS-yes += common/swapyv12buffer.c
VP9_COMMON_SRCS-$(CONFIG_POSTPROC_VISUALIZER) += common/textblit.c
VP9_COMMON_SRCS-yes += common/treecoder.c
VP9_COMMON_SRCS-$(CONFIG_IMPLICIT_SEGMENTATION) += common/implicit_segmentation.c

VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/idct_x86.h
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/subpixel_x86.h
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/loopfilter_x86.h
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/postproc_x86.h
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/x86_systemdependent.c
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/vp8_asm_stubs.c
VP9_COMMON_SRCS-$(ARCH_X86)$(ARCH_X86_64) += common/x86/loopfilter_x86.c
VP9_COMMON_SRCS-$(CONFIG_POSTPROC) += common/postproc.h
VP9_COMMON_SRCS-$(CONFIG_POSTPROC) += common/postproc.c
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/idctllm_mmx.asm
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/iwalsh_mmx.asm
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/recon_mmx.asm
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/subpixel_mmx.asm
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/loopfilter_mmx.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/idctllm_sse2.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/recon_sse2.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/recon_wrapper_sse2.c
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/subpixel_sse2.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/loopfilter_sse2.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/iwalsh_sse2.asm
VP9_COMMON_SRCS-$(HAVE_SSSE3) += common/x86/subpixel_8t_ssse3.asm
VP9_COMMON_SRCS-$(HAVE_SSSE3) += common/x86/subpixel_ssse3.asm
ifeq ($(CONFIG_POSTPROC),yes)
VP9_COMMON_SRCS-$(HAVE_MMX) += common/x86/postproc_mmx.asm
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/postproc_sse2.asm
endif

# common (c)
ifeq ($(CONFIG_CSM),yes)
VP9_COMMON_SRCS-yes += common/maskingmv.c
VP9_COMMON_SRCS-$(HAVE_SSE3) += common/x86/mask_sse3.asm
endif

VP9_COMMON_SRCS-$(HAVE_SSE4_1) += common/x86/filter_sse4.c
ifeq ($(HAVE_SSE4_1),yes)
$(call xform_obj_path_o_d,vp9/common/x86/filter_sse4.c): CFLAGS += -msse4
endif

VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/filter_sse2.c
VP9_COMMON_SRCS-$(HAVE_SSE2) += common/x86/sadmxn_x86.c
ifeq ($(HAVE_SSE2),yes)
$(call xform_obj_path_o_d,vp9/common/x86/filter_sse2.c): CFLAGS += -msse2
$(call xform_obj_path_o_d,vp9/common/x86/loopfilter_x86.c): CFLAGS += -msse2
$(call xform_obj_path_o_d,vp9/common/x86/sadmxn_x86.c): CFLAGS += -msse2
endif

VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/arm_systemdependent.c
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/bilinearfilter_arm.c
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/bilinearfilter_arm.h
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/filter_arm.c
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/idct_arm.h
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/loopfilter_arm.c
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/loopfilter_arm.h
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/recon_arm.h
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/reconintra_arm.c
VP9_COMMON_SRCS-$(ARCH_ARM)  += common/arm/subpixel_arm.h

# common (armv6)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/bilinearfilter_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem8x4_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem8x8_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/copymem16x16_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/dc_only_idct_add_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/iwalsh_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/filter_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/idct_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/loopfilter_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/recon_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/simpleloopfilter_v6$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV6)  += common/arm/armv6/sixtappredict8x4_v6$(ASM)

# common (neon)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict4x4_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict8x4_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict8x8_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/bilinearpredict16x16_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem8x4_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem8x8_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/copymem16x16_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/dc_only_idct_add_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/iwalsh_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfilter_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfiltersimplehorizontaledge_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/loopfiltersimpleverticaledge_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/mbloopfilter_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/recon2b_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/recon4b_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/reconb_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/shortidct4x4llm_1_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/shortidct4x4llm_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict4x4_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict8x4_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict8x8_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/sixtappredict16x16_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/recon16x16mb_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/buildintrapredictorsmby_neon$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/save_neon_reg$(ASM)
VP9_COMMON_SRCS-$(HAVE_ARMV7)  += common/arm/neon/recon_neon.c


$(eval $(call asm_offsets_template,\
         vp9_asm_com_offsets.asm, $(VP9_PREFIX)common/asm_com_offsets.c))

$(eval $(call rtcd_h_template,vp9_rtcd,vp9/common/rtcd_defs.sh))
