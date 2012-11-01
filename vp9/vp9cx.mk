##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


include $(SRC_PATH_BARE)/$(VP9_PREFIX)vp9_common.mk

VP9_CX_EXPORTS += exports_enc

VP9_CX_SRCS-yes += $(VP9_COMMON_SRCS-yes)
VP9_CX_SRCS-no  += $(VP9_COMMON_SRCS-no)
VP9_CX_SRCS_REMOVE-yes += $(VP9_COMMON_SRCS_REMOVE-yes)
VP9_CX_SRCS_REMOVE-no  += $(VP9_COMMON_SRCS_REMOVE-no)

ifeq ($(ARCH_ARM),yes)
  include $(SRC_PATH_BARE)/$(VP9_PREFIX)vp9cx_arm.mk
endif

VP9_CX_SRCS-yes += vp9_cx_iface.c

# encoder
#INCLUDES += algo/vpx_common/vpx_mem/include
#INCLUDES += common
#INCLUDES += common
#INCLUDES += common
#INCLUDES += algo/vpx_ref/cpu_id/include
#INCLUDES += common
#INCLUDES += encoder

VP9_CX_SRCS-yes += encoder/asm_enc_offsets.c
VP9_CX_SRCS-yes += encoder/bitstream.c
VP9_CX_SRCS-yes += encoder/boolhuff.c
VP9_CX_SRCS-yes += encoder/dct.c
VP9_CX_SRCS-yes += encoder/encodeframe.c
VP9_CX_SRCS-yes += encoder/encodeintra.c
VP9_CX_SRCS-yes += encoder/encodemb.c
VP9_CX_SRCS-yes += encoder/encodemv.c
VP9_CX_SRCS-yes += encoder/firstpass.c
VP9_CX_SRCS-yes += encoder/generic/csystemdependent.c
VP9_CX_SRCS-yes += encoder/block.h
VP9_CX_SRCS-yes += encoder/boolhuff.h
VP9_CX_SRCS-yes += encoder/bitstream.h
VP9_CX_SRCS-yes += encoder/encodeintra.h
VP9_CX_SRCS-yes += encoder/encodemb.h
VP9_CX_SRCS-yes += encoder/encodemv.h
VP9_CX_SRCS-yes += encoder/firstpass.h
VP9_CX_SRCS-yes += encoder/lookahead.c
VP9_CX_SRCS-yes += encoder/lookahead.h
VP9_CX_SRCS-yes += encoder/mcomp.h
VP9_CX_SRCS-yes += encoder/modecosts.h
VP9_CX_SRCS-yes += encoder/onyx_int.h
VP9_CX_SRCS-yes += encoder/psnr.h
VP9_CX_SRCS-yes += encoder/quantize.h
VP9_CX_SRCS-yes += encoder/ratectrl.h
VP9_CX_SRCS-yes += encoder/rdopt.h
VP9_CX_SRCS-yes += encoder/tokenize.h
VP9_CX_SRCS-yes += encoder/treewriter.h
VP9_CX_SRCS-yes += encoder/variance.h
VP9_CX_SRCS-yes += encoder/mcomp.c
VP9_CX_SRCS-yes += encoder/modecosts.c
VP9_CX_SRCS-yes += encoder/onyx_if.c
VP9_CX_SRCS-yes += encoder/picklpf.c
VP9_CX_SRCS-yes += encoder/psnr.c
VP9_CX_SRCS-yes += encoder/quantize.c
VP9_CX_SRCS-yes += encoder/ratectrl.c
VP9_CX_SRCS-yes += encoder/rdopt.c
VP9_CX_SRCS-yes += encoder/sad_c.c
VP9_CX_SRCS-yes += encoder/satd_c.c
VP9_CX_SRCS-yes += encoder/segmentation.c
VP9_CX_SRCS-yes += encoder/segmentation.h
VP9_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/ssim.c
VP9_CX_SRCS-yes += encoder/tokenize.c
VP9_CX_SRCS-yes += encoder/treewriter.c
VP9_CX_SRCS-yes += encoder/variance_c.c
ifeq ($(CONFIG_POSTPROC),yes)
VP9_CX_SRCS-$(CONFIG_INTERNAL_STATS) += common/postproc.h
VP9_CX_SRCS-$(CONFIG_INTERNAL_STATS) += common/postproc.c
endif
VP9_CX_SRCS-yes += encoder/temporal_filter.c
VP9_CX_SRCS-yes += encoder/temporal_filter.h
VP9_CX_SRCS-yes += encoder/mbgraph.c
VP9_CX_SRCS-yes += encoder/mbgraph.h


VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/mcomp_x86.h
VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/quantize_x86.h
VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/temporal_filter_x86.h
VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/x86_csystemdependent.c
VP9_CX_SRCS-$(HAVE_MMX) += encoder/x86/variance_mmx.c
VP9_CX_SRCS-$(HAVE_MMX) += encoder/x86/variance_impl_mmx.asm
VP9_CX_SRCS-$(HAVE_MMX) += encoder/x86/sad_mmx.asm
VP9_CX_SRCS-$(HAVE_MMX) += encoder/x86/dct_mmx.asm
VP9_CX_SRCS-$(HAVE_MMX) += encoder/x86/subtract_mmx.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/dct_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/variance_sse2.c
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/variance_impl_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/sad_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/fwalsh_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/quantize_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/subtract_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE2) += encoder/x86/temporal_filter_apply_sse2.asm
VP9_CX_SRCS-$(HAVE_SSE3) += encoder/x86/sad_sse3.asm
VP9_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/sad_ssse3.asm
VP9_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/variance_ssse3.c
VP9_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/variance_impl_ssse3.asm
VP9_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/quantize_ssse3.asm
VP9_CX_SRCS-$(HAVE_SSE4_1) += encoder/x86/sad_sse4.asm
VP9_CX_SRCS-$(HAVE_SSE4_1) += encoder/x86/quantize_sse4.asm
VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/quantize_mmx.asm
VP9_CX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += encoder/x86/encodeopt.asm
VP9_CX_SRCS-$(ARCH_X86_64) += encoder/x86/ssim_opt.asm


VP9_CX_SRCS-yes := $(filter-out $(VP9_CX_SRCS_REMOVE-yes),$(VP9_CX_SRCS-yes))
