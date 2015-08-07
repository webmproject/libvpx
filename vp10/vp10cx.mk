##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

VP10_CX_EXPORTS += exports_enc

VP10_CX_SRCS-yes += $(VP10_COMMON_SRCS-yes)
VP10_CX_SRCS-no  += $(VP10_COMMON_SRCS-no)
VP10_CX_SRCS_REMOVE-yes += $(VP10_COMMON_SRCS_REMOVE-yes)
VP10_CX_SRCS_REMOVE-no  += $(VP10_COMMON_SRCS_REMOVE-no)

VP10_CX_SRCS-yes += vp10_cx_iface.c

VP10_CX_SRCS-yes += encoder/vp9_avg.c
VP10_CX_SRCS-yes += encoder/vp9_bitstream.c
VP10_CX_SRCS-yes += encoder/vp9_context_tree.c
VP10_CX_SRCS-yes += encoder/vp9_context_tree.h
VP10_CX_SRCS-yes += encoder/vp9_cost.h
VP10_CX_SRCS-yes += encoder/vp9_cost.c
VP10_CX_SRCS-yes += encoder/vp9_dct.c
VP10_CX_SRCS-$(CONFIG_VP9_TEMPORAL_DENOISING) += encoder/vp9_denoiser.c
VP10_CX_SRCS-$(CONFIG_VP9_TEMPORAL_DENOISING) += encoder/vp9_denoiser.h
VP10_CX_SRCS-yes += encoder/vp9_encodeframe.c
VP10_CX_SRCS-yes += encoder/vp9_encodeframe.h
VP10_CX_SRCS-yes += encoder/vp9_encodemb.c
VP10_CX_SRCS-yes += encoder/vp9_encodemv.c
VP10_CX_SRCS-yes += encoder/vp9_ethread.h
VP10_CX_SRCS-yes += encoder/vp9_ethread.c
VP10_CX_SRCS-yes += encoder/vp9_extend.c
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/vp9_fastssim.c
VP10_CX_SRCS-yes += encoder/vp9_firstpass.c
VP10_CX_SRCS-yes += encoder/vp9_block.h
VP10_CX_SRCS-yes += encoder/vp9_bitstream.h
VP10_CX_SRCS-yes += encoder/vp9_encodemb.h
VP10_CX_SRCS-yes += encoder/vp9_encodemv.h
VP10_CX_SRCS-yes += encoder/vp9_extend.h
VP10_CX_SRCS-yes += encoder/vp9_firstpass.h
VP10_CX_SRCS-yes += encoder/vp9_lookahead.c
VP10_CX_SRCS-yes += encoder/vp9_lookahead.h
VP10_CX_SRCS-yes += encoder/vp9_mcomp.h
VP10_CX_SRCS-yes += encoder/vp9_encoder.h
VP10_CX_SRCS-yes += encoder/vp9_quantize.h
VP10_CX_SRCS-yes += encoder/vp9_ratectrl.h
VP10_CX_SRCS-yes += encoder/vp9_rd.h
VP10_CX_SRCS-yes += encoder/vp9_rdopt.h
VP10_CX_SRCS-yes += encoder/vp9_pickmode.h
VP10_CX_SRCS-yes += encoder/vp9_svc_layercontext.h
VP10_CX_SRCS-yes += encoder/vp9_tokenize.h
VP10_CX_SRCS-yes += encoder/vp9_treewriter.h
VP10_CX_SRCS-yes += encoder/vp9_mcomp.c
VP10_CX_SRCS-yes += encoder/vp9_encoder.c
VP10_CX_SRCS-yes += encoder/vp9_picklpf.c
VP10_CX_SRCS-yes += encoder/vp9_picklpf.h
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/vp9_psnrhvs.c
VP10_CX_SRCS-yes += encoder/vp9_quantize.c
VP10_CX_SRCS-yes += encoder/vp9_ratectrl.c
VP10_CX_SRCS-yes += encoder/vp9_rd.c
VP10_CX_SRCS-yes += encoder/vp9_rdopt.c
VP10_CX_SRCS-yes += encoder/vp9_pickmode.c
VP10_CX_SRCS-yes += encoder/vp9_segmentation.c
VP10_CX_SRCS-yes += encoder/vp9_segmentation.h
VP10_CX_SRCS-yes += encoder/vp9_speed_features.c
VP10_CX_SRCS-yes += encoder/vp9_speed_features.h
VP10_CX_SRCS-yes += encoder/vp9_subexp.c
VP10_CX_SRCS-yes += encoder/vp9_subexp.h
VP10_CX_SRCS-yes += encoder/vp9_svc_layercontext.c
VP10_CX_SRCS-yes += encoder/vp9_resize.c
VP10_CX_SRCS-yes += encoder/vp9_resize.h
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/vp9_ssim.c
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/vp9_ssim.h
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += encoder/vp9_blockiness.c

VP10_CX_SRCS-yes += encoder/vp9_tokenize.c
VP10_CX_SRCS-yes += encoder/vp9_treewriter.c
VP10_CX_SRCS-yes += encoder/vp9_aq_variance.c
VP10_CX_SRCS-yes += encoder/vp9_aq_variance.h
VP10_CX_SRCS-yes += encoder/vp9_aq_cyclicrefresh.c
VP10_CX_SRCS-yes += encoder/vp9_aq_cyclicrefresh.h
VP10_CX_SRCS-yes += encoder/vp9_aq_complexity.c
VP10_CX_SRCS-yes += encoder/vp9_aq_complexity.h
VP10_CX_SRCS-yes += encoder/vp9_skin_detection.c
VP10_CX_SRCS-yes += encoder/vp9_skin_detection.h
ifeq ($(CONFIG_VP9_POSTPROC),yes)
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += common/vp9_postproc.h
VP10_CX_SRCS-$(CONFIG_INTERNAL_STATS) += common/vp9_postproc.c
endif
VP10_CX_SRCS-yes += encoder/vp9_temporal_filter.c
VP10_CX_SRCS-yes += encoder/vp9_temporal_filter.h
VP10_CX_SRCS-yes += encoder/vp9_mbgraph.c
VP10_CX_SRCS-yes += encoder/vp9_mbgraph.h

VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_avg_intrin_sse2.c
VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_temporal_filter_apply_sse2.asm
VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_quantize_sse2.c
ifeq ($(CONFIG_VP9_HIGHBITDEPTH),yes)
VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_highbd_block_error_intrin_sse2.c
endif

ifeq ($(CONFIG_USE_X86INC),yes)
VP10_CX_SRCS-$(HAVE_MMX) += encoder/x86/vp9_dct_mmx.asm
VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_error_sse2.asm
endif

ifeq ($(ARCH_X86_64),yes)
ifeq ($(CONFIG_USE_X86INC),yes)
VP10_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/vp9_quantize_ssse3_x86_64.asm
VP10_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/vp9_dct_ssse3_x86_64.asm
endif
endif
VP10_CX_SRCS-$(ARCH_X86_64) += encoder/x86/vp9_ssim_opt_x86_64.asm

VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_dct_sse2.c
VP10_CX_SRCS-$(HAVE_SSSE3) += encoder/x86/vp9_dct_ssse3.c

ifeq ($(CONFIG_VP9_TEMPORAL_DENOISING),yes)
VP10_CX_SRCS-$(HAVE_SSE2) += encoder/x86/vp9_denoiser_sse2.c
endif

VP10_CX_SRCS-$(HAVE_AVX2) += encoder/x86/vp9_error_intrin_avx2.c

ifneq ($(CONFIG_VP9_HIGHBITDEPTH),yes)
VP10_CX_SRCS-$(HAVE_NEON) += encoder/arm/neon/vp9_dct_neon.c
VP10_CX_SRCS-$(HAVE_NEON) += encoder/arm/neon/vp9_error_neon.c
endif
VP10_CX_SRCS-$(HAVE_NEON) += encoder/arm/neon/vp9_avg_neon.c
VP10_CX_SRCS-$(HAVE_NEON) += encoder/arm/neon/vp9_quantize_neon.c

VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_avg_msa.c
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_error_msa.c
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_fdct4x4_msa.c
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_fdct8x8_msa.c
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_fdct16x16_msa.c
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_fdct_msa.h
VP10_CX_SRCS-$(HAVE_MSA) += encoder/mips/msa/vp9_temporal_filter_msa.c

VP10_CX_SRCS-yes := $(filter-out $(VP10_CX_SRCS_REMOVE-yes),$(VP10_CX_SRCS-yes))
