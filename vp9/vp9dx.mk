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

VP9_DX_EXPORTS += exports_dec

VP9_DX_SRCS-yes += $(VP9_COMMON_SRCS-yes)
VP9_DX_SRCS-no  += $(VP9_COMMON_SRCS-no)
VP9_DX_SRCS_REMOVE-yes += $(VP9_COMMON_SRCS_REMOVE-yes)
VP9_DX_SRCS_REMOVE-no  += $(VP9_COMMON_SRCS_REMOVE-no)

ifeq ($(ARCH_ARM),yes)
  include $(SRC_PATH_BARE)/$(VP9_PREFIX)vp9dx_arm.mk
endif

VP9_DX_SRCS-yes += vp9_dx_iface.c

# common
#define ARM
#define DISABLE_THREAD

#INCLUDES += algo/vpx_common/vpx_mem/include
#INCLUDES += common
#INCLUDES += common
#INCLUDES += common
#INCLUDES += common
#INCLUDES += decoder



# decoder
#define ARM
#define DISABLE_THREAD

#INCLUDES += algo/vpx_common/vpx_mem/include
#INCLUDES += common
#INCLUDES += common
#INCLUDES += common
#INCLUDES += common
#INCLUDES += decoder

VP9_DX_SRCS-yes += decoder/asm_dec_offsets.c
VP9_DX_SRCS-yes += decoder/dboolhuff.c
VP9_DX_SRCS-yes += decoder/decodemv.c
VP9_DX_SRCS-yes += decoder/decodframe.c
VP9_DX_SRCS-yes += decoder/dequantize.c
VP9_DX_SRCS-yes += decoder/detokenize.c
VP9_DX_SRCS-yes += decoder/dboolhuff.h
VP9_DX_SRCS-yes += decoder/decodemv.h
VP9_DX_SRCS-yes += decoder/dequantize.h
VP9_DX_SRCS-yes += decoder/detokenize.h
VP9_DX_SRCS-yes += decoder/onyxd_int.h
VP9_DX_SRCS-yes += decoder/treereader.h
VP9_DX_SRCS-yes += decoder/onyxd_if.c
VP9_DX_SRCS-yes += decoder/idct_blk.c

VP9_DX_SRCS-yes := $(filter-out $(VP9_DX_SRCS_REMOVE-yes),$(VP9_DX_SRCS-yes))

VP9_DX_SRCS-$(ARCH_X86)$(ARCH_X86_64) += decoder/x86/x86_dsystemdependent.c
VP9_DX_SRCS-$(HAVE_MMX) += decoder/x86/dequantize_mmx.asm
VP9_DX_SRCS-$(HAVE_MMX) += decoder/x86/idct_blk_mmx.c
VP9_DX_SRCS-$(HAVE_SSE2) += decoder/x86/idct_blk_sse2.c
