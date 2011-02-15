##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


#VP8_DX_SRCS list is modified according to different platforms.

VP8_DX_SRCS-$(ARCH_ARM)  += decoder/arm/arm_dsystemdependent.c
VP8_CX_SRCS-$(ARCH_ARM)  += decoder/asm_dec_offsets.c

VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/dequantize_arm.c

#File list for armv6
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequant_dc_idct_v6$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequant_idct_v6$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequantize_v6$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/idct_blk_v6.c

#File list for neon
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/idct_dequant_dc_full_2x_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/idct_dequant_dc_0_2x_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/dequant_idct_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/idct_dequant_full_2x_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/idct_dequant_0_2x_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/dequantizeb_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/idct_blk_neon.c
