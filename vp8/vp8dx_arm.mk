##
##  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


#VP8_DX_SRCS list is modified according to different platforms.

#File list for arm
# decoder
#VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/decodframe_arm.c
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/dequantize_arm.c
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/dsystemdependent.c

#VP8_DX_SRCS_REMOVE-$(HAVE_ARMV6)  += decoder/decodframe.c
VP8_DX_SRCS_REMOVE-$(HAVE_ARMV6)  += decoder/dequantize.c
VP8_DX_SRCS_REMOVE-$(HAVE_ARMV6)  += decoder/generic/dsystemdependent.c

#File list for armv6
# decoder
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequantdcidct_v6$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequantidct_v6$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV6)  += decoder/arm/armv6/dequantize_v6$(ASM)

#File list for neon
# decoder
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/dequantdcidct_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/dequantidct_neon$(ASM)
VP8_DX_SRCS-$(HAVE_ARMV7)  += decoder/arm/neon/dequantizeb_neon$(ASM)


#for new token test
ifeq ($(ARCH_ARM),yes)
VP8_DX_SRCS-$(CONFIG_NEW_TOKENS)  += decoder/arm/detokenize_arm_sjl.c
VP8_DX_SRCS-$(CONFIG_NEW_TOKENS)  += decoder/arm/detokenize_arm_v6$(ASM)
VP8_DX_SRCS-$(CONFIG_NEW_TOKENS)  += decoder/onyxd_if_sjl.c

VP8_DX_SRCS_REMOVE-$(CONFIG_NEW_TOKENS)  += decoder/arm/detokenize_arm.c
VP8_DX_SRCS_REMOVE-$(CONFIG_NEW_TOKENS)  += decoder/onyxd_if.c
endif
