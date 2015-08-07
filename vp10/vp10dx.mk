##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

VP10_DX_EXPORTS += exports_dec

VP10_DX_SRCS-yes += $(VP10_COMMON_SRCS-yes)
VP10_DX_SRCS-no  += $(VP10_COMMON_SRCS-no)
VP10_DX_SRCS_REMOVE-yes += $(VP10_COMMON_SRCS_REMOVE-yes)
VP10_DX_SRCS_REMOVE-no  += $(VP10_COMMON_SRCS_REMOVE-no)

VP10_DX_SRCS-yes += vp10_dx_iface.c

VP10_DX_SRCS-yes += decoder/vp9_decodemv.c
VP10_DX_SRCS-yes += decoder/vp9_decodeframe.c
VP10_DX_SRCS-yes += decoder/vp9_decodeframe.h
VP10_DX_SRCS-yes += decoder/vp9_detokenize.c
VP10_DX_SRCS-yes += decoder/vp9_decodemv.h
VP10_DX_SRCS-yes += decoder/vp9_detokenize.h
VP10_DX_SRCS-yes += decoder/vp9_dthread.c
VP10_DX_SRCS-yes += decoder/vp9_dthread.h
VP10_DX_SRCS-yes += decoder/vp9_decoder.c
VP10_DX_SRCS-yes += decoder/vp9_decoder.h
VP10_DX_SRCS-yes += decoder/vp9_dsubexp.c
VP10_DX_SRCS-yes += decoder/vp9_dsubexp.h

VP10_DX_SRCS-yes := $(filter-out $(VP10_DX_SRCS_REMOVE-yes),$(VP10_DX_SRCS-yes))
