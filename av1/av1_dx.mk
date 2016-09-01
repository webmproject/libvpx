##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

AV1_DX_EXPORTS += exports_dec

AV1_DX_SRCS-yes += $(AV1_COMMON_SRCS-yes)
AV1_DX_SRCS-no  += $(AV1_COMMON_SRCS-no)
AV1_DX_SRCS_REMOVE-yes += $(AV1_COMMON_SRCS_REMOVE-yes)
AV1_DX_SRCS_REMOVE-no  += $(AV1_COMMON_SRCS_REMOVE-no)

AV1_DX_SRCS-yes += av1_dx_iface.c

AV1_DX_SRCS-yes += decoder/decodemv.c
AV1_DX_SRCS-yes += decoder/decodeframe.c
AV1_DX_SRCS-yes += decoder/decodeframe.h
AV1_DX_SRCS-yes += decoder/detokenize.c
AV1_DX_SRCS-yes += decoder/decodemv.h
AV1_DX_SRCS-yes += decoder/detokenize.h
AV1_DX_SRCS-yes += decoder/dthread.c
AV1_DX_SRCS-yes += decoder/dthread.h
AV1_DX_SRCS-yes += decoder/decoder.c
AV1_DX_SRCS-yes += decoder/decoder.h
AV1_DX_SRCS-yes += decoder/dsubexp.c
AV1_DX_SRCS-yes += decoder/dsubexp.h
AV1_DX_SRCS-yes += decoder/bitreader.h

AV1_DX_SRCS-yes := $(filter-out $(AV1_DX_SRCS_REMOVE-yes),$(AV1_DX_SRCS-yes))
