##
##  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


API_EXPORTS += exports

API_SRCS-$(CONFIG_AV1_ENCODER) += aom.h
API_SRCS-$(CONFIG_AV1_ENCODER) += aomcx.h
API_DOC_SRCS-$(CONFIG_AV1_ENCODER) += aom.h
API_DOC_SRCS-$(CONFIG_AV1_ENCODER) += aomcx.h

API_SRCS-$(CONFIG_AV1_DECODER) += aom.h
API_SRCS-$(CONFIG_AV1_DECODER) += aomdx.h
API_DOC_SRCS-$(CONFIG_AV1_DECODER) += aom.h
API_DOC_SRCS-$(CONFIG_AV1_DECODER) += aomdx.h

API_DOC_SRCS-yes += aom_codec.h
API_DOC_SRCS-yes += aom_decoder.h
API_DOC_SRCS-yes += aom_encoder.h
API_DOC_SRCS-yes += aom_frame_buffer.h
API_DOC_SRCS-yes += aom_image.h

API_SRCS-yes += src/aom_decoder.c
API_SRCS-yes += aom_decoder.h
API_SRCS-yes += src/aom_encoder.c
API_SRCS-yes += aom_encoder.h
API_SRCS-yes += internal/aom_codec_internal.h
API_SRCS-yes += src/aom_codec.c
API_SRCS-yes += src/aom_image.c
API_SRCS-yes += aom_codec.h
API_SRCS-yes += aom_codec.mk
API_SRCS-yes += aom_frame_buffer.h
API_SRCS-yes += aom_image.h
API_SRCS-yes += aom_integer.h
