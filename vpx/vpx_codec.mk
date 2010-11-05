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

API_SRCS-$(CONFIG_DECODERS) += src/vpx_decoder.c
API_SRCS-$(CONFIG_DECODERS) += src/vpx_decoder_compat.c
API_SRCS-$(CONFIG_DECODERS) += vpx_decoder.h
API_SRCS-$(CONFIG_DECODERS) += vpx_decoder_compat.h
API_SRCS-$(CONFIG_ENCODERS) += src/vpx_encoder.c
API_SRCS-$(CONFIG_ENCODERS) += vpx_encoder.h
API_SRCS-yes                += internal/vpx_codec_internal.h
API_SRCS-yes                += src/vpx_codec.c
API_SRCS-yes                += src/vpx_image.c
API_SRCS-yes                += vpx_codec.h
API_SRCS-yes                += vpx_codec.mk
API_SRCS-yes                += vpx_codec_impl_bottom.h
API_SRCS-yes                += vpx_codec_impl_top.h
API_SRCS-yes                += vpx_image.h
