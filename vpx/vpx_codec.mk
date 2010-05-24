##
##  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license and patent
##  grant that can be found in the LICENSE file in the root of the source
##  tree. All contributing project authors may be found in the AUTHORS
##  file in the root of the source tree.
##


API_EXPORTS += exports

API_SRCS-yes += internal/vpx_codec_internal.h
API_SRCS-yes += vpx_codec.h
API_SRCS-yes += vpx_codec.mk
API_SRCS-yes += vpx_codec_impl_top.h
API_SRCS-yes += vpx_codec_impl_bottom.h
API_SRCS-yes += vpx_decoder.h
API_SRCS-yes += vpx_decoder_compat.h
API_SRCS-yes += vpx_encoder.h
API_SRCS-yes += vpx_image.h
API_SRCS-yes += src/vpx_codec.c
API_SRCS-yes += src/vpx_decoder.c
API_SRCS-yes += src/vpx_decoder_compat.c
API_SRCS-yes += src/vpx_image.c
API_SRCS-yes += src/vpx_encoder.c
