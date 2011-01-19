/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vp8.h"

/*!\defgroup vp8_decoder WebM VP8 Decoder
 * \ingroup vp8
 *
 * @{
 */
/*!\file vp8dx.h
 * \brief Provides definitions for using the VP8 algorithm within the vpx Decoder
 *        interface.
 */
#ifndef VP8DX_H
#define VP8DX_H
#include "vpx/vpx_codec_impl_top.h"

/*!\brief Algorithm interface for VP8
 *
 * This interface provides the capability to decode raw VP8 streams, as would
 * be found in AVI files and other non-Flash uses.
 */
extern vpx_codec_iface_t  vpx_codec_vp8_dx_algo;
extern vpx_codec_iface_t* vpx_codec_vp8_dx(void);

/* Include controls common to both the encoder and decoder */
#include "vp8.h"


/*!\brief VP8 decoder control functions
 *
 * The set of macros define the control functions of VP8 decoder interface
 */
enum vp8d_dec_control_id
{
    VP8_DECODER_CTRL_ID_START   = 256,
    VP8D_GET_LAST_REF_UPDATES,              /**< control function to get info on which reference frames were updated
                                            by the last decode */
    VP8D_GET_FRAME_CORRUPTED,               /**< check if the indicated frame is corrupted */
    VP8_DECODER_CTRL_ID_MAX
} ;


/*!\brief VP8 encoder control function parameter type
 *
 * Defines the data types that VP8E control functions take. Note that
 * additional common controls are defined in vp8.h
 *
 */


VPX_CTRL_USE_TYPE(VP8D_GET_LAST_REF_UPDATES,   int *)
VPX_CTRL_USE_TYPE(VP8D_GET_FRAME_CORRUPTED,    int *)


/*! @} - end defgroup vp8_decoder */


#include "vpx/vpx_codec_impl_bottom.h"
#endif
