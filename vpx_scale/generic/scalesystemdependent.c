/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "vpx_scale/vpxscale.h"


void (*vp8_yv12_extend_frame_borders_ptr)(YV12_BUFFER_CONFIG *ybf);
extern void vp8_yv12_extend_frame_borders(YV12_BUFFER_CONFIG *ybf);

void (*vp8_yv12_copy_frame_yonly_ptr)(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
extern void vp8_yv12_copy_frame_yonly(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);

void (*vp8_yv12_copy_frame_ptr)(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);
extern void vp8_yv12_copy_frame(YV12_BUFFER_CONFIG *src_ybc, YV12_BUFFER_CONFIG *dst_ybc);

/****************************************************************************
*  Imports
*****************************************************************************/

/****************************************************************************
 *
 *  ROUTINE       : vp8_scale_machine_specific_config
 *
 *  INPUTS        : UINT32 Version : Codec version number.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Checks for machine specifc features such as MMX support
 *                  sets appropriate flags and function pointers.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void vp8_scale_machine_specific_config() {
  vp8_yv12_extend_frame_borders_ptr      = vp8_yv12_extend_frame_borders;
  vp8_yv12_copy_frame_yonly_ptr          = vp8_yv12_copy_frame_yonly;
  vp8_yv12_copy_frame_ptr           = vp8_yv12_copy_frame;

}
