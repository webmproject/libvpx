/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     preproc_if.h
*
*   Description  :     Pre-processor interface header file.
*
****************************************************************************/

#ifndef __PREPROC_IF_H
#define __PREPROC_IF_H

/****************************************************************************
*  Header Files
****************************************************************************/
#include "type_aliases.h"

/****************************************************************************
*  Types
****************************************************************************/

typedef struct
{
    UINT8 *Yuv0ptr;
    UINT8 *Yuv1ptr;

    UINT8   *frag_info;              // blocks coded : passed in
    UINT32   frag_info_element_size;   // size of each element
    UINT32   frag_info_coded_mask;     // mask to get at whether fragment is coded

    UINT32 *region_index;            // Gives pixel index for top left of each block
    UINT32 video_frame_height;
    UINT32 video_frame_width;
    UINT8 hfrag_pixels;
    UINT8 vfrag_pixels;

} SCAN_CONFIG_DATA;

typedef enum
{
    SCP_FILTER_ON_OFF,
    SCP_SET_SRF_OFFSET,
    SCP_SET_EBO_ON_OFF,
    SCP_SET_VCAP_LEVEL_OFFSET,
    SCP_SET_SHOW_LOCAL

} SCP_SETTINGS;

typedef struct PP_INSTANCE *x_pp_inst;

/****************************************************************************
*  Module statics
****************************************************************************/
/* Controls whether Early break out is on or off in default case */
#define EARLY_BREAKOUT_DEFAULT  TRUE

/****************************************************************************
*  Functions
****************************************************************************/
extern  void set_scan_param(x_pp_inst ppi, UINT32 param_id, INT32 param_value);
extern  UINT32 yuvanalyse_frame(x_pp_inst ppi, UINT32 *KFIndicator);
extern  x_pp_inst create_pp_instance(void);
extern  void delete_pp_instance(x_pp_inst *);
extern  BOOL scan_yuvinit(x_pp_inst,  SCAN_CONFIG_DATA *scan_config_ptr);

#endif
