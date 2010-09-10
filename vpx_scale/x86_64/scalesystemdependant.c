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
*   Module Title :     system_dependant.c
*
*   Description  :     Miscellaneous system dependant functions
*
****************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/
#include "vpx_scale/vpxscale.h"
#include "cpuidlib.h"

/****************************************************************************
*  Imports
*****************************************************************************/
extern void register_generic_scalers(void);
extern void register_mmxscalers(void);

/****************************************************************************
 *
 *  ROUTINE       : post_proc_machine_specific_config
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
void
vp8_scale_machine_specific_config(void)
{
    int wmt_enabled = 1;

    if (wmt_enabled)
    {
        register_mmxscalers();
    }
    else
    {
        register_generic_scalers();
    }
}
