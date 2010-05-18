/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
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

/****************************************************************************
*  Imports
*****************************************************************************/
extern int register_generic_scalers(void);
extern int de_register_generic_scalers(void);

/****************************************************************************
 *
 *  ROUTINE       : vp8_scale_machine_specific_config
 *
 *  INPUTS        : UINT32 Version : Codec version number.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : int
 *
 *  FUNCTION      : Checks for machine specifc features such as MMX support
 *                  sets appropriate flags and function pointers.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
int
vp8_scale_machine_specific_config()
{
    return register_generic_scalers();
}

/****************************************************************************
 *
 *  ROUTINE       : vp8_scale_machine_specific_config
 *
 *  INPUTS        : UINT32 Version : Codec version number.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : int
 *
 *  FUNCTION      : Resets the funtion pointers and deallocates memory.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
int
scale_machine_specific_de_config()
{
    return de_register_generic_scalers();
}
