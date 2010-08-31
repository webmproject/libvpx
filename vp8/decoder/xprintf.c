/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


/****************************************************************************
*
*   Module Title :     xprintf.cpp
*
*   Description  :     Display a printf style message on the current video frame.
*
****************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#ifdef _WIN32_WCE
#include <windows.h>
#endif
#include "xprintf.h"

/****************************************************************************
 *
 *  ROUTINE       : xprintf
 *
 *  INPUTS        : const PB_INSTANCE *ppbi : Pointer to decoder instance.
 *                  long n_pixel             : Offset into buffer to write text.
 *                  const char *format      : Format string for print.
 *                  ...                     : Variable length argument list.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : int: Size (in bytes) of the formatted text.
 *
 *  FUNCTION      : Display a printf style message on the current video frame.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
int onyx_xprintf(unsigned char *ppbuffer, long n_pixel, long n_size, long n_stride, const char *format, ...)
{
    BOOL b_rc;
    va_list arglist;
    HFONT hfont, hfonto;

    int rc = 0;
    char sz_formatted[256] = "";
    unsigned char *p_dest = &ppbuffer[n_pixel];

#ifdef _WIN32_WCE
    //  Set up temporary bitmap
    HDC hdc_memory   = NULL;
    HBITMAP hbm_temp = NULL;
    HBITMAP hbm_orig = NULL;

    RECT rect;

    //  Copy bitmap to video frame
    long x;
    long y;

    //  Format text
    va_start(arglist, format);
    _vsnprintf(sz_formatted, sizeof(sz_formatted), format, arglist);
    va_end(arglist);

    rect.left   = 0;
    rect.top    = 0;
    rect.right  = 8 * strlen(sz_formatted);
    rect.bottom = 8;

    hdc_memory = create_compatible_dc(NULL);

    if (hdc_memory == NULL)
        goto Exit;

    hbm_temp = create_bitmap(rect.right, rect.bottom, 1, 1, NULL);

    if (hbm_temp == NULL)
        goto Exit;

    hbm_orig = (HBITMAP)(select_object(hdc_memory, hbm_temp));

    if (!hbm_orig)
        goto Exit;

    //  Write text into bitmap
    //  font?
    hfont = create_font(8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, VARIABLE_PITCH | FF_SWISS, "");

    if (hfont == NULL)
        goto Exit;

    hfonto = (HFONT)(select_object(hdc_memory, hbm_temp));

    if (!hfonto)
        goto Exit;

    select_object(hdc_memory, hfont);
    set_text_color(hdc_memory, 1);
    set_bk_color(hdc_memory, 0);
    set_bk_mode(hdc_memory, TRANSPARENT);

    b_rc = bit_blt(hdc_memory, rect.left, rect.top, rect.right, rect.bottom, hdc_memory, rect.left, rect.top, BLACKNESS);

    if (!b_rc)
        goto Exit;

    b_rc = ext_text_out(hdc_memory, 0, 0, ETO_CLIPPED, &rect, sz_formatted, strlen(sz_formatted), NULL);

    if (!b_rc)
        goto Exit;

    for (y = rect.top; y < rect.bottom; ++y)
    {
        for (x = rect.left; x < rect.right; ++x)
        {
            if (get_pixel(hdc_memory, x, rect.bottom - 1 - y))
                p_dest[x] = 255;
        }

        p_dest += n_stride;
    }

    rc = strlen(sz_formatted);

Exit:

    if (hbm_temp != NULL)
    {
        if (hbm_orig != NULL)
        {
            select_object(hdc_memory, hbm_orig);
        }

        delete_object(hbm_temp);
    }

    if (hfont != NULL)
    {
        if (hfonto != NULL)
            select_object(hdc_memory, hfonto);

        delete_object(hfont);
    }

    if (hdc_memory != NULL)
        delete_dc(hdc_memory);

    hdc_memory = 0;

#endif

    return rc;
}
