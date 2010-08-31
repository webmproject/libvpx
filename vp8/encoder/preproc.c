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
*   Module Title :     preproc.c
*
*   Description  :     Simple pre-processor.
*
****************************************************************************/

/****************************************************************************
*  Header Files
****************************************************************************/

#include "memory.h"
#include "preproc7.h"
#include "vpx_mem/vpx_mem.h"

/****************************************************************************
*  Macros
****************************************************************************/
#define FRAMECOUNT 7
#define ROUNDUP32(X) ( ( ( (unsigned long) X ) + 31 )&( 0xFFFFFFE0 ) )

/****************************************************************************
*  Imports
****************************************************************************/
extern void vp8_get_processor_flags(int *mmx_enabled, int *xmm_enabled, int *wmt_enabled);

/****************************************************************************
*  Exported Global Variables
****************************************************************************/
void (*temp_filter)(pre_proc_instance *ppi, unsigned char *s, unsigned char *d, int bytes, int strength);
void temp_filter_mmx
(
    pre_proc_instance *ppi,
    unsigned char *s,
    unsigned char *d,
    int bytes,
    int strength
);
void temp_filter_wmt
(
    pre_proc_instance *ppi,
    unsigned char *s,
    unsigned char *d,
    int bytes,
    int strength
);

/****************************************************************************
 *
 *  ROUTINE       : temp_filter_c
 *
 *  INPUTS        : pre_proc_instance *ppi : Pointer to pre-processor instance.
 *                  unsigned char *s     : Pointer to source frame.
 *                  unsigned char *d     : Pointer to destination frame.
 *                  int bytes            : Number of bytes to filter.
 *                  int strength         : Strength of filter to apply.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Performs a closesness adjusted temporarl blur
 *
 *  SPECIAL NOTES : Destination frame can be same as source frame.
 *
 ****************************************************************************/
void temp_filter_c
(
    pre_proc_instance *ppi,
    unsigned char *s,
    unsigned char *d,
    int bytes,
    int strength
)
{
    int byte = 0;
    unsigned char *frameptr = ppi->frame_buffer;

    if (ppi->frame == 0)
    {
        do
        {
            int frame = 0;

            do
            {
                *frameptr = s[byte];
                ++frameptr;
                ++frame;
            }
            while (frame < FRAMECOUNT);

            d[byte] = s[byte];

            ++byte;
        }
        while (byte < bytes);
    }
    else
    {
        int modifier;
        int offset = (ppi->frame % FRAMECOUNT);

        do
        {
            int accumulator = 0;
            int count = 0;
            int frame = 0;

            frameptr[offset] = s[byte];

            do
            {
                int pixel_value = *frameptr;

                modifier   = s[byte];
                modifier  -= pixel_value;
                modifier  *= modifier;
                modifier >>= strength;
                modifier  *= 3;

                if (modifier > 16)
                    modifier = 16;

                modifier = 16 - modifier;

                accumulator += modifier * pixel_value;

                count += modifier;

                frameptr++;

                ++frame;
            }
            while (frame < FRAMECOUNT);

            accumulator += (count >> 1);
            accumulator *= ppi->fixed_divide[count];
            accumulator >>= 16;

            d[byte] = accumulator;

            ++byte;
        }
        while (byte < bytes);
    }

    ++ppi->frame;
}
/****************************************************************************
 *
 *  ROUTINE       : delete_pre_proc
 *
 *  INPUTS        : pre_proc_instance *ppi : Pointer to pre-processor instance.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void
 *
 *  FUNCTION      : Deletes a pre-processing instance.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void delete_pre_proc(pre_proc_instance *ppi)
{
    if (ppi->frame_buffer_alloc)
        vpx_free(ppi->frame_buffer_alloc);

    ppi->frame_buffer_alloc = 0;
    ppi->frame_buffer      = 0;

    if (ppi->fixed_divide_alloc)
        vpx_free(ppi->fixed_divide_alloc);

    ppi->fixed_divide_alloc = 0;
    ppi->fixed_divide      = 0;
}

/****************************************************************************
 *
 *  ROUTINE       : init_pre_proc
 *
 *  INPUTS        : pre_proc_instance *ppi : Pointer to pre-processor instance.
 *                  int frame_size        : Number of bytes in one frame.
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : int: 1 if successful, 0 if failed.
 *
 *  FUNCTION      : Initializes prepprocessor instance.
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
int init_pre_proc7(pre_proc_instance *ppi, int frame_size)
{
    int i;
    int mmx_enabled;
    int xmm_enabled;
    int wmt_enabled;

    vp8_get_processor_flags(&mmx_enabled, &xmm_enabled, &wmt_enabled);

    if (wmt_enabled)
        temp_filter = temp_filter_wmt;
    else if (mmx_enabled)
        temp_filter = temp_filter_mmx;
    else
        temp_filter = temp_filter_c;


    delete_pre_proc(ppi);

    ppi->frame_buffer_alloc = vpx_malloc(32 + frame_size * FRAMECOUNT * sizeof(unsigned char));

    if (!ppi->frame_buffer_alloc)
    {
        delete_pre_proc(ppi);
        return 0;
    }

    ppi->frame_buffer = (unsigned char *) ROUNDUP32(ppi->frame_buffer_alloc);

    ppi->fixed_divide_alloc = vpx_malloc(32 + 255 * sizeof(unsigned int));

    if (!ppi->fixed_divide_alloc)
    {
        delete_pre_proc(ppi);
        return 0;
    }

    ppi->fixed_divide = (unsigned int *) ROUNDUP32(ppi->fixed_divide_alloc);

    for (i = 1; i < 255; i++)
        ppi->fixed_divide[i] = 0x10000 / i;

    return 1;
}
