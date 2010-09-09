/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef VPXBLIT_H_INCL
#define VPXBLIT_H_INCL
/*==============================================================================
                              Includes
==============================================================================*/

/*==============================================================================
                              Defines
==============================================================================*/


#ifdef VPX_BIG_ENDIAN
#define BYTE_ZERO(X)    ((X & 0xFF000000) >> (24 - 2)   )
#define BYTE_ONE(X)     ((X & 0x00FF0000) >> (16 - 2)   )
#define BYTE_TWO(X)     ((X & 0x0000FF00) >> (8 - 2)    )
#define BYTE_THREE(X)   ((X & 0x000000FF) << (0 + 2)    )

#define BYTE_ZERO_UV(X) ((X & 0x0000FF00) >> (8 - 2)    )
#define BYTE_ONE_UV(X)  ((X & 0x000000FF) << (0 + 2)    )

#define REREFERENCE(X) (*((int *) &(X)))

#else

#define BYTE_THREE(X)   ((X & 0xFF000000) >> (24 - 2)   )
#define BYTE_TWO(X)     ((X & 0x00FF0000) >> (16 - 2)   )
#define BYTE_ONE(X)     ((X & 0x0000FF00) >> (8 - 2)    )
#define BYTE_ZERO(X)    ((X & 0x000000FF) << (0 + 2)    )

#define BYTE_ONE_UV(X) ((X & 0x0000FF00) >> (8 - 2) )
#define BYTE_ZERO_UV(X)     ((X & 0x000000FF) << (0 + 2)    )

#define REREFERENCE(X) (*((int *) &(X)))

#endif


/*==============================================================================
                            Type Definitions
==============================================================================*/
typedef struct  // YUV buffer configuration structure
{
    int   y_width;
    int   y_height;
    int   y_stride;

    int   uv_width;
    int   uv_height;
    int   uv_stride;

    char *y_buffer;
    char *u_buffer;
    char *v_buffer;

    char *uv_start;
    int   uv_dst_area;
    int   uv_used_area;

} VPX_BLIT_CONFIG;

typedef struct tx86_params
{
    unsigned int pushed_registers[6];
    unsigned int return_address;
    unsigned int dst;
    unsigned int scrn_pitch;
    VPX_BLIT_CONFIG *buff_config;
} x86_params;

/*=============================================================================
                                Enums
==============================================================================*/


/*==============================================================================
                              Structures
==============================================================================*/

/*==============================================================================
                             Constants
==============================================================================*/


/*==============================================================================
                               Variables
==============================================================================*/




/*==============================================================================
                            Function Protoypes/MICROS
==============================================================================*/
int vpx_get_size_of_pixel(unsigned int bd);
void *vpx_get_blitter(unsigned int bd);
void vpx_set_blit(void);
void vpx_destroy_blit(void);



#endif //VPXBLIT_H_INCL
