/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef CODEC_COMMON_INTERFACE_H
#define CODEC_COMMON_INTERFACE_H

#define __export
#define _export
#define dll_export   __declspec( dllexport )
#define dll_import   __declspec( dllimport )

// Playback ERROR Codes.
#define NO_DECODER_ERROR            0
#define REMOTE_DECODER_ERROR        -1

#define DFR_BAD_DCT_COEFF           -100
#define DFR_ZERO_LENGTH_FRAME       -101
#define DFR_FRAME_SIZE_INVALID      -102
#define DFR_OUTPUT_BUFFER_OVERFLOW  -103
#define DFR_INVALID_FRAME_HEADER    -104
#define FR_INVALID_MODE_TOKEN       -110
#define ETR_ALLOCATION_ERROR        -200
#define ETR_INVALID_ROOT_PTR        -201
#define SYNCH_ERROR                 -400
#define BUFFER_UNDERFLOW_ERROR      -500
#define PB_IB_OVERFLOW_ERROR        -501

// External error triggers
#define PB_HEADER_CHECKSUM_ERROR    -601
#define PB_DATA_CHECKSUM_ERROR      -602

// DCT Error Codes
#define DDCT_EXPANSION_ERROR        -700
#define DDCT_INVALID_TOKEN_ERROR    -701

// exception_errors
#define GEN_EXCEPTIONS              -800
#define EX_UNQUAL_ERROR             -801

// Unrecoverable error codes
#define FATAL_PLAYBACK_ERROR        -1000
#define GEN_ERROR_CREATING_CDC      -1001
#define GEN_THREAD_CREATION_ERROR   -1002
#define DFR_CREATE_BMP_FAILED       -1003

// YUV buffer configuration structure
typedef struct
{
    int     y_width;
    int     y_height;
    int     y_stride;

    int     uv_width;
    int     uv_height;
    int     uv_stride;

    unsigned char   *y_buffer;
    unsigned char   *u_buffer;
    unsigned char   *v_buffer;

} YUV_BUFFER_CONFIG;
typedef enum
{
    C_SET_KEY_FRAME,
    C_SET_FIXED_Q,
    C_SET_FIRSTPASS_FILE,
    C_SET_EXPERIMENTAL_MIN,
    C_SET_EXPERIMENTAL_MAX = C_SET_EXPERIMENTAL_MIN + 255,
    C_SET_CHECKPROTECT,
    C_SET_TESTMODE,
    C_SET_INTERNAL_SIZE,
    C_SET_RECOVERY_FRAME,
    C_SET_REFERENCEFRAME,
    C_SET_GOLDENFRAME

#ifndef VP50_COMP_INTERFACE
    // Specialist test facilities.
//    C_VCAP_PARAMS,              // DO NOT USE FOR NOW WITH VFW CODEC
#endif

} C_SETTING;

typedef unsigned long C_SET_VALUE;


#endif
