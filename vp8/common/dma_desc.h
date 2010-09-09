/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#ifndef _dma_desc_h
#define _dma_desc_h

#if defined(__cplusplus)
extern "C" {
#endif


#define NDSIZE_LG   0x00000900  // Next Descriptor Size
#define NDSIZE_SM   0x00000800  // Next Descriptor Size
#define NDSIZE_7    0x00000700  // Next Descriptor Size
#define NDSIZE_6    0x00000600  // Next Descriptor Size
#define NDSIZE_5    0x00000500  // Next Descriptor Size
#define NDSIZE_4    0x00000400  // Next Descriptor Size
#define NDSIZE_3    0x00000300  // Next Descriptor Size
#define NDSIZE_2    0x00000200  // Next Descriptor Size
#define NDSIZE_1    0x00000100  // Next Descriptor Size

#define FLOW_STOP       0x0000
#define FLOW_AUTO       0x1000
#define FLOW_DESC_AR    0x4000
#define FLOW_DESC_SM    0x6000
#define FLOW_DESC_LG    0x7000

    typedef struct
    {
        unsigned int ndp;
        //unsigned short ndpl;
        //unsigned short ndph;
        unsigned int sa;
        //unsigned short sal;
        //unsigned short sah;

        unsigned short dmacfg;
        unsigned short xcnt;
        unsigned short xmod;
        unsigned short ycnt;
        unsigned short ymod;

    } LARGE_DESC;

    typedef struct
    {
        unsigned short ndpl;
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
        unsigned short xcnt;
        unsigned short xmod;
        unsigned short ycnt;
        unsigned short ymod;
    } SMALL_DESC;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
        unsigned short xcnt;
        unsigned short xmod;
        unsigned short ycnt;
        unsigned short ymod;
    } ARRAY_DESC_7;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
        unsigned short xcnt;
        unsigned short xmod;
        unsigned short ycnt;
    } ARRAY_DESC_6;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
        unsigned short xcnt;
        unsigned short xmod;
    } ARRAY_DESC_5;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
        unsigned short xcnt;
    } ARRAY_DESC_4;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
        unsigned short dmacfg;
    } ARRAY_DESC_3;

    typedef struct
    {
        unsigned short sal;
        unsigned short sah;
    } ARRAY_DESC_2;

    typedef struct
    {
        unsigned short sal;
    } ARRAY_DESC_1;

#if defined(__cplusplus)
}
#endif

#endif //_dma_desc_h
