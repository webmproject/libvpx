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
*   Module Title :     scaleopt.cpp
*
*   Description  :     Optimized scaling functions
*
****************************************************************************/
#include "pragmas.h"



/****************************************************************************
*  Module Statics
****************************************************************************/
__declspec(align(16)) const static unsigned short one_fifth[]  = { 51, 51, 51, 51 };
__declspec(align(16)) const static unsigned short two_fifths[] = { 102, 102, 102, 102 };
__declspec(align(16)) const static unsigned short three_fifths[] = { 154, 154, 154, 154 };
__declspec(align(16)) const static unsigned short four_fifths[] = { 205, 205, 205, 205 };
__declspec(align(16)) const static unsigned short round_values[] = { 128, 128, 128, 128 };
__declspec(align(16)) const static unsigned short four_ones[] = { 1, 1, 1, 1};
__declspec(align(16)) const static unsigned short const45_2[] = {205, 154, 102,  51 };
__declspec(align(16)) const static unsigned short const45_1[] = { 51, 102, 154, 205 };
__declspec(align(16)) const static unsigned char  mask45[] = { 0, 0, 0, 0, 0, 0, 255, 0};
__declspec(align(16)) const static unsigned short const35_2[] = { 154,  51, 205, 102 };
__declspec(align(16)) const static unsigned short const35_1[] = { 102, 205,  51, 154 };



#include "vpx_scale/vpxscale.h"
#include "vpx_mem/vpx_mem.h"

/****************************************************************************
*
*  ROUTINE       : horizontal_line_3_5_scale_mmx
*
*  INPUTS        : const unsigned char *source :
*                  unsigned int source_width    :
*                  unsigned char *dest         :
*                  unsigned int dest_width      :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 3 to 5 up-scaling of a horizontal line of pixels.
*
*  SPECIAL NOTES : None.
*
****************************************************************************/
static
void horizontal_line_3_5_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    (void) dest_width;

    __asm
    {

        push        rbx

        mov         rsi,    source
        mov         rdi,    dest

        mov         ecx,    source_width
        lea         rdx,    [rsi+rcx-3];

        movq        mm5,    const35_1       // mm5 = 66 xx cd xx 33 xx 9a xx
        movq        mm6,    const35_2       // mm6 = 9a xx 33 xx cd xx 66 xx

        movq        mm4,    round_values     // mm4 = 80 xx 80 xx 80 xx 80 xx
        pxor        mm7,    mm7             // clear mm7

        horiz_line_3_5_loop:

        mov         eax,    DWORD PTR [rsi] // eax = 00 01 02 03
        mov         ebx,    eax

        and         ebx,    0xffff00        // ebx = xx 01 02 xx
        mov         ecx,    eax             // ecx = 00 01 02 03

        and         eax,    0xffff0000      // eax = xx xx 02 03
        xor         ecx,    eax             // ecx = 00 01 xx xx

        shr         ebx,    8               // ebx = 01 02 xx xx
        or          eax,    ebx             // eax = 01 02 02 03

        shl         ebx,    16              // ebx = xx xx 01 02
        movd        mm1,    eax             // mm1 = 01 02 02 03 xx xx xx xx

        or          ebx,    ecx             // ebx = 00 01 01 02
        punpcklbw   mm1,    mm7             // mm1 = 01 xx 02 xx 02 xx 03 xx

        movd        mm0,    ebx             // mm0 = 00 01 01 02
        pmullw      mm1,    mm6             //

        punpcklbw   mm0,    mm7             // mm0 = 00 xx 01 xx 01 xx 02 xx
        pmullw      mm0,    mm5             //

        mov         [rdi],  ebx             // writeoutput 00 xx xx xx
        add         rsi,    3

        add         rdi,    5
        paddw       mm0,    mm1

        paddw       mm0,    mm4
        psrlw       mm0,    8

        cmp         rsi,    rdx
        packuswb    mm0,    mm7

        movd        DWORD Ptr [rdi-4], mm0
        jl          horiz_line_3_5_loop

//Exit:
        mov         eax,    DWORD PTR [rsi] // eax = 00 01 02 03
        mov         ebx,    eax

        and         ebx,    0xffff00        // ebx = xx 01 02 xx
        mov         ecx,    eax             // ecx = 00 01 02 03

        and         eax,    0xffff0000      // eax = xx xx 02 03
        xor         ecx,    eax             // ecx = 00 01 xx xx

        shr         ebx,    8               // ebx = 01 02 xx xx
        or          eax,    ebx             // eax = 01 02 02 03

        shl         eax,    8               // eax = xx 01 02 02
        and         eax,    0xffff0000      // eax = xx xx 02 02

        or          eax,    ebx             // eax = 01 02 02 02

        shl         ebx,    16              // ebx = xx xx 01 02
        movd        mm1,    eax             // mm1 = 01 02 02 02 xx xx xx xx

        or          ebx,    ecx             // ebx = 00 01 01 02
        punpcklbw   mm1,    mm7             // mm1 = 01 xx 02 xx 02 xx 02 xx

        movd        mm0,    ebx             // mm0 = 00 01 01 02
        pmullw      mm1,    mm6             //

        punpcklbw   mm0,    mm7             // mm0 = 00 xx 01 xx 01 xx 02 xx
        pmullw      mm0,    mm5             //

        mov         [rdi],  ebx             // writeoutput 00 xx xx xx
        paddw       mm0,    mm1

        paddw       mm0,    mm4
        psrlw       mm0,    8

        packuswb    mm0,    mm7
        movd        DWORD Ptr [rdi+1], mm0

        pop rbx

    }

}


/****************************************************************************
*
*  ROUTINE       : horizontal_line_4_5_scale_mmx
*
*  INPUTS        : const unsigned char *source :
*                  unsigned int source_width    :
*                  unsigned char *dest         :
*                  unsigned int dest_width      :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 4 to 5 up-scaling of a horizontal line of pixels.
*
*  SPECIAL NOTES : None.
*
****************************************************************************/
static
void horizontal_line_4_5_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    (void)dest_width;

    __asm
    {

        mov         rsi,    source
        mov         rdi,    dest

        mov         ecx,    source_width
        lea         rdx,    [rsi+rcx-8];

        movq        mm5,    const45_1       // mm5 = 33 xx 66 xx 9a xx cd xx
        movq        mm6,    const45_2       // mm6 = cd xx 9a xx 66 xx 33 xx

        movq        mm4,    round_values     // mm4 = 80 xx 80 xx 80 xx 80 xx
        pxor        mm7,    mm7             // clear mm7

        horiz_line_4_5_loop:

        movq        mm0,    QWORD PTR [rsi]           // mm0 = 00 01 02 03 04 05 06 07
        movq        mm1,    QWORD PTR [rsi+1];        // mm1 = 01 02 03 04 05 06 07 08

        movq        mm2,    mm0             // mm2 = 00 01 02 03 04 05 06 07
        movq        mm3,    mm1             // mm3 = 01 02 03 04 05 06 07 08

        movd        DWORD PTR [rdi],  mm0             // write output 00 xx xx xx
        punpcklbw   mm0,    mm7             // mm0 = 00 xx 01 xx 02 xx 03 xx

        punpcklbw   mm1,    mm7             // mm1 = 01 xx 02 xx 03 xx 04 xx
        pmullw      mm0,    mm5             // 00* 51 01*102 02*154 03*205

        pmullw      mm1,    mm6             // 01*205 02*154 03*102 04* 51
        punpckhbw   mm2,    mm7             // mm2 = 04 xx 05 xx 06 xx 07 xx

        movd        DWORD PTR [rdi+5], mm2            // write ouput 05 xx xx xx
        pmullw      mm2,    mm5             // 04* 51 05*102 06*154 07*205

        punpckhbw   mm3,    mm7             // mm3 = 05 xx 06 xx 07 xx 08 xx
        pmullw      mm3,    mm6             // 05*205 06*154 07*102 08* 51

        paddw       mm0,    mm1             // added round values
        paddw       mm0,    mm4

        psrlw       mm0,    8               // output: 01 xx 02 xx 03 xx 04 xx
        packuswb    mm0,    mm7

        movd        DWORD PTR [rdi+1], mm0  // write output 01 02 03 04
        add         rdi,    10

        add         rsi,    8
        paddw       mm2,    mm3             //

        paddw       mm2,    mm4             // added round values
        cmp         rsi,    rdx

        psrlw       mm2,    8
        packuswb    mm2,    mm7

        movd        DWORD PTR [rdi-4], mm2 // writeoutput 06 07 08 09
        jl         horiz_line_4_5_loop

//Exit:
        movq        mm0,    [rsi]           // mm0 = 00 01 02 03 04 05 06 07
        movq        mm1,    mm0             // mm1 = 00 01 02 03 04 05 06 07

        movq        mm2,    mm0             // mm2 = 00 01 02 03 04 05 06 07
        psrlq       mm1,    8               // mm1 = 01 02 03 04 05 06 07 00

        movq        mm3,    mask45          // mm3 = 00 00 00 00 00 00 ff 00
        pand        mm3,    mm1             // mm3 = 00 00 00 00 00 00 07 00

        psllq       mm3,    8               // mm3 = 00 00 00 00 00 00 00 07
        por         mm1,    mm3             // mm1 = 01 02 03 04 05 06 07 07

        movq        mm3,    mm1

        movd        DWORD PTR [rdi],  mm0   // write output 00 xx xx xx
        punpcklbw   mm0,    mm7             // mm0 = 00 xx 01 xx 02 xx 03 xx

        punpcklbw   mm1,    mm7             // mm1 = 01 xx 02 xx 03 xx 04 xx
        pmullw      mm0,    mm5             // 00* 51 01*102 02*154 03*205

        pmullw      mm1,    mm6             // 01*205 02*154 03*102 04* 51
        punpckhbw   mm2,    mm7             // mm2 = 04 xx 05 xx 06 xx 07 xx

        movd        DWORD PTR [rdi+5], mm2  // write ouput 05 xx xx xx
        pmullw      mm2,    mm5             // 04* 51 05*102 06*154 07*205

        punpckhbw   mm3,    mm7             // mm3 = 05 xx 06 xx 07 xx 08 xx
        pmullw      mm3,    mm6             // 05*205 06*154 07*102 07* 51

        paddw       mm0,    mm1             // added round values
        paddw       mm0,    mm4

        psrlw       mm0,    8               // output: 01 xx 02 xx 03 xx 04 xx
        packuswb    mm0,    mm7             // 01 02 03 04 xx xx xx xx

        movd        DWORD PTR [rdi+1], mm0  // write output 01 02 03 04
        paddw       mm2,    mm3             //

        paddw       mm2,    mm4             // added round values
        psrlw       mm2,    8

        packuswb    mm2,    mm7
        movd        DWORD PTR [rdi+6], mm2  // writeoutput 06 07 08 09


    }
}

/****************************************************************************
*
*  ROUTINE       : vertical_band_4_5_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 4 to 5 up-scaling of a 4 pixel high band of pixels.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has a "C" only
*                  version.
*
****************************************************************************/
static
void vertical_band_4_5_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {

        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        lea         rdi,    [rsi+rcx*2]             // tow lines below
        add         rdi,    rcx                     // three lines below

        pxor        mm7,    mm7                     // clear out mm7
        mov         edx,    dest_width               // Loop counter

        vs_4_5_loop:

        movq        mm0,    QWORD ptr [rsi]         // src[0];
        movq        mm1,    QWORD ptr [rsi+rcx]     // src[1];

        movq        mm2,    mm0                     // Make a copy
        punpcklbw   mm0,    mm7                     // unpack low to word

        movq        mm5,    one_fifth
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm0,    mm5                     // a * 1/5

        movq        mm3,    mm1                     // make a copy
        punpcklbw   mm1,    mm7                     // unpack low to word

        pmullw      mm2,    mm5                     // a * 1/5
        movq        mm6,    four_fifths               // constan

        movq        mm4,    mm1                     // copy of low b
        pmullw      mm4,    mm6                     // b * 4/5

        punpckhbw   mm3,    mm7                     // unpack high to word
        movq        mm5,    mm3                     // copy of high b

        pmullw      mm5,    mm6                     // b * 4/5
        paddw       mm0,    mm4                     // a * 1/5 + b * 4/5

        paddw       mm2,    mm5                     // a * 1/5 + b * 4/5
        paddw       mm0,    round_values             // + 128

        paddw       mm2,    round_values             // + 128
        psrlw       mm0,    8

        psrlw       mm2,    8
        packuswb    mm0,    mm2                     // des [1]

        movq        QWORD ptr [rsi+rcx], mm0        // write des[1]
        movq        mm0,    [rsi+rcx*2]             // mm0 = src[2]

        // mm1, mm3 --- Src[1]
        // mm0 --- Src[2]
        // mm7 for unpacking

        movq        mm5,    two_fifths
        movq        mm2,    mm0                     // make a copy

        pmullw      mm1,    mm5                     // b * 2/5
        movq        mm6,    three_fifths


        punpcklbw   mm0,    mm7                     // unpack low to word
        pmullw      mm3,    mm5                     // b * 2/5

        movq        mm4,    mm0                     // make copy of c
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm4,    mm6                     // c * 3/5
        movq        mm5,    mm2

        pmullw      mm5,    mm6                     // c * 3/5
        paddw       mm1,    mm4                     // b * 2/5 + c * 3/5

        paddw       mm3,    mm5                     // b * 2/5 + c * 3/5
        paddw       mm1,    round_values             // + 128

        paddw       mm3,    round_values             // + 128
        psrlw       mm1,    8

        psrlw       mm3,    8
        packuswb    mm1,    mm3                     // des[2]

        movq        QWORD ptr [rsi+rcx*2], mm1      // write des[2]
        movq        mm1,    [rdi]                   // mm1=Src[3];

        // mm0, mm2 --- Src[2]
        // mm1 --- Src[3]
        // mm6 --- 3/5
        // mm7 for unpacking

        pmullw      mm0,    mm6                     // c * 3/5
        movq        mm5,    two_fifths               // mm5 = 2/5

        movq        mm3,    mm1                     // make a copy
        pmullw      mm2,    mm6                     // c * 3/5

        punpcklbw   mm1,    mm7                     // unpack low
        movq        mm4,    mm1                     // make a copy

        punpckhbw   mm3,    mm7                     // unpack high
        pmullw      mm4,    mm5                     // d * 2/5

        movq        mm6,    mm3                     // make a copy
        pmullw      mm6,    mm5                     // d * 2/5

        paddw       mm0,    mm4                     // c * 3/5 + d * 2/5
        paddw       mm2,    mm6                     // c * 3/5 + d * 2/5

        paddw       mm0,    round_values             // + 128
        paddw       mm2,    round_values             // + 128

        psrlw       mm0,    8
        psrlw       mm2,    8

        packuswb    mm0,    mm2                     // des[3]
        movq        QWORD ptr [rdi], mm0            // write des[3]

        //  mm1, mm3 --- Src[3]
        //  mm7 -- cleared for unpacking

        movq        mm0,    [rdi+rcx*2]             // mm0, Src[0] of the next group

        movq        mm5,    four_fifths              // mm5 = 4/5
        pmullw      mm1,    mm5                     // d * 4/5

        movq        mm6,    one_fifth                // mm6 = 1/5
        movq        mm2,    mm0                     // make a copy

        pmullw      mm3,    mm5                     // d * 4/5
        punpcklbw   mm0,    mm7                     // unpack low

        pmullw      mm0,    mm6                     // an * 1/5
        punpckhbw   mm2,    mm7                     // unpack high

        paddw       mm1,    mm0                     // d * 4/5 + an * 1/5
        pmullw      mm2,    mm6                     // an * 1/5

        paddw       mm3,    mm2                     // d * 4/5 + an * 1/5
        paddw       mm1,    round_values             // + 128

        paddw       mm3,    round_values             // + 128
        psrlw       mm1,    8

        psrlw       mm3,    8
        packuswb    mm1,    mm3                     // des[4]

        movq        QWORD ptr [rdi+rcx], mm1        // write des[4]

        add         rdi,    8
        add         rsi,    8

        sub         rdx,    8
        jg          vs_4_5_loop
    }
}

/****************************************************************************
*
*  ROUTINE       : last_vertical_band_4_5_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : None
*
*  FUNCTION      : 4 to 5 up-scaling of the last 4-pixel high band in an image.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has an "C" only
*                  version.
*
****************************************************************************/
static
void last_vertical_band_4_5_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {
        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        lea         rdi,    [rsi+rcx*2]             // tow lines below
        add         rdi,    rcx                     // three lines below

        pxor        mm7,    mm7                     // clear out mm7
        mov         edx,    dest_width               // Loop counter

        last_vs_4_5_loop:

        movq        mm0,    QWORD ptr [rsi]         // src[0];
        movq        mm1,    QWORD ptr [rsi+rcx]     // src[1];

        movq        mm2,    mm0                     // Make a copy
        punpcklbw   mm0,    mm7                     // unpack low to word

        movq        mm5,    one_fifth
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm0,    mm5                     // a * 1/5

        movq        mm3,    mm1                     // make a copy
        punpcklbw   mm1,    mm7                     // unpack low to word

        pmullw      mm2,    mm5                     // a * 1/5
        movq        mm6,    four_fifths               // constan

        movq        mm4,    mm1                     // copy of low b
        pmullw      mm4,    mm6                     // b * 4/5

        punpckhbw   mm3,    mm7                     // unpack high to word
        movq        mm5,    mm3                     // copy of high b

        pmullw      mm5,    mm6                     // b * 4/5
        paddw       mm0,    mm4                     // a * 1/5 + b * 4/5

        paddw       mm2,    mm5                     // a * 1/5 + b * 4/5
        paddw       mm0,    round_values             // + 128

        paddw       mm2,    round_values             // + 128
        psrlw       mm0,    8

        psrlw       mm2,    8
        packuswb    mm0,    mm2                     // des [1]

        movq        QWORD ptr [rsi+rcx], mm0        // write des[1]
        movq        mm0,    [rsi+rcx*2]             // mm0 = src[2]

        // mm1, mm3 --- Src[1]
        // mm0 --- Src[2]
        // mm7 for unpacking

        movq        mm5,    two_fifths
        movq        mm2,    mm0                     // make a copy

        pmullw      mm1,    mm5                     // b * 2/5
        movq        mm6,    three_fifths


        punpcklbw   mm0,    mm7                     // unpack low to word
        pmullw      mm3,    mm5                     // b * 2/5

        movq        mm4,    mm0                     // make copy of c
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm4,    mm6                     // c * 3/5
        movq        mm5,    mm2

        pmullw      mm5,    mm6                     // c * 3/5
        paddw       mm1,    mm4                     // b * 2/5 + c * 3/5

        paddw       mm3,    mm5                     // b * 2/5 + c * 3/5
        paddw       mm1,    round_values             // + 128

        paddw       mm3,    round_values             // + 128
        psrlw       mm1,    8

        psrlw       mm3,    8
        packuswb    mm1,    mm3                     // des[2]

        movq        QWORD ptr [rsi+rcx*2], mm1      // write des[2]
        movq        mm1,    [rdi]                   // mm1=Src[3];

        movq        QWORD ptr [rdi+rcx], mm1        // write des[4];

        // mm0, mm2 --- Src[2]
        // mm1 --- Src[3]
        // mm6 --- 3/5
        // mm7 for unpacking

        pmullw      mm0,    mm6                     // c * 3/5
        movq        mm5,    two_fifths               // mm5 = 2/5

        movq        mm3,    mm1                     // make a copy
        pmullw      mm2,    mm6                     // c * 3/5

        punpcklbw   mm1,    mm7                     // unpack low
        movq        mm4,    mm1                     // make a copy

        punpckhbw   mm3,    mm7                     // unpack high
        pmullw      mm4,    mm5                     // d * 2/5

        movq        mm6,    mm3                     // make a copy
        pmullw      mm6,    mm5                     // d * 2/5

        paddw       mm0,    mm4                     // c * 3/5 + d * 2/5
        paddw       mm2,    mm6                     // c * 3/5 + d * 2/5

        paddw       mm0,    round_values             // + 128
        paddw       mm2,    round_values             // + 128

        psrlw       mm0,    8
        psrlw       mm2,    8

        packuswb    mm0,    mm2                     // des[3]
        movq        QWORD ptr [rdi], mm0            // write des[3]

        //  mm1, mm3 --- Src[3]
        //  mm7 -- cleared for unpacking
        add         rdi,    8
        add         rsi,    8

        sub         rdx,    8
        jg          last_vs_4_5_loop
    }
}

/****************************************************************************
*
*  ROUTINE       : vertical_band_3_5_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 3 to 5 up-scaling of a 3-pixel high band of pixels.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has an "C" only
*                  version.
*
****************************************************************************/
static
void vertical_band_3_5_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {
        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        lea         rdi,    [rsi+rcx*2]             // two lines below
        add         rdi,    rcx                     // three lines below

        pxor        mm7,    mm7                     // clear out mm7
        mov         edx,    dest_width               // Loop counter

        vs_3_5_loop:

        movq        mm0,    QWORD ptr [rsi]         // src[0];
        movq        mm1,    QWORD ptr [rsi+rcx]     // src[1];

        movq        mm2,    mm0                     // Make a copy
        punpcklbw   mm0,    mm7                     // unpack low to word

        movq        mm5,    two_fifths               // mm5 = 2/5
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm0,    mm5                     // a * 2/5

        movq        mm3,    mm1                     // make a copy
        punpcklbw   mm1,    mm7                     // unpack low to word

        pmullw      mm2,    mm5                     // a * 2/5
        movq        mm6,    three_fifths             // mm6 = 3/5

        movq        mm4,    mm1                     // copy of low b
        pmullw      mm4,    mm6                     // b * 3/5

        punpckhbw   mm3,    mm7                     // unpack high to word
        movq        mm5,    mm3                     // copy of high b

        pmullw      mm5,    mm6                     // b * 3/5
        paddw       mm0,    mm4                     // a * 2/5 + b * 3/5

        paddw       mm2,    mm5                     // a * 2/5 + b * 3/5
        paddw       mm0,    round_values             // + 128

        paddw       mm2,    round_values             // + 128
        psrlw       mm0,    8

        psrlw       mm2,    8
        packuswb    mm0,    mm2                     // des [1]

        movq        QWORD ptr [rsi+rcx], mm0        // write des[1]
        movq        mm0,    [rsi+rcx*2]             // mm0 = src[2]

        // mm1, mm3 --- Src[1]
        // mm0 --- Src[2]
        // mm7 for unpacking

        movq        mm4,    mm1                     // b low
        pmullw      mm1,    four_fifths              // b * 4/5 low

        movq        mm5,    mm3                     // b high
        pmullw      mm3,    four_fifths              // b * 4/5 high

        movq        mm2,    mm0                     // c
        pmullw      mm4,    one_fifth                // b * 1/5

        punpcklbw   mm0,    mm7                     // c low
        pmullw      mm5,    one_fifth                // b * 1/5

        movq        mm6,    mm0                     // make copy of c low
        punpckhbw   mm2,    mm7                     // c high

        pmullw      mm6,    one_fifth                // c * 1/5 low
        movq        mm7,    mm2                     // make copy of c high

        pmullw      mm7,    one_fifth                // c * 1/5 high
        paddw       mm1,    mm6                     // b * 4/5 + c * 1/5 low

        paddw       mm3,    mm7                     // b * 4/5 + c * 1/5 high
        movq        mm6,    mm0                     // make copy of c low

        pmullw      mm6,    four_fifths              // c * 4/5 low
        movq        mm7,    mm2                     // make copy of c high

        pmullw      mm7,    four_fifths              // c * 4/5 high

        paddw       mm4,    mm6                     // b * 1/5 + c * 4/5 low
        paddw       mm5,    mm7                     // b * 1/5 + c * 4/5 high

        paddw       mm1,    round_values             // + 128
        paddw       mm3,    round_values             // + 128

        psrlw       mm1,    8
        psrlw       mm3,    8

        packuswb    mm1,    mm3                     // des[2]
        movq        QWORD ptr [rsi+rcx*2], mm1      // write des[2]

        paddw       mm4,    round_values             // + 128
        paddw       mm5,    round_values             // + 128

        psrlw       mm4,    8
        psrlw       mm5,    8

        packuswb    mm4,    mm5                     // des[3]
        movq        QWORD ptr [rdi], mm4            // write des[3]

        //  mm0, mm2 --- Src[3]

        pxor        mm7,    mm7                     // clear mm7 for unpacking
        movq        mm1,    [rdi+rcx*2]             // mm1 = Src[0] of the next group

        movq        mm5,    three_fifths             // mm5 = 3/5
        pmullw      mm0,    mm5                     // d * 3/5

        movq        mm6,    two_fifths                // mm6 = 2/5
        movq        mm3,    mm1                     // make a copy

        pmullw      mm2,    mm5                     // d * 3/5
        punpcklbw   mm1,    mm7                     // unpack low

        pmullw      mm1,    mm6                     // an * 2/5
        punpckhbw   mm3,    mm7                     // unpack high

        paddw       mm0,    mm1                     // d * 3/5 + an * 2/5
        pmullw      mm3,    mm6                     // an * 2/5

        paddw       mm2,    mm3                     // d * 3/5 + an * 2/5
        paddw       mm0,    round_values             // + 128

        paddw       mm2,    round_values             // + 128
        psrlw       mm0,    8

        psrlw       mm2,    8
        packuswb    mm0,    mm2                     // des[4]

        movq        QWORD ptr [rdi+rcx], mm0        // write des[4]

        add         rdi,    8
        add         rsi,    8

        sub         rdx,    8
        jg          vs_3_5_loop
    }
}

/****************************************************************************
*
*  ROUTINE       : last_vertical_band_3_5_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 3 to 5 up-scaling of a 3-pixel high band of pixels.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has an "C" only
*                  version.
*
****************************************************************************/
static
void last_vertical_band_3_5_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {
        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        lea         rdi,    [rsi+rcx*2]             // tow lines below
        add         rdi,    rcx                     // three lines below

        pxor        mm7,    mm7                     // clear out mm7
        mov         edx,    dest_width               // Loop counter


        last_vs_3_5_loop:

        movq        mm0,    QWORD ptr [rsi]         // src[0];
        movq        mm1,    QWORD ptr [rsi+rcx]     // src[1];

        movq        mm2,    mm0                     // Make a copy
        punpcklbw   mm0,    mm7                     // unpack low to word

        movq        mm5,    two_fifths               // mm5 = 2/5
        punpckhbw   mm2,    mm7                     // unpack high to word

        pmullw      mm0,    mm5                     // a * 2/5

        movq        mm3,    mm1                     // make a copy
        punpcklbw   mm1,    mm7                     // unpack low to word

        pmullw      mm2,    mm5                     // a * 2/5
        movq        mm6,    three_fifths             // mm6 = 3/5

        movq        mm4,    mm1                     // copy of low b
        pmullw      mm4,    mm6                     // b * 3/5

        punpckhbw   mm3,    mm7                     // unpack high to word
        movq        mm5,    mm3                     // copy of high b

        pmullw      mm5,    mm6                     // b * 3/5
        paddw       mm0,    mm4                     // a * 2/5 + b * 3/5

        paddw       mm2,    mm5                     // a * 2/5 + b * 3/5
        paddw       mm0,    round_values             // + 128

        paddw       mm2,    round_values             // + 128
        psrlw       mm0,    8

        psrlw       mm2,    8
        packuswb    mm0,    mm2                     // des [1]

        movq        QWORD ptr [rsi+rcx], mm0        // write des[1]
        movq        mm0,    [rsi+rcx*2]             // mm0 = src[2]



        // mm1, mm3 --- Src[1]
        // mm0 --- Src[2]
        // mm7 for unpacking

        movq        mm4,    mm1                     // b low
        pmullw      mm1,    four_fifths              // b * 4/5 low

        movq        QWORD ptr [rdi+rcx], mm0        // write des[4]

        movq        mm5,    mm3                     // b high
        pmullw      mm3,    four_fifths              // b * 4/5 high

        movq        mm2,    mm0                     // c
        pmullw      mm4,    one_fifth                // b * 1/5

        punpcklbw   mm0,    mm7                     // c low
        pmullw      mm5,    one_fifth                // b * 1/5

        movq        mm6,    mm0                     // make copy of c low
        punpckhbw   mm2,    mm7                     // c high

        pmullw      mm6,    one_fifth                // c * 1/5 low
        movq        mm7,    mm2                     // make copy of c high

        pmullw      mm7,    one_fifth                // c * 1/5 high
        paddw       mm1,    mm6                     // b * 4/5 + c * 1/5 low

        paddw       mm3,    mm7                     // b * 4/5 + c * 1/5 high
        movq        mm6,    mm0                     // make copy of c low

        pmullw      mm6,    four_fifths              // c * 4/5 low
        movq        mm7,    mm2                     // make copy of c high

        pmullw      mm7,    four_fifths              // c * 4/5 high

        paddw       mm4,    mm6                     // b * 1/5 + c * 4/5 low
        paddw       mm5,    mm7                     // b * 1/5 + c * 4/5 high

        paddw       mm1,    round_values             // + 128
        paddw       mm3,    round_values             // + 128

        psrlw       mm1,    8
        psrlw       mm3,    8

        packuswb    mm1,    mm3                     // des[2]
        movq        QWORD ptr [rsi+rcx*2], mm1      // write des[2]

        paddw       mm4,    round_values             // + 128
        paddw       mm5,    round_values             // + 128

        psrlw       mm4,    8
        psrlw       mm5,    8

        packuswb    mm4,    mm5                     // des[3]
        movq        QWORD ptr [rdi], mm4            // write des[3]

        //  mm0, mm2 --- Src[3]

        add         rdi,    8
        add         rsi,    8

        sub         rdx,    8
        jg          last_vs_3_5_loop
    }
}

/****************************************************************************
*
*  ROUTINE       : vertical_band_1_2_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 1 to 2 up-scaling of a band of pixels.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has an "C" only
*                  version.
*
****************************************************************************/
static
void vertical_band_1_2_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {

        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        pxor        mm7,    mm7                     // clear out mm7
        mov         edx,    dest_width               // Loop counter

        vs_1_2_loop:

        movq        mm0,    [rsi]                   // get Src[0]
        movq        mm1,    [rsi + rcx * 2]         // get Src[1]

        movq        mm2,    mm0                     // make copy before unpack
        movq        mm3,    mm1                     // make copy before unpack

        punpcklbw   mm0,    mm7                     // low Src[0]
        movq        mm6,    four_ones                // mm6= 1, 1, 1, 1

        punpcklbw   mm1,    mm7                     // low Src[1]
        paddw       mm0,    mm1                     // low (a + b)

        punpckhbw   mm2,    mm7                     // high Src[0]
        paddw       mm0,    mm6                     // low (a + b + 1)

        punpckhbw   mm3,    mm7
        paddw       mm2,    mm3                     // high (a + b )

        psraw       mm0,    1                       // low (a + b +1 )/2
        paddw       mm2,    mm6                     // high (a + b + 1)

        psraw       mm2,    1                       // high (a + b + 1)/2
        packuswb    mm0,    mm2                     // pack results

        movq        [rsi+rcx], mm0                  // write out eight bytes
        add         rsi,    8

        sub         rdx,    8
        jg          vs_1_2_loop
    }

}

/****************************************************************************
*
*  ROUTINE       : last_vertical_band_1_2_scale_mmx
*
*  INPUTS        : unsigned char *dest    :
*                  unsigned int dest_pitch :
*                  unsigned int dest_width :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 1 to 2 up-scaling of band of pixels.
*
*  SPECIAL NOTES : The routine uses the first line of the band below
*                  the current band. The function also has an "C" only
*                  version.
*
****************************************************************************/
static
void last_vertical_band_1_2_scale_mmx
(
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {
        mov         rsi,    dest                    // Get the source and destination pointer
        mov         ecx,    dest_pitch               // Get the pitch size

        mov         edx,    dest_width               // Loop counter

        last_vs_1_2_loop:

        movq        mm0,    [rsi]                   // get Src[0]
        movq        [rsi+rcx], mm0                  // write out eight bytes

        add         rsi,    8
        sub         rdx,    8

        jg          last_vs_1_2_loop
    }
}

/****************************************************************************
*
*  ROUTINE       : horizontal_line_1_2_scale
*
*  INPUTS        : const unsigned char *source :
*                  unsigned int source_width    :
*                  unsigned char *dest         :
*                  unsigned int dest_width      :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 1 to 2 up-scaling of a horizontal line of pixels.
*
*  SPECIAL NOTES : None.
*
****************************************************************************/
static
void horizontal_line_1_2_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    (void) dest_width;

    __asm
    {
        mov         rsi,    source
        mov         rdi,    dest

        pxor        mm7,    mm7
        movq        mm6,    four_ones

        mov         ecx,    source_width

        hs_1_2_loop:

        movq        mm0,    [rsi]
        movq        mm1,    [rsi+1]

        movq        mm2,    mm0
        movq        mm3,    mm1

        movq        mm4,    mm0
        punpcklbw   mm0,    mm7

        punpcklbw   mm1,    mm7
        paddw       mm0,    mm1

        paddw       mm0,    mm6
        punpckhbw   mm2,    mm7

        punpckhbw   mm3,    mm7
        paddw       mm2,    mm3

        paddw       mm2,    mm6
        psraw       mm0,    1

        psraw       mm2,    1
        packuswb    mm0,    mm2

        movq        mm2,    mm4
        punpcklbw   mm2,    mm0

        movq        [rdi],  mm2
        punpckhbw   mm4,    mm0

        movq        [rdi+8], mm4
        add         rsi,    8

        add         rdi,    16
        sub         rcx,    8

        cmp         rcx,    8
        jg          hs_1_2_loop

// last eight pixel

        movq        mm0,    [rsi]
        movq        mm1,    mm0

        movq        mm2,    mm0
        movq        mm3,    mm1

        psrlq       mm1,    8
        psrlq       mm3,    56

        psllq       mm3,    56
        por         mm1,    mm3

        movq        mm3,    mm1
        movq        mm4,    mm0

        punpcklbw   mm0,    mm7
        punpcklbw   mm1,    mm7

        paddw       mm0,    mm1
        paddw       mm0,    mm6

        punpckhbw   mm2,    mm7
        punpckhbw   mm3,    mm7

        paddw       mm2,    mm3
        paddw       mm2,    mm6

        psraw       mm0,    1
        psraw       mm2,    1

        packuswb    mm0,    mm2
        movq        mm2,    mm4

        punpcklbw   mm2,    mm0
        movq        [rdi],  mm2

        punpckhbw   mm4,    mm0
        movq        [rdi+8], mm4
    }
}





__declspec(align(16)) const static unsigned short const54_2[] = {  0,  64, 128, 192 };
__declspec(align(16)) const static unsigned short const54_1[] = {256, 192, 128,  64 };


/****************************************************************************
*
*  ROUTINE       : horizontal_line_5_4_scale_mmx
*
*  INPUTS        : const unsigned char *source : Pointer to source data.
*                  unsigned int source_width    : Stride of source.
*                  unsigned char *dest         : Pointer to destination data.
*                  unsigned int dest_width      : Stride of destination (NOT USED).
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : Copies horizontal line of pixels from source to
*                  destination scaling up by 4 to 5.
*
*  SPECIAL NOTES : None.
*
****************************************************************************/
static
void horizontal_line_5_4_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    /*
    unsigned i;
    unsigned int a, b, c, d, e;
    unsigned char *des = dest;
    const unsigned char *src = source;

    (void) dest_width;

    for ( i=0; i<source_width; i+=5 )
    {
        a = src[0];
        b = src[1];
        c = src[2];
        d = src[3];
        e = src[4];

        des[0] = a;
        des[1] = ((b*192 + c* 64 + 128)>>8);
        des[2] = ((c*128 + d*128 + 128)>>8);
        des[3] = ((d* 64 + e*192 + 128)>>8);

        src += 5;
        des += 4;
    }
    */
    __asm
    {

        mov         rsi,        source              ;
        mov         rdi,        dest                ;

        mov         ecx,        source_width         ;
        movq        mm5,        const54_1           ;

        pxor        mm7,        mm7                 ;
        movq        mm6,        const54_2           ;

        movq        mm4,        round_values         ;
        lea         rdx,        [rsi+rcx]           ;
        horizontal_line_5_4_loop:

        movq        mm0,        QWORD PTR  [rsi]    ;
        00 01 02 03 04 05 06 07
        movq        mm1,        mm0                 ;
        00 01 02 03 04 05 06 07

        psrlq       mm0,        8                   ;
        01 02 03 04 05 06 07 xx
        punpcklbw   mm1,        mm7                 ;
        xx 00 xx 01 xx 02 xx 03

        punpcklbw   mm0,        mm7                 ;
        xx 01 xx 02 xx 03 xx 04
        pmullw      mm1,        mm5

        pmullw      mm0,        mm6
        add         rsi,        5

        add         rdi,        4
        paddw       mm1,        mm0

        paddw       mm1,        mm4
        psrlw       mm1,        8

        cmp         rsi,        rdx
        packuswb    mm1,        mm7

        movd        DWORD PTR [rdi-4], mm1

        jl          horizontal_line_5_4_loop

    }

}
__declspec(align(16)) const static unsigned short one_fourths[]   = {  64,  64,  64, 64  };
__declspec(align(16)) const static unsigned short two_fourths[]   = { 128, 128, 128, 128 };
__declspec(align(16)) const static unsigned short three_fourths[] = { 192, 192, 192, 192 };

static
void vertical_band_5_4_scale_mmx
(
    unsigned char *source,
    unsigned int src_pitch,
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{

    __asm
    {

        mov         rsi,    source                    // Get the source and destination pointer
        mov         ecx,    src_pitch               // Get the pitch size

        mov         rdi,    dest                    // tow lines below
        pxor        mm7,    mm7                     // clear out mm7

        mov         edx,    dest_pitch               // Loop counter
        mov         ebx,    dest_width

        vs_5_4_loop:

        movd        mm0,    DWORD ptr [rsi]         // src[0];
        movd        mm1,    DWORD ptr [rsi+rcx]     // src[1];

        movd        mm2,    DWORD ptr [rsi+rcx*2]
        lea         rax,    [rsi+rcx*2]             //

        punpcklbw   mm1,    mm7
        punpcklbw   mm2,    mm7

        movq        mm3,    mm2
        pmullw      mm1,    three_fourths

        pmullw      mm2,    one_fourths
        movd        mm4,    [rax+rcx]

        pmullw      mm3,    two_fourths
        punpcklbw   mm4,    mm7

        movq        mm5,    mm4
        pmullw      mm4,    two_fourths

        paddw       mm1,    mm2
        movd        mm6,    [rax+rcx*2]

        pmullw      mm5,    one_fourths
        paddw       mm1,    round_values;

        paddw       mm3,    mm4
        psrlw       mm1,    8

        punpcklbw   mm6,    mm7
        paddw       mm3,    round_values

        pmullw      mm6,    three_fourths
        psrlw       mm3,    8

        packuswb    mm1,    mm7
        packuswb    mm3,    mm7

        movd        DWORD PTR [rdi], mm0
        movd        DWORD PTR [rdi+rdx], mm1


        paddw       mm5,    mm6
        movd        DWORD PTR [rdi+rdx*2], mm3

        lea         rax,    [rdi+rdx*2]
        paddw       mm5,    round_values

        psrlw       mm5,    8
        add         rdi,    4

        packuswb    mm5,    mm7
        movd        DWORD PTR [rax+rdx], mm5

        add         rsi,    4
        sub         rbx,    4

        jg         vs_5_4_loop
    }
}


__declspec(align(16)) const static unsigned short const53_1[] = {  0,  85, 171, 0 };
__declspec(align(16)) const static unsigned short const53_2[] = {256, 171,  85, 0 };


static
void horizontal_line_5_3_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    __asm
    {

        mov         rsi,        source              ;
        mov         rdi,        dest                ;

        mov         ecx,        source_width         ;
        movq        mm5,        const53_1           ;

        pxor        mm7,        mm7                 ;
        movq        mm6,        const53_2           ;

        movq        mm4,        round_values         ;
        lea         rdx,        [rsi+rcx-5]         ;
        horizontal_line_5_3_loop:

        movq        mm0,        QWORD PTR  [rsi]    ;
        00 01 02 03 04 05 06 07
        movq        mm1,        mm0                 ;
        00 01 02 03 04 05 06 07

        psllw       mm0,        8                   ;
        xx 00 xx 02 xx 04 xx 06
        psrlw       mm1,        8                   ;
        01 xx 03 xx 05 xx 07 xx

        psrlw       mm0,        8                   ;
        00 xx 02 xx 04 xx 06 xx
        psllq       mm1,        16                  ;
        xx xx 01 xx 03 xx 05 xx

        pmullw      mm0,        mm6

        pmullw      mm1,        mm5
        add         rsi,        5

        add         rdi,        3
        paddw       mm1,        mm0

        paddw       mm1,        mm4
        psrlw       mm1,        8

        cmp         rsi,        rdx
        packuswb    mm1,        mm7

        movd        DWORD PTR [rdi-3], mm1
        jl          horizontal_line_5_3_loop

//exit condition
        movq        mm0,        QWORD PTR  [rsi]    ;
        00 01 02 03 04 05 06 07
        movq        mm1,        mm0                 ;
        00 01 02 03 04 05 06 07

        psllw       mm0,        8                   ;
        xx 00 xx 02 xx 04 xx 06
        psrlw       mm1,        8                   ;
        01 xx 03 xx 05 xx 07 xx

        psrlw       mm0,        8                   ;
        00 xx 02 xx 04 xx 06 xx
        psllq       mm1,        16                  ;
        xx xx 01 xx 03 xx 05 xx

        pmullw      mm0,        mm6

        pmullw      mm1,        mm5
        paddw       mm1,        mm0

        paddw       mm1,        mm4
        psrlw       mm1,        8

        packuswb    mm1,        mm7
        movd        rax,        mm1

        mov         rdx,        rax
        shr         rdx,        16

        mov         WORD PTR[rdi],   ax
        mov         BYTE PTR[rdi+2], dl

    }

}

__declspec(align(16)) const static unsigned short one_thirds[] = {  85,  85,  85,  85 };
__declspec(align(16)) const static unsigned short two_thirds[] = { 171, 171, 171, 171 };

static
void vertical_band_5_3_scale_mmx
(
    unsigned char *source,
    unsigned int src_pitch,
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{

    __asm
    {

        mov         rsi,    source                    // Get the source and destination pointer
        mov         ecx,    src_pitch               // Get the pitch size

        mov         rdi,    dest                    // tow lines below
        pxor        mm7,    mm7                     // clear out mm7

        mov         edx,    dest_pitch               // Loop counter
        movq        mm5,    one_thirds

        movq        mm6,    two_thirds
        mov         ebx,    dest_width;

        vs_5_3_loop:

        movd        mm0,    DWORD ptr [rsi]         // src[0];
        movd        mm1,    DWORD ptr [rsi+rcx]     // src[1];

        movd        mm2,    DWORD ptr [rsi+rcx*2]
        lea         rax,    [rsi+rcx*2]             //

        punpcklbw   mm1,    mm7
        punpcklbw   mm2,    mm7

        pmullw      mm1,    mm5
        pmullw      mm2,    mm6

        movd        mm3,    DWORD ptr [rax+rcx]
        movd        mm4,    DWORD ptr [rax+rcx*2]

        punpcklbw   mm3,    mm7
        punpcklbw   mm4,    mm7

        pmullw      mm3,    mm6
        pmullw      mm4,    mm5


        movd        DWORD PTR [rdi], mm0
        paddw       mm1,    mm2

        paddw       mm1,    round_values
        psrlw       mm1,    8

        packuswb    mm1,    mm7
        paddw       mm3,    mm4

        paddw       mm3,    round_values
        movd        DWORD PTR [rdi+rdx], mm1

        psrlw       mm3,    8
        packuswb    mm3,    mm7

        movd        DWORD PTR [rdi+rdx*2], mm3


        add         rdi,    4
        add         rsi,    4

        sub         rbx,    4
        jg          vs_5_3_loop
    }
}




/****************************************************************************
*
*  ROUTINE       : horizontal_line_2_1_scale
*
*  INPUTS        : const unsigned char *source :
*                  unsigned int source_width    :
*                  unsigned char *dest         :
*                  unsigned int dest_width      :
*
*  OUTPUTS       : None.
*
*  RETURNS       : void
*
*  FUNCTION      : 1 to 2 up-scaling of a horizontal line of pixels.
*
*  SPECIAL NOTES : None.
*
****************************************************************************/
static
void horizontal_line_2_1_scale_mmx
(
    const unsigned char *source,
    unsigned int source_width,
    unsigned char *dest,
    unsigned int dest_width
)
{
    (void) dest_width;

    __asm
    {
        mov         rsi,    source
        mov         rdi,    dest

        pxor        mm7,    mm7
        mov         ecx,    dest_width

        xor         rdx,    rdx
        hs_2_1_loop:

        movq        mm0,    [rsi+rdx*2]
        psllw       mm0,    8

        psrlw       mm0,    8
        packuswb    mm0,    mm7

        movd        DWORD Ptr [rdi+rdx], mm0;
        add         rdx,    4

        cmp         rdx,    rcx
        jl          hs_2_1_loop

    }
}



static
void vertical_band_2_1_scale_mmx
(
    unsigned char *source,
    unsigned int src_pitch,
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width)
{
    vpx_memcpy(dest, source, dest_width);
}


__declspec(align(16)) const static unsigned short three_sixteenths[] = {  48,  48,  48,  48 };
__declspec(align(16)) const static unsigned short ten_sixteenths[]   = { 160, 160, 160, 160 };

static
void vertical_band_2_1_scale_i_mmx
(
    unsigned char *source,
    unsigned int src_pitch,
    unsigned char *dest,
    unsigned int dest_pitch,
    unsigned int dest_width
)
{
    __asm
    {
        mov         rsi,        source
        mov         rdi,        dest

        mov         eax,        src_pitch
        mov         edx,        dest_width

        pxor        mm7,        mm7
        sub         rsi,        rax             //back one line


        lea         rcx,        [rsi+rdx];
        movq        mm6,        round_values;

        movq        mm5,        three_sixteenths;
        movq        mm4,        ten_sixteenths;

        vs_2_1_i_loop:
        movd        mm0,        [rsi]           //
        movd        mm1,        [rsi+rax]       //

        movd        mm2,        [rsi+rax*2]     //
        punpcklbw   mm0,        mm7

        pmullw      mm0,        mm5
        punpcklbw   mm1,        mm7

        pmullw      mm1,        mm4
        punpcklbw   mm2,        mm7

        pmullw      mm2,        mm5
        paddw       mm0,        round_values

        paddw       mm1,        mm2
        paddw       mm0,        mm1

        psrlw       mm0,        8
        packuswb    mm0,        mm7

        movd        DWORD PTR [rdi],        mm0
        add         rsi,        4

        add         rdi,        4;
        cmp         rsi,        rcx
        jl          vs_2_1_i_loop

    }
}



void
register_mmxscalers(void)
{
    vp8_horizontal_line_1_2_scale        = horizontal_line_1_2_scale_mmx;
    vp8_horizontal_line_3_5_scale        = horizontal_line_3_5_scale_mmx;
    vp8_horizontal_line_4_5_scale        = horizontal_line_4_5_scale_mmx;
    vp8_vertical_band_1_2_scale          = vertical_band_1_2_scale_mmx;
    vp8_last_vertical_band_1_2_scale      = last_vertical_band_1_2_scale_mmx;
    vp8_vertical_band_3_5_scale          = vertical_band_3_5_scale_mmx;
    vp8_last_vertical_band_3_5_scale      = last_vertical_band_3_5_scale_mmx;
    vp8_vertical_band_4_5_scale          = vertical_band_4_5_scale_mmx;
    vp8_last_vertical_band_4_5_scale      = last_vertical_band_4_5_scale_mmx;

    vp8_vertical_band_5_4_scale           = vertical_band_5_4_scale_mmx;
    vp8_vertical_band_5_3_scale           = vertical_band_5_3_scale_mmx;
    vp8_vertical_band_2_1_scale           = vertical_band_2_1_scale_mmx;
    vp8_vertical_band_2_1_scale_i         = vertical_band_2_1_scale_i_mmx;
    vp8_horizontal_line_2_1_scale         = horizontal_line_2_1_scale_mmx;
    vp8_horizontal_line_5_3_scale         = horizontal_line_5_3_scale_mmx;
    vp8_horizontal_line_5_4_scale         = horizontal_line_5_4_scale_mmx;
}
