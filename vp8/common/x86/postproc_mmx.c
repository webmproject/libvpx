/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include <math.h>
#include <stdlib.h>
#include "vpx_scale/yv12config.h"
#include "pragmas.h"

#define VP8_FILTER_WEIGHT 128
#define VP8_FILTER_SHIFT  7



/* static constants */
__declspec(align(16))
const static short  Blur[48] =
{

    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    64, 64, 64, 64, 64, 64, 64, 64,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    0,  0,  0,  0,  0,  0,  0,  0,

};
#define RD  __declspec(align(16)) __int64 rd  = 0x0040004000400040;
#define R4D2 __declspec(align(16)) __int64 rd42[2] = {0x0004000400040004,0x0004000400040004};

#ifndef RELOCATEABLE
const static RD;
const static R4D2;
#endif


/* external references */
extern double vp8_gaussian(double sigma, double mu, double x);
extern short vp8_rv[];
extern int vp8_q2mbl(int x) ;



void vp8_post_proc_down_and_across_mmx
(
    unsigned char *src_ptr,
    unsigned char *dst_ptr,
    int src_pixels_per_line,
    int dst_pixels_per_line,
    int rows,
    int cols,
    int flimit
)
{
#ifdef RELOCATEABLE
    RD
    R4D2
#endif

    __asm
    {
        push        ebx
        lea         ebx, Blur
        movd        mm2, flimit
        punpcklwd   mm2, mm2
        punpckldq   mm2, mm2

        mov         esi,        src_ptr
        mov         edi,        dst_ptr

        mov         ecx, DWORD PTR rows
        mov         eax, src_pixels_per_line ;
        destination pitch?
        pxor        mm0, mm0              ;
        mm0 = 00000000

        nextrow:

        xor         edx,        edx       ;

        clear out edx for use as loop counter
        nextcol:

        pxor        mm7, mm7              ;

    mm7 = 00000000
    movq        mm6, [ebx + 32 ]      ;
        mm6 = kernel 2 taps
        movq        mm3, [esi]            ;
        mm4 = r0 p0..p7
        punpcklbw   mm3, mm0              ;
        mm3 = p0..p3
        movq        mm1, mm3              ;
        mm1 = p0..p3
        pmullw      mm3, mm6              ;
        mm3 *= kernel 2 modifiers

        movq        mm6, [ebx + 48]       ;
        mm6 = kernel 3 taps
        movq        mm5, [esi + eax]      ;
        mm4 = r1 p0..p7
        punpcklbw   mm5, mm0              ;
        mm5 = r1 p0..p3
        pmullw      mm6, mm5              ;
        mm6 *= p0..p3 * kernel 3 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm6

        ;
        thresholding
        movq        mm7, mm1              ;
        mm7 = r0 p0..p3
        psubusw     mm7, mm5              ;
        mm7 = r0 p0..p3 - r1 p0..p3
        psubusw     mm5, mm1              ;
        mm5 = r1 p0..p3 - r0 p0..p3
        paddusw     mm7, mm5              ;
        mm7 = abs(r0 p0..p3 - r1 p0..p3)
        pcmpgtw     mm7, mm2

        movq        mm6, [ebx + 64 ]      ;
        mm6 = kernel 4 modifiers
        movq        mm5, [esi + 2*eax]    ;
        mm4 = r2 p0..p7
        punpcklbw   mm5, mm0              ;
        mm5 = r2 p0..p3
        pmullw      mm6, mm5              ;
        mm5 *= kernel 4 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = r0 p0..p3
        psubusw     mm6, mm5              ;
        mm6 = r0 p0..p3 - r2 p0..p3
        psubusw     mm5, mm1              ;
        mm5 = r2 p0..p3 - r2 p0..p3
        paddusw     mm6, mm5              ;
        mm6 = abs(r0 p0..p3 - r2 p0..p3)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds


        neg         eax
        movq        mm6, [ebx ]           ;
        kernel 0 taps
        movq        mm5, [esi+2*eax]      ;
        mm4 = r-2 p0..p7
        punpcklbw   mm5, mm0              ;
        mm5 = r-2 p0..p3
        pmullw      mm6, mm5              ;
        mm5 *= kernel 0 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = r0 p0..p3
        psubusw     mm6, mm5              ;
        mm6 = p0..p3 - r-2 p0..p3
        psubusw     mm5, mm1              ;
        mm5 = r-2 p0..p3 - p0..p3
        paddusw     mm6, mm5              ;
        mm6 = abs(r0 p0..p3 - r-2 p0..p3)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds

        movq        mm6, [ebx + 16]       ;
        kernel 1 taps
        movq        mm4, [esi+eax]        ;
        mm4 = r-1 p0..p7
        punpcklbw   mm4, mm0              ;
        mm4 = r-1 p0..p3
        pmullw      mm6, mm4              ;
        mm4 *= kernel 1 modifiers.
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = r0 p0..p3
        psubusw     mm6, mm4              ;
        mm6 = p0..p3 - r-2 p0..p3
        psubusw     mm4, mm1              ;
        mm5 = r-1 p0..p3 - p0..p3
        paddusw     mm6, mm4              ;
        mm6 = abs(r0 p0..p3 - r-1 p0..p3)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds


        paddusw     mm3, rd               ;
        mm3 += round value
        psraw       mm3, VP8_FILTER_SHIFT     ;
        mm3 /= 128

        pand        mm1, mm7              ;
        mm1 select vals > thresh from source
        pandn       mm7, mm3              ;
        mm7 select vals < thresh from blurred result
        paddusw     mm1, mm7              ;
        combination

        packuswb    mm1, mm0              ;
        pack to bytes

        movd        [edi], mm1            ;
        neg         eax                   ;
        pitch is positive


        add         esi, 4
        add         edi, 4
        add         edx, 4

        cmp         edx, cols
        jl          nextcol
        // done with the all cols, start the across filtering in place
        sub         esi, edx
        sub         edi, edx


        push        eax
        xor         edx,    edx
        mov         eax,    [edi-4];

        acrossnextcol:
        pxor        mm7, mm7              ;
        mm7 = 00000000
        movq        mm6, [ebx + 32 ]      ;
        movq        mm4, [edi+edx]        ;
        mm4 = p0..p7
        movq        mm3, mm4              ;
        mm3 = p0..p7
        punpcklbw   mm3, mm0              ;
        mm3 = p0..p3
        movq        mm1, mm3              ;
        mm1 = p0..p3
        pmullw      mm3, mm6              ;
        mm3 *= kernel 2 modifiers

        movq        mm6, [ebx + 48]
        psrlq       mm4, 8                ;
        mm4 = p1..p7
        movq        mm5, mm4              ;
        mm5 = p1..p7
        punpcklbw   mm5, mm0              ;
        mm5 = p1..p4
        pmullw      mm6, mm5              ;
        mm6 *= p1..p4 * kernel 3 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm6

        ;
        thresholding
        movq        mm7, mm1              ;
        mm7 = p0..p3
        psubusw     mm7, mm5              ;
        mm7 = p0..p3 - p1..p4
        psubusw     mm5, mm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     mm7, mm5              ;
        mm7 = abs(p0..p3 - p1..p4)
        pcmpgtw     mm7, mm2

        movq        mm6, [ebx + 64 ]
        psrlq       mm4, 8                ;
        mm4 = p2..p7
        movq        mm5, mm4              ;
        mm5 = p2..p7
        punpcklbw   mm5, mm0              ;
        mm5 = p2..p5
        pmullw      mm6, mm5              ;
        mm5 *= kernel 4 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = p0..p3
        psubusw     mm6, mm5              ;
        mm6 = p0..p3 - p1..p4
        psubusw     mm5, mm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     mm6, mm5              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds


        movq        mm6, [ebx ]
        movq        mm4, [edi+edx-2]      ;
        mm4 = p-2..p5
        movq        mm5, mm4              ;
        mm5 = p-2..p5
        punpcklbw   mm5, mm0              ;
        mm5 = p-2..p1
        pmullw      mm6, mm5              ;
        mm5 *= kernel 0 modifiers
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = p0..p3
        psubusw     mm6, mm5              ;
        mm6 = p0..p3 - p1..p4
        psubusw     mm5, mm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     mm6, mm5              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds

        movq        mm6, [ebx + 16]
        psrlq       mm4, 8                ;
        mm4 = p-1..p5
        punpcklbw   mm4, mm0              ;
        mm4 = p-1..p2
        pmullw      mm6, mm4              ;
        mm4 *= kernel 1 modifiers.
        paddusw     mm3, mm6              ;
        mm3 += mm5

        ;
        thresholding
        movq        mm6, mm1              ;
        mm6 = p0..p3
        psubusw     mm6, mm4              ;
        mm6 = p0..p3 - p1..p4
        psubusw     mm4, mm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     mm6, mm4              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     mm6, mm2
        por         mm7, mm6              ;
        accumulate thresholds

        paddusw     mm3, rd               ;
        mm3 += round value
        psraw       mm3, VP8_FILTER_SHIFT     ;
        mm3 /= 128

        pand        mm1, mm7              ;
        mm1 select vals > thresh from source
        pandn       mm7, mm3              ;
        mm7 select vals < thresh from blurred result
        paddusw     mm1, mm7              ;
        combination

        packuswb    mm1, mm0              ;
        pack to bytes
        mov         DWORD PTR [edi+edx-4],  eax   ;
        store previous four bytes
        movd        eax,    mm1

        add         edx, 4
        cmp         edx, cols
        jl          acrossnextcol;

        mov         DWORD PTR [edi+edx-4],  eax
        pop         eax

        // done with this rwo
        add         esi, eax               ;
        next line
        mov         eax, dst_pixels_per_line ;
        destination pitch?
        add         edi, eax               ;
        next destination
        mov         eax, src_pixels_per_line ;
        destination pitch?

        dec         ecx                   ;
        decrement count
        jnz         nextrow               ;
        next row
        pop         ebx

    }
}



void vp8_post_proc_down_and_across_xmm
(
    unsigned char *src_ptr,
    unsigned char *dst_ptr,
    int src_pixels_per_line,
    int dst_pixels_per_line,
    int rows,
    int cols,
    int flimit
)
{
#ifdef RELOCATEABLE
    R4D2
#endif

    __asm
    {
        movd        xmm2,       flimit
        punpcklwd   xmm2,       xmm2
        punpckldq   xmm2,       xmm2
        punpcklqdq  xmm2,       xmm2

        mov         esi,        src_ptr
        mov         edi,        dst_ptr

        mov         ecx,        DWORD PTR rows
        mov         eax,        src_pixels_per_line ;
        destination pitch?
        pxor        xmm0,       xmm0              ;
        mm0 = 00000000

        nextrow:

        xor         edx,        edx       ;

        clear out edx for use as loop counter
        nextcol:
        movq        xmm3,       QWORD PTR [esi]         ;

        mm4 = r0 p0..p7
        punpcklbw   xmm3,       xmm0                    ;
        mm3 = p0..p3
        movdqa      xmm1,       xmm3                    ;
        mm1 = p0..p3
        psllw       xmm3,       2                       ;

        movq        xmm5,       QWORD PTR [esi + eax]   ;
        mm4 = r1 p0..p7
        punpcklbw   xmm5,       xmm0                    ;
        mm5 = r1 p0..p3
        paddusw     xmm3,       xmm5                    ;
        mm3 += mm6

        ;
        thresholding
        movdqa      xmm7,       xmm1                    ;
        mm7 = r0 p0..p3
        psubusw     xmm7,       xmm5                    ;
        mm7 = r0 p0..p3 - r1 p0..p3
        psubusw     xmm5,       xmm1                    ;
        mm5 = r1 p0..p3 - r0 p0..p3
        paddusw     xmm7,       xmm5                    ;
        mm7 = abs(r0 p0..p3 - r1 p0..p3)
        pcmpgtw     xmm7,       xmm2

        movq        xmm5,       QWORD PTR [esi + 2*eax] ;
        mm4 = r2 p0..p7
        punpcklbw   xmm5,       xmm0                    ;
        mm5 = r2 p0..p3
        paddusw     xmm3,       xmm5                    ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1                    ;
        mm6 = r0 p0..p3
        psubusw     xmm6,       xmm5                    ;
        mm6 = r0 p0..p3 - r2 p0..p3
        psubusw     xmm5,       xmm1                    ;
        mm5 = r2 p0..p3 - r2 p0..p3
        paddusw     xmm6,       xmm5                    ;
        mm6 = abs(r0 p0..p3 - r2 p0..p3)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6                    ;
        accumulate thresholds


        neg         eax
        movq        xmm5,       QWORD PTR [esi+2*eax]   ;
        mm4 = r-2 p0..p7
        punpcklbw   xmm5,       xmm0                    ;
        mm5 = r-2 p0..p3
        paddusw     xmm3,       xmm5                    ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1                    ;
        mm6 = r0 p0..p3
        psubusw     xmm6,       xmm5                    ;
        mm6 = p0..p3 - r-2 p0..p3
        psubusw     xmm5,       xmm1                    ;
        mm5 = r-2 p0..p3 - p0..p3
        paddusw     xmm6,       xmm5                    ;
        mm6 = abs(r0 p0..p3 - r-2 p0..p3)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6                    ;
        accumulate thresholds

        movq        xmm4,       QWORD PTR [esi+eax]     ;
        mm4 = r-1 p0..p7
        punpcklbw   xmm4,       xmm0                    ;
        mm4 = r-1 p0..p3
        paddusw     xmm3,       xmm4                    ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1                    ;
        mm6 = r0 p0..p3
        psubusw     xmm6,       xmm4                    ;
        mm6 = p0..p3 - r-2 p0..p3
        psubusw     xmm4,       xmm1                    ;
        mm5 = r-1 p0..p3 - p0..p3
        paddusw     xmm6,       xmm4                    ;
        mm6 = abs(r0 p0..p3 - r-1 p0..p3)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6                    ;
        accumulate thresholds


        paddusw     xmm3,       rd42                    ;
        mm3 += round value
        psraw       xmm3,       3                       ;
        mm3 /= 8

        pand        xmm1,       xmm7                    ;
        mm1 select vals > thresh from source
        pandn       xmm7,       xmm3                    ;
        mm7 select vals < thresh from blurred result
        paddusw     xmm1,       xmm7                    ;
        combination

        packuswb    xmm1,       xmm0                    ;
        pack to bytes
        movq        QWORD PTR [edi], xmm1             ;

        neg         eax                   ;
        pitch is positive
        add         esi,        8
        add         edi,        8

        add         edx,        8
        cmp         edx,        cols

        jl          nextcol

        // done with the all cols, start the across filtering in place
        sub         esi,        edx
        sub         edi,        edx

        xor         edx,        edx
        movq        mm0,        QWORD PTR [edi-8];

        acrossnextcol:
        movq        xmm7,       QWORD PTR [edi +edx -2]
        movd        xmm4,       DWORD PTR [edi +edx +6]

        pslldq      xmm4,       8
        por         xmm4,       xmm7

        movdqa      xmm3,       xmm4
        psrldq      xmm3,       2
        punpcklbw   xmm3,       xmm0              ;
        mm3 = p0..p3
        movdqa      xmm1,       xmm3              ;
        mm1 = p0..p3
        psllw       xmm3,       2


        movdqa      xmm5,       xmm4
        psrldq      xmm5,       3
        punpcklbw   xmm5,       xmm0              ;
        mm5 = p1..p4
        paddusw     xmm3,       xmm5              ;
        mm3 += mm6

        ;
        thresholding
        movdqa      xmm7,       xmm1              ;
        mm7 = p0..p3
        psubusw     xmm7,       xmm5              ;
        mm7 = p0..p3 - p1..p4
        psubusw     xmm5,       xmm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     xmm7,       xmm5              ;
        mm7 = abs(p0..p3 - p1..p4)
        pcmpgtw     xmm7,       xmm2

        movdqa      xmm5,       xmm4
        psrldq      xmm5,       4
        punpcklbw   xmm5,       xmm0              ;
        mm5 = p2..p5
        paddusw     xmm3,       xmm5              ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1              ;
        mm6 = p0..p3
        psubusw     xmm6,       xmm5              ;
        mm6 = p0..p3 - p1..p4
        psubusw     xmm5,       xmm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     xmm6,       xmm5              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6              ;
        accumulate thresholds


        movdqa      xmm5,       xmm4              ;
        mm5 = p-2..p5
        punpcklbw   xmm5,       xmm0              ;
        mm5 = p-2..p1
        paddusw     xmm3,       xmm5              ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1              ;
        mm6 = p0..p3
        psubusw     xmm6,       xmm5              ;
        mm6 = p0..p3 - p1..p4
        psubusw     xmm5,       xmm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     xmm6,       xmm5              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6              ;
        accumulate thresholds

        psrldq      xmm4,       1                   ;
        mm4 = p-1..p5
        punpcklbw   xmm4,       xmm0              ;
        mm4 = p-1..p2
        paddusw     xmm3,       xmm4              ;
        mm3 += mm5

        ;
        thresholding
        movdqa      xmm6,       xmm1              ;
        mm6 = p0..p3
        psubusw     xmm6,       xmm4              ;
        mm6 = p0..p3 - p1..p4
        psubusw     xmm4,       xmm1              ;
        mm5 = p1..p4 - p0..p3
        paddusw     xmm6,       xmm4              ;
        mm6 = abs(p0..p3 - p1..p4)
        pcmpgtw     xmm6,       xmm2
        por         xmm7,       xmm6              ;
        accumulate thresholds

        paddusw     xmm3,       rd42              ;
        mm3 += round value
        psraw       xmm3,       3                 ;
        mm3 /= 8

        pand        xmm1,       xmm7              ;
        mm1 select vals > thresh from source
        pandn       xmm7,       xmm3              ;
        mm7 select vals < thresh from blurred result
        paddusw     xmm1,       xmm7              ;
        combination

        packuswb    xmm1,       xmm0              ;
        pack to bytes
        movq        QWORD PTR [edi+edx-8],  mm0   ;
        store previous four bytes
        movdq2q     mm0,        xmm1

        add         edx,        8
        cmp         edx,        cols
        jl          acrossnextcol;

        // last 8 pixels
        movq        QWORD PTR [edi+edx-8],  mm0

        // done with this rwo
        add         esi, eax               ;
        next line
        mov         eax, dst_pixels_per_line ;
        destination pitch?
        add         edi, eax               ;
        next destination
        mov         eax, src_pixels_per_line ;
        destination pitch?

        dec         ecx                   ;
        decrement count
        jnz         nextrow               ;
        next row
    }
}


void vp8_mbpost_proc_down_mmx(unsigned char *dst, int pitch, int rows, int cols, int flimit)
{
    int c, i;
    __declspec(align(16))
    int flimit2[2];
    __declspec(align(16))
    unsigned char d[16][8];

    flimit = vp8_q2mbl(flimit);

    for (i = 0; i < 2; i++)
        flimit2[i] = flimit;

    rows += 8;

    for (c = 0; c < cols; c += 4)
    {
        unsigned char *s = &dst[c];

        __asm
        {
            mov         esi,        s           ;
            pxor        mm0,        mm0     ;

            mov         eax,        pitch       ;
            neg         eax                                     // eax = -pitch

            lea         esi,        [esi + eax*8];              // edi = s[-pitch*8]
            neg         eax


            pxor        mm5,        mm5
            pxor        mm6,        mm6     ;

            pxor        mm7,        mm7     ;
            mov         edi,        esi

            mov         ecx,        15          ;

            loop_initvar:
            movd        mm1,        DWORD PTR [edi];
            punpcklbw   mm1,        mm0     ;

            paddw       mm5,        mm1     ;
            pmullw      mm1,        mm1     ;

            movq        mm2,        mm1     ;
            punpcklwd   mm1,        mm0     ;

            punpckhwd   mm2,        mm0     ;
            paddd       mm6,        mm1     ;

            paddd       mm7,        mm2     ;
            lea         edi,        [edi+eax]   ;

            dec         ecx
            jne         loop_initvar
            //save the var and sum
            xor         edx,        edx
            loop_row:
            movd        mm1,        DWORD PTR [esi]     // [s-pitch*8]
            movd        mm2,        DWORD PTR [edi]     // [s+pitch*7]

            punpcklbw   mm1,        mm0
            punpcklbw   mm2,        mm0

            paddw       mm5,        mm2
            psubw       mm5,        mm1

            pmullw      mm2,        mm2
            movq        mm4,        mm2

            punpcklwd   mm2,        mm0
            punpckhwd   mm4,        mm0

            paddd       mm6,        mm2
            paddd       mm7,        mm4

            pmullw      mm1,        mm1
            movq        mm2,        mm1

            punpcklwd   mm1,        mm0
            psubd       mm6,        mm1

            punpckhwd   mm2,        mm0
            psubd       mm7,        mm2


            movq        mm3,        mm6
            pslld       mm3,        4

            psubd       mm3,        mm6
            movq        mm1,        mm5

            movq        mm4,        mm5
            pmullw      mm1,        mm1

            pmulhw      mm4,        mm4
            movq        mm2,        mm1

            punpcklwd   mm1,        mm4
            punpckhwd   mm2,        mm4

            movq        mm4,        mm7
            pslld       mm4,        4

            psubd       mm4,        mm7

            psubd       mm3,        mm1
            psubd       mm4,        mm2

            psubd       mm3,        flimit2
            psubd       mm4,        flimit2

            psrad       mm3,        31
            psrad       mm4,        31

            packssdw    mm3,        mm4
            packsswb    mm3,        mm0

            movd        mm1,        DWORD PTR [esi+eax*8]

            movq        mm2,        mm1
            punpcklbw   mm1,        mm0

            paddw       mm1,        mm5
            mov         ecx,        edx

            and         ecx,        127
            movq        mm4,        vp8_rv[ecx*2]

            paddw       mm1,        mm4
            //paddw     xmm1,       eight8s
            psraw       mm1,        4

            packuswb    mm1,        mm0
            pand        mm1,        mm3

            pandn       mm3,        mm2
            por         mm1,        mm3

            and         ecx,        15
            movd        DWORD PTR  d[ecx*4], mm1

            mov         ecx,        edx
            sub         ecx,        8

            and         ecx,        15
            movd        mm1,        DWORD PTR d[ecx*4]

            movd        [esi],      mm1
            lea         esi,        [esi+eax]

            lea         edi,        [edi+eax]
            add         edx,        1

            cmp         edx,        rows
            jl          loop_row

        }

    }
}

void vp8_mbpost_proc_down_xmm(unsigned char *dst, int pitch, int rows, int cols, int flimit)
{
    int c, i;
    __declspec(align(16))
    int flimit4[4];
    __declspec(align(16))
    unsigned char d[16][8];

    flimit = vp8_q2mbl(flimit);

    for (i = 0; i < 4; i++)
        flimit4[i] = flimit;

    rows += 8;

    for (c = 0; c < cols; c += 8)
    {
        unsigned char *s = &dst[c];

        __asm
        {
            mov         esi,        s           ;
            pxor        xmm0,       xmm0        ;

            mov         eax,        pitch       ;
            neg         eax                                     // eax = -pitch

            lea         esi,        [esi + eax*8];              // edi = s[-pitch*8]
            neg         eax


            pxor        xmm5,       xmm5
            pxor        xmm6,       xmm6        ;

            pxor        xmm7,       xmm7        ;
            mov         edi,        esi

            mov         ecx,        15          ;

            loop_initvar:
            movq        xmm1,       QWORD PTR [edi];
            punpcklbw   xmm1,       xmm0        ;

            paddw       xmm5,       xmm1        ;
            pmullw      xmm1,       xmm1        ;

            movdqa      xmm2,       xmm1        ;
            punpcklwd   xmm1,       xmm0        ;

            punpckhwd   xmm2,       xmm0        ;
            paddd       xmm6,       xmm1        ;

            paddd       xmm7,       xmm2        ;
            lea         edi,        [edi+eax]   ;

            dec         ecx
            jne         loop_initvar
            //save the var and sum
            xor         edx,        edx
            loop_row:
            movq        xmm1,       QWORD PTR [esi]     // [s-pitch*8]
            movq        xmm2,       QWORD PTR [edi]     // [s+pitch*7]

            punpcklbw   xmm1,       xmm0
            punpcklbw   xmm2,       xmm0

            paddw       xmm5,       xmm2
            psubw       xmm5,       xmm1

            pmullw      xmm2,       xmm2
            movdqa      xmm4,       xmm2

            punpcklwd   xmm2,       xmm0
            punpckhwd   xmm4,       xmm0

            paddd       xmm6,       xmm2
            paddd       xmm7,       xmm4

            pmullw      xmm1,       xmm1
            movdqa      xmm2,       xmm1

            punpcklwd   xmm1,       xmm0
            psubd       xmm6,       xmm1

            punpckhwd   xmm2,       xmm0
            psubd       xmm7,       xmm2


            movdqa      xmm3,       xmm6
            pslld       xmm3,       4

            psubd       xmm3,       xmm6
            movdqa      xmm1,       xmm5

            movdqa      xmm4,       xmm5
            pmullw      xmm1,       xmm1

            pmulhw      xmm4,       xmm4
            movdqa      xmm2,       xmm1

            punpcklwd   xmm1,       xmm4
            punpckhwd   xmm2,       xmm4

            movdqa      xmm4,       xmm7
            pslld       xmm4,       4

            psubd       xmm4,       xmm7

            psubd       xmm3,       xmm1
            psubd       xmm4,       xmm2

            psubd       xmm3,       flimit4
            psubd       xmm4,       flimit4

            psrad       xmm3,       31
            psrad       xmm4,       31

            packssdw    xmm3,       xmm4
            packsswb    xmm3,       xmm0

            movq        xmm1,       QWORD PTR [esi+eax*8]

            movq        xmm2,       xmm1
            punpcklbw   xmm1,       xmm0

            paddw       xmm1,       xmm5
            mov         ecx,        edx

            and         ecx,        127
            movdqu      xmm4,       vp8_rv[ecx*2]

            paddw       xmm1,       xmm4
            //paddw     xmm1,       eight8s
            psraw       xmm1,       4

            packuswb    xmm1,       xmm0
            pand        xmm1,       xmm3

            pandn       xmm3,       xmm2
            por         xmm1,       xmm3

            and         ecx,        15
            movq        QWORD PTR  d[ecx*8], xmm1

            mov         ecx,        edx
            sub         ecx,        8

            and         ecx,        15
            movq        mm0,        d[ecx*8]

            movq        [esi],      mm0
            lea         esi,        [esi+eax]

            lea         edi,        [edi+eax]
            add         edx,        1

            cmp         edx,        rows
            jl          loop_row

        }

    }
}
#if 0
/****************************************************************************
 *
 *  ROUTINE       : plane_add_noise_wmt
 *
 *  INPUTS        : unsigned char *Start    starting address of buffer to add gaussian
 *                                  noise to
 *                  unsigned int Width    width of plane
 *                  unsigned int Height   height of plane
 *                  int  Pitch    distance between subsequent lines of frame
 *                  int  q        quantizer used to determine amount of noise
 *                                  to add
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void.
 *
 *  FUNCTION      : adds gaussian noise to a plane of pixels
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void vp8_plane_add_noise_wmt(unsigned char *Start, unsigned int Width, unsigned int Height, int Pitch, int q, int a)
{
    unsigned int i;

    __declspec(align(16)) unsigned char blackclamp[16];
    __declspec(align(16)) unsigned char whiteclamp[16];
    __declspec(align(16)) unsigned char bothclamp[16];
    char char_dist[300];
    char Rand[2048];
    double sigma;
//    return;
    __asm emms
    sigma = a + .5 + .6 * (63 - q) / 63.0;

    // set up a lookup table of 256 entries that matches
    // a gaussian distribution with sigma determined by q.
    //
    {
        double i;
        int next, j;

        next = 0;

        for (i = -32; i < 32; i++)
        {
            double g = 256 * vp8_gaussian(sigma, 0, 1.0 * i);
            int a = (int)(g + .5);

            if (a)
            {
                for (j = 0; j < a; j++)
                {
                    char_dist[next+j] = (char) i;
                }

                next = next + j;
            }

        }

        for (next = next; next < 256; next++)
            char_dist[next] = 0;

    }

    for (i = 0; i < 2048; i++)
    {
        Rand[i] = char_dist[rand() & 0xff];
    }

    for (i = 0; i < 16; i++)
    {
        blackclamp[i] = -char_dist[0];
        whiteclamp[i] = -char_dist[0];
        bothclamp[i] = -2 * char_dist[0];
    }

    for (i = 0; i < Height; i++)
    {
        unsigned char *Pos = Start + i * Pitch;
        char  *Ref = Rand + (rand() & 0xff);

        __asm
        {
            mov ecx, [Width]
            mov esi, Pos
            mov edi, Ref
            xor         eax, eax

            nextset:
            movdqu      xmm1, [esi+eax]        // get the source

            psubusb     xmm1, blackclamp       // clamp both sides so we don't outrange adding noise
            paddusb     xmm1, bothclamp
            psubusb     xmm1, whiteclamp

            movdqu      xmm2, [edi+eax]        // get the noise for this line
            paddb       xmm1, xmm2             // add it in
            movdqu      [esi+eax], xmm1        // store the result

            add         eax, 16                // move to the next line

            cmp         eax, ecx
            jl          nextset


        }

    }
}
#endif
__declspec(align(16))
static const int four8s[4] = { 8, 8, 8, 8};
void vp8_mbpost_proc_across_ip_xmm(unsigned char *src, int pitch, int rows, int cols, int flimit)
{
    int r, i;
    __declspec(align(16))
    int flimit4[4];
    unsigned char *s = src;
    int sumsq;
    int sum;


    flimit = vp8_q2mbl(flimit);
    flimit4[0] =
        flimit4[1] =
            flimit4[2] =
                flimit4[3] = flimit;

    for (r = 0; r < rows; r++)
    {


        sumsq = 0;
        sum = 0;

        for (i = -8; i <= 6; i++)
        {
            sumsq += s[i] * s[i];
            sum   += s[i];
        }

        __asm
        {
            mov         eax,    sumsq
            movd        xmm7,   eax

            mov         eax,    sum
            movd        xmm6,   eax

            mov         esi,    s
            xor         ecx,    ecx

            mov         edx,    cols
            add         edx,    8
            pxor        mm0,    mm0
            pxor        mm1,    mm1

            pxor        xmm0,   xmm0
            nextcol4:

            movd        xmm1,   DWORD PTR [esi+ecx-8]   // -8 -7 -6 -5
            movd        xmm2,   DWORD PTR [esi+ecx+7]   // +7 +8 +9 +10

            punpcklbw   xmm1,   xmm0                    // expanding
            punpcklbw   xmm2,   xmm0                    // expanding

            punpcklwd   xmm1,   xmm0                    // expanding to dwords
            punpcklwd   xmm2,   xmm0                    // expanding to dwords

            psubd       xmm2,   xmm1                    // 7--8   8--7   9--6 10--5
            paddd       xmm1,   xmm1                    // -8*2   -7*2   -6*2 -5*2

            paddd       xmm1,   xmm2                    // 7+-8   8+-7   9+-6 10+-5
            pmaddwd     xmm1,   xmm2                    // squared of 7+-8   8+-7   9+-6 10+-5

            paddd       xmm6,   xmm2
            paddd       xmm7,   xmm1

            pshufd      xmm6,   xmm6,   0               // duplicate the last ones
            pshufd      xmm7,   xmm7,   0               // duplicate the last ones

            psrldq      xmm1,       4                   // 8--7   9--6 10--5  0000
            psrldq      xmm2,       4                   // 8--7   9--6 10--5  0000

            pshufd      xmm3,   xmm1,   3               // 0000  8--7   8--7   8--7 squared
            pshufd      xmm4,   xmm2,   3               // 0000  8--7   8--7   8--7 squared

            paddd       xmm6,   xmm4
            paddd       xmm7,   xmm3

            pshufd      xmm3,   xmm1,   01011111b       // 0000  0000   9--6   9--6 squared
            pshufd      xmm4,   xmm2,   01011111b       // 0000  0000   9--6   9--6 squared

            paddd       xmm7,   xmm3
            paddd       xmm6,   xmm4

            pshufd      xmm3,   xmm1,   10111111b       // 0000  0000   8--7   8--7 squared
            pshufd      xmm4,   xmm2,   10111111b       // 0000  0000   8--7   8--7 squared

            paddd       xmm7,   xmm3
            paddd       xmm6,   xmm4

            movdqa      xmm3,   xmm6
            pmaddwd     xmm3,   xmm3

            movdqa      xmm5,   xmm7
            pslld       xmm5,   4

            psubd       xmm5,   xmm7
            psubd       xmm5,   xmm3

            psubd       xmm5,   flimit4
            psrad       xmm5,   31

            packssdw    xmm5,   xmm0
            packsswb    xmm5,   xmm0

            movd        xmm1,   DWORD PTR [esi+ecx]
            movq        xmm2,   xmm1

            punpcklbw   xmm1,   xmm0
            punpcklwd   xmm1,   xmm0

            paddd       xmm1,   xmm6
            paddd       xmm1,   four8s

            psrad       xmm1,   4
            packssdw    xmm1,   xmm0

            packuswb    xmm1,   xmm0
            pand        xmm1,   xmm5

            pandn       xmm5,   xmm2
            por         xmm5,   xmm1

            movd        [esi+ecx-8],  mm0
            movq        mm0,    mm1

            movdq2q     mm1,    xmm5
            psrldq      xmm7,   12

            psrldq      xmm6,   12
            add         ecx,    4

            cmp         ecx,    edx
            jl          nextcol4

        }
        s += pitch;
    }
}

#if 0

/****************************************************************************
 *
 *  ROUTINE       : plane_add_noise_mmx
 *
 *  INPUTS        : unsigned char *Start    starting address of buffer to add gaussian
 *                                  noise to
 *                  unsigned int Width    width of plane
 *                  unsigned int Height   height of plane
 *                  int  Pitch    distance between subsequent lines of frame
 *                  int  q        quantizer used to determine amount of noise
 *                                  to add
 *
 *  OUTPUTS       : None.
 *
 *  RETURNS       : void.
 *
 *  FUNCTION      : adds gaussian noise to a plane of pixels
 *
 *  SPECIAL NOTES : None.
 *
 ****************************************************************************/
void vp8_plane_add_noise_mmx(unsigned char *Start, unsigned int Width, unsigned int Height, int Pitch, int q, int a)
{
    unsigned int i;
    int Pitch4 = Pitch * 4;
    const int noise_amount = 2;
    const int noise_adder = 2 * noise_amount + 1;

    __declspec(align(16)) unsigned char blackclamp[16];
    __declspec(align(16)) unsigned char whiteclamp[16];
    __declspec(align(16)) unsigned char bothclamp[16];

    char char_dist[300];
    char Rand[2048];

    double sigma;
    __asm emms
    sigma = a + .5 + .6 * (63 - q) / 63.0;

    // set up a lookup table of 256 entries that matches
    // a gaussian distribution with sigma determined by q.
    //
    {
        double i, sum = 0;
        int next, j;

        next = 0;

        for (i = -32; i < 32; i++)
        {
            int a = (int)(.5 + 256 * vp8_gaussian(sigma, 0, i));

            if (a)
            {
                for (j = 0; j < a; j++)
                {
                    char_dist[next+j] = (char) i;
                }

                next = next + j;
            }

        }

        for (next = next; next < 256; next++)
            char_dist[next] = 0;

    }

    for (i = 0; i < 2048; i++)
    {
        Rand[i] = char_dist[rand() & 0xff];
    }

    for (i = 0; i < 16; i++)
    {
        blackclamp[i] = -char_dist[0];
        whiteclamp[i] = -char_dist[0];
        bothclamp[i] = -2 * char_dist[0];
    }

    for (i = 0; i < Height; i++)
    {
        unsigned char *Pos = Start + i * Pitch;
        char  *Ref = Rand + (rand() & 0xff);

        __asm
        {
            mov ecx, [Width]
            mov esi, Pos
            mov edi, Ref
            xor         eax, eax

            nextset:
            movq        mm1, [esi+eax]        // get the source

            psubusb     mm1, blackclamp       // clamp both sides so we don't outrange adding noise
            paddusb     mm1, bothclamp
            psubusb     mm1, whiteclamp

            movq        mm2, [edi+eax]        // get the noise for this line
            paddb       mm1, mm2             // add it in
            movq        [esi+eax], mm1        // store the result

            add         eax, 8                // move to the next line

            cmp         eax, ecx
            jl          nextset


        }

    }
}
#else
extern char an[8][64][3072];
extern int cd[8][64];

void vp8_plane_add_noise_mmx(unsigned char *Start, unsigned int Width, unsigned int Height, int Pitch, int q, int a)
{
    unsigned int i;
    __declspec(align(16)) unsigned char blackclamp[16];
    __declspec(align(16)) unsigned char whiteclamp[16];
    __declspec(align(16)) unsigned char bothclamp[16];


    __asm emms

    for (i = 0; i < 16; i++)
    {
        blackclamp[i] = -cd[a][q];
        whiteclamp[i] = -cd[a][q];
        bothclamp[i] = -2 * cd[a][q];
    }

    for (i = 0; i < Height; i++)
    {
        unsigned char *Pos = Start + i * Pitch;
        char  *Ref = an[a][q] + (rand() & 0xff);

        __asm
        {
            mov ecx, [Width]
            mov esi, Pos
            mov edi, Ref
            xor         eax, eax

            nextset:
            movq        mm1, [esi+eax]        // get the source

            psubusb     mm1, blackclamp       // clamp both sides so we don't outrange adding noise
            paddusb     mm1, bothclamp
            psubusb     mm1, whiteclamp

            movq        mm2, [edi+eax]        // get the noise for this line
            paddb       mm1, mm2             // add it in
            movq        [esi+eax], mm1        // store the result

            add         eax, 8                // move to the next line

            cmp         eax, ecx
            jl          nextset
        }
    }
}


void vp8_plane_add_noise_wmt(unsigned char *Start, unsigned int Width, unsigned int Height, int Pitch, int q, int a)
{
    unsigned int i;

    __declspec(align(16)) unsigned char blackclamp[16];
    __declspec(align(16)) unsigned char whiteclamp[16];
    __declspec(align(16)) unsigned char bothclamp[16];

    __asm emms

    for (i = 0; i < 16; i++)
    {
        blackclamp[i] = -cd[a][q];
        whiteclamp[i] = -cd[a][q];
        bothclamp[i] = -2 * cd[a][q];
    }

    for (i = 0; i < Height; i++)
    {
        unsigned char *Pos = Start + i * Pitch;
        char *Ref = an[a][q] + (rand() & 0xff);

        __asm
        {
            mov ecx,    [Width]
            mov esi,    Pos
            mov edi,    Ref
            xor         eax, eax

            nextset:
            movdqu      xmm1, [esi+eax]        // get the source

            psubusb     xmm1, blackclamp       // clamp both sides so we don't outrange adding noise
            paddusb     xmm1, bothclamp
            psubusb     xmm1, whiteclamp

            movdqu      xmm2, [edi+eax]        // get the noise for this line
            paddb       xmm1, xmm2             // add it in
            movdqu      [esi+eax], xmm1        // store the result

            add         eax, 16                // move to the next line

            cmp         eax, ecx
            jl          nextset
        }
    }
}

#endif
