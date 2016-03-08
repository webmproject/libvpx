;
;  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%include "third_party/x86inc/x86inc.asm"

SECTION_RODATA
align 16
tfe:
    times 8 db 0xfe
t80:
    times 8 db 0x80
t3:
    times 8 db 0x03
t4:
    times 8 db 0x04
ones:
    times 4 dw 0x0001

SECTION .text

%define stkreg rsp

%define t0                  0
%define t1            t0 + 16
%define p1            t1 + 16
%define p0            p1 + 16
%define q0            p0 + 16
%define q1            q0 + 16
%define lstacksize    q1 + 16

%define goffsetq _limitq

;void vpx_lpf_horizontal_4_mmx(unsigned char *src_ptr, int  src_pixel_step,
;                              const char *blimit, const char *limit,
;                              const char *thresh);
INIT_MMX mmx
cglobal lpf_horizontal_4, 5, 6, 8, 0 - lstacksize, \
                                s, p, _blimit, _limit, _thresh, s1
    movq                  m7, [_limitq]
    GET_GOT         goffsetq
%if GET_GOT_DEFINED=1
    add rsp, gprsize                          ; restore stack
%endif
    lea                  s1q, [sq + pq]       ; s1q points to row +1

    ; calculate breakout conditions
    movq                  m2, [s1q + 2 * pq]  ; q3
    movq                  m1, [ sq + 2 * pq]  ; q2
    movq                  m6, m1              ; q2
    psubusb               m1, m2              ; q2-=q3
    psubusb               m2, m6              ; q3-=q2
    por                   m1, m2              ; abs(q3-q2)
    psubusb               m1, m7
    movq                  m4, [sq + pq]       ; q1
    movq                  m3, m4              ; q1
    psubusb               m4, m6              ; q1-=q2
    psubusb               m6, m3              ; q2-=q1
    por                   m4, m6              ; abs(q2-q1)
    psubusb               m4, m7
    por                   m1, m4
    movq                  m4, [sq]            ; q0
    movq                  m0, m4              ; q0
    psubusb               m4, m3              ; q0-=q1
    psubusb               m3, m0              ; q1-=q0
    por                   m4, m3              ; abs(q0-q1)
    movq       [stkreg + t0], m4              ; save to t0
    psubusb               m4, m7
    por                   m1, m4
    neg                   pq                  ; negate pitch to deal with
                                              ; above border
    movq                  m2, [ sq + 4 * pq]  ; p3
    movq                  m4, [s1q + 4 * pq]  ; p2
    movq                  m5, m4              ; p2
    psubusb               m4, m2              ; p2-=p3
    psubusb               m2, m5              ; p3-=p2
    por                   m4, m2              ; abs(p3 - p2)
    psubusb               m4, m7
    por                   m1, m4
    movq                  m4, [sq + 2 * pq]   ; p1
    movq                  m3, m4              ; p1
    psubusb               m4, m5              ; p1-=p2
    psubusb               m5, m3              ; p2-=p1
    por                   m4, m5              ; abs(p2 - p1)
    psubusb               m4, m7
    por                   m1, m4
    movq                  m2, m3              ; p1
    movq                  m4, [sq + pq]       ; p0
    movq                  m5, m4              ; p0
    psubusb               m4, m3              ; p0-=p1
    psubusb               m3, m5              ; p1-=p0
    por                   m4, m3              ; abs(p1 - p0)
    movq       [stkreg + t1], m4              ; save to t1
    psubusb               m4, m7
    por                   m1, m4
    movq                  m3, [s1q]           ; q1
    movq                  m4, m3              ; q1
    psubusb               m3, m2              ; q1-=p1
    psubusb               m2, m4              ; p1-=q1
    por                   m2, m3              ; abs(p1-q1)
    pand                  m2, [GLOBAL(tfe)]   ; set lsb of each byte to zero
    psrlw                 m2, 1               ; abs(p1-q1)/2
    movq                  m6, m5              ; p0
    movq                  m3, [sq]            ; q0
    psubusb               m5, m3              ; p0-=q0
    psubusb               m3, m6              ; q0-=p0
    por                   m5, m3              ; abs(p0 - q0)
    paddusb               m5, m5              ; abs(p0-q0)*2
    paddusb               m5, m2              ; abs (p0 - q0) * 2 + abs(p1-q1)/2
    movq                  m7, [_blimitq]            ; blimit
    psubusb               m5, m7              ; abs (p0 - q0) * 2 +
                                              ; abs(p1-q1)/2  > blimit
    por                   m1, m5
    pxor                  m5, m5
    pcmpeqb               m1, m5              ; mask m1

    ; calculate high edge variance
    movq                  m7, [_threshq]
    movq                  m4, [stkreg + t0]   ; get abs (q1 - q0)
    psubusb               m4, m7
    movq                  m3, [stkreg + t1]   ; get abs (p1 - p0)
    psubusb               m3, m7
    paddb                 m4, m3              ; abs(q1 - q0) > thresh ||
                                              ; abs(p1 - p0) > thresh
    pcmpeqb               m4, m5
    pcmpeqb               m5, m5
    movq                  m3, [GLOBAL(t80)]
    pxor                  m4, m5

    ; start work on filters
    movq                  m2, [sq + 2 * pq]   ; p1
    movq                  m7, [s1q]           ; q1
    pxor                  m2, m3              ; p1 converted to signed values
    pxor                  m7, m3              ; q1 converted to signed values
    psubsb                m2, m7              ; p1 - q1
    pand                  m2, m4              ; high var mask (hvm)(p1 - q1)
    pxor                  m6, m3              ; p0 converted to signed values
    pxor                  m0, m3              ; q0 converted to signed values
    movq                  m3, m0              ; q0
    psubsb                m0, m6              ; q0 - p0
    paddsb                m2, m0              ; 1 * (q0 - p0) + hvm(p1 - q1)
    paddsb                m2, m0              ; 2 * (q0 - p0) + hvm(p1 - q1)
    paddsb                m2, m0              ; 3 * (q0 - p0) + hvm(p1 - q1)
    pand                  m1, m2              ; mask filter values we don't
                                              ; care about
    movq                  m2, m1
    paddsb                m1, [GLOBAL(t4)]    ; 3* (q0 - p0) + hvm(p1 - q1) + 4
    paddsb                m2, [GLOBAL(t3)]    ; 3* (q0 - p0) + hvm(p1 - q1) + 3

    pxor                  m0, m0
    pxor                  m5, m5
    punpcklbw             m0, m2
    punpckhbw             m5, m2
    psraw                 m0, 11
    psraw                 m5, 11
    packsswb              m0, m5
    movq                  m2, m0              ; (3* (q0 - p0) + hvm(p1 - q1)
                                              ; + 3) >> 3;
    pxor                  m0, m0
    movq                  m5, m1              ; abcdefgh
    punpcklbw             m0, m1              ; e0f0g0h0
    psraw                 m0, 11              ; sign extended shift right by 3
    pxor                  m1, m1
    punpckhbw             m1, m5              ; a0b0c0d0
    psraw                 m1, 11              ; sign extended shift right by 3
    movq                  m5, m0              ; save results

    packsswb              m0, m1              ; (3* (q0 - p0) + hvm(p1 - q1)
                                              ; + 4) >>3
    paddsw                m5, [GLOBAL(ones)]
    paddsw                m1, [GLOBAL(ones)]
    psraw                 m5, 1
    psraw                 m1, 1
    packsswb              m5, m1              ; (3* (q0 - p0) + hvm(p1 - q1)
                                              ; + 4) >>4
    movq                  m1, [GLOBAL(t80)]
    pandn                 m4, m5              ; high edge variance additive
    paddsb                m6, m2              ; p0+= p0 add
    pxor                  m6, m1              ; unoffset
    movq           [sq + pq], m6              ; write back
    movq                  m6, [sq + 2 * pq]   ; p1
    pxor                  m6, m1              ; reoffset
    paddsb                m6, m4              ; p1+= p1 add
    pxor                  m6, m1              ; unoffset
    movq       [sq + 2 * pq], m6              ; write back
    psubsb                m3, m0              ; q0-= q0 add
    pxor                  m3, m1              ; unoffset
    movq                [sq], m3              ; write back
    psubsb                m7, m4              ; q1-= q1 add
    pxor                  m7, m1              ; unoffset
    movq               [s1q], m7              ; write back
    RET

;void vpx_lpf_vertical_4_mmx(unsigned char *src_ptr, int  src_pixel_step,
;                            const char *blimit, const char *limit,
;                            const char *thresh);
INIT_MMX mmx
cglobal lpf_vertical_4, 5, 6, 8, 0 - lstacksize, \
                              s, p, _blimit, _limit, _thresh, s1
    lea                   sq, [sq + pq * 4 - 4]
    lea                  s1q, [sq + pq]       ; s1q points to row +1
    ;transpose
    movq                  m6, [ sq + 2 * pq]  ; 67 66 65 64 63 62 61 60
    movq                  m7, m6              ; 77 76 75 74 73 72 71 70
    punpckhbw             m7, [s1q + 2 * pq]  ; 77 67 76 66 75 65 74 64
    punpcklbw             m6, [s1q + 2 * pq]  ; 73 63 72 62 71 61 70 60
    movq                  m4, [sq]            ; 47 46 45 44 43 42 41 40
    movq                  m5, m4              ; 47 46 45 44 43 42 41 40
    punpckhbw             m5, [sq + pq]       ; 57 47 56 46 55 45 54 44
    punpcklbw             m4, [sq + pq]       ; 53 43 52 42 51 41 50 40
    movq                  m3, m5              ; 57 47 56 46 55 45 54 44
    punpckhwd             m5, m7              ; 77 67 57 47 76 66 56 46
    punpcklwd             m3, m7              ; 75 65 55 45 74 64 54 44
    movq                  m2, m4              ; 53 43 52 42 51 41 50 40
    punpckhwd             m4, m6              ; 73 63 53 43 72 62 52 42
    punpcklwd             m2, m6              ; 71 61 51 41 70 60 50 40
    neg                   pq
    movq                  m6, [ sq + pq * 2]  ; 27 26 25 24 23 22 21 20
    movq                  m1, m6              ; 27 26 25 24 23 22 21 20
    punpckhbw             m6, [ sq + pq    ]  ; 37 27 36 36 35 25 34 24
    punpcklbw             m1, [ sq + pq    ]  ; 33 23 32 22 31 21 30 20
    movq                  m7, [ sq + pq * 4]; ; 07 06 05 04 03 02 01 00
    punpckhbw             m7, [s1q + pq * 4]  ; 17 07 16 06 15 05 14 04
    movq                  m0, m7              ; 17 07 16 06 15 05 14 04
    punpckhwd             m7, m6              ; 37 27 17 07 36 26 16 06
    punpcklwd             m0, m6              ; 35 25 15 05 34 24 14 04
    movq                  m6, m7              ; 37 27 17 07 36 26 16 06
    punpckhdq             m7, m5              ; 77 67 57 47 37 27 17 07  = q3
    punpckldq             m6, m5              ; 76 66 56 46 36 26 16 06  = q2
    movq                  m5, m6              ; 76 66 56 46 36 26 16 06
    psubusb               m5, m7              ; q2-q3
    psubusb               m7, m6              ; q3-q2
    por                   m7, m5;             ; m7=abs (q3-q2)
    movq                  m5, m0              ; 35 25 15 05 34 24 14 04
    punpckhdq             m5, m3              ; 75 65 55 45 35 25 15 05 = q1
    punpckldq             m0, m3              ; 74 64 54 44 34 24 15 04 = q0
    movq                  m3, m5              ; 75 65 55 45 35 25 15 05 = q1
    psubusb               m3, m6              ; q1-q2
    psubusb               m6, m5              ; q2-q1
    por                   m6, m3              ; m6=abs(q2-q1)

    movq       [stkreg + q1], m5              ; save q1
    movq       [stkreg + q0], m0              ; save q0

    movq                  m3, [ sq + pq * 4]  ; 07 06 05 04 03 02 01 00
    punpcklbw             m3, [s1q + pq * 4]  ; 13 03 12 02 11 01 10 00
    movq                  m0, m3              ; 13 03 12 02 11 01 10 00
    punpcklwd             m0, m1              ; 31 21 11 01 30 20 10 00
    punpckhwd             m3, m1              ; 33 23 13 03 32 22 12 02
    movq                  m1, m0              ; 31 21 11 01 30 20 10 00
    punpckldq             m0, m2              ; 70 60 50 40 30 20 10 00  =p3
    punpckhdq             m1, m2              ; 71 61 51 41 31 21 11 01  =p2
    movq                  m2, m1              ; 71 61 51 41 31 21 11 01  =p2
    psubusb               m2, m0              ; p2-p3
    psubusb               m0, m1              ; p3-p2
    por                   m0, m2              ; m0=abs(p3-p2)
    movq                  m2, m3              ; 33 23 13 03 32 22 12 02
    punpckldq             m2, m4              ; 72 62 52 42 32 22 12 02 = p1
    punpckhdq             m3, m4              ; 73 63 53 43 33 23 13 03 = p0

    movq       [stkreg + p0], m3              ; save p0
    movq       [stkreg + p1], m2              ; save p1
    movq                  m5, m2              ; m5 = p1
    psubusb               m2, m1              ; p1-p2
    psubusb               m1, m5              ; p2-p1
    por                   m1, m2              ; m1=abs(p2-p1)
    movq                  m4, [_limitq]
    GET_GOT         goffsetq
%if GET_GOT_DEFINED=1
    add rsp, gprsize                          ; restore stack
%endif
    psubusb               m7, m4
    psubusb               m0, m4
    psubusb               m1, m4
    psubusb               m6, m4
    por                   m7, m6
    por                   m0, m1
    por                   m0, m7              ; abs(q3-q2) > limit ||
                                              ; abs(p3-p2) > limit ||
                                              ; abs(p2-p1) > limit ||
                                              ; abs(q2-q1) > limit
    movq                  m1, m5              ; p1
    movq                  m7, m3              ; m3=m7=p0
    psubusb               m7, m5              ; p0 - p1
    psubusb               m5, m3              ; p1 - p0
    por                   m5, m7              ; abs(p1-p0)
    movq       [stkreg + t0], m5              ; save abs(p1-p0)
    psubusb               m5, m4
    por                   m0, m5              ; m0=mask
    movq                  m5, [stkreg + q0]   ; m5=q0
    movq                  m7, [stkreg + q1]   ; m7=q1
    movq                  m6, m5              ; m6=q0
    movq                  m2, m7              ; q1
    psubusb               m5, m7              ; q0-q1
    psubusb               m7, m6              ; q1-q0
    por                   m7, m5              ; abs(q1-q0)
    movq       [stkreg + t1], m7              ; save abs(q1-q0)
    psubusb               m7, m4
    por                   m0, m7              ; mask
    movq                  m5, m2              ; q1
    psubusb               m5, m1              ; q1-=p1
    psubusb               m1, m2              ; p1-=q1
    por                   m5, m1              ; abs(p1-q1)
    pand                  m5, [GLOBAL(tfe)]   ; set lsb of each byte to zero
    psrlw                 m5, 1               ; abs(p1-q1)/2
    movq                  m4, [_blimitq]
    movq                  m1, m3              ; m1=m3=p0
    movq                  m7, m6              ; m7=m6=q0
    psubusb               m1, m7              ; p0-q0
    psubusb               m7, m3              ; q0-p0
    por                   m1, m7              ; abs(q0-p0)
    paddusb               m1, m1              ; abs(q0-p0)*2
    paddusb               m1, m5              ; abs(p0 - q0)*2 + abs(p1-q1)/2
    psubusb               m1, m4              ; abs(p0 - q0)*2 + abs(p1-q1)/2
                                              ; > blimit
    por                   m1, m0;             ; mask
    pxor                  m0, m0
    pcmpeqb               m1, m0

    ; calculate high edge variance
    movq                  m7, [_threshq]
    movq                  m4, [stkreg + t0]   ; get abs (q1 - q0)
    psubusb               m4, m7
    movq                  m3, [stkreg + t1]   ; get abs (p1 - p0)
    psubusb               m3, m7
    por                   m4, m3              ; abs(q1 - q0) > thresh ||
                                              ; abs(p1 - p0) > thresh
    pcmpeqb               m4, m0
    pcmpeqb               m0, m0
    movq                  m3, [GLOBAL(t80)]
    pxor                  m4, m0

    ; start work on filters
    movq                  m2, [stkreg + p1]
    movq                  m7, [stkreg + q1]
    movq                  m6, [stkreg + p0]
    movq                  m0, [stkreg + q0]
    pxor                  m2, m3
    pxor                  m7, m3
    psubsb                m2, m7              ; p1 - q1
    pand                  m2, m4              ; high var mask (hvm)(p1 - q1)
    pxor                  m6, m3
    pxor                  m0, m3
    movq                  m3, m0              ; q0
    psubsb                m0, m6              ; q0 - p0
    paddsb                m2, m0              ; 1 * (q0 - p0) + hvm(p1 - q1)
    paddsb                m2, m0              ; 2 * (q0 - p0) + hvm(p1 - q1)
    paddsb                m2, m0              ; 3 * (q0 - p0) + hvm(p1 - q1)
    pand                  m1, m2              ; mask filter values we don't
                                              ; care about
    movq                  m2, m1
    paddsb                m1, [GLOBAL(t4)]    ; 3*(q0 - p0) + hvm(p1 - q1) + 4
    paddsb                m2, [GLOBAL(t3)]    ; 3*(q0 - p0) + hvm(p1 - q1) + 3
    pxor                  m0, m0
    pxor                  m5, m5
    punpcklbw             m0, m2
    punpckhbw             m5, m2
    psraw                 m0, 11
    psraw                 m5, 11
    packsswb              m0, m5
    movq                  m2, m0              ; (3*(q0 - p0) + hvm(p1 - q1)
                                              ; + 3) >> 3;
    pxor                  m0, m0
    movq                  m5, m1              ; abcdefgh
    punpcklbw             m0, m1              ; e0f0g0h0
    psraw                 m0, 11              ; sign extended shift right by 3
    pxor                  m1, m1
    punpckhbw             m1, m5              ; a0b0c0d0
    psraw                 m1, 11              ; sign extended shift right by 3
    movq                  m5, m0              ; save results
    packsswb              m0, m1              ; (3*(q0 - p0) + hvm(p1 - q1)
                                              ; + 4) >>3
    paddsw                m5, [GLOBAL(ones)]
    paddsw                m1, [GLOBAL(ones)]
    psraw                 m5, 1
    psraw                 m1, 1
    packsswb              m5, m1              ; (3* (q0 - p0) + hvm(p1 - q1)
                                              ; + 4) >>4
    pandn                 m4, m5              ; high edge variance additive
    movq                  m5, [GLOBAL(t80)]
    paddsb                m6, m2              ; p0+= p0 add
    pxor                  m6, m5              ; unoffset
    ; m6=p0
    movq                  m1, [stkreg + p1]
    pxor                  m1, m5              ; reoffset
    paddsb                m1, m4              ; p1+= p1 add
    pxor                  m1, m5              ; unoffset
    ; m6 = p0 m1 = p1
    psubsb                m3, m0              ; q0-= q0 add
    pxor                  m3, m5              ; unoffset
    ; m3 = q0
    psubsb                m7, m4              ; q1-= q1 add
    pxor                  m7, m5              ; unoffset
    ; m7 = q1
    ; transpose and write back
    ; m1 =    72 62 52 42 32 22 12 02
    ; m6 =    73 63 53 43 33 23 13 03
    ; m3 =    74 64 54 44 34 24 14 04
    ; m7 =    75 65 55 45 35 25 15 05
    movq                  m2, m1              ; 72 62 52 42 32 22 12 02
    punpcklbw             m2, m6              ; 33 32 23 22 13 12 03 02
    movq                  m4, m3              ; 74 64 54 44 34 24 14 04
    punpckhbw             m1, m6              ; 73 72 63 62 53 52 43 42
    punpcklbw             m4, m7              ; 35 34 25 24 15 14 05 04
    punpckhbw             m3, m7              ; 75 74 65 64 55 54 45 44
    movq                  m6, m2              ; 33 32 23 22 13 12 03 02
    punpcklwd             m2, m4              ; 15 14 13 12 05 04 03 02
    punpckhwd             m6, m4              ; 35 34 33 32 25 24 23 22
    movq                  m5, m1              ; 73 72 63 62 53 52 43 42
    punpcklwd             m1, m3              ; 55 54 53 52 45 44 43 42
    punpckhwd             m5, m3              ; 75 74 73 72 65 64 63 62

    ; m2 = 15 14 13 12 05 04 03 02
    ; m6 = 35 34 33 32 25 24 23 22
    ; m5 = 55 54 53 52 45 44 43 42
    ; m1 = 75 74 73 72 65 64 63 62
    movd   [sq + pq * 4 + 2], m2
    psrlq                 m2, 32
    movd  [s1q + pq * 4 + 2], m2
    movd   [sq + pq * 2 + 2], m6
    psrlq                 m6, 32
    movd       [sq + pq + 2], m6
    movd            [sq + 2], m1
    psrlq                 m1, 32
    movd           [s1q + 2], m1
    neg                   pq
    movd      [s1q + pq + 2], m5
    psrlq                 m5, 32
    movd  [s1q + pq * 2 + 2], m5
    RET
