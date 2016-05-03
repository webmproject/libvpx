;
;  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
;
;  Use of this source code is governed by a BSD-style license
;  that can be found in the LICENSE file in the root of the source
;  tree. An additional intellectual property rights grant can be found
;  in the file PATENTS.  All contributing project authors may
;  be found in the AUTHORS file in the root of the source tree.
;

%include "vpx_ports/x86_abi_support.asm"

;void vpx_plane_add_noise_mmx (unsigned char *Start, unsigned char *noise,
;                            unsigned char blackclamp[16],
;                            unsigned char whiteclamp[16],
;                            unsigned char bothclamp[16],
;                            unsigned int Width, unsigned int Height, int Pitch)
global sym(vpx_plane_add_noise_mmx) PRIVATE
sym(vpx_plane_add_noise_mmx):
    push        rbp
    mov         rbp, rsp
    SHADOW_ARGS_TO_STACK 8
    GET_GOT     rbx
    push        rsi
    push        rdi
    ; end prolog

.addnoise_loop:
    call sym(LIBVPX_RAND) WRT_PLT
    mov     rcx, arg(1) ;noise
    and     rax, 0xff
    add     rcx, rax

    ; we rely on the fact that the clamping vectors are stored contiguously
    ; in black/white/both order. Note that we have to reload this here because
    ; rdx could be trashed by rand()
    mov     rdx, arg(2) ; blackclamp


            mov     rdi, rcx
            movsxd  rcx, dword arg(5) ;[Width]
            mov     rsi, arg(0) ;Pos
            xor         rax,rax

.addnoise_nextset:
            movq        mm1,[rsi+rax]         ; get the source

            psubusb     mm1, [rdx]    ;blackclamp        ; clamp both sides so we don't outrange adding noise
            paddusb     mm1, [rdx+32] ;bothclamp
            psubusb     mm1, [rdx+16] ;whiteclamp

            movq        mm2,[rdi+rax]         ; get the noise for this line
            paddb       mm1,mm2              ; add it in
            movq        [rsi+rax],mm1         ; store the result

            add         rax,8                 ; move to the next line

            cmp         rax, rcx
            jl          .addnoise_nextset

    movsxd  rax, dword arg(7) ; Pitch
    add     arg(0), rax ; Start += Pitch
    sub     dword arg(6), 1   ; Height -= 1
    jg      .addnoise_loop

    ; begin epilog
    pop rdi
    pop rsi
    RESTORE_GOT
    UNSHADOW_ARGS
    pop         rbp
    ret


SECTION_RODATA
align 16
Blur:
    times 16 dw 16
    times  8 dw 64
    times 16 dw 16
    times  8 dw  0

rd:
    times 4 dw 0x40
