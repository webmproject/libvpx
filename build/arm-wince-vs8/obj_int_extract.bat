@echo off
REM   Copyright (c) 2010 The WebM project authors. All Rights Reserved.
REM
REM   Use of this source code is governed by a BSD-style license
REM   that can be found in the LICENSE file in the root of the source
REM   tree. An additional intellectual property rights grant can be found
REM   in the file PATENTS.  All contributing project authors may
REM   be found in the AUTHORS file in the root of the source tree.
echo on


cl /I ".\\" /I "..\vp6_decoder_sdk" /I "..\vp6_decoder_sdk\vpx_ports" /D "NDEBUG" /D "_WIN32_WCE=0x420" /D "UNDER_CE" /D "WIN32_PLATFORM_PSPC" /D "WINCE" /D "_LIB" /D "ARM" /D "_ARM_" /D "_UNICODE" /D "UNICODE" /FD /EHsc /MT /GS- /fp:fast /GR- /Fo"Pocket_PC_2003__ARMV4_\%1/" /Fd"Pocket_PC_2003__ARMV4_\%1/vc80.pdb" /W3 /nologo /c /TC ..\vp6_decoder_sdk\vp6_decoder\algo\common\arm\dec_asm_offsets_arm.c
obj_int_extract.exe rvds "Pocket_PC_2003__ARMV4_\%1/dec_asm_offsets_arm.obj"
