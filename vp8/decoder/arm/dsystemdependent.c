/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_ports/config.h"
#include "blockd.h"
#include "pragmas.h"
#include "postproc.h"
#include "dboolhuff.h"
#include "dequantize.h"
#include "onyxd_int.h"

void vp8_dmachine_specific_config(VP8D_COMP *pbi)
{
#if CONFIG_RUNTIME_CPU_DETECT
    pbi->mb.rtcd         = &pbi->common.rtcd;
#if HAVE_ARMV7
    pbi->dequant.block   = vp8_dequantize_b_neon;
    pbi->dequant.idct    = vp8_dequant_idct_neon;
    pbi->dequant.idct_dc = vp8_dequant_dc_idct_neon;
    pbi->dboolhuff.start = vp8dx_start_decode_c;
    pbi->dboolhuff.fill  = vp8dx_bool_decoder_fill_c;
    pbi->dboolhuff.debool = vp8dx_decode_bool_c;
    pbi->dboolhuff.devalue = vp8dx_decode_value_c;

#elif HAVE_ARMV6
    pbi->dequant.block   = vp8_dequantize_b_v6;
    pbi->dequant.idct    = vp8_dequant_idct_v6;
    pbi->dequant.idct_dc = vp8_dequant_dc_idct_v6;
    pbi->dboolhuff.start = vp8dx_start_decode_c;
    pbi->dboolhuff.fill  = vp8dx_bool_decoder_fill_c;
    pbi->dboolhuff.debool = vp8dx_decode_bool_c;
    pbi->dboolhuff.devalue = vp8dx_decode_value_c;
#endif
#endif
}
