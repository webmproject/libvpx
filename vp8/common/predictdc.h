/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#ifndef __PREDICTDC_H
#define __PREDICTDC_H

void uvvp8_predict_dc(short *lastdc, short *thisdc, short quant, short *cons);
void vp8_predict_dc(short *lastdc, short *thisdc, short quant, short *cons);

#endif
