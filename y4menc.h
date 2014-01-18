/*
 *  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef Y4MENC_H_
#define Y4MENC_H_

#include <stdio.h>

#include "./tools_common.h"

#include "vpx/vpx_decoder.h"

void y4m_write_file_header(FILE *file, int width, int height,
                           const struct VpxRational *framerate,
                           vpx_img_fmt_t fmt);

void y4m_write_frame_header(FILE *file);


#endif  // Y4MENC_H_
