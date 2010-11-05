/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include <stdio.h>
#include "tools_common.h"
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

FILE* set_binary_mode(FILE *stream)
{
    (void)stream;
#ifdef _WIN32
    _setmode(_fileno(stream), _O_BINARY);
#endif
    return stream;
}
