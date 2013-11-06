/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#ifndef TOOLS_COMMON_H_
#define TOOLS_COMMON_H_

#include <stdio.h>

#define VP8_FOURCC (0x30385056)
#define VP9_FOURCC (0x30395056)
#define VP8_FOURCC_MASK (0x00385056)
#define VP9_FOURCC_MASK (0x00395056)

/* Sets a stdio stream into binary mode */
FILE *set_binary_mode(FILE *stream);

void die(const char *fmt, ...);
void fatal(const char *fmt, ...);
void warn(const char *fmt, ...);

/* The tool including this file must define usage_exit() */
void usage_exit();

#endif  // TOOLS_COMMON_H_
