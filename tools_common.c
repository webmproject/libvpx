/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "tools_common.h"

#include <stdarg.h>
#include <stdlib.h>

#if defined(_WIN32) || defined(__OS2__)
#include <io.h>
#include <fcntl.h>

#ifdef __OS2__
#define _setmode    setmode
#define _fileno     fileno
#define _O_BINARY   O_BINARY
#endif
#endif

#define LOG_ERROR(label) do {\
  const char *l = label;\
  va_list ap;\
  va_start(ap, fmt);\
  if (l)\
    fprintf(stderr, "%s: ", l);\
  vfprintf(stderr, fmt, ap);\
  fprintf(stderr, "\n");\
  va_end(ap);\
} while (0)


FILE *set_binary_mode(FILE *stream) {
  (void)stream;
#if defined(_WIN32) || defined(__OS2__)
  _setmode(_fileno(stream), _O_BINARY);
#endif
  return stream;
}

void die(const char *fmt, ...) {
  LOG_ERROR(NULL);
  usage_exit();
}

void fatal(const char *fmt, ...) {
  LOG_ERROR("Fatal");
  exit(EXIT_FAILURE);
}

void warn(const char *fmt, ...) {
  LOG_ERROR("Warning");
}
