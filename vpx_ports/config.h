/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


/* This file uses some preprocessor magic to expand the value of HAVE_CONFIG_H,
 * as defined by the build system, so that different projects can use the file
 * name for config.h that suits them.
 */
#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)
#include QUOTE(HAVE_CONFIG_H)
#undef QUOTE
#undef QUOTE_
