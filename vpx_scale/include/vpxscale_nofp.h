/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#if defined(__S60_V20__) || defined(__SYMBIAN32__) && !defined(__WINS__)
#include "symbian\vpxscale_nofp.h"
#else
#include "generic\vpxscale_nofp.h"
#endif
