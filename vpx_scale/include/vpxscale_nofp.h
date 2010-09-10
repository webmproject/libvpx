/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#if defined(__S60_V20__) || defined(__SYMBIAN32__) && !defined(__WINS__)
#include "symbian\vpxscale_nofp.h"
#else
#include "generic\vpxscale_nofp.h"
#endif
