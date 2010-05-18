/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */


#if !defined(_mac_specs_h)
#define _mac_specs_h


#if defined(__cplusplus)
extern "C" {
#endif

    extern unsigned int vp8_read_tsc();

    extern unsigned int vp8_get_processor_freq();

    extern unsigned int vpx_has_altivec();

#if defined(__cplusplus)
}
#endif


#endif
