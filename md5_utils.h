/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license and patent
 *  grant that can be found in the LICENSE file in the root of the source
 *  tree. All contributing project authors may be found in the AUTHORS
 *  file in the root of the source tree.
 */

/*
Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
rights reserved.

License to copy and use this software is granted provided that it
is identified as the "RSA Data Security, Inc. MD5 Message-Digest
Algorithm" in all material mentioning or referencing this software
or this function.

License is also granted to make and use derivative works provided
that such works are identified as "derived from the RSA Data
Security, Inc. MD5 Message-Digest Algorithm" in all material
mentioning or referencing the derived work.

RSA Data Security, Inc. makes no representations concerning either
the merchantability of this software or the suitability of this
software for any particular purpose. It is provided "as is"
without express or implied warranty of any kind.

These notices must be retained in any copies of any part of this
documentation and/or software.
*/
#ifndef HAVE_VPX_PORTS
#define HAVE_VPX_PORTS 0
#endif
#if HAVE_VPX_PORTS
#include "vpx_ports/vpx_integer.h"
#else
#include "vpx_integer.h"
#endif

/* MD5 context. */
typedef struct
{
    uint32_t state[4];        /* state (ABCD) */
    uint32_t count[2];        /* number of bits, modulo 2^64 (lsb first) */
    uint8_t  buffer[64];      /* input buffer */
} md5_ctx_t;

void md5_init(md5_ctx_t *ctx);
void md5_update(md5_ctx_t *ctx, const uint8_t *buf, unsigned int len);
void md5_finalize(md5_ctx_t *ctx, uint8_t md5[16]);
