/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#include "vpx_config.h"
#include "vpx_rtcd.h"
#include "vp8/common/subpixel.h"
#include "vp8/common/loopfilter.h"
#include "vp8/common/onyxc_int.h"

#if CONFIG_MULTITHREAD
#if HAVE_UNISTD_H
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);
#endif
#endif

extern void vp8_arch_x86_common_init(VP8_COMMON *ctx);
extern void vp8_arch_arm_common_init(VP8_COMMON *ctx);

#if CONFIG_MULTITHREAD
static int get_cpu_count()
{
    int core_count = 16;

#if HAVE_UNISTD_H
#if defined(_SC_NPROCESSORS_ONLN)
    core_count = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    core_count = sysconf(_SC_NPROC_ONLN);
#endif
#elif defined(_WIN32)
    {
        PGNSI pGNSI;
        SYSTEM_INFO sysinfo;

        /* Call GetNativeSystemInfo if supported or
         * GetSystemInfo otherwise. */

        pGNSI = (PGNSI) GetProcAddress(
                GetModuleHandle(TEXT("kernel32.dll")), "GetNativeSystemInfo");
        if (pGNSI != NULL)
            pGNSI(&sysinfo);
        else
            GetSystemInfo(&sysinfo);

        core_count = sysinfo.dwNumberOfProcessors;
    }
#else
    /* other platforms */
#endif

    return core_count > 0 ? core_count : 1;
}
#endif

void vp8_machine_specific_config(VP8_COMMON *ctx)
{
#if CONFIG_RUNTIME_CPU_DETECT
    VP8_COMMON_RTCD *rtcd = &ctx->rtcd;

    rtcd->subpix.sixtap16x16   = vp8_sixtap_predict16x16_c;
    rtcd->subpix.sixtap8x8     = vp8_sixtap_predict8x8_c;
    rtcd->subpix.sixtap8x4     = vp8_sixtap_predict8x4_c;
    rtcd->subpix.sixtap4x4     = vp8_sixtap_predict_c;
    rtcd->subpix.bilinear16x16 = vp8_bilinear_predict16x16_c;
    rtcd->subpix.bilinear8x8   = vp8_bilinear_predict8x8_c;
    rtcd->subpix.bilinear8x4   = vp8_bilinear_predict8x4_c;
    rtcd->subpix.bilinear4x4   = vp8_bilinear_predict4x4_c;
#endif

#if ARCH_X86 || ARCH_X86_64
    vp8_arch_x86_common_init(ctx);
#endif

#if ARCH_ARM
    vp8_arch_arm_common_init(ctx);
#endif

#if CONFIG_MULTITHREAD
    ctx->processor_core_count = get_cpu_count();
#endif /* CONFIG_MULTITHREAD */

    vpx_rtcd();
}
