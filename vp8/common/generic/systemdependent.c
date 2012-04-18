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
#if ARCH_ARM
#include "vpx_ports/arm.h"
#elif ARCH_X86 || ARCH_X86_64
#include "vpx_ports/x86.h"
#endif
#include "vp8/common/onyxc_int.h"

#if CONFIG_MULTITHREAD
#if HAVE_UNISTD_H && !defined(__OS2__)
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);
#elif defined(__OS2__)
#define INCL_DOS
#define INCL_DOSSPINLOCK
#include <os2.h>
#endif
#endif

#if CONFIG_MULTITHREAD
static int get_cpu_count()
{
    int core_count = 16;

#if HAVE_UNISTD_H && !defined(__OS2__)
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
#elif defined(__OS2__)
    {
        ULONG proc_id;
        ULONG status;

        core_count = 0;
        for (proc_id = 1; ; proc_id++)
        {
            if (DosGetProcessorStatus(proc_id, &status))
                break;

            if (status == PROC_ONLINE)
                core_count++;
        }
    }
#else
    /* other platforms */
#endif

    return core_count > 0 ? core_count : 1;
}
#endif


#if HAVE_PTHREAD_H
#include <pthread.h>
static void once(void (*func)(void))
{
    static pthread_once_t lock = PTHREAD_ONCE_INIT;
    pthread_once(&lock, func);
}


#elif defined(_WIN32)
static void once(void (*func)(void))
{
    /* Using a static initializer here rather than InitializeCriticalSection()
     * since there's no race-free context in which to execute it. Protecting
     * it with an atomic op like InterlockedCompareExchangePointer introduces
     * an x86 dependency, and InitOnceExecuteOnce requires Vista.
     */
    static CRITICAL_SECTION lock = {(void *)-1, -1, 0, 0, 0, 0};
    static int done;

    EnterCriticalSection(&lock);

    if (!done)
    {
        func();
        done = 1;
    }

    LeaveCriticalSection(&lock);
}


#else
/* No-op version that performs no synchronization. vpx_rtcd() is idempotent,
 * so as long as your platform provides atomic loads/stores of pointers
 * no synchronization is strictly necessary.
 */

static void once(void (*func)(void))
{
    static int done;

    if(!done)
    {
        func();
        done = 1;
    }
}
#endif


void vp8_machine_specific_config(VP8_COMMON *ctx)
{
#if CONFIG_MULTITHREAD
    ctx->processor_core_count = get_cpu_count();
#endif /* CONFIG_MULTITHREAD */

#if ARCH_ARM
    ctx->cpu_caps = arm_cpu_caps();
#elif ARCH_X86 || ARCH_X86_64
    ctx->cpu_caps = x86_simd_caps();
#endif

    once(vpx_rtcd);
}
