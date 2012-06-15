/*
 *  Copyright (c) 2011 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */
#include "vpx_config.h"
#define RTCD_C
#include "vpx_rtcd.h"

#if CONFIG_MULTITHREAD && HAVE_PTHREAD_H
#include <pthread.h>
static void once(void (*func)(void))
{
    static pthread_once_t lock = PTHREAD_ONCE_INIT;
    pthread_once(&lock, func);
}


#elif CONFIG_MULTITHREAD && defined(_WIN32)
#include <windows.h>
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


void vpx_rtcd()
{
    once(setup_rtcd_internal);
}
