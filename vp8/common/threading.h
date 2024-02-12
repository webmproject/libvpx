/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_VP8_COMMON_THREADING_H_
#define VPX_VP8_COMMON_THREADING_H_

#include "./vpx_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_OS_SUPPORT && CONFIG_MULTITHREAD

/* Thread management macros */
#if defined(_WIN32) && !HAVE_PTHREAD_H
/* Win32 */
#include <process.h>
#include <windows.h>
#if defined(__GNUC__) && \
    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))
#define THREAD_FUNCTION \
  __attribute__((force_align_arg_pointer)) unsigned int __stdcall
#else
#define THREAD_FUNCTION unsigned int __stdcall
#endif
#define THREAD_FUNCTION_RETURN DWORD
#define THREAD_SPECIFIC_INDEX DWORD
#define pthread_t HANDLE
#define pthread_detach(thread) \
  if (thread != NULL) CloseHandle(thread)
#define thread_sleep(nms) Sleep(nms)
#define pthread_cancel(thread) terminate_thread(thread, 0)
#define ts_key_create(ts_key, destructor) \
  { ts_key = TlsAlloc(); };
#define pthread_getspecific(ts_key) TlsGetValue(ts_key)
#define pthread_setspecific(ts_key, value) TlsSetValue(ts_key, (void *)value)
#define pthread_self() GetCurrentThreadId()

#else
#ifdef __APPLE__
#include <mach/mach_init.h>
#include <mach/semaphore.h>
#include <mach/task.h>
#include <time.h>
#include <unistd.h>

#else
#include <semaphore.h>
#endif

#include <pthread.h>
/* pthreads */
/* Nearly everything is already defined */
#define THREAD_FUNCTION void *
#define THREAD_FUNCTION_RETURN void *
#define THREAD_SPECIFIC_INDEX pthread_key_t
#define ts_key_create(ts_key, destructor) \
  pthread_key_create(&(ts_key), destructor);
#endif

/* Synchronization macros: Win32 and Pthreads */
#if defined(_WIN32) && !HAVE_PTHREAD_H
#define sem_t HANDLE
#define pause(voidpara) __asm PAUSE
#define sem_init(sem, sem_attr1, sem_init_value) \
  (int)((*sem = CreateSemaphore(NULL, 0, 32768, NULL)) == NULL)
#define sem_wait(sem) \
  (int)(WAIT_OBJECT_0 != WaitForSingleObject(*sem, INFINITE))
#define sem_post(sem) ReleaseSemaphore(*sem, 1, NULL)
#define sem_destroy(sem) \
  if (*sem) ((int)(CloseHandle(*sem)) == TRUE)
#define thread_sleep(nms) Sleep(nms)

#else

#ifdef __APPLE__
#define sem_t semaphore_t
#define sem_init(X, Y, Z) \
  semaphore_create(mach_task_self(), X, SYNC_POLICY_FIFO, Z)
#define sem_wait(sem) (semaphore_wait(*sem))
#define sem_post(sem) semaphore_signal(*sem)
#define sem_destroy(sem) semaphore_destroy(mach_task_self(), *sem)
#else
#include <unistd.h>
#include <sched.h>
#endif /* __APPLE__ */
/* Not Windows. Assume pthreads */

/* thread_sleep implementation: yield unless Linux/Unix. */
#if defined(__unix__) || defined(__APPLE__)
#define thread_sleep(nms)
/* {struct timespec ts;ts.tv_sec=0;
    ts.tv_nsec = 1000*nms;nanosleep(&ts, NULL);} */
#else
#define thread_sleep(nms) sched_yield();
#endif /* __unix__ || __APPLE__ */

#endif

#if VPX_ARCH_X86 || VPX_ARCH_X86_64
#include "vpx_ports/x86.h"
#else
#define x86_pause_hint()
#endif

#include "vpx_util/vpx_atomics.h"

static INLINE void vp8_atomic_spin_wait(
    int mb_col, const vpx_atomic_int *last_row_current_mb_col,
    const int nsync) {
  while (mb_col > (vpx_atomic_load_acquire(last_row_current_mb_col) - nsync)) {
    x86_pause_hint();
    thread_sleep(0);
  }
}

#endif /* CONFIG_OS_SUPPORT && CONFIG_MULTITHREAD */

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_VP8_COMMON_THREADING_H_
