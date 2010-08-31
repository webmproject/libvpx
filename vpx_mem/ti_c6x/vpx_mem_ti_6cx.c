/*
 *  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#define __VPX_MEM_C__

#include "..\include\vpx_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "..\include\vpx_mem_intrnl.h"

void *vpx_mem_alloc(int id, size_t size, size_t align)
{
#if defined CHIP_DM642 || defined __uClinux__
    void *mem = (void *)mem_alloc(id, size, align);

    if (!mem)
    {
        _P(fprintf(stderr,
                   "\n"
                   "*********************************************************\n"
                   "WARNING: mem_alloc returned 0 for id=%p size=%u align=%u.\n"
                   "*********************************************************\n",
                   mem, size, align));
        // should no longer need this.  Softier says it's fixed. 2005-01-21 tjf
        //#if defined __uClinux__
        //while(1)usleep(1000000);
        //#endif
    }

#if defined __uClinux__
    else if (mem == (void *)0xFFFFFFFF)
    {
        // out of memory/error
        mem = (void *)0;

        _P(fprintf(stderr,
                   "\n"
                   "******************************************************\n"
                   "ERROR: mem_alloc id=%p size=%u align=%u OUT OF MEMORY.\n"
                   "******************************************************\n",
                   mem, size, align));
    }

#endif  // __uClinux__

    return mem;
#else
    (void)id;
    (void)size;
    (void)align;
    return (void *)0;
#endif
}

void vpx_mem_free(int id, void *mem, size_t size)
{
#if defined CHIP_DM642 || defined __uClinux__

    if (!mem)
    {
        _P(fprintf(stderr,
                   "\n"
                   "**************************************\n"
                   "WARNING: 0 being free'd id=%p size=%u.\n"
                   "**************************************\n",
                   id, size));

        // should no longer need this.  Softier says it's fixed. 2005-01-21 tjf
        //#if defined __uClinux__
        //while(1)usleep(1000000);
        //#endif
    }

    mem_free(id, mem, size);
#else
    (void)id;
    (void)mem;
    (void)size;
#endif
}

#if CONFIG_MEM_TRACKER
void *xvpx_mem_alloc(int id, size_t size, size_t align, char *file, int line)
{
    void *mem = vpx_mem_alloc(id, size, align);

    vpx_memory_tracker_add((size_t)mem, size, file, line, 0);

    return mem;
}

void xvpx_mem_free(int id, void *mem, size_t size, char *file, int line)
{
    if (vpx_memory_tracker_remove((size_t)mem) == -2)
    {
#if REMOVE_PRINTFS
        (void)file;
        (void)line;
#endif
        _P(fprintf(stderr, "[vpx_mem][xvpx_mem_free] addr: %p (id=%p size=%u) "
                   "not found in list; freed from file:%s"
                   " line:%d\n", mem, id, size, file, line));
    }

    vpx_mem_free(id, mem, size);
}
#endif /*CONFIG_MEM_TRACKER*/
