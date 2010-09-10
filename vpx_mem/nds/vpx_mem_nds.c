/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */


#define __VPX_MEM_C__
#include "vpx_mem.h"
#include <nitro.h>
#include "vpx_mem_intrnl.h"

// Allocate memory from the Arena specified by id.  Align it to
//  the value specified by align.
void *vpx_mem_nds_alloc(osarena_id id, osheap_handle handle, size_t size, size_t align)
{
    void *addr,
         * x = NULL;

    addr = os_alloc_from_heap((osarena_id) id, handle,
                              size + align - 1 + ADDRESS_STORAGE_SIZE);

    if (addr)
    {
        x = align_addr((unsigned char *)addr + ADDRESS_STORAGE_SIZE, (int)align);

        // save the actual malloc address
        ((size_t *)x)[-1] = (size_t)addr;
    }

    return x;
}

// Free them memory allocated by vpx_mem_nds_alloc
void vpx_mem_nds_free(osarena_id id, osheap_handle handle, void *mem)
{
    if (mem)
    {
        void *addr = (void *)(((size_t *)mem)[-1]);
        os_free_to_heap(id, handle, addr);
    }
}

int vpx_nds_alloc_heap(osarena_id id, u32 size)
{
    osheap_handle    arena_handle;
    void           *nstart;
    void           *heap_start;

    nstart = os_init_alloc(id, os_get_arena_lo(id), os_get_arena_hi(id), 1);
    os_set_arena_lo(id, nstart);

    heap_start = os_alloc_from_arena_lo(id, size, 32);
    arena_handle = os_create_heap(id, heap_start, (void *)((u32)heap_start + size));

    if (os_check_heap(id, arena_handle) == -1)
        return -1;      //ERROR: DTCM heap is not consistent

    (void)os_set_current_heap(id, arena_handle);

    return arena_handle;
}
