MEM_SRCS-yes += vpx_mem.mk
MEM_SRCS-yes += vpx_mem.c
MEM_SRCS-yes += vpx_mem.h
MEM_SRCS-yes += include/vpx_mem_intrnl.h

MEM_SRCS-$(CONFIG_MEM_TRACKER) += vpx_mem_tracker.c
MEM_SRCS-$(CONFIG_MEM_TRACKER) += include/vpx_mem_tracker.h
