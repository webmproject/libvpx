all : filtering horiz_filter

CFLAGS += -mssse3

filtering : filtering.c
horiz_filter : horiz_filter.c

.PHONY : clean
clean :
	$(RM) filtering horiz_filter

