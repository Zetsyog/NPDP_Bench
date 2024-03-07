BENCHMARKS=counting counting_im

include ../common.mk

CFLAGS=-march=native -O3 -fopenmp
LDLIBS=-lm -DPOLYBENCH_TIME
EXTRA_FLAGS=

counting: counting.c mem.h
	$(CC) $(CFLAGS) -o counting counting.c $(LDLIBS) -I ../include/ $(EXTRA_FLAGS)

counting_im: counting_im.c mem.h
	$(CC) $(CFLAGS) -o counting_im counting_im.c $(LDLIBS) -I ../include/ $(EXTRA_FLAGS)

