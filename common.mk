CC=gcc
POLYCC=${PLUTO}
# -ffp-contract=off: without it, gcc fuses a*b+c into FMA instructions wherever
# it can, which changes floating-point results (mcc, mea) depending on how the
# loops are written, so original and pluto checksums would not match.
CFLAGS=-march=native -O3 -fopenmp -ffp-contract=off
CPPFLAGS=-I ../include/
LDLIBS=-lm
PLUTO_FLAGS=--tile --parallel
HEADERS=$(wildcard ../include/*.h) $(wildcard *.h)


all: $(PROG) $(PROG)_pluto

$(PROG): $(PROG).c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(PROG) $(PROG).c $(LDLIBS) $(EXTRA_FLAGS)

$(PROG)_pluto: $(PROG).pluto.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $(PROG)_pluto $(PROG).pluto.c $(LDLIBS) $(EXTRA_FLAGS)

$(PROG).pluto.c: $(PROG).c
ifeq ($(strip $(POLYCC)),)
	$(error PLUTO is not set: export PLUTO=/path/to/polycc to build $(PROG)_pluto)
endif
	$(POLYCC) -o $@ $< $(PLUTO_FLAGS) $(PLUTO_EXTRA_FLAGS) $(PLUTO_QUIET)

clean:
	rm -f $(PROG) $(PROG)_pluto  *.pluto.c *.cloog

.PHONY: all clean
