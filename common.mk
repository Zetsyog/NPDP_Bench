CC=gcc
POLYCC=${PLUTO}
CFLAGS=-march=native -O3 -fopenmp
LDLIBS=-lm -I ../include/
PLUTO_FLAGS=--tile --parallel


all: $(PROG) $(PROG)_pluto

$(PROG): $(PROG).c
	$(CC) $(CFLAGS) -o $(PROG) $(PROG).c $(LDLIBS) $(EXTRA_FLAGS)

$(PROG)_pluto: $(PROG).pluto.c
	$(CC) $(CFLAGS) -o $(PROG)_pluto $(PROG).pluto.c $(LDLIBS) $(EXTRA_FLAGS)

$(PROG).pluto.c: $(PROG).c
	$(POLYCC) -o $@ $< $(PLUTO_FLAGS) $(PLUTO_EXTRA_FLAGS)

clean:
	rm -f $(PROG) $(PROG)_pluto  *.pluto.c *.cloog