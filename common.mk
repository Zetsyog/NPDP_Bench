CC=gcc
POLYCC=${POLYCC}

all: $(BENCHMARKS)

clean:
	rm -f $(BENCHMARKS) *.cloog