# Top-level targets. Each benchmark directory also has its own Makefile
# (see common.mk) to build <bench> and <bench>_pluto.

BENCHES = counting counting_im knuth mcc mcm mea nussinov nw sw sw3d triang zuker

.PHONY: all test golden clean $(BENCHES)

# Build the original (non-pluto) version of every benchmark.
all: $(BENCHES)

$(BENCHES):
	$(MAKE) -C $@ $@

# Run the test suite (warnings, sanitizers, golden checksums, known answers, valgrind).
test:
	bash tests/run_tests.sh

# Regenerate tests/golden/ (only when a change of results is intended).
golden:
	bash tests/update_golden.sh

clean:
	for b in $(BENCHES); do $(MAKE) -C $$b clean; done
