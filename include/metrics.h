// Kernel timing.
//
//   NPDP_TIME              print the kernel execution time (seconds) on stdout
//   NPDP_FLUSH_CACHE       flush the caches right before starting the timer,
//                          so that every run starts cold
//   NPDP_FLUSH_CACHE_SIZE  bytes written/read to flush (default 64 MiB; use
//                          more than the last level cache)
//
// Usage in a benchmark:
//   NPDP_TIMER_START;
//   kernel_xxx(N, ...);
//   NPDP_TIMER_STOP;
//   NPDP_TIMER_PRINT;
#ifndef NPDP_METRICS_H
#define NPDP_METRICS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef NPDP_FLUSH_CACHE
#ifndef NPDP_FLUSH_CACHE_SIZE
#define NPDP_FLUSH_CACHE_SIZE (64UL * 1024 * 1024)
#endif
static inline void npdp_flush_cache(void) {
	size_t n = NPDP_FLUSH_CACHE_SIZE / sizeof(double);
	double *buf = calloc(n, sizeof(double));
	if (buf == NULL) {
		return;
	}
	for (size_t i = 0; i < n; i++) {
		buf[i] = (double)i;
	}
	volatile double sum = 0.0;
	for (size_t i = 0; i < n; i++) {
		sum += buf[i];
	}
	(void)sum;
	free(buf);
}
#define NPDP_FLUSH() npdp_flush_cache()
#else
#define NPDP_FLUSH() ((void)0)
#endif // NPDP_FLUSH_CACHE

#ifdef NPDP_TIME
static double npdp_t_start, npdp_t_end;

// Monotonic wall-clock time in seconds (not affected by clock adjustments).
static inline double npdp_clock(void) {
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

#define NPDP_TIMER_START                                                                                                                                       \
	do {                                                                                                                                                       \
		NPDP_FLUSH();                                                                                                                                          \
		npdp_t_start = npdp_clock();                                                                                                                           \
	} while (0)
#define NPDP_TIMER_STOP (npdp_t_end = npdp_clock())
#define NPDP_TIMER_PRINT printf("%0.6lf\n", npdp_t_end - npdp_t_start)
#else
#define NPDP_TIMER_START NPDP_FLUSH()
#define NPDP_TIMER_STOP ((void)0)
#define NPDP_TIMER_PRINT ((void)0)
#endif // NPDP_TIME

#endif // NPDP_METRICS_H
