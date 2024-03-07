#ifndef METRICS_H
#define METRICS_H

#include <stdio.h>
#include <sys/time.h>

#define metrics_start_timer
#define metrics_stop_timer
#define metrics_print_timer

#ifdef METRICS_TIME
#undef metrics_start_timer
#undef metrics_stop_timer
#undef metrics_print_timer
#define metrics_start_timer _metrics_timer_start();
#define metrics_stop_timer _metrics_timer_stop();
#define metrics_print_timer _metrics_timer_print();
#endif // METRICS_TIME

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef RSEED
#define RSEED time(NULL)
#endif

double metric_t_start, metric_t_end;

double rtclock() {
	struct timezone Tzp;
	struct timeval Tp;
	int stat;
	stat = gettimeofday(&Tp, &Tzp);
	if (stat != 0)
		fprintf(stderr, "Error return from gettimeofday: %d", stat);

	return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void _metrics_timer_start() { metric_t_start = rtclock(); }

void _metrics_timer_stop() { metric_t_end = rtclock(); }

void _metrics_timer_print() { printf("%0.6lf\n", metric_t_end - metric_t_start); }

#endif // METRICS_H