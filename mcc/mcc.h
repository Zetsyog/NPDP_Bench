#ifndef NPDP_MCC_H
#define NPDP_MCC_H

// Energy weight of a base pair: -2, -1, 0, 1, 2
#ifndef NPDP_EBP
#define NPDP_EBP 0
#endif

// 'Normalized' temperature: 1, 2, 3, 4, 5
#ifndef NPDP_RT
#define NPDP_RT 1
#endif

// Minimum loop length 0-5. Passed to the kernel as the `l` parameter.
#ifndef NPDP_MIN_LOOP
#define NPDP_MIN_LOOP 0
#endif

#endif // NPDP_MCC_H
