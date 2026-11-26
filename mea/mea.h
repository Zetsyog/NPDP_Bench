#ifndef MEA_H
#define MEA_H

#ifndef Ebp
#define Ebp 0 // Energy weight of base pair  -2, -1, 0, 1, 2
#endif		  // Ebp

#ifndef RT
#define RT 1 // 'Normalized' temperature 1,2,3,4,5
#endif		 // RT

// Minimum loop length 0-5. Passed to the kernel as the `l` parameter.
#ifndef MIN_LOOP
#define MIN_LOOP 0
#endif // MIN_LOOP

#ifndef delta
#define delta 1 // Base pair weighting  1-5
#endif			// delta

#endif