#ifndef RND_H
#define RND_H

#include "cfg_decoding.h"

#define PI			3.1415926536
#define N 			624
#define M 			397
#define MATRIX_A 	0x9908b0dfUL	/* constant vector a */
#define UPPER_MASK 0x80000000UL 	/* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL		/* least significant r bits */

extern void init_genrand(unsigned long s);
extern unsigned long genrand_int32(void);
extern float gaussrand();

#endif
