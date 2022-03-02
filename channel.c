#include <stdio.h>
#include <stdlib.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "rnd.h"
#include "math.h"

float eb2n0 = 15;

float awgn_gen(float eb2n0_snr)
{
	float val = 0;

	#define E_B		1.0
	float w = 0, r = 0;

	w = (float)rand() / (float)RAND_MAX;
	if (w == 1.0)
	{
		w = 0.999999;
	}

#if 1
	r = gaussrand();
	val = r * (sqrt((E_B / ((float)MESSAGE_LEN / (float)CODEWORD_LEN)) / pow(10, eb2n0_snr / 10) / 2));
#else
	r = sqrt(2.0 * log(1.0 / (1.0 - w)));

	/*convert it to modulation constellation*/
	val = r * (float)cos(2 * PI * w);

	//val = val * (sqrt(E_B / pow(10, eb2n0 / 10) / 2));

	val = val * sqrt((E_B / ((float)MESSAGE_LEN / (float)CODEWORD_LEN)) / pow(10.0, (float)(eb2n0) / 10.0) / 2);
	
#endif
	DEBUG_NOTICE("val: %f %f %f\n", r, w, val);
	DEBUG_NOTICE("gauss_val: %f\n", val);

#if 0
#if (1 == OUTPUT_LOG)
	FILE *frc;
	frc = fopen("gauss.txt", "a+");
	fprintf(frc, "%f\n", val);
	fclose(frc);
	frc = NULL;
	FILE *frc1;
	frc1 = fopen("val.txt", "a+");
	fprintf(frc1, "%f\n", r);
	fclose(frc1);
	frc1 = NULL;
#endif	
#endif
	return val;
}
