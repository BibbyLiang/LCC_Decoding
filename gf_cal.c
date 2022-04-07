#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"

#if (3 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x3},
	{0x4, 0x6},
	{0x5, 0x7},
	{0x6, 0x5}
};
#endif

#if (4 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x8},
	{0x4, 0x3},
	{0x5, 0x6},
	{0x6, 0xc},
	{0x7, 0xb},
	{0x8, 0x5},
	{0x9, 0xa},
	{0xa, 0x7},
	{0xb, 0xe},
	{0xc, 0xf},
	{0xd, 0xd},
	{0xe, 0x9}
};
#endif

#if (6 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0},
	{0, 1},
	{1, 2},
	{2, 4},
	{3, 8},
	{4, 16},
	{5, 32},
	{6, 3},
	{7, 6},
	{8, 12},
	{9, 24},
	{10, 48},
	{11, 35},
	{12, 5},
	{13, 10},
	{14, 20},
	{15, 40},
	{16, 19},
	{17, 38},
	{18, 15},
	{19, 30},
	{20, 60},
	{21, 59},
	{22, 53},
	{23, 41},
	{24, 17},
	{25, 34},
	{26, 7},
	{27, 14},
	{28, 28},
	{29, 56},
	{30, 51},
	{31, 37},
	{32, 9},
	{33, 18},
	{34, 36},
	{35, 11},
	{36, 22},
	{37, 44},
	{38, 27},
	{39, 54},
	{40, 47},
	{41, 29},
	{42, 58},
	{43, 55},
	{44, 45},
	{45, 25},
	{46, 50},
	{47, 39},
	{48, 13},
	{49, 26},
	{50, 52},
	{51, 43},
	{52, 21},
	{53, 42},
	{54, 23},
	{55, 46},
	{56, 31},
	{57, 62},
	{58, 63},
	{59, 61},
	{60, 57},
	{61, 49},
	{62, 33}
};
#endif

#if (8 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0},
	{0, 1},
	{1, 2},
	{2, 4},
	{3, 8},
	{4, 16},
	{5, 32},
	{6, 64},
	{7, 128},
	{8, 29},
	{9, 58},
	{10, 116},
	{11, 232},
	{12, 205},
	{13, 135},
	{14, 19},
	{15, 38},
	{16, 76},
	{17, 152},
	{18, 45},
	{19, 90},
	{20, 180},
	{21, 117},
	{22, 234},
	{23, 201},
	{24, 143},
	{25, 3},
	{26, 6},
	{27, 12},
	{28, 24},
	{29, 48},
	{30, 96},
	{31, 192},
	{32, 157},
	{33, 39},
	{34, 78},
	{35, 156},
	{36, 37},
	{37, 74},
	{38, 148},
	{39, 53},
	{40, 106},
	{41, 212},
	{42, 181},
	{43, 119},
	{44, 238},
	{45, 193},
	{46, 159},
	{47, 35},
	{48, 70},
	{49, 140},
	{50, 5},
	{51, 10},
	{52, 20},
	{53, 40},
	{54, 80},
	{55, 160},
	{56, 93},
	{57, 186},
	{58, 105},
	{59, 210},
	{60, 185},
	{61, 111},
	{62, 222},
	{63, 161},
	{64, 95},
	{65, 190},
	{66, 97},
	{67, 194},
	{68, 153},
	{69, 47},
	{70, 94},
	{71, 188},
	{72, 101},
	{73, 202},
	{74, 137},
	{75, 15},
	{76, 30},
	{77, 60},
	{78, 120},
	{79, 240},
	{80, 253},
	{81, 231},
	{82, 211},
	{83, 187},
	{84, 107},
	{85, 214},
	{86, 177},
	{87, 127},
	{88, 254},
	{89, 225},
	{90, 223},
	{91, 163},
	{92, 91},
	{93, 182},
	{94, 113},
	{95, 226},
	{96, 217},
	{97, 175},
	{98, 67},
	{99, 134},
	{100, 17},
	{101, 34},
	{102, 68},
	{103, 136},
	{104, 13},
	{105, 26},
	{106, 52},
	{107, 104},
	{108, 208},
	{109, 189},
	{110, 103},
	{111, 206},
	{112, 129},
	{113, 31},
	{114, 62},
	{115, 124},
	{116, 248},
	{117, 237},
	{118, 199},
	{119, 147},
	{120, 59},
	{121, 118},
	{122, 236},
	{123, 197},
	{124, 151},
	{125, 51},
	{126, 102},
	{127, 204},
	{128, 133},
	{129, 23},
	{130, 46},
	{131, 92},
	{132, 184},
	{133, 109},
	{134, 218},
	{135, 169},
	{136, 79},
	{137, 158},
	{138, 33},
	{139, 66},
	{140, 132},
	{141, 21},
	{142, 42},
	{143, 84},
	{144, 168},
	{145, 77},
	{146, 154},
	{147, 41},
	{148, 82},
	{149, 164},
	{150, 85},
	{151, 170},
	{152, 73},
	{153, 146},
	{154, 57},
	{155, 114},
	{156, 228},
	{157, 213},
	{158, 183},
	{159, 115},
	{160, 230},
	{161, 209},
	{162, 191},
	{163, 99},
	{164, 198},
	{165, 145},
	{166, 63},
	{167, 126},
	{168, 252},
	{169, 229},
	{170, 215},
	{171, 179},
	{172, 123},
	{173, 246},
	{174, 241},
	{175, 255},
	{176, 227},
	{177, 219},
	{178, 171},
	{179, 75},
	{180, 150},
	{181, 49},
	{182, 98},
	{183, 196},
	{184, 149},
	{185, 55},
	{186, 110},
	{187, 220},
	{188, 165},
	{189, 87},
	{190, 174},
	{191, 65},
	{192, 130},
	{193, 25},
	{194, 50},
	{195, 100},
	{196, 200},
	{197, 141},
	{198, 7},
	{199, 14},
	{200, 28},
	{201, 56},
	{202, 112},
	{203, 224},
	{204, 221},
	{205, 167},
	{206, 83},
	{207, 166},
	{208, 81},
	{209, 162},
	{210, 89},
	{211, 178},
	{212, 121},
	{213, 242},
	{214, 249},
	{215, 239},
	{216, 195},
	{217, 155},
	{218, 43},
	{219, 86},
	{220, 172},
	{221, 69},
	{222, 138},
	{223, 9},
	{224, 18},
	{225, 36},
	{226, 72},
	{227, 144},
	{228, 61},
	{229, 122},
	{230, 244},
	{231, 245},
	{232, 247},
	{233, 243},
	{234, 251},
	{235, 235},
	{236, 203},
	{237, 139},
	{238, 11},
	{239, 22},
	{240, 44},
	{241, 88},
	{242, 176},
	{243, 125},
	{244, 250},
	{245, 233},
	{246, 207},
	{247, 131},
	{248, 27},
	{249, 54},
	{250, 108},
	{251, 216},
	{252, 173},
	{253, 71},
	{254, 142}
};
#endif

#if (1 == GF_CAL_COUNT)
unsigned char cnt_switch = 0;

long long add_cnt = 0;
long long mul_cnt = 0;
long long div_cnt = 0;
long long real_cbm_cnt = 0;
long long real_mul_ff_cnt = 0;
long long pow_cnt = 0;
long long add_cnt_prev = 0;
long long mul_cnt_prev = 0;
long long div_cnt_prev = 0;
long long real_cbm_cnt_prev = 0;
long long real_mul_ff_cnt_prev = 0;
long long pow_cnt_prev = 0;
/*redundant in fact*/
long long err_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long add_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long mul_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long div_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long real_cbm_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long real_mul_ff_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
long long pow_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
#endif

#if (1 == CFG_PARTIALLY_PARALLEL)
long long storage_cnt = 0;
long long *gf_add_in_batch;
long long *gf_mul_in_batch;
long long *running_time;

long long cost_parallel = 0;
float speed_up = 0;
float efficiency = 0;

long long add_cnt_prev_batch = 0;
long long mul_cnt_prev_batch = 0;
long long div_cnt_prev_batch = 0;

long long *gf_add_in_batch_best;
long long *gf_mul_in_batch_best;
long long *gf_add_in_batch_worst;
long long *gf_mul_in_batch_worst;
#endif

unsigned char gf_pow2poly(unsigned char val_in_pow)
{
	unsigned char val_in_poly = 0;
	if(0xFF == val_in_pow)
	{
		return power_polynomial_table[0][1];
	}
	else
	{
		return power_polynomial_table[val_in_pow + 1][1];
	}
}

unsigned char gf_poly2pow(unsigned char val_in_poly)
{
	unsigned char i = 0;
	unsigned char val_in_pow = 0;

	for(i = 0; i < GF_FIELD; i++)
	{
		if(power_polynomial_table[i][1] == val_in_poly)
		{
			val_in_pow = power_polynomial_table[i][0];
			break;
		}
	}

	return val_in_pow;
}

unsigned char gf_location(unsigned char val)
{
	unsigned char val_location = 0xFF;

	if(0xFF == val)
	{
		val_location = power_polynomial_table[0][0];
	}
	else
	{
		val_location = power_polynomial_table[val + 1][0];
	}

	return val_location;
}

unsigned char gf_add(unsigned char a, unsigned char b)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		add_cnt++;
	}
#endif

	if(0xFF == a)
	{
		return b;
	}
	if(0xFF == b)
	{
		return a;
	}

	unsigned char i = 0;
	unsigned char sum_in_pow = 0;
	
	unsigned char sum_in_poly = gf_pow2poly(a) ^ gf_pow2poly(b);
#if 0	
	if(26 == b)
	{
		DEBUG_NOTICE("sum_in_poly: %x = %x + %x | %x ^ %x\n", sum_in_poly, a, b, gf_pow2poly(a), gf_pow2poly(b));
	}
#endif	
	sum_in_pow = gf_poly2pow(sum_in_poly);

	return sum_in_pow;
}

unsigned char gf_multp(unsigned char a, unsigned char b)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		mul_cnt++;
	}
#endif

	if((0xFF == a) || (0xFF == b))
	{
		return 0xFF;
	}

	if(0x0 == a)
	{
		return b;
	}
	
	if(0x0 == b)
	{
		return a;
	}

	unsigned char product_in_pow = (a + b) % (GF_FIELD - 1);

	return product_in_pow;
}

unsigned char gf_div(unsigned char a, unsigned char b)
{	
	if(0xFF == a)
	{
		return 0xFF;
	}
	if(0xFF == b)
	{
		DEBUG_NOTICE("div err.\n");
		return 0xFF;
	}
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		div_cnt++;
	}
#endif
	//DEBUG_NOTICE("div: %x %x\n", a, b);
	unsigned char quotient_in_pow = 0;
	if(a >= b)
	{
		quotient_in_pow = (a - b) % (GF_FIELD - 1);
	}
	else
	{
		quotient_in_pow = ((b / (GF_FIELD - 1) + 1) * (GF_FIELD - 1) + a - b) % (GF_FIELD - 1);
	}

	return quotient_in_pow;
}

unsigned char gf_mod_single_term(unsigned char a, unsigned char b)
{
	unsigned char remainder_in_pow = 0xFF;
	if(a >= b)
	{
		/*remove small term*/
		remainder_in_pow = 0xFF;
	}
	else
	{
		/*keep big term*/
		remainder_in_pow = 0;
	}

	return remainder_in_pow;
}

unsigned char gf_degree(unsigned char* a, unsigned char len_a)
{
	unsigned char i = 0;

	for(i = len_a - 1; i >= 0; i--)
	{
		if(0xFF != a[i])
		{
			break;
		}
	}

	return i;
}

unsigned char gf_div_q_r(unsigned char* dividend, unsigned char len_dividend,
							unsigned char* divisor, unsigned char len_divisor,
							unsigned char* quotien, unsigned char len_quotien,
							unsigned char* remainder, unsigned char len_remainder)
{
	unsigned char i = 0, j = 0, k = 0;
	unsigned char locator = 0, factor = 0, locator_rmd = 0, factor_rmd = 0;
	unsigned char dividend_tmp[len_dividend], remainder_tmp[len_remainder];
	memset(dividend_tmp, 0xFF, sizeof(unsigned char) * len_dividend);
	memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);

	if(gf_degree(divisor, len_divisor) > gf_degree(dividend, len_dividend))
	{
		for(i = 0; i < len_remainder; i++)
		{
			remainder[i] = dividend[i];
		}
		//DEBUG_NOTICE("quotien is zero: %d %d\n", gf_degree(divisor, len_divisor), gf_degree(dividend, len_dividend));

		return 0;
	}

	memcpy(dividend_tmp, dividend, sizeof(unsigned char) * len_dividend);
	for(i = 0; i < len_dividend; i++)
	{
		locator = gf_degree(dividend_tmp, len_dividend) - gf_degree(divisor, len_divisor);
		factor = gf_div(dividend_tmp[gf_degree(dividend_tmp, len_dividend)], divisor[gf_degree(divisor, len_divisor)]);

		quotien[locator] = factor;
		//DEBUG_NOTICE("quotien: %x %d %x\n", quotien[locator], locator, factor);

		for(j = 0; j < len_divisor; j++)
		{
			factor_rmd = gf_multp(factor, divisor[gf_degree(divisor, len_divisor) - j]);
			locator_rmd = locator + gf_degree(divisor, len_divisor) - j;
			remainder_tmp[locator_rmd] = gf_add(dividend_tmp[locator_rmd], factor_rmd);
		}
#if 0
		DEBUG_NOTICE("remainder_tmp:\n");
		for(k = 0; k < len_remainder; k++)
		{
			DEBUG_NOTICE("%x ", remainder_tmp[k]);
		}
		DEBUG_NOTICE("\n");
#endif
		if(gf_degree(divisor, len_divisor) > gf_degree(remainder_tmp, len_remainder))
		{
			for(k = 0; k < len_remainder; k++)
			{
				remainder[k] = remainder_tmp[k];
			}
			break;
		}
		else
		{
			memcpy(dividend_tmp, remainder_tmp, sizeof(unsigned char) * len_dividend);
			memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);
		}
	}

	return 0;
}

unsigned char gf_multp_poly(unsigned char* a, unsigned char len_a,
								unsigned char* b, unsigned char len_b,
								unsigned char* product, unsigned char len_product)
{
	unsigned char i = 0, j = 0, idx = 0;

	for(i = 0; i < len_a; i++)
	{
		for(j = 0; j < len_b; j++)
		{
			idx = i + j;
			if(len_product <= idx)
			{
				//DEBUG_NOTICE("product len err: %d\n", idx);
				continue;
			}
			product[idx] = gf_add(product[idx], gf_multp(a[i], b[j]));
		}
	}
}

int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 		   unsigned char* b, unsigned char len_b,
				 		   unsigned char* product, unsigned char len_product)
{
	unsigned int i = 0, j = 0, idx = len_product - 1;
	unsigned char reg[len_a - 1];
	memset(reg, 0xFF, sizeof(unsigned char) * (len_a - 1));
	unsigned char pd_tmp = 0xFF;

	/*high -> low*/
	for(i = 0; i < len_b; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), *(b + (len_b - 1 - i)));
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
		//DEBUG_NOTICE("unrel_group_seq: %x\n", a[0]);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		*(b + (len_b - 1 - i)));
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
#if 0
			DEBUG_NOTICE("idx: %d\n", len_a - 1 - j);
			DEBUG_NOTICE("%x %x %x %x %x\n", reg[len_a - 1 - j], reg[len_a - 1 - j - 1], pd_tmp,
								 	   *(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
#endif
		}
		reg[0] = gf_multp(*(a + 0), *(b + (len_b - 1 - i)));
		//DEBUG_NOTICE("%x %x %x\n", reg[0], *(a + 0), *(b + (len_b - 1 - i)));
	}

	for(i = 0; i < len_a; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), 0xFF);
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		0xFF);
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j - 1), 0xFF);
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
		}
		reg[0] = gf_multp(*(a + 0), 0xFF);
	}
	
	return 0;
}

long long real_combine(long long n, long long k)
{
#if (1 == GF_CAL_COUNT)	
	if(1 == cnt_switch)
	{
		real_cbm_cnt++;
	}
#endif
	long long combine_num = 0;

#if 0//it is useless when values are too large
	int i = 0;
	long tmp_n = 1, tmp_k = 1, tmp_n_k = 1;

	for(i = 1; i < (n + 1); i++)
	{
		tmp_n = tmp_n * i;
	}

	for(i = 1; i < (k + 1); i++)
	{
		tmp_k = tmp_k * i;
	}

	for(i = 1; i < (n - k + 1); i++)
	{
		tmp_n_k = tmp_n_k * i;
	}

	combine_num = tmp_n / tmp_k / tmp_n_k;
#else//fast calculation for finite field
	if(k == (n & k))
	{
		combine_num = 1;
	}
	else
	{
		combine_num = 2;
	}
#endif

	return combine_num;
}

unsigned char gf_real_mutp_ff(long long n, unsigned char ff)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		real_mul_ff_cnt++;
	}
#endif
	unsigned char val = 0xFF;

	if(0 != (n % 2))
	{
		val = ff;
	}
	else
	{
		val = 0xFF;
	}

	return val;
}

unsigned char gf_pow_cal(unsigned char ff, long long n)
{
#if (1 == GF_CAL_COUNT)
	if(1 == cnt_switch)
	{
		pow_cnt++;
	}
#endif	
	unsigned char val = 0xFF;
	if(0xFF == ff)
	{
		if(0 != n)
		{
			return 0xFF;
		}
		else
		{
			return 0x0;
		}
	}

	if(0 <= n)
	{
		val = (ff * n) % (GF_FIELD - 1);
	}
	else
	{
		val = (power_polynomial_table[-n + 1][0]) % (GF_FIELD - 1);
		val = (ff * val) % (GF_FIELD - 1);
	}

	return val;
}

unsigned char phase_trans(unsigned char phase)
{
	unsigned char val = 0;

	if(0 == phase)
	{
		val = 1;
	}
	else
	{
		val = 0;
	}

	return val;
}

#if (1 == GF_CAL_COUNT)
int gf_count_hist(long long err_cnt)
{
	if(err_num < (CODEWORD_LEN - MESSAGE_LEN - 1))
	{
		err_hist[err_num]++;
		add_cnt_hist[err_num] = add_cnt_hist[err_num] + add_cnt - add_cnt_prev;
		mul_cnt_hist[err_num] = mul_cnt_hist[err_num] + mul_cnt - mul_cnt_prev;;
		div_cnt_hist[err_num] = div_cnt_hist[err_num] + div_cnt - div_cnt_prev;;
		real_cbm_cnt_hist[err_num] = real_cbm_cnt_hist[err_num] + real_cbm_cnt - real_cbm_cnt_prev;;
		real_mul_ff_cnt_hist[err_num] = real_mul_ff_cnt_hist[err_num] + real_mul_ff_cnt - real_mul_ff_cnt_prev;;
		pow_cnt_hist[err_num] = pow_cnt_hist[err_num] + pow_cnt - pow_cnt_prev;;
	}
	else
	{
		err_hist[CODEWORD_LEN - MESSAGE_LEN - 2]++;
		add_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = add_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + add_cnt - add_cnt_prev;
		mul_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = mul_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + mul_cnt - mul_cnt_prev;;
		div_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = div_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + div_cnt - div_cnt_prev;;
		real_cbm_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = real_cbm_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + real_cbm_cnt - real_cbm_cnt_prev;;
		real_mul_ff_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = real_mul_ff_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + real_mul_ff_cnt - real_mul_ff_cnt_prev;;
		pow_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] = pow_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 2] + pow_cnt - pow_cnt_prev;;
	}

	add_cnt_prev = add_cnt;
	mul_cnt_prev = mul_cnt;
	div_cnt_prev = div_cnt;
	real_cbm_cnt_prev = real_cbm_cnt;
	real_mul_ff_cnt_prev = real_mul_ff_cnt;
	pow_cnt_prev = pow_cnt;

	return 0;
}

void gf_count_reset()
{
	add_cnt = 0;
	mul_cnt = 0;
	div_cnt = 0;
	real_cbm_cnt = 0;
	real_mul_ff_cnt = 0;
	pow_cnt = 0;
	add_cnt_prev = 0;
	mul_cnt_prev = 0;
	div_cnt_prev = 0;
	real_cbm_cnt_prev = 0;
	real_mul_ff_cnt_prev = 0;
	pow_cnt_prev = 0;

	return;
}
#endif

void gf_count_switch(unsigned char count_switch)
{
#if (1 == GF_CAL_COUNT)	
	cnt_switch = count_switch;
#endif
}

void BubbleSort4(float *A, int len, long long *A_idx)
{
    int low = 0, high = len - 1;
    int i = 0;
    float tmp = 0;

    while(low < high)
    {
        for(i = low; i < high; i++)  // 正向冒泡,找到最大者
        {
            if(A[i] > A[i + 1])
            {
                tmp = A[i];
                A[i] = A[i + 1];
                A[i + 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i + 1];
                A_idx[i + 1] = tmp;
            }
        }
        high--;            // 修改high值, 前移一位 
        for(i = high; i > low; i--)     // 反向冒泡,找到最小者 
        {
            if(A[i] < A[i - 1])
            {
                tmp = A[i];
                A[i] = A[i - 1];
                A[i - 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i - 1];
                A_idx[i - 1] = tmp;
            }
        }
        low++;            // 修改low值,后移一位
    }
}

void BubbleSort5(float *A, int len, long long *A_idx)
{
    int low = 0, high = len - 1;
    int i = 0;
    float tmp = 0;

    while(low < high)
    {
        for(i = low; i < high; i++)  // 正向冒泡,找到最小者
        {
            if(A[i] < A[i + 1])
            {
                tmp = A[i];
                A[i] = A[i + 1];
                A[i + 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i + 1];
                A_idx[i + 1] = tmp;
            }
        }
        high--;            // 修改high值, 前移一位 
        for(i = high; i > low; i--)     // 反向冒泡,找到最大者 
        {
            if(A[i] > A[i - 1])
            {
                tmp = A[i];
                A[i] = A[i - 1];
                A[i - 1] = tmp;

                tmp = A_idx[i];
                A_idx[i] = A_idx[i - 1];
                A_idx[i - 1] = tmp;
            }
        }
        low++;            // 修改low值,后移一位
    }
}

#if (1 == CFG_PARTIALLY_PARALLEL)
int gf_partially_parallel_init()
{
	gf_add_in_batch = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	gf_mul_in_batch = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	running_time = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	
	gf_add_in_batch_best = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	gf_mul_in_batch_best = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	gf_add_in_batch_worst = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	gf_mul_in_batch_worst = (long long*)malloc(sizeof(long long) * PARALLEL_BATCH_NUM);
	
	memset(gf_add_in_batch, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	memset(gf_mul_in_batch, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	memset(running_time, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	
	memset(gf_add_in_batch_best, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	memset(gf_mul_in_batch_best, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	memset(gf_add_in_batch_worst, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	memset(gf_mul_in_batch_worst, 0, sizeof(long long) * PARALLEL_BATCH_NUM);
	
	long long batch_size = pow_val / PARALLEL_BATCH_NUM;
	cmplx_per_round_add = (long long*)malloc(sizeof(long long) * batch_size);
	latency_per_round_add = (long long*)malloc(sizeof(long long) * batch_size);
	memset(cmplx_per_round_add, 0, sizeof(long long) * batch_size);
	memset(latency_per_round_add, 0, sizeof(long long) * batch_size);
	cmplx_per_round_mul = (long long*)malloc(sizeof(long long) * batch_size);
	latency_per_round_mul = (long long*)malloc(sizeof(long long) * batch_size);
	memset(cmplx_per_round_mul, 0, sizeof(long long) * batch_size);
	memset(latency_per_round_mul, 0, sizeof(long long) * batch_size);
	
	return 0;
}

int gf_partially_parallel_exit()
{
	free(gf_add_in_batch);
	gf_add_in_batch = NULL;
	free(gf_mul_in_batch);
	gf_mul_in_batch = NULL;
	free(running_time);
	running_time = NULL;
	
	free(gf_add_in_batch_best);
	gf_add_in_batch_best = NULL;
	free(gf_mul_in_batch_best);
	gf_mul_in_batch_best = NULL;
	free(gf_add_in_batch_worst);
	gf_add_in_batch_worst = NULL;
	free(gf_mul_in_batch_worst);
	gf_mul_in_batch_worst = NULL;
	
	free(cmplx_per_round_add);
	cmplx_per_round_add = NULL;
	free(latency_per_round_add);
	latency_per_round_add = NULL;
	free(cmplx_per_round_mul);
	cmplx_per_round_mul = NULL;
	free(latency_per_round_mul);
	latency_per_round_mul = NULL;
	
	return 0;
}

int gf_partially_parallel_cnt(long long batch_idx)
{
	storage_cnt = 0;

	long long add_inc = add_cnt - add_cnt_prev_batch;
	long long mul_inc = (mul_cnt + div_cnt)
						- (mul_cnt_prev_batch + div_cnt_prev_batch);

	gf_add_in_batch[batch_idx] = gf_add_in_batch[batch_idx] + add_inc;
	gf_mul_in_batch[batch_idx] = gf_mul_in_batch[batch_idx] + mul_inc;
	running_time[batch_idx] = 0;

	if((0 == gf_add_in_batch_best[batch_idx])
	   && (0 == gf_mul_in_batch_best[batch_idx]))
	{
		gf_add_in_batch_best[batch_idx] = add_inc;
		gf_mul_in_batch_best[batch_idx] = mul_inc;
	}

	if((add_inc < gf_add_in_batch_best[batch_idx])
		&& (0 != add_inc))
	{
		gf_add_in_batch_best[batch_idx] = add_inc;
	}
	if((mul_inc < gf_mul_in_batch_best[batch_idx])
		&& (0 != mul_inc))
	{
		gf_mul_in_batch_best[batch_idx] = mul_inc;
	}

	if(add_inc > gf_add_in_batch_worst[batch_idx])
	{
		gf_add_in_batch_worst[batch_idx] = add_inc;
	}
	if(mul_inc > gf_mul_in_batch_worst[batch_idx])
	{
		gf_mul_in_batch_worst[batch_idx] = mul_inc;
	}

	add_cnt_prev_batch = add_cnt;
	mul_cnt_prev_batch = mul_cnt;
	div_cnt_prev_batch = div_cnt;

	cost_parallel = 0;
	speed_up = 0;
	efficiency = 0;

	return 0;
}
#endif
