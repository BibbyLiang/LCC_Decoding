#ifndef GF_CAL_H
#define GF_CAL_H

#include "cfg_decoding.h"

#if (3 == GF_Q)
#define GF_FIELD        8
#define MESSAGE_LEN     5
#endif
#if (4 == GF_Q)
#define GF_FIELD        16
#define MESSAGE_LEN     7
#endif
#if (6 == GF_Q)
#define GF_FIELD        64
#define MESSAGE_LEN     21
#endif
#if (8 == GF_Q)
#define GF_FIELD        256
#define MESSAGE_LEN     239
#endif
#define CODEWORD_LEN    (GF_FIELD - 1)   

extern unsigned char power_polynomial_table[GF_FIELD][2];

#if (1 == GF_CAL_COUNT)
extern unsigned char cnt_switch;
extern long long add_cnt;
extern long long mul_cnt;
extern long long div_cnt;
extern long long real_cbm_cnt;
extern long long real_mul_ff_cnt;
extern long long pow_cnt;
extern long long add_cnt_prev;
extern long long mul_cnt_prev;
extern long long div_cnt_prev;
extern long long real_cbm_cnt_prev;
extern long long real_mul_ff_cnt_prev;
extern long long pow_cnt_prev;
extern long long err_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long add_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long mul_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long div_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long real_cbm_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long real_mul_ff_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
extern long long pow_cnt_hist[CODEWORD_LEN - MESSAGE_LEN - 1];
#endif

extern unsigned char gf_pow2poly(unsigned char val_in_pow);
extern unsigned char gf_poly2pow(unsigned char val_in_poly);
extern unsigned char gf_location(unsigned char val);
extern unsigned char gf_add(unsigned char a, unsigned char b);
extern unsigned char gf_multp(unsigned char a, unsigned char b);
extern unsigned char gf_div(unsigned char a, unsigned char b);
extern unsigned char gf_mod_single_term(unsigned char a, unsigned char b);
extern unsigned char gf_degree(unsigned char* a, unsigned char len_a);
extern unsigned char gf_div_q_r(unsigned char* dividend, unsigned char len_dividend,
								   unsigned char* divisor, unsigned char len_divisor,
								   unsigned char* quotien, unsigned char len_quotien,
								   unsigned char* remainder, unsigned char len_remainder);
extern unsigned char gf_multp_poly(unsigned char* a, unsigned char len_a,
									   unsigned char* b, unsigned char len_b,
									   unsigned char* product, unsigned char len_product);

extern int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 				  unsigned char* b, unsigned char len_b,
				 				  unsigned char* product, unsigned char len_product);
extern long long real_combine(long long n, long long k);
extern unsigned char gf_real_mutp_ff(long long n, unsigned char ff);
unsigned char gf_pow_cal(unsigned char ff, long long n);
extern unsigned char phase_trans(unsigned char phase);
#if (1 == GF_CAL_COUNT)
extern int gf_count_hist(long long err_cnt);
#endif
#endif
