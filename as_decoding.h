#ifndef AS_DECODING_H
#define AS_DECODING_H

#include "cfg_decoding.h"

//#define K_M						S_MUL * (S_MUL - 1) * MESSAGE_LEN
#define LAYER_NUM				(MESSAGE_LEN + 1)
//#define LAYER_NUM				2
/*approximate and sufficient size allocation*/
#define C_SIZE					(S_MUL * (S_MUL + 1) * CODEWORD_LEN / 2)
//#define TERM_SIZE_Y				(((((sqrt(8 * C_SIZE / (MESSAGE_LEN - 1) + 1)) + 1) / 2) - 1) << 1)
//#define TERM_SIZE_Y				(8 * (S_MUL + 1))
//#define TERM_SIZE_X				((C_SIZE / ((TERM_SIZE_Y >> 0) + 1) + (TERM_SIZE_Y >> 0) / (MESSAGE_LEN - 1) / 2) << 1)
#define TERM_SIZE				((5 * (S_MUL + 1) * CODEWORD_LEN) / 6)//REAL_SIZE = TERM_SIZE^2
#define TERM_SIZE_Y				15
#define TERM_SIZE_X				TERM_SIZE
#if (0 == RECUR_RR)
#define POLY_NUM				TERM_SIZE_X * TERM_SIZE_Y * LAYER_NUM// notice there may be dangerous! POLY_NUM = ROOT_CNT_LAST_LAYER * GF_FIELD
#define ROOT_SIZE				POLY_NUM
#else
#define POLY_NUM				1// notice there may be dangerous! POLY_NUM = ROOT_CNT_LAST_LAYER * GF_FIELD
#define ROOT_SIZE				1
#endif

extern unsigned char received_polynomial[CODEWORD_LEN];
extern unsigned char output_polynomial[CODEWORD_LEN];
extern unsigned char decoded_codeword[CODEWORD_LEN];
extern unsigned char decoded_message[MESSAGE_LEN];

extern long long err_num;
extern unsigned char decoding_ok_flag;
extern long long weight_stored;
extern long long hamm_distance_debug;
extern long long rr_err;
extern long long max_dx, max_dy;
extern long long term_size_x, term_size_y;
extern unsigned char ***g_term_c_p;

extern int as_decoding();
extern int g_term_malloc();
extern int g_term_destroy();
extern int dfs_rr_recur(unsigned char *g_c_q,
					unsigned char *g_c_0_y,
					long long l);
extern long long hamm_distance_code_cal(unsigned char *a,
									  					  unsigned char *b,
									  					  long long len);
extern long long hamm_distance_bit_cal(unsigned char *a,
									  			  unsigned char *b,
									  			  long long len);

#endif
