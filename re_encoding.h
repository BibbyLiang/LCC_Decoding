#ifndef RE_ENCODING_H
#define RE_ENCODING_H

#include "cfg_decoding.h"

extern float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern float chnl_rel_matrix_tmp[CODEWORD_LEN + 1][CODEWORD_LEN];
extern unsigned char mul_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern unsigned char beta_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern unsigned char rel_group_seq[MESSAGE_LEN];
extern unsigned char unrel_group_seq[CODEWORD_LEN - MESSAGE_LEN];
extern unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
extern unsigned char tao[CODEWORD_LEN - MESSAGE_LEN + 1];
extern unsigned char omega[CODEWORD_LEN - MESSAGE_LEN];
extern unsigned char sigma[((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)];
extern unsigned char erasure_polynomial[CODEWORD_LEN];
extern unsigned char phi[CODEWORD_LEN];
extern unsigned char v[MESSAGE_LEN + 1];
#if (1 == RE_ENCODING)
extern unsigned char re_encoded_codeword[CODEWORD_LEN];
#endif
extern unsigned char rel_group_seq[MESSAGE_LEN];
extern unsigned char unrel_group_seq[CODEWORD_LEN - MESSAGE_LEN];
extern long long chnl_rel_order_idx[CODEWORD_LEN];
extern long long chnl_rel_max_id[CODEWORD_LEN];
extern long long chnl_rel_scd_id[CODEWORD_LEN];

extern int chnl_rel_cal(float **input_seq, long long input_len);
extern int mul_assign();
extern int l_cal(unsigned char locator_j, unsigned char *L);
extern int tao_cal();
extern unsigned char poly_eva(unsigned char *poly, long long poly_len, unsigned char input_val);
extern int bm_re_encoding(unsigned char *msg_phi, unsigned char *tmp_cw);
extern int re_encoding();
extern int recover_codeword();
extern int chnl_rel_seq_order();

#endif
