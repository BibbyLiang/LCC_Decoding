#ifndef ENCODING_H
#define ENCODING_H

#include "cfg_decoding.h"

extern unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1];
extern unsigned char message_polynomial[MESSAGE_LEN];
extern unsigned char error_polynomial[CODEWORD_LEN];
extern unsigned char encoded_polynomial[CODEWORD_LEN];

extern unsigned char evaluation_encoding();
extern unsigned char systematic_encoding();
extern unsigned char evaluation_encoding_v2(unsigned char *message,
											       unsigned char *codeword_output);
extern unsigned char systematic_encoding_v2(unsigned char *message,
									               unsigned char *codeword_output);
extern int gen_poly_trans();
extern void encoding2();
#endif