#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define GF_Q			 8
#define GF_CAL_COUNT	 0

#define SYS_ENC			 1

#define TEST_MODE		 0

#define EARLY_TERMINATION		1
#define EARLY_TERMINATION_NUM	100
#define OUTPUT_LOG				1

#define RELEX_ORDER				1
#define SIMPLE_ASD				0
#define RE_ENCODING				1
#define RECUR_RR				1
#define DYNAMIC_MEM				1
#define REDUNDANT_SIZE			0
#define DYNAMIC_TERM			1
#if (1 == DYNAMIC_TERM)
#define DYNAMIC_TERM_ITER		1
#define DYNAMIC_TERM_X			9/8
#define DYNAMIC_TERM_Y			10/9
#endif
#if (1 == RE_ENCODING)
#define SYN_LEN		(MESSAGE_LEN / 2 * 2)//need to be checked
#endif

#define S_MUL					1
#define LEX_TABLE_EXPAND_SIZE	4
#define YITA					1

#if (1 == RE_ENCODING)
#define Y_WEIGHT				(-1)
#else
#define Y_WEIGHT				(MESSAGE_LEN - 1)
#endif

/*RR_MODE for Re-encoding:
0: Normal RR
1: BMA-Erasure_Decoding
2: Fast RR for m=1
3: Normal RR for Systematic Encoding, tmp for m=1
4: RR for Evaluating Encoding, tmp for m=1
5: RR with IDFT for Evaluating Encoding, tmp for m=1*/
#define BMA_RR				1
#define FAST_RR_M1			2
#define CONV_RE_ENC_SYS		3
#define CONV_RE_ENC			4
#define CONV_RE_ENC_IDFT	5
#if (0 == RE_ENCODING)
#define CFG_RR_MODE			0
#else
#define CFG_RR_MODE			FAST_RR_M1
#endif
#endif
