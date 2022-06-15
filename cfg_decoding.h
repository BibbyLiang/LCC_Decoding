#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define GF_Q			 8
#define GF_CAL_COUNT	 1

#define SYS_ENC			 1

#define TEST_MODE		 0

#define EARLY_TERMINATION		0
#define EARLY_TERMINATION_NUM	100
#define OUTPUT_LOG				0
#define DEBUG_LOG				0

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
#define SYN_LEN		(MESSAGE_LEN * 2 / 20)//need to be checked
#endif
#define CFG_FAST_SKIP_TST_VCT	0

#define CFG_PARTIALLY_PARALLEL	1
#define PARALLEL_BATCH_NUM		1//BATCH_NUM * BATCH_SIZE = TOTAL_NUM
#define CFG_PARALLEL_FAST_SKIP	0
#define CFG_STORE_PARALEL		0
#define CFG_STORE_LEN			0
#define CFG_GROUP_SCHEME		1

#define CFG_Q0_FAST_SKIP		1//44
#define CFG_PRG_DECODING		0
#define CFG_IMD_STORE			1

#define CFG_ADAPTIVE_PARALLEL	0
#define CFG_ADAPTIVE_SIZE		8
#define CFG_CHNL_REL_THRD		50
#define CFG_ADAPTIVE_DELTA		8
#define CFG_GROUP_METHOD		1

#define CFG_RE_ARRANGE_TST_VCT	0//not finish now, need to be checked and move to adaptive mode

#define CFG_Y_ROOTS_SKIP		1

#define BF_INTERPOLATION		0

#define S_MUL					1
#define LEX_TABLE_EXPAND_SIZE	4
#define YITA					0

#define SHORTEN_LEN				7

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
#define CFG_RR_MODE			CONV_RE_ENC_SYS
#endif

/*there may be some problems. and it should be used with BF_INTERPOLATION*/
#define TERM_SIZE_DBG		1//33
#define FAST_RR_M1_DBG		1

#endif
