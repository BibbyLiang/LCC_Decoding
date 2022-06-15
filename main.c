#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"
#include "re_encoding.h"
#include "encoding.h"
#include "mod.h"
#include "rnd.h"
#include "channel.h"
#include "time.h"

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

#if (1 == GF_CAL_COUNT)
	memset(err_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(add_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(mul_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(div_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(real_cbm_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(real_mul_ff_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(pow_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
#endif
	DEBUG_SYS("init_simulation OK\n");

	g_term_malloc();

	mod_init();

#if (1 == SYS_ENC)
	gen_poly_trans();
#endif

#if (1 == CFG_PARTIALLY_PARALLEL)
	gf_partially_parallel_init();
#endif

	return;
}

void main()
{
	long long i = 0, j = 0, k = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;
	long long iter = 0;
	long long bit_err = 0, symbol_err = 0, symbol_err_prev = 0, frame_err = 0;
	long long uncoded_bit_err = 0, uncoded_symbol_err = 0, uncoded_frame_err = 0;
	unsigned char frame_err_flag = 0, uncoded_frame_err_flag = 0;
	unsigned char tmp_mes = 0, tmp_dec = 0;
	long long hamm_err = 0;
#if (0 == RECUR_RR)	
	long long rr_err_cnt = 0;
#endif
#if (1 == GF_CAL_COUNT)	
	long long symbol_err_this_frame = 0;
#endif
#if (1 == CFG_PARTIALLY_PARALLEL)
	long long batch_size = 0;
	double latency_add = 0, latency_mul = 0;
#endif	
	float avg_round = 0;

	clock_t start, stop;
	float runtime;

	/*input simulation parameters*/
	float eb2n0_start = 2.5, eb2n0_stop = 2.5, eb2n0_step = 1, eb2n0 = 5;
	long long iter_cnt = 1, monitor_cnt = 1;
#if (0 == TEST_MODE)
#if 1
	printf("Please Input Eb/N0 Start: ");
	scanf("%f", &eb2n0_start);
	printf("Please Input Eb/N0 Stop: ");
	scanf("%f", &eb2n0_stop);
	printf("Please Input Eb/N0 Step: ");
	scanf("%f", &eb2n0_step);
	printf("Please Input Simulation Times: ");
	scanf("%ld", &iter_cnt);
	printf("Please Monitor Times: ");
	scanf("%ld", &monitor_cnt);
#endif	
#endif

#if (1 == OUTPUT_LOG)
	/*init file log*/
	char log_name[255];
	sprintf(log_name, "n_%d-k_%d-m_%d-yita_%d-snr_%f_%f_%f-cnt_%ld_%ld.txt",
					  CODEWORD_LEN - SHORTEN_LEN,
					  MESSAGE_LEN - SHORTEN_LEN,
					  S_MUL,
					  YITA,
					  eb2n0_start,
					  eb2n0_step,
					  eb2n0_stop,
					  iter_cnt,
					  monitor_cnt);
	FILE *frc;
#endif

#if (1 == DEBUG_LOG)
	frc_debug = fopen("debug_aaa_log.txt", "a+");
	//fclose(frc_debug);
	//frc_debug = NULL;
#endif

	init_simulation();

	start = clock();

	/*for SNR*/
	for(eb2n0 = eb2n0_start; eb2n0 <= eb2n0_stop; eb2n0 = eb2n0 + eb2n0_step)
	{
		/*init counts for certain SNR*/
		gf_count_reset();
		bit_err = 0;
		symbol_err = 0;
		symbol_err_prev = 0;
		frame_err = 0;
		uncoded_bit_err = 0;
		uncoded_symbol_err = 0;
		uncoded_frame_err = 0;
		hamm_err = 0;
		time_measure = 0;
		memset(skip_hist, 0, sizeof(long long) * pow_val);
		memset(pgd_hist, 0, sizeof(long long) * pow_val);
#if (1 == CFG_ADAPTIVE_PARALLEL)
		memset(round_hist, 0, sizeof(long long) * (pow_val + 1));
#else
		memset(round_hist, 0, sizeof(long long) * (pow_val / PARALLEL_BATCH_NUM + 1));
#endif
		memset(latency_per_round_add, 0, sizeof(long long) * batch_size);
		memset(latency_per_round_mul, 0, sizeof(long long) * batch_size);

		/*for every frame*/
		for(iter = 0; iter < iter_cnt; iter++)
		{
#if (1 == DEBUG_LOG)
			printf("iter: %ld\n", iter);
			fclose(frc_debug);
			frc_debug = NULL;
			frc_debug = fopen("debug_aaa_log.txt", "w+");
#endif

			/*clear some counts for this frame*/
			decoding_ok_flag = 0;
			err_num = 0;
			memset(recv_rel, 0.0, sizeof(float) * CODEWORD_LEN);
			symbol_err_prev = uncoded_symbol_err;

			/*generate messages*/
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				j = genrand_int32() % GF_FIELD;
#if (1 == TEST_MODE)
				message_polynomial[i] = 0x0;
#else
				message_polynomial[i] = power_polynomial_table[j][0];
#endif
			}
#if (0 == TEST_MODE)			
			shorten_msg(MESSAGE_LEN, SHORTEN_LEN);
#endif			
			//memset(message_polynomial, 0x0, sizeof(unsigned char) * MESSAGE_LEN);

#if (1 == TEST_MODE)
			/*generate messages for test mode*/
			message_polynomial[0] = 0x2;
			message_polynomial[1] = 0x2;
			message_polynomial[2] = 0x6;
			message_polynomial[3] = 0x2;
			message_polynomial[4] = 0x3;
			//message_polynomial[5] = 0xC;
			//message_polynomial[6] = 0xA;
			//memset(message_polynomial, 0x0, sizeof(unsigned char) * MESSAGE_LEN);
#endif

			/*encoding*/
#if (1 == SYS_ENC)
			//systematic_encoding();
			encoding2();
#else
			evaluation_encoding();
#endif

			/*modulation*/
			bpsk_mod(encoded_polynomial,
					 CODEWORD_LEN,
					 recv_seq,
					 symbol_num);

			/*transmission over channel*/
			DEBUG_IMPOTANT("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				if(0 == (i % GF_Q))
				{
					DEBUG_NOTICE("---------------\n");
				}
				recv_seq[i][0] = recv_seq[i][0] + awgn_gen(eb2n0);
				recv_seq[i][1] = recv_seq[i][1] + awgn_gen(eb2n0);
				DEBUG_NOTICE("%f %f\n", recv_seq[i][0], recv_seq[i][1]);
			}
			DEBUG_IMPOTANT("\n");

			/*demodulation*/
			bpsk_demod((float **)recv_seq,
					 symbol_num,
					 received_polynomial,
					 CODEWORD_LEN);
#if (0 == TEST_MODE)
			shorten_trans(CODEWORD_LEN, SHORTEN_LEN);
#endif

			/*add errors for test mode*/
#if (1 == TEST_MODE)//test
#if 0
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
			}
#else
			//memcpy(received_polynomial, encoded_polynomial, sizeof(unsigned char) * CODEWORD_LEN);
			//received_polynomial[0] = gf_add(encoded_polynomial[0], 0x1);
			//received_polynomial[3] = gf_add(encoded_polynomial[3], 0x1);
			//received_polynomial[4] = gf_add(encoded_polynomial[4], 0x2);
			//received_polynomial[5] = gf_add(encoded_polynomial[5], 0x0);
			//received_polynomial[9] = gf_add(encoded_polynomial[9], 0x2);
			//received_polynomial[0] = 0x5;
			//received_polynomial[1] = 0x4;
			//received_polynomial[2] = 0x4;
			//received_polynomial[3] = 0x0;
			//received_polynomial[4] = 0x1;
			//received_polynomial[5] = 0x1;
			//received_polynomial[6] = 0x0;
#endif
#endif

			/*count uncoded errors*/
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_NOTICE("%x %x\n", received_polynomial[i], encoded_polynomial[i]);
				if(received_polynomial[i] != encoded_polynomial[i])
				{
					uncoded_symbol_err = uncoded_symbol_err + 1;
					err_num = err_num + 1;
					if(0 == uncoded_frame_err_flag)
					{
						uncoded_frame_err_flag = 1;
						uncoded_frame_err = uncoded_frame_err + 1;
					}
					for(j = 0; j < GF_Q; j++)
					{
						tmp_mes = ((received_polynomial[i] >> j) & 0x1);
						tmp_dec = ((encoded_polynomial[i] >> j) & 0x1);
						if(tmp_mes != tmp_dec)
						{
							uncoded_bit_err = uncoded_bit_err + 1;
						}
					}
#if (1 == GF_CAL_COUNT)
					symbol_err_this_frame++;
#endif
				}
			}

			/*channel reliability calculation*/
#if (1 == RE_ENCODING)
			chnl_rel_cal(recv_seq, symbol_num);
#endif

			/*nultiplicity assignment*/
			mul_assign();
			gf_count_switch(1);
			/*re-encoding transform*/

			re_encoding();

			/*GS decoding*/
			as_decoding();
			gf_count_switch(0);
#if (1 == GF_CAL_COUNT)
			/*count gf field calculating complexity*/
			gf_count_hist(symbol_err_this_frame);
			symbol_err_this_frame = 0;
#endif
			uncoded_frame_err_flag = 0;

			/*count decoded errors*/
#if (1 == SYS_ENC)			
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_NOTICE("%x %x\n", decoded_message[i], message_polynomial[i]);
				if(decoded_message[i] != message_polynomial[i])
				{
					symbol_err = symbol_err + 1;
					if(0 == frame_err_flag)
					{
						frame_err_flag = 1;
						frame_err = frame_err + 1;
					}
					for(j = 0; j < GF_Q; j++)
					{
						tmp_mes = ((message_polynomial[i] >> j) & 0x1);
						tmp_dec = ((decoded_message[i] >> j) & 0x1);
						if(tmp_mes != tmp_dec)
						{
							bit_err = bit_err + 1;
						}
					}
				}
			}
#else
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_NOTICE("%x %x\n", decoded_codeword[i], encoded_polynomial[i]);
				if(decoded_codeword[i] != encoded_polynomial[i])
				{
					symbol_err = symbol_err + 1;
					if(0 == frame_err_flag)
					{
						frame_err_flag = 1;
						frame_err = frame_err + 1;
					}
					for(j = 0; j < GF_Q; j++)
					{
						tmp_mes = ((encoded_polynomial[i] >> j) & 0x1);
						tmp_dec = ((decoded_codeword[i] >> j) & 0x1);
						if(tmp_mes != tmp_dec)
						{
							bit_err = bit_err + 1;
						}
					}
				}
			}
#endif

			if((((CODEWORD_LEN - MESSAGE_LEN) / 2) >= (uncoded_symbol_err - symbol_err_prev))
			   && (1 == frame_err_flag))
			{
				DEBUG_SYS("Radius Err. for Decoding 1: %d %d %ld\n",
						  frame_err_flag,
						  (uncoded_symbol_err - symbol_err_prev),
						  best_tst_vct_diff);
				for(i = 0; i < CODEWORD_LEN; i++)
				{
					if((encoded_polynomial[i] != received_polynomial[i])
						|| (encoded_polynomial[i] != tst_vct_debug[i]))
					{
						DEBUG_SYS("Radius Check: %ld | %d %d %d\n",
								  i,
								  encoded_polynomial[i],
								  received_polynomial[i],
								  tst_vct_debug[i]);
					}
				}
#if (1 == OUTPUT_LOG)					
				frc = fopen(log_name, "a+");
				fprintf(frc, "Radius Err. for Decoding 1: %d %d %ld\n",
							 frame_err_flag,
							 (uncoded_symbol_err - symbol_err_prev),
							 best_tst_vct_diff);
				for(i = 0; i < CODEWORD_LEN; i++)
				{
					if((encoded_polynomial[i] != received_polynomial[i])
						|| (encoded_polynomial[i] != tst_vct_debug[i]))
					{
						fprintf(frc, "Radius Check: %ld | %d %d %d\n",
								i,
								encoded_polynomial[i],
								received_polynomial[i],
								tst_vct_debug[i]);
					}
				}			 
				fclose(frc);
				frc = NULL;
#endif		  
			}
			
			if((((CODEWORD_LEN - MESSAGE_LEN) / 2) >= best_tst_vct_diff)
			   && (1 == frame_err_flag))
			{
				DEBUG_SYS("Radius Err. for Decoding 2: %d %d\n",
						  frame_err_flag,
						  best_tst_vct_diff);
#if (1 == DEBUG_LOG)
				printf("Radius Err. for Decoding 2: %d %d\n",
						frame_err_flag,
						best_tst_vct_diff);
				fclose(frc_debug);
				frc_debug = NULL;	  	
				while(1);	
#endif				
						  
#if (1 == OUTPUT_LOG)					
				frc = fopen(log_name, "a+");
				fprintf(frc, "Radius Err. for Decoding 2: %d %d\n",
						  	 frame_err_flag,
						  	 best_tst_vct_diff);
				fclose(frc);
				frc = NULL;
#endif
			}

			/*print program error*/
			if(((1 == decoding_ok_flag)
				&& (1 == frame_err_flag))
				|| ((((CODEWORD_LEN - MESSAGE_LEN) / 2) >= (uncoded_symbol_err - symbol_err_prev))
					&& (1 == frame_err_flag)))
			{
				if(0 == hamm_distance_debug)
				{
					hamm_err = hamm_err + 1;

					for(i = 0; i < MESSAGE_LEN; i++)
					{
						DEBUG_SYS("Strage Err: %x %x\n",
							      decoded_message[i],
							      message_polynomial[i]);
					}
				}
				else
				{
					DEBUG_SYS("Prog. Err. for Decoding\n");
#if (1 == OUTPUT_LOG)					
					frc = fopen(log_name, "a+");
					fprintf(frc, "Prog. Err. for Decoding\n");
#endif					

					DEBUG_SYS("Para.: %ld %ld %ld %ld %ld %ld\n",
							  (S_MUL * (CODEWORD_LEN - err_num)),
							  weight_stored,
							  hamm_distance_debug,
							  uncoded_symbol_err,
							  symbol_err_prev,
							  uncoded_symbol_err - symbol_err_prev);
#if (1 == OUTPUT_LOG)
					fprintf(frc, "Para.: %ld %ld %ld %ld %ld %ld\n",
							(S_MUL * (CODEWORD_LEN - err_num)),
							weight_stored,
							hamm_distance_debug,
							uncoded_symbol_err,
							symbol_err_prev,
							uncoded_symbol_err - symbol_err_prev);
#endif

#if 0

#if (1 == RE_ENCODING)

					for(i = 0; i < MESSAGE_LEN; i++)
					{
						DEBUG_SYS("message: %x\n", message_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "message: %x\n", message_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("encoded: %x\n", encoded_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "encoded: %x\n", encoded_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("recv: %x\n", received_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "recv: %x\n", received_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
#endif
					}

#endif
					for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
					{
						DEBUG_SYS("unrel_group_seq: %d\n", unrel_group_seq[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "unrel_group_seq: %d\n", unrel_group_seq[i]);
#endif
					}
					
#endif					
					
#if (1 == OUTPUT_LOG)
					fclose(frc);
					frc = NULL;
#endif
				}
			}
			frame_err_flag = 0;
#if (0 == RECUR_RR)
			if(1 == rr_err)
			{
				rr_err_cnt = rr_err_cnt + 1;
			}
#endif

			/*output ongoing results*/
			if(0 == (iter % monitor_cnt))
			{
				stop = clock();
				runtime = (stop - start) / 1000.0000;
				
				DEBUG_SYS("---------------------\n");
				DEBUG_SYS("Time: %fs\n", runtime);
				DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
				DEBUG_SYS("Frame: %ld\n", iter + 1);
				DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
				DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
				DEBUG_SYS("Frame Error: %ld\n", frame_err);
				DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
				DEBUG_SYS("Bit Error: %ld\n", bit_err);
				//DEBUG_SYS("Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)				
				DEBUG_SYS("RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)
				DEBUG_SYS("Add/Mul Cnt: %f %f\n", (double)add_cnt / (double)(iter + 1), (double)(mul_cnt + div_cnt) / (double)(iter + 1));
				//DEBUG_SYS("Add Cnt: %ld\n", add_cnt);
				//DEBUG_SYS("Mul Cnt: %ld\n", mul_cnt);
				//DEBUG_SYS("Div Cnt: %ld\n", div_cnt);
				//DEBUG_SYS("RCB Cnt: %ld\n", real_cbm_cnt);
				//DEBUG_SYS("RMF Cnt: %ld\n", real_mul_ff_cnt);
				//DEBUG_SYS("Pow Cnt: %ld\n", pow_cnt);
#endif				
				DEBUG_SYS("Max DX: %ld\n", max_dx);
				DEBUG_SYS("Max DY: %ld\n", max_dy);
				DEBUG_SYS("TERM_SIZE: %ld %ld\n", term_size_x, term_size_y);

#if (1 == CFG_PARTIALLY_PARALLEL)
				for(k = 0; k < PARALLEL_BATCH_NUM; k++)
				{
					//DEBUG_SYS("batch_add: %ld %ld\n", k, gf_add_in_batch[k]);
					//DEBUG_SYS("batch_mul: %ld %ld\n", k, gf_mul_in_batch[k]);
					DEBUG_INFO("batch_add: %ld %f %ld %d\n",
					          k,
					          (double)gf_add_in_batch[k] / (double)(iter + 1),
					          gf_add_in_batch_best[k],
					          gf_add_in_batch_worst[k]);
					DEBUG_INFO("batch_mul: %ld %f %ld %ld\n",
					          k,
					          (double)gf_mul_in_batch[k] / (double)(iter + 1),
					          gf_mul_in_batch_best[k],
					          gf_mul_in_batch_worst[k]);
				}
				
				batch_size = pow_val / PARALLEL_BATCH_NUM;
				latency_add = 0, latency_mul = 0;
				for(k = 0; k < batch_size; k++)
				{
					DEBUG_INFO("round_add: %ld %f %f\n",
					          k,
					          (double)cmplx_per_round_add[k] / (double)(iter + 1),
					          (double)latency_per_round_add[k] / (double)(iter + 1));
					DEBUG_INFO("round_mul: %ld %f %f\n",
					          k,
					          (double)cmplx_per_round_mul[k] / (double)(iter + 1),
					          (double)latency_per_round_mul[k] / (double)(iter + 1));
					          
					latency_add = latency_add + (double)latency_per_round_add[k] / (double)(iter + 1);
					latency_mul = latency_mul + (double)latency_per_round_mul[k] / (double)(iter + 1);
				}
				DEBUG_SYS("Latency: %f %f\n", latency_add, latency_mul);
				avg_round = 0;
#if (1 == CFG_ADAPTIVE_PARALLEL)
				for(i = 0; i < (pow_val + 1); i++)
#else
				for(i = 0; i < (pow_val / PARALLEL_BATCH_NUM + 1); i++)
#endif				
				{
					avg_round = avg_round + (float)(round_hist[i] * (i));
				}
				DEBUG_SYS("avg_round: %f\n", avg_round / (float)(iter + 1));
#endif		

#if (1 == OUTPUT_LOG)
				frc = fopen(log_name, "a+");
				fprintf(frc, "---------------------\n");
				fprintf(frc, "Time: %fs\n", runtime);
				fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
				fprintf(frc, "Frame: %ld\n", iter + 1);
				fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
				fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
				fprintf(frc, "Frame Error: %ld\n", frame_err);
				fprintf(frc, "Symbol Error: %ld\n", symbol_err);
				fprintf(frc, "Bit Error: %ld\n", bit_err);
				//fprintf(frc, "Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)				
				fprintf(frc, "RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)
				fprintf(frc, "Add/Mul Cnt: %f %f\n", (double)add_cnt / (double)(iter + 1), (double)(mul_cnt + div_cnt) / (double)(iter + 1));
				//fprintf(frc, "Add Cnt: %ld\n", add_cnt);
				//fprintf(frc, "Mul Cnt: %ld\n", mul_cnt);
				//fprintf(frc, "Div Cnt: %ld\n", div_cnt);
				//fprintf(frc, "RCB Cnt: %ld\n", real_cbm_cnt);
				//fprintf(frc, "RMF Cnt: %ld\n", real_mul_ff_cnt);
				//fprintf(frc, "Pow Cnt: %ld\n", pow_cnt);
#endif
				fprintf(frc, "TERM_SIZE: %ld %ld\n", term_size_x, term_size_y);
				fprintf(frc, "Latency: %f %f\n", latency_add, latency_mul);
				fprintf(frc, "avg_round: %f\n", avg_round / (float)(iter + 1));
			    fclose(frc);
				frc = NULL;
#endif

#if (1 == GF_CAL_COUNT)
				for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN - 1); k++)
				{
					if(0 != err_hist[k])
					{
						//DEBUG_SYS("-------------------\n");
						//DEBUG_SYS("Err Hist: %ld %ld\n", k, err_hist[k]);
						//DEBUG_SYS("Add/Mul Hist: %ld %f %f\n", k, (float)(add_cnt_hist[k] / err_hist[k]), (float)((mul_cnt_hist[k] + div_cnt_hist[k]) / err_hist[k]));
						//DEBUG_SYS("Mul Hist: %ld %f\n", k, (float)(mul_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("Div Hist: %ld %f\n", k, (float)(div_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("RCB Hist: %ld %f\n", k, (float)(real_cbm_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("RMF Hist: %ld %f\n", k, (float)(real_mul_ff_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("Pow Hist: %ld %f\n", k, (float)(pow_cnt_hist[k] / err_hist[k]));

#if (1 == OUTPUT_LOG)
						frc = fopen(log_name, "a+");
						//fprintf(frc, "-------------------\n");
						//fprintf(frc, "Err Hist: %ld %ld\n", k, err_hist[k]);
						//fprintf(frc, "Add/Mul Hist: %ld %f %f\n", k, (float)(add_cnt_hist[k] / err_hist[k]), (float)((mul_cnt_hist[k] + div_cnt_hist[k]) / err_hist[k]));
						//fprintf(frc, "Mul Hist: %ld %f\n", k, (float)(mul_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "Div Hist: %ld %f\n", k, (float)(div_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "RCB Hist: %ld %f\n", k, (float)(real_cbm_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "RMF Hist: %ld %f\n", k, (float)(real_mul_ff_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "Pow Hist: %ld %f\n", k, (float)(pow_cnt_hist[k] / err_hist[k]));
						fclose(frc);
						frc = NULL;
#endif						
					}
				}
#endif				

#if (1 == CFG_IMD_STORE)
				DEBUG_INFO("intp_cnt: %ld\n", intp_cnt);
#endif				

			}

/*early termination for simulation, reduce iteration times*/
/*more than 10 errors are found, and 10% simulation times have been excuted*/
#if (1 == EARLY_TERMINATION)
			if((EARLY_TERMINATION_NUM <= (frame_err - hamm_err))
				&& ((iter_cnt / 100000) < iter))
			{
				/*simulation times are enough, go to next Eb/N0 point*/
				break;
			}
#endif

#if (1 == DEBUG_LOG)		
			fclose(frc_debug);
			frc_debug = NULL;
#endif

		}

		stop = clock();
		runtime = (stop - start) / 1000.0000;

		/*output final results*/
		DEBUG_SYS("*********************************\n");
		DEBUG_SYS("Time: %fs\n", runtime);
		DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
		DEBUG_SYS("Frame: %ld\n", iter_cnt + 1);
		DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
		DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
		DEBUG_SYS("Frame Error: %ld\n", frame_err);
		DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
		DEBUG_SYS("Bit Error: %ld\n", bit_err);
		//DEBUG_SYS("Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)		
		DEBUG_SYS("RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)		
		DEBUG_SYS("Add/Mul Cnt: %f %f\n", (double)add_cnt / (double)(iter + 1), (double)(mul_cnt + div_cnt) / (double)(iter + 1));
		//DEBUG_SYS("Mul Cnt: %ld\n", mul_cnt);
		//DEBUG_SYS("Div Cnt: %ld\n", div_cnt);
		//DEBUG_SYS("RCB Cnt: %ld\n", real_cbm_cnt);
		//DEBUG_SYS("RMF Cnt: %ld\n", real_mul_ff_cnt);
		//DEBUG_SYS("Pow Cnt: %ld\n", pow_cnt);
#endif

#if (1 == CFG_PARTIALLY_PARALLEL)
		for(k = 0; k < PARALLEL_BATCH_NUM; k++)
		{
			//DEBUG_SYS("batch_add: %ld %ld\n", k, gf_add_in_batch[k]);
			//DEBUG_SYS("batch_mul: %ld %ld\n", k, gf_mul_in_batch[k]);
			DEBUG_INFO("batch_add: %ld %f %ld %d\n",
			          k,
			          (double)gf_add_in_batch[k] / (double)(iter + 1),
			          gf_add_in_batch_best[k],
			          gf_add_in_batch_worst[k]);
			DEBUG_INFO("batch_mul: %ld %f %ld %ld\n",
			          k,
			          (double)gf_mul_in_batch[k] / (double)(iter + 1),
			          gf_mul_in_batch_best[k],
			          gf_mul_in_batch_worst[k]);
		}
		
		batch_size = pow_val / PARALLEL_BATCH_NUM;
		latency_add = 0, latency_mul = 0;
		for(k = 0; k < batch_size; k++)
		{
			DEBUG_INFO("round_add: %ld %f %f\n",
			          k,
			          (double)cmplx_per_round_add[k] / (double)(iter + 1),
			          (double)latency_per_round_add[k] / (double)(iter + 1));
			DEBUG_INFO("round_mul: %ld %f %f\n",
			          k,
			          (double)cmplx_per_round_mul[k] / (double)(iter + 1),
			          (double)latency_per_round_mul[k] / (double)(iter + 1));
			          
			latency_add = latency_add + (double)latency_per_round_add[k] / (double)(iter + 1);
			latency_mul = latency_mul + (double)latency_per_round_mul[k] / (double)(iter + 1);
		}
		DEBUG_SYS("Latency: %f %f\n", latency_add, latency_mul);
#endif		

		for(i = 0; i < pow_val; i++)
		{
			if(0 != skip_hist[i])
			{
				DEBUG_SYS("skip_hist: %ld %ld\n", i, skip_hist[i]);
			}
		}
		for(i = 0; i < pow_val; i++)
		{
			if(0 != pgd_hist[i])
			{
				DEBUG_SYS("pgd_hist: %ld %ld\n", i, pgd_hist[i]);
			}
		}
		avg_round = 0;
#if (1 == CFG_ADAPTIVE_PARALLEL)
		for(i = 0; i < (pow_val + 1); i++)
		{
#else
		for(i = 0; i < (pow_val / PARALLEL_BATCH_NUM + 1); i++)
		{
#endif		
			if(0 != round_hist[i])
			{
				DEBUG_SYS("round_hist: %ld %ld\n", i, round_hist[i]);
			}
			avg_round = avg_round + (float)(round_hist[i] * (i));
		}
	
		DEBUG_SYS("avg_round: %f\n", avg_round / (float)(iter + 1));
		DEBUG_SYS("time_measure: %f\n", time_measure / (float)(iter + 1));

		DEBUG_SYS("Uncoded Results: %.10lf %.10lf %.10lf\n", 
			    (double)uncoded_frame_err / (double)(iter + 1),
			    (double)uncoded_symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
			    (double)uncoded_bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q);
		DEBUG_SYS("Decoded Results: %.10lf %.10lf %.10lf %.10lf\n", 
				(double)frame_err / (double)(iter + 1),
				(double)symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
				(double)bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q,
				(double)(frame_err - hamm_err) / (double)(iter + 1));

#if (1 == OUTPUT_LOG)
		frc = fopen(log_name, "a+");
		fprintf(frc, "*********************************\n");
		fprintf(frc, "Time: %fs\n", runtime);
		fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
		fprintf(frc, "Frame: %ld\n", iter_cnt + 1);
		fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
		fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
		fprintf(frc, "Frame Error: %ld\n", frame_err);
		fprintf(frc, "Symbol Error: %ld\n", symbol_err);
		fprintf(frc, "Bit Error: %ld\n", bit_err);
		//fprintf(frc, "Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)		
		fprintf(frc, "RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)
		fprintf(frc, "Add/Mul Cnt: %f %f\n", (double)add_cnt / (double)(iter + 1), (double)(mul_cnt + div_cnt) / (double)(iter + 1));
		//fprintf(frc, "Add Cnt: %ld\n", add_cnt);
		//fprintf(frc, "Mul Cnt: %ld\n", mul_cnt);
		//fprintf(frc, "Div Cnt: %ld\n", div_cnt);
		//fprintf(frc, "RCB Cnt: %ld\n", real_cbm_cnt);
		//fprintf(frc, "RMF Cnt: %ld\n", real_mul_ff_cnt);
		//fprintf(frc, "Pow Cnt: %ld\n", pow_cnt);
#endif		

		fprintf(frc, "Latency: %f %f\n", latency_add, latency_mul);
		for(i = 0; i < pow_val; i++)
		{
			fprintf(frc, "skip_hist: %ld %ld\n", i, skip_hist[i]);
		}
		for(i = 0; i < pow_val; i++)
		{
			fprintf(frc, "pgd_hist: %ld %ld\n", i, pgd_hist[i]);
		}
#if (1 == CFG_ADAPTIVE_PARALLEL)
		for(i = 0; i < (pow_val + 1); i++)
#else
		for(i = 0; i < (pow_val / PARALLEL_BATCH_NUM + 1); i++)
#endif		
		{
			fprintf(frc, "round_hist: %ld %ld\n", i, round_hist[i]);
		}
		fprintf(frc, "avg_round: %f\n", avg_round / (float)(iter + 1));
		fprintf(frc, "time_measure: %f\n", time_measure / (float)(iter + 1));
		fprintf(frc, "Uncoded Results: %.10lf %.10lf %.10lf\n", 
			    (double)uncoded_frame_err / (double)(iter + 1),
			    (double)uncoded_symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
			    (double)uncoded_bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q);
		fprintf(frc, "Decoded Results: %.10lf %.10lf %.10lf %.10lf\n", 
				(double)frame_err / (double)(iter + 1),
				(double)symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
				(double)bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q,
				(double)(frame_err - hamm_err) / (double)(iter + 1));
	    fclose(frc);
		frc = NULL;
#endif

	}

#if 0
	for (i = 0; i < symbol_num; i++)
	{
  		free(mod_seq[i]);
		mod_seq[i] = NULL;
  	}
	free(mod_seq);
	mod_seq = NULL;
#endif

	/*clear and exit*/
	mod_exit();

#if (0 == TEST_MODE)
	g_term_destroy();
#endif

#if (1 == CFG_PARTIALLY_PARALLEL)
	gf_partially_parallel_exit();
#endif

#if (1 == DEBUG_LOG)
	fclose(frc_debug);
	free(frc_debug);
	frc_debug = NULL;
#endif

	return;
}
