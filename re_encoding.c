#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "rnd.h"
#include "channel.h"
#include "mod.h"
#include "as_decoding.h"
#include "re_encoding.h"

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_matrix_tmp[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_order[CODEWORD_LEN];
long long chnl_rel_order_idx[CODEWORD_LEN];
long long chnl_rel_max_id[CODEWORD_LEN];
long long chnl_rel_scd_id[CODEWORD_LEN];

unsigned char mul_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
unsigned char beta_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
unsigned char rel_group_seq[MESSAGE_LEN];
unsigned char unrel_group_seq[CODEWORD_LEN - MESSAGE_LEN];
unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
unsigned char tao[CODEWORD_LEN - MESSAGE_LEN + 1];
unsigned char omega[CODEWORD_LEN - MESSAGE_LEN];
unsigned char sigma[((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)];
unsigned char erasure_polynomial[CODEWORD_LEN];
unsigned char phi[CODEWORD_LEN];
unsigned char v[MESSAGE_LEN + 1];
unsigned char g_v_val[CODEWORD_LEN];
#if (1 == RE_ENCODING)
unsigned char re_encoded_codeword[CODEWORD_LEN];
#endif
unsigned char H_msg[MESSAGE_LEN + 1];

void find_max_val(float matrix[][CODEWORD_LEN], long long col,
					 unsigned char* m_ptr, unsigned char* n_ptr)
{
	long long i = 0, j = 0;
	float max_val = 0;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(matrix[i][j] >= max_val)
			{
				max_val = matrix[i][j];
				*m_ptr = i;
				*n_ptr = j;
			}
		}
	}

	//DEBUG_INFO("i: %d, j: %d, val: %f\n", *m_ptr, *n_ptr, max_val);

	return;
}

int chnl_rel_seq_order()
{
	long long i = 0, j = 0;
	float tmp = 0, max_val = 0, scd_val = 0;
	long long max_idx = 0, scd_idx = 0;
	float chnl_rel[CODEWORD_LEN];

	memset(chnl_rel_order, 0, sizeof(float) * CODEWORD_LEN);
	memset(chnl_rel, 0, sizeof(float) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		chnl_rel_order_idx[i] = i;
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		max_val = 0;
		scd_val = 0;
	
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if(max_val < chnl_rel_matrix[j][i])
			{
				max_val = chnl_rel_matrix[j][i];
				max_idx = j;
			}

			chnl_rel_max_id[i] = max_idx;
		}
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if((scd_val < chnl_rel_matrix[j][i])
				&& (max_val > chnl_rel_matrix[j][i]))
			{
				scd_val = chnl_rel_matrix[j][i];
				scd_idx = j;
			}

			chnl_rel_scd_id[i] = scd_idx;
		}

		if(0 == scd_val)
		{
			scd_val = 0.00001;
		}
		chnl_rel[i] = max_val / scd_val;
		DEBUG_NOTICE("chnl_rel_order: %f %f %f %ld %ld\n",
		             max_val,
		             scd_val,
		             chnl_rel[i],
		             chnl_rel_max_id[i],
		             chnl_rel_scd_id[i]);
	}

	memcpy(chnl_rel_order, chnl_rel, sizeof(float) * CODEWORD_LEN);

	BubbleSort4(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}

	return 0;
}

int chnl_rel_cal(float **input_seq,
				    long long input_len)
{
	long long i = 0, j = 0, k = 0;
	float n0, d0 = 0, d1 = 0, temp = 0;
	unsigned char tmp_bit = 0;
	/*for BPSK*/
	float map[input_len][2];

	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		memset(chnl_rel_matrix[i], 0, sizeof(float) * CODEWORD_LEN);
	}


	for(i = 0; i < CODEWORD_LEN + 1; i++)
	{
		memcpy(chnl_rel_matrix_tmp[i], chnl_rel_matrix[i], sizeof(float) * CODEWORD_LEN);
	}

#if 0
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < (GF_Q / BITS_PER_SYMBOL_BPSK); j++)
		{
			DEBUG_NOTICE("recv_seq: %f\n", fabs(recv_seq[i * (GF_Q / BITS_PER_SYMBOL_BPSK) + j][0]));
			recv_rel[i] = recv_rel[i] + fabs(recv_seq[i * (GF_Q / BITS_PER_SYMBOL_BPSK) + j][0]);
		}
		DEBUG_NOTICE("recv_rel: %f\n", recv_rel[i]);
	}
#endif

	n0 = 1 / ((float)MESSAGE_LEN / (float)CODEWORD_LEN) / (pow(10, eb2n0 / 10) / 2);

	for(i = 0; i < input_len; i++)
	{
		d0 = (input_seq[i][0] - (1.0)) * (input_seq[i][0] - (1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		d1 = (input_seq[i][0] - (-1.0)) * (input_seq[i][0] - (-1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));

#if 1
		d0 = d0 / 1e4;
		d1 = d1 / 1e4;
#endif

		/*for BPSK*/
		map[i][0] = 1 / (PI * n0) * exp((-d0) / n0);
        map[i][1] = 1 / (PI * n0) * exp((-d1) / n0);
        DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        		   i,
        		   n0,
        		   d0,
        		   d1,
        		   map[i][0],
        		   map[i][1]);
        if((0 == map[i][0])
        	|| (0 == map[i][1]))
        {
        	DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        			  i,
        			  n0,
        			  map[i][0],
        			  map[i][0],
        			  d0,
        			  d1);
        }
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			chnl_rel_matrix[i][j] = 1;

			for(k = 0; k < GF_Q; k++)
			{
				/*for BPSK*/
				tmp_bit = ((power_polynomial_table[i][1] >> k) & 0x1);
				DEBUG_INFO("chnl_rel_matrix: %d %d | %d | %d | %f %f\n",
		        		   i,
		        		   j,
		        		   k,
		        		   tmp_bit,
		        		   chnl_rel_matrix[i][j],
		        		   map[j * GF_Q + k][tmp_bit]);
				chnl_rel_matrix[i][j] = chnl_rel_matrix[i][j]
									   * map[j * GF_Q + k][tmp_bit];
			}
		}
	}

	for (i = 0; i < CODEWORD_LEN; i++)
    {
		DEBUG_NOTICE("recv: %x\n", received_polynomial[i]);

        temp = 0;
        for (j = 0; j < GF_FIELD; j++)
        {
            temp = temp + chnl_rel_matrix[j][i];
        }
        
        /*special process*/
#if 0        
        if(0 == temp)
        {
        	if(0xFF == received_polynomial[i])
        	{
        		chnl_rel_matrix[0][i] = 1;
        	}
        	else
        	{
        		chnl_rel_matrix[received_polynomial[i] + 1][i] = 1;
        	}
        	continue;
        }
#endif        

        for (j = 0; j < GF_FIELD; j++)
        {
            chnl_rel_matrix[j][i] = chnl_rel_matrix[j][i] / temp;
            if(chnl_rel_matrix[j][i] != chnl_rel_matrix[j][i])
            {
            	DEBUG_SYS("NAN: %d %d | %f %f | %f %f %f %f %f %f %f %f\n",
            			   i,
            	           j,
            	           temp,
            	           chnl_rel_matrix[j][i],
            	           input_seq[i * GF_Q][0],
            	           input_seq[i * GF_Q + 1][0],
            	           input_seq[i * GF_Q + 2][0],
            	           input_seq[i * GF_Q + 3][0],
            	           input_seq[i * GF_Q + 4][0],
            	           input_seq[i * GF_Q + 5][0],
            	           input_seq[i * GF_Q + 6][0],
            	           input_seq[i * GF_Q + 7][0]);
            }
			DEBUG_NOTICE("chnl_rel: %ld %ld | %f\n", i, power_polynomial_table[j][0], chnl_rel_matrix[j][i]);
        }
    }

	chnl_rel_seq_order();

	return 0;
}

int mul_assign()
{
	long long i = 0, j = 0;
	
#if 0	
	long long s = 0;
	unsigned char *m_ptr = (unsigned char*)malloc(sizeof(unsigned char));
	unsigned char *n_ptr = (unsigned char*)malloc(sizeof(unsigned char));

	for(i = 0; i < CODEWORD_LEN + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			mul_matrix[i][j] = 0;
		}
	}

	for(s = 0; s < S_MUL; s++)
	{
		find_max_val(chnl_rel_matrix_tmp, CODEWORD_LEN + 1,
					  m_ptr, n_ptr);
		i = *m_ptr;
		j = *n_ptr;
		chnl_rel_matrix_tmp[i][j] = chnl_rel_matrix[i][j] / (mul_matrix[i][j] + 2);
		mul_matrix[i][j] = mul_matrix[i][j] + 1;
	}

	DEBUG_IMPOTANT("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			DEBUG_IMPOTANT("%d ", mul_matrix[j][i]);
		}
		DEBUG_IMPOTANT("\n");
	}

	free(m_ptr);
	free(n_ptr);
#else//GS algorithm
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			mul_matrix[j][i] = 0;

			if(power_polynomial_table[j][0] == received_polynomial[i])
			{
				mul_matrix[j][i] = S_MUL;

#if (1 == TEST_MODE)
				if(received_polynomial[i] != encoded_polynomial[i])
				{
					//mul_matrix[j][i] = S_MUL / 2;
				}
#endif
#if 0//for ASD test
				if(encoded_polynomial[i] != received_polynomial[i])
				{
					mul_matrix[j][i] = 1;
				}
#endif
#if (1 == SIMPLE_ASD)//for ASD test, a simple multiplicity assignment strategy
				if(GF_Q > recv_rel[i])
				{
					mul_matrix[j][i] = S_MUL / 2;
				}
				if((GF_Q / 2) > recv_rel[i])
				{
					mul_matrix[j][i] = S_MUL / 4;
				}
#endif

			}
		}
	}

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			DEBUG_IMPOTANT("%d ", mul_matrix[j][i]);
		}
		DEBUG_IMPOTANT("\n");
	}
#endif	
#endif
	
	return 0;
}

unsigned char syndrome_cal(unsigned char *recv, unsigned char *synd,
								long long cw_len, long long msg_len)
{
	long long i = 0, j = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;

	for(i = 0; i < (cw_len - msg_len); i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < cw_len; j++)
		{
			if(0xFF != recv[j])
			{
				tmp = gf_multp(recv[j], ((i + 1) * j) % (GF_FIELD - 1));
				tmp_sum = gf_add(tmp, tmp_sum);
			}
#if 0			
			if(0xFF == tmp_sum)
			{
				DEBUG_INFO("s_tmp: %d %d\n", recv[j], power_polynomial_table[0][1]);
			}
			else
			{
				DEBUG_INFO("s_tmp: %d %d %d %d %d %d\n",
						   i,
						   j,
						   recv[j],
						   power_polynomial_table[tmp + 1][1],
					       power_polynomial_table[recv[j] + 1][1],
					       power_polynomial_table[tmp_sum + 1][1]);
			}
#endif			
		}
		synd[i] = tmp_sum;
#if 0		
		if(0xFF == synd[i])
		{
			DEBUG_INFO("%d %d\n", i, power_polynomial_table[0][1]);
		}
		else
		{
			DEBUG_INFO("%d %d\n", i, power_polynomial_table[synd[i] + 1][1]);
		}
#endif		
	}
#if (1 == CFG_DEBUG_INFO)	
	DEBUG_INFO("Syndrome Polynomial:\n");
	for(i = 0; i < cw_len - msg_len; i++)
	{
		if(0xFF == synd[i])
		{
			DEBUG_INFO("%d %d\n", i, power_polynomial_table[0][0]);
		}
		else
		{
			DEBUG_INFO("%d %d\n", i, power_polynomial_table[synd[i] + 1][0]);
		}
	}
	DEBUG_INFO("\n");
#endif	
}

int rel_group()
{
	long long i = 0, j = 0;
#if 0
	float rel_thrd = 1.0;
	long long rel_cnt = 0, rel_flag = 0;

	while(MESSAGE_LEN > rel_cnt)
	{
		rel_flag = 0;
		memset(rel_group_seq, 0, sizeof(unsigned char) * MESSAGE_LEN);
		memset(unrel_group_seq, 0, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			rel_flag = 0;
			for(i = 0; i < (CODEWORD_LEN + 1); i++)
			{
				if(rel_thrd < chnl_rel_matrix[i][j])
				{
					rel_flag = 1;
					DEBUG_INFO("rel_val: %d %d %f %f %d\n", i, j, rel_thrd, chnl_rel_matrix[i][j], rel_cnt);
					break;
				}
			}

			if(1 == rel_flag)
			{
				if(MESSAGE_LEN > rel_cnt)
				{
					rel_group_seq[rel_cnt] = j;
				}
				rel_cnt = rel_cnt + 1;
			}
		}

		if(MESSAGE_LEN < rel_cnt)
		{
			rel_thrd = rel_thrd + 0.01;
		}
		else
		{
			rel_thrd = rel_thrd - 0.01;
		}
	}

	rel_cnt = 0;
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		rel_flag = 0;
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			if(j == rel_group_seq[i])
			{
				rel_flag = 1;
				break;
			}
		}
		if(0 == rel_flag)
		{
			unrel_group_seq[rel_cnt] = j;
			rel_cnt = rel_cnt + 1;
		}
	}
#endif
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		rel_group_seq[i] = chnl_rel_order_idx[CODEWORD_LEN - 1 - i];
	}
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		unrel_group_seq[i] = chnl_rel_order_idx[(CODEWORD_LEN - MESSAGE_LEN) - 1 - i];
	}

#if 0//(1 == TEST_MODE)//just for GF(8) test
	rel_group_seq[0] = 4;
	rel_group_seq[1] = 2;
	rel_group_seq[2] = 3;
	rel_group_seq[3] = 5;
	rel_group_seq[4] = 6;
	unrel_group_seq[0] = 1;
	unrel_group_seq[1] = 0;
	chnl_rel_order_idx[0] = unrel_group_seq[1];
	chnl_rel_order_idx[1] = unrel_group_seq[0];
	chnl_rel_order_idx[2] = rel_group_seq[6];
	chnl_rel_order_idx[3] = rel_group_seq[5];
	chnl_rel_order_idx[4] = rel_group_seq[3];
	chnl_rel_order_idx[5] = rel_group_seq[2];
	chnl_rel_order_idx[6] = rel_group_seq[4];
	//unrel_group_seq[2] = 4;
	//unrel_group_seq[3] = 5;
	//unrel_group_seq[4] = 6;
#endif

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("rel_seq: ");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_INFO("%d ", rel_group_seq[i]);
	}
	DEBUG_INFO("\n");
	DEBUG_INFO("unrel_seq: ");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		DEBUG_INFO("%d ", unrel_group_seq[i]);
	}
	DEBUG_INFO("\n");
#endif

#if 0
	DEBUG_SYS("rel_seq: ");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_SYS("%d ", rel_group_seq[i]);
	}
	DEBUG_SYS("\n");
#endif

#if 0
	/*check errors in reliable group*/
	long long err_in_rel = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			if((i == rel_group_seq[j])
				&& (received_polynomial[i] != encoded_polynomial[i]))
			{
				err_in_rel++;
				break;
			}
		}
	}
	DEBUG_SYS("err_in_rel: %ld\n", err_in_rel);
#endif

	return 0;
}

int v_cal()
{
	long long i = 0, j = 0;
	memset(v, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));
	v[0] = rel_group_seq[0];
	v[1] = 0;
	unsigned char v_tmp_1[MESSAGE_LEN + 1], v_tmp_2[MESSAGE_LEN + 1];

	for(i = 0; i < (MESSAGE_LEN - 1); i++)
	{
		memset(v_tmp_1, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));
		memset(v_tmp_2, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));

		for(j = 0; j < (MESSAGE_LEN + 1); j++)
		{
			if(0 != j)//increase degree because of x term
			{
				v_tmp_1[j] = v[j - 1];
				DEBUG_NOTICE("v_tmp_1: %d | %x\n", j, v_tmp_1[j]);
			}

			if(0xFF != v[j])
			{
				v_tmp_2[j] = gf_multp(v[j], rel_group_seq[i + 1]);//calculation of a_i term
				DEBUG_NOTICE("v_tmp_2: %d | %x=%x*%x\n", j, v_tmp_2[j], v[j], rel_group_seq[i]);
			}
		}

		for(j = 0; j < (MESSAGE_LEN + 1); j++)
		{
			v[j] = gf_add(v_tmp_1[j], v_tmp_2[j]);//add 2 parts
			DEBUG_NOTICE("v_tmp: %d | %x\n", j, v[j]);
		}
	}

	for(i = 0; i < (MESSAGE_LEN + 1); i++)
	{
		DEBUG_NOTICE("v: %d | %x\n", i, v[i]);
	}
	
	return 0;
}

int g_v_val_cal()
{
	long long i = 0, j = 0;
	unsigned char skip_v_flag = 0;
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		skip_v_flag = 0;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			if(i == rel_group_seq[j])
			{
				skip_v_flag = 1;
				break;
			}
		}
		if(1 == skip_v_flag)
		{
			g_v_val[i] = 0xFF;
			continue;
		}
	
		g_v_val[i] = poly_eva(v,
				 			 (MESSAGE_LEN + 1),
				 			 power_polynomial_table[i + 1][0]);
		DEBUG_NOTICE("g_v_val: %d | %x\n",
			         i,
			         g_v_val[i]);		 			 
	}
	
	return 0;
}

int l_cal(unsigned char locator_j, unsigned char *L)
{
	long long i = 0, j = 0;
	unsigned char v_tmp_1[MESSAGE_LEN + 1], v_tmp_2[MESSAGE_LEN + 1];
	unsigned char tmp_div = 0, tmp_inv = 0xFF;

	memset(L, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	L[0] = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(rel_group_seq[i] == locator_j)
		{
			continue;
		}
		
		memset(v_tmp_1, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));
		memset(v_tmp_2, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));

		for(j = 0; j < (MESSAGE_LEN + 0); j++)
		{
			v_tmp_1[j + 1] = L[j];//increase degree because of x term
			DEBUG_NOTICE("v_tmp_1: %d | %x\n", j + 1, v_tmp_1[j + 1]);

			if(0x0 == L[j])
			{
				v_tmp_2[j] = rel_group_seq[i];
				continue;
			}
			if(0xFF != L[j])
			{
				v_tmp_2[j] = gf_multp(L[j], rel_group_seq[i]);//calculation of a_i term
			}
			DEBUG_NOTICE("v_tmp_2: %d | %x=%x*%x\n", j, v_tmp_2[j], L[j], rel_group_seq[i]);
		}

		for(j = 0; j < (MESSAGE_LEN + 1); j++)
		{
			if(0xFF == v_tmp_1[j])
			{
				L[j]= v_tmp_2[j];
				continue;
			}
			if(0xFF == v_tmp_2[j])
			{
				L[j]= v_tmp_1[j];
				continue;
			}
			L[j] = gf_add(v_tmp_1[j], v_tmp_2[j]);//add 2 parts
			DEBUG_NOTICE("v_tmp: %d | %x\n", j, L[j]);
		}
	}

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(rel_group_seq[i] == locator_j)
		{
			continue;
		}

		tmp_inv = gf_div(0x0, gf_add(locator_j, rel_group_seq[i]));
		DEBUG_NOTICE("tmp_inv: %d | %x\n", i, tmp_inv);
		tmp_div = gf_multp(tmp_div, tmp_inv);
		DEBUG_NOTICE("tmp_div: %d | %x\n", i, tmp_div);
	}

	for(i = 0; i < (MESSAGE_LEN + 1); i++)
	{
		if(0x0 == L[i])
		{
			L[i] = tmp_div;
			continue;
		}
		if(0x0 == tmp_div)
		{
			continue;
		}
	
		L[i] = gf_multp(L[i], tmp_div);
		DEBUG_NOTICE("L: %d | %x\n", i, L[i]);
	}
	
	return 0;
}

int tao_cal(unsigned char *erasure_group)
{
	long long i = 0, j = 0;
	memset(tao, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	unsigned char a[2] = {0, 0};
	unsigned char reg[CODEWORD_LEN - MESSAGE_LEN + 1];/*n-k+1, R' has (n-k) terms*/
	memset(reg, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
#if 0
	reg[0] = unrel_group_seq[0];
	reg[1] = 0;
#endif
	unsigned char b[CODEWORD_LEN - MESSAGE_LEN];
	memset(b, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
#if 0
	for(i = 1; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		a[0] = unrel_group_seq[i];
		a[1] = 0;
		memcpy(b, reg, sizeof(unsigned char) * (1 + i));
		gf_multp_poly_hw(a, 2,
						  b, (1 + i),
						  reg, (2 + i));
	}
#else
	reg[0] = 0;
	reg[1] = erasure_group[0];
	for(i = 1; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		a[0] = 0;
		a[1] = erasure_group[i];
		memcpy(b, reg, sizeof(unsigned char) * (1 + i));
		gf_multp_poly_hw(a, 2,
						  b, (1 + i),
						  reg, (2 + i));
	}
#endif
	memcpy(tao, reg, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("tao: ");
	for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN + 1); j++)
	{
		DEBUG_INFO("%d ", reg[j]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int sigma_cal()
{
	long long i = 0, j = 0;
	memset(sigma, 0xFF, sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

	unsigned char tmp[(CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1];
	memset(tmp, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	gf_multp_poly_hw(syndrome, (CODEWORD_LEN - MESSAGE_LEN),
					  tao, 		  (CODEWORD_LEN - MESSAGE_LEN + 1),
					  tmp, 		  ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	memcpy(sigma, tmp + (CODEWORD_LEN - MESSAGE_LEN), sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

#if (1 == CFG_DEBUG_NOTICE)
	DEBUG_NOTICE("sigma: ");
	for(i = 0; i < (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)); i++)
	{
		DEBUG_NOTICE("%d ", sigma[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	memcpy(omega, tmp, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

#if (1 == CFG_DEBUG_NOTICE)
	DEBUG_NOTICE("omega: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%d ", omega[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	return 0;
}

unsigned char poly_eva(unsigned char *poly, long long poly_len, unsigned char input_val)
{
	long long i = 0;
	unsigned char poly_val = 0xFF, tmp_product = 0;

	for(i = 0; i < poly_len; i++)
	{
		if((0xFF == (*(poly + i)))
			|| (0xFF == input_val))
		{
			continue;
		}
	
		tmp_product = gf_multp(*(poly + i), (input_val * i) % (GF_FIELD - 1));
		//DEBUG_NOTICE("poly_eva: %d %x %x %x %x %x\n", i, *(poly + i), input_val, (input_val * i) % (GF_FIELD - 1), tmp_product, poly_val);
		poly_val = gf_add(tmp_product, poly_val);
	}
	
	return poly_val;
}

int phi_cal()
{
	long long i = 0, k = 0;
	unsigned char find_flag = 0;
	unsigned char tao_dev[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1]; 
	unsigned char tmp = 0xFF, locator = 0xFF;

	memset(phi, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(erasure_polynomial, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(unrel_group_seq[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			erasure_polynomial[i] = 0xFF;
		}
		if(0xFF == erasure_polynomial[i])
		{
			DEBUG_NOTICE("erasure_polynomial: %d %x %x %d\n", i, received_polynomial[i], erasure_polynomial[i], power_polynomial_table[0][1]);
		}
		else
		{
			DEBUG_NOTICE("erasure_polynomial: %d %x %x %d\n", i, received_polynomial[i], erasure_polynomial[i], power_polynomial_table[erasure_polynomial[i] + 1][1]);
		}
	}
	find_flag = 0;

	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	syndrome_cal(erasure_polynomial,
				  syndrome,
				  CODEWORD_LEN,
				  MESSAGE_LEN);
	tao_cal(unrel_group_seq);
	sigma_cal();

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		if(0 != (i % 2))
		{
			tao_dev[i] = 0xFF;
		}
		else
		{
			tao_dev[i] = tao[i + 1];
		}
	}

#if (1 == CFG_DEBUG_NOTICE)	
	DEBUG_NOTICE("tao_dev: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%d ", tao_dev[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(unrel_group_seq[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			locator = power_polynomial_table[i + 1][0];
			locator = ((GF_FIELD - 1) - locator) % (GF_FIELD - 1);
			//DEBUG_NOTICE("%x: %x\n", i, locator);
			tmp = poly_eva(omega, CODEWORD_LEN - MESSAGE_LEN, locator);
			//DEBUG_NOTICE("%x ", tmp);
			//tmp = gf_div(tmp, (locator * (CODEWORD_LEN - MESSAGE_LEN + 1)) % (GF_FIELD - 1));
			tmp = gf_div(tmp, poly_eva(tao_dev, CODEWORD_LEN - MESSAGE_LEN, locator));
			//DEBUG_NOTICE("%x ", tmp);
			//phi[i] = gf_add(tmp, received_polynomial[i]);
			phi[i] = tmp;
			//DEBUG_NOTICE("\n");
		}
		else
		{
			phi[i] = received_polynomial[i];
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("phi:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0xFF == phi[i])
		{
			DEBUG_INFO("%x  %x  %x  %x  %d\n",
						encoded_polynomial[i],
						received_polynomial[i],
						phi[i],
						gf_add(received_polynomial[i], phi[i]),
						power_polynomial_table[0][1]);
		}
		else
		{
			DEBUG_INFO("%x  %x  %x  %x  %d\n",
						encoded_polynomial[i],
						received_polynomial[i],
						phi[i],
						gf_add(received_polynomial[i], phi[i]),
						power_polynomial_table[phi[i] + 1][1]);
		}
	}
	DEBUG_INFO("\n");
#if 0
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(phi, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif	
#endif

	return 0;
}

unsigned char coordinate_trans(unsigned char locator, unsigned char r, unsigned char r_hd)
{
	long long i = 0, j = 0;
	unsigned char coordinate = 0xFF;
	unsigned char locator_product = 0, tmp_sum = 0xFF;
#if 0
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		if(locator == unrel_group_seq[i])
		{
			continue;
		}
		tmp_sum = gf_add(locator, unrel_group_seq[i]);
		locator_product = gf_multp(locator_product, tmp_sum);
	}
	//locator_product = gf_multp(0, locator);
	locator_product = gf_multp(locator_product, locator);
	tmp_sum = gf_add(r_hd, r);
	locator_product = gf_multp(tmp_sum, locator_product);

	tmp_sum = poly_eva(sigma, ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN), locator);

	coordinate = gf_add(locator_product, tmp_sum);
#else
	if(r_hd == r)
	{
		return 0xFF;
	}

#if 0
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		tmp_sum = gf_add(locator, rel_group_seq[i]);
		locator_product = gf_multp(locator_product, tmp_sum);
		DEBUG_NOTICE("v_cal: (%x+%x)*%x=%x\n",
					 locator,
					 rel_group_seq[i],
					 tmp_sum,
					 locator_product);
	}
#endif

	locator_product = poly_eva(v, MESSAGE_LEN + 1, locator);

	//locator_product = 0;
	coordinate = gf_div(gf_add(r_hd, r), locator_product);
	//coordinate = gf_add(r_hd, r);
#if 0	
	DEBUG_NOTICE("locator_product: %x | %x + %x = %x | %x | %x\n",
				 locator,
				 r_hd,
				 r,
				 gf_add(r_hd, r),
				 locator_product,
				 coordinate);
#endif				 
#endif
	return coordinate;
}

int erasure_decoding(unsigned char *r_seq, unsigned char *erasure_group)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char find_flag = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	//syndrome_cal(r_seq, syndrome, CODEWORD_LEN, MESSAGE_LEN);

	unsigned char tao_dev[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1]; 
	unsigned char tmp = 0xFF, locator = 0xFF;
	memset(phi, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(erasure_polynomial, r_seq, sizeof(unsigned char) * CODEWORD_LEN);
	//for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(erasure_group[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			erasure_polynomial[i] = 0xFF;
			DEBUG_NOTICE("erasure_polynomial: %d %x %x\n", i, r_seq[i], erasure_polynomial[i]);
		}
	}
	find_flag = 0;
	
	syndrome_cal(erasure_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
	tao_cal(erasure_group);
	sigma_cal();

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		if(0 != (i % 2))
		{
			tao_dev[i] = 0xFF;
		}
		else
		{
			tao_dev[i] = tao[i + 1];
		}
	}

#if (1 == CFG_DEBUG_NOTICE)	
	DEBUG_NOTICE("tao_dev: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%x ", tao_dev[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(erasure_group[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			locator = power_polynomial_table[i + 1][0];
			locator = ((GF_FIELD - 1) - locator) % (GF_FIELD - 1);
			DEBUG_NOTICE("%x: %x\n", i, locator);
			tmp = poly_eva(omega, CODEWORD_LEN - MESSAGE_LEN, locator);
			DEBUG_NOTICE("%x ", tmp);
			//tmp = gf_div(tmp, (locator * (CODEWORD_LEN - MESSAGE_LEN + 1)) % (GF_FIELD - 1));
			DEBUG_NOTICE("%x ", tmp);
			tmp = gf_div(tmp, poly_eva(tao_dev, CODEWORD_LEN - MESSAGE_LEN, locator));
			DEBUG_NOTICE("%x ", tmp);
			//phi[i] = gf_add(tmp, received_polynomial[i]);
			phi[i] = tmp;
			DEBUG_NOTICE("\n");
		}
		else
		{
			DEBUG_NOTICE("cal: %d | %x\n", i, r_seq[i]);
			phi[i] = r_seq[i];
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("recover_erasure_codeword: ");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_INFO("%x ", phi[i]);
	}
	DEBUG_INFO("\n");
#if 0
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(phi, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif	
#endif

	return 0;
}

int erasure_decoding_use(unsigned char *r_seq, unsigned char *erasure_group)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char find_flag = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	//syndrome_cal(r_seq, syndrome, CODEWORD_LEN, MESSAGE_LEN);

	unsigned char tao_dev[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1]; 
	unsigned char tmp = 0xFF, locator = 0xFF;
	//memset(phi, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(erasure_polynomial, r_seq, sizeof(unsigned char) * CODEWORD_LEN);
	//for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(erasure_group[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			erasure_polynomial[i] = 0xFF;
			DEBUG_NOTICE("erasure_polynomial: %d %x %x\n", i, r_seq[i], erasure_polynomial[i]);
		}
	}
	find_flag = 0;
	
	syndrome_cal(erasure_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
	tao_cal(erasure_group);
	sigma_cal();

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		if(0 != (i % 2))
		{
			tao_dev[i] = 0xFF;
		}
		else
		{
			tao_dev[i] = tao[i + 1];
		}
	}

#if (1 == CFG_DEBUG_NOTICE)	
	DEBUG_NOTICE("tao_dev: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%x ", tao_dev[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(erasure_group[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			locator = power_polynomial_table[i + 1][0];
			locator = ((GF_FIELD - 1) - locator) % (GF_FIELD - 1);
			DEBUG_NOTICE("%x: %x\n", i, locator);
			tmp = poly_eva(omega, CODEWORD_LEN - MESSAGE_LEN, locator);
			DEBUG_NOTICE("%x ", tmp);
			//tmp = gf_div(tmp, (locator * (CODEWORD_LEN - MESSAGE_LEN + 1)) % (GF_FIELD - 1));
			DEBUG_NOTICE("%x ", tmp);
			tmp = gf_div(tmp, poly_eva(tao_dev, CODEWORD_LEN - MESSAGE_LEN, locator));
			DEBUG_NOTICE("%x ", tmp);
			//phi[i] = gf_add(tmp, received_polynomial[i]);
			r_seq[i] = tmp;
			DEBUG_NOTICE("\n");
		}
		else
		{
			DEBUG_NOTICE("cal: %d | %x\n", i, r_seq[i]);
			r_seq[i] = r_seq[i];
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("recover_erasure_codeword: ");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_INFO("%x ", r_seq[i]);
	}
	DEBUG_INFO("\n");
#if 0
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(phi, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif	
#endif

	return 0;
}

int erasure_decoding_lag(unsigned char *r_seq)
{
	long long i = 0, j = 0;

	unsigned char L[MESSAGE_LEN][MESSAGE_LEN + 1];
	memset(H_msg, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		l_cal(rel_group_seq[i], L[i]);
	}
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < (MESSAGE_LEN + 1); j++)
		{
			if(0xFF != L[i][j])
			{
				L[i][j] = gf_multp(L[i][j], r_seq[rel_group_seq[i]]);
			}
		}
	}
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < (MESSAGE_LEN + 1); j++)
		{
			H_msg[j] = gf_add(H_msg[j], L[i][j]);
		}
	}
	for(i = 0; i < (MESSAGE_LEN + 1); i++)
	{
		DEBUG_NOTICE("H: %d | %x\n", i, H_msg[i]);
	}
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
#if 1	
		phi[i] = poly_eva(H_msg,
		                  (MESSAGE_LEN + 1),
		                  power_polynomial_table[i + 1][0]);
		DEBUG_INFO("%x ", phi[i]);
#endif		
	}
}

int syndrome_re_encoding(unsigned char *synd)
{
	long long i = 0, j = 0;
	unsigned char tmp_f = 0xFF, tmp_v = 0xFF;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_INFO("decoded_message: %ld %x\n", i, decoded_message[i]);
	}
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_INFO("rel_group_seq: %ld\n", rel_group_seq[i]);
	}

	for(i = 0; i < MESSAGE_LEN; i++)//locator
	{
		tmp_f = 0xFF;
		tmp_v = 0xFF;

		for(j = 0; j < MESSAGE_LEN; j++)//msg
		{
			tmp_f = gf_add(tmp_f,
						    gf_multp(decoded_message[j],
						    		gf_pow_cal(i, (((i + 1) * j) % (GF_FIELD - 1)))));

			//tmp_v = gf_multp(tmp_v, gf_add(rel_group_seq[i], rel_group_seq[j]));
			tmp_v = 0x0;
			DEBUG_INFO("tmp: %ld | %x %x %x %x %x %x\n",
				       j,
				       decoded_message[j],
				       ((i + 1) * j) % (GF_FIELD - 1),
				       tmp_f,
				       rel_group_seq[i],
				       rel_group_seq[j],
				       tmp_v);
		}

		//synd[i] = gf_div(tmp_f, tmp_v);
		synd[i] = decoded_message[i];
		DEBUG_INFO("synd: %ld | %x %x %x\n", i, tmp_f, tmp_v, synd[i]);
	}
	
	return 0;
}

#if (1 == RE_ENCODING)
int bm_re_encoding(unsigned char *msg_phi, unsigned char *tmp_cw)
{
	long long i = 0, j = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF, lambda_root = 0;
	unsigned char syndrome_tmp[SYN_LEN];
	unsigned char syndrome[SYN_LEN + 1];
	unsigned char lambda[MESSAGE_LEN];
	unsigned char omega[MESSAGE_LEN];
	unsigned char err_location[CODEWORD_LEN];
	unsigned char lambda_dev[MESSAGE_LEN];
	unsigned char v_tmp[MESSAGE_LEN + 1];
	unsigned char err_mag[CODEWORD_LEN];

	memset(syndrome_tmp, 0xFF, sizeof(unsigned char) * SYN_LEN);
	memset(syndrome, 0xFF, sizeof(unsigned char) * (SYN_LEN + 1));
	memset(lambda, 0x0, sizeof(unsigned char) * MESSAGE_LEN);
	memset(omega, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	memset(err_location, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memset(lambda_dev, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	memset(v_tmp, 0xFF, sizeof(unsigned char) * (MESSAGE_LEN + 1));
	memset(err_mag, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

	unsigned char sigma = 0xFF;
	unsigned char B[MESSAGE_LEN];
	unsigned char lambda_tmp[MESSAGE_LEN];
	unsigned char L = 0;

	memset(B, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	B[1] = 0;
	memset(lambda_tmp, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	lambda_tmp[0] = 0;
	memcpy(lambda, lambda_tmp, sizeof(unsigned char) * MESSAGE_LEN);

	DEBUG_NOTICE("-------------------------------------------------\n");
	DEBUG_NOTICE("BM Decoding:\n");
	DEBUG_NOTICE("-------------------------------------------------\n");

	/*compute the syndrome polynomial from received symbols*/
	//syndrome_re_encoding(syndrome_tmp);
#if 0
	syndrome_tmp[0] = 0x1;
	syndrome_tmp[1] = 0x3;
	syndrome_tmp[2] = 0x5;
	syndrome_tmp[3] = 0x0;
	syndrome_tmp[4] = 0x2;
#endif

	memcpy(syndrome_tmp, msg_phi, sizeof(unsigned char) * SYN_LEN);
	memcpy(syndrome + 1, syndrome_tmp, sizeof(unsigned char) * SYN_LEN);

	for(i = 1; i < (SYN_LEN + 1); i++)
	{
		tmp = 0xFF;
		for(j = 1; j <= L; j++)
		{
			DEBUG_NOTICE("tmp: %d %d %x %x %x\n", i, j, tmp, lambda[j], syndrome[i - j]);
			tmp = gf_add(tmp, gf_multp(syndrome[i - j], lambda[j]));
		}
		sigma = gf_add(syndrome[i], tmp);
		DEBUG_NOTICE("sigma: %x %x %x\n", sigma, syndrome[i], tmp);

		for(j = 0; j < SYN_LEN; j++)
		{
			tmp = gf_multp(B[j], sigma);
			lambda[j] = gf_add(lambda[j], tmp);
		}
		for(j = 0; j < SYN_LEN; j++)
		{
			DEBUG_NOTICE("lambda: %d %x\n", j, lambda[j]);
		}

		if(0xFF != sigma)
		{
			if(2 * L < (i))
			{
				L = (i) - L;
				for(j = 0; j < SYN_LEN; j++)
				{
					B[j] = gf_div(lambda_tmp[j], sigma);
					DEBUG_NOTICE("lambda_tmp: %x %x %x\n", lambda_tmp[j], sigma, B[j]);
				}
			}
		}

		for(j = (SYN_LEN - 1); j >= 1; j--)
		{
			//DEBUG_NOTICE("B: %d %x %x\n", j, B[j], B[j - 1]);
			B[j] = B[j - 1];
		}
		B[0] = 0xFF;

		for(j = 0; j < SYN_LEN; j++)
		{
			DEBUG_NOTICE("B: %d %x\n", j, B[j]);
		}

		memcpy(lambda_tmp, lambda, sizeof(unsigned char) * SYN_LEN);
	}

	for(i = 0; i < SYN_LEN; i++)
	{
		for(j = 0; j <= i; j++)
		{
			tmp = gf_multp(lambda[j], syndrome_tmp[i - j]);
			DEBUG_NOTICE("omega: %d %d | %x %x %x %x\n",
						 i,
						 j,
						 lambda[j],
						 syndrome_tmp[i - j],
						 tmp,
						 omega[i]);
			omega[i] = gf_add(omega[i], tmp);
		}
	}

	DEBUG_NOTICE("Omega Polynomial:\n");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("%x ", omega[i]);
	}
	DEBUG_NOTICE("\n");

	long long max_lambda_degree = MESSAGE_LEN;
	DEBUG_NOTICE("Lambda Polynomial:\n");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("%x ", lambda[i]);
		if(0xFF != lambda[i])
		{
			if(max_lambda_degree < i)
			{
				max_lambda_degree = i;
			}
		}
	}
	DEBUG_NOTICE("\n");
	DEBUG_NOTICE("max_lambda_degree: %ld\n", max_lambda_degree);

	unsigned char skip_root_flag = 0, root_found_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)//search for root
	{
		/*It is useless to check roots in unreliable group.*/
		skip_root_flag = 0;
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN); j++)
		{
			if(i == unrel_group_seq[j])
			{
				skip_root_flag = 1;
				break;
			}
		}
		if(1 == skip_root_flag)
		{
			continue;
		}

		/*Although even more errors in reliable group may exist, they cannot be corrected.*/
		if(root_found_cnt >= (SYN_LEN / 2))
		{
			break;
		}

		lambda_root = 0xFF;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			if(j > max_lambda_degree)
			{
				break;
			}
			if(0xFF == lambda[j])
			{
				continue;
			}
			//DEBUG_NOTICE("%d %d: %x %x %x %x\n", i, j, lambda_root, lambda[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda[j], (j * power_polynomial_table[i + 1][0])));
			lambda_root = gf_add(lambda_root, gf_multp(lambda[j], (j * power_polynomial_table[i + 1][0]) % (GF_FIELD - 1)));
		}
		DEBUG_NOTICE("%d: %x %x\n", i, lambda_root, power_polynomial_table[i + 1][0]);
		if(0xFF == lambda_root)
		{
			err_location[i] = i;
			root_found_cnt++;
		}
	}

	DEBUG_NOTICE("Error Location:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("%x ", err_location[i]);
		if(0xFF != err_location[i])
		{
			DEBUG_NOTICE("err_location: %d %d\n", i, err_location[i]);
		}
	}
	DEBUG_NOTICE("\n");

	/*derivative of lambda*/
	for(i = 0; i < (MESSAGE_LEN - 1); i++)
	{
		if(0 != ((i + 1) % 2))
		{
			lambda_dev[i] = lambda[i + 1];
		}
	}
	DEBUG_NOTICE("Lambda Derivative Polynomial:\n");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("%x ", lambda_dev[i]);
	}
	DEBUG_NOTICE("\n");

	for(i = 0; i < (MESSAGE_LEN + 1); i++)
	{
		if(0 != ((i + 1) % 2))
		{
			v_tmp[i] = v[i + 1];
		}
		else
		{
			v_tmp[i] = 0xFF;
		}
	}
	DEBUG_NOTICE("v_tmp Derivative Polynomial:\n");
	for(i = 0; i < (MESSAGE_LEN + 1); i++)
	{
		DEBUG_NOTICE("%x ", v_tmp[i]);
	}
	DEBUG_NOTICE("\n");

	/*magnitude of error*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0xFF != err_location[i])
		{
			DEBUG_NOTICE("err_location: %d | %x %x\n", i, err_location[i], err_mag[i]);
			tmp = 0xFF;

			tmp = poly_eva(omega, MESSAGE_LEN, err_location[i]);
			tmp = gf_multp(tmp,
				           poly_eva(v_tmp, MESSAGE_LEN, err_location[i]));
			tmp = gf_div(tmp,
				         poly_eva(lambda_dev, MESSAGE_LEN, err_location[i]));
			
			err_mag[i] = tmp;
		}
	}

	DEBUG_NOTICE("Error Magnitude:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("%x ", err_mag[i]);
	}
	DEBUG_NOTICE("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0xFF != err_mag[i])
		{
			tmp_cw[i] = gf_add(received_polynomial[i], err_mag[i]);
		}
		else
		{
			tmp_cw[i] = received_polynomial[i];
		}
	}

	return 0;
}
#endif

int recover_codeword()
{
	DEBUG_NOTICE("recover\n");
	long long i = 0, j = 0;
	long long cnt = 0;

	/*test*/
#if 0
	evaluation_encoding_v2(decoded_message, decoded_codeword);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		decoded_codeword[i] = gf_add(decoded_codeword[i], phi[i]);
		DEBUG_NOTICE("decoded_codeword: %d %x\n", i, decoded_codeword[i]);
	}
#endif

	//syndrome_re_encoding();//test
	//bm_re_encoding(decoded_message);

#if 0
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			if(i == rel_group_seq[j])
			{
				decoded_codeword[i] = decoded_message[cnt];
				cnt = cnt + 1;
			}
		}
		if(MESSAGE_LEN <= cnt)
		{
			break;
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		decoded_codeword[i] = gf_add(decoded_codeword[i], received_polynomial[i]);
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			if(i == unrel_group_seq[j])
			{
				decoded_codeword[i] = 0xFF;
			}
		}
		DEBUG_NOTICE("recover_decoded_codeword: %x\n", decoded_codeword[i]);
	}

	erasure_decoding(decoded_codeword);

	memcpy(decoded_codeword, phi, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(decoded_message, decoded_codeword + (CODEWORD_LEN - MESSAGE_LEN), sizeof(unsigned char) * MESSAGE_LEN);
#endif

	//memcpy(decoded_codeword, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN); j++)
		{
			if(i == unrel_group_seq[j])
			{
				decoded_codeword[i] = 0xFF;
				cnt++;
			}
		}

		if(cnt >= (CODEWORD_LEN - MESSAGE_LEN))
		{
			break;
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("decoded_codeword: %x\n", decoded_codeword[i]);
	}
	erasure_decoding(decoded_codeword, unrel_group_seq);
	memcpy(decoded_codeword, phi, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(decoded_message, decoded_codeword + (CODEWORD_LEN - MESSAGE_LEN), sizeof(unsigned char) * MESSAGE_LEN);

	return 0;
}

int re_encoding()
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long tmp = 0;
	unsigned char find_flag = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

#if 0
	syndrome_cal(received_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
#endif				  
#if (1 == RE_ENCODING)//use re-encoding
#if 0
	DEBUG_INFO("syn_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(syndrome, (CODEWORD_LEN - MESSAGE_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(encoded_polynomial, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif

	rel_group();
#if 0
	tao_cal();
	sigma_cal();
#endif

	v_cal();

#if 0//test
	unsigned char L[MESSAGE_LEN];
	l_cal(0x0, L);
#endif

	phi_cal();

	//erasure_decoding(erasure_polynomial);
#if 0
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			find_flag = 0;
			for(k = 0; k < MESSAGE_LEN; k++)/*it is inefficient in C, but it is easy in verilog*/
			{
				if(rel_group_seq[k] == j)
				{
					find_flag = 1;
					break;
				}
			}
#if 1			
			//if((0 == mul_matrix[i][j]) || (1 == find_flag))
			if(1 == find_flag)
			{
				beta_matrix[i][j] = 0xFF;
			}
			else
			{
				//DEBUG_INFO("%d %d %x\n", i, j, mul_matrix[i][j]);
				beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], received_polynomial[j]);
			}
#else
			beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], received_polynomial[j]);
#endif
		}
	}
#else

	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			/*for LCC, it is not necessary to calculate every beta*/
#if 1
			if((0 == mul_matrix[i][j])
				&& (i != chnl_rel_max_id[j])
				&& (i != chnl_rel_scd_id[j]))
			{
				beta_matrix[i][j] = 0xFF;
				continue;
			}
#endif			
		
			find_flag = 0;
			for(k = 0; k < MESSAGE_LEN; k++)/*it is inefficient in C, but it is easy in verilog*/
			{
				if(rel_group_seq[k] == j)
				{
					find_flag = 1;
					break;
				}
			}

			if(1 == find_flag)
			{
				beta_matrix[i][j] = 0xFF;
			}
			else
			{
				beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], phi[j]);
			}

			if(power_polynomial_table[i][0] == received_polynomial[j])
			{
				re_encoded_codeword[i] = beta_matrix[i][j];
				DEBUG_INFO("re_encoded_codeword: %d %d | %x %x | %x\n",
					       i,
					       j,
					       received_polynomial[j],
					       power_polynomial_table[i][0],
					       re_encoded_codeword[i]);
			}

			DEBUG_NOTICE("beta_matrix: %d %d %x %x %x %x\n",
				         i,
				         j,
				         power_polynomial_table[j + 1][0],
				         power_polynomial_table[i][0],
				         phi[j],
				         beta_matrix[i][j]);
		}
	}

#endif
#else//do not use re-encoding
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)
		{
			if(0 != mul_matrix[i][j])
			{
				beta_matrix[i][j] = power_polynomial_table[i][0];
			}
			else
			{
				beta_matrix[i][j] = 0xFF;
			}
		}
	}
#endif

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("beta:\n");
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)
		{
			DEBUG_IMPOTANT("%x ", beta_matrix[i][j]);
		}
		DEBUG_IMPOTANT("\n");
	}
	DEBUG_IMPOTANT("\n");
#endif

	return 0;
}
