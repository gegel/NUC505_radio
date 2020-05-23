/* vim: set tabstop=4:softtabstop=0:shiftwidth=4:noexpandtab */

 /*
  * ===================================================================
  *  TS 26.104
  *  REL-5 V5.4.0 2004-03
  *  REL-6 V6.1.0 2004-03
  *  3GPP AMR Floating-point Speech Codec
  * ===================================================================
  *
  */

// double

/*
 * sp_dec.c
 *
 *
 * Project:
 *    AMR Floating-Point Codec
 *
 * Contains:
 *    This module contains all the functions needed decoding AMR
 *    encoder parameters to 16-bit speech samples
 *
 */
/*
 * include files
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "platform.h"
//#include <memory.h>
#include <math.h>
#include "sp_dec.h"
#include "rom_dec.h"
//#include <ophtools.h>

int dec_mem=0;

void mtrap(int line);
void memadd(int len, unsigned char id);

/*
 * Declare structure types
 */
enum DTXStateType {
	SPEECH = 0, DTX, DTX_MUTE
};

/*
 * Decoder memory structure
 */


 //typedef struct {
	/* history vector of past synthesis speech energy */
 //	int32_t frameEnergyHist[L_ENERGYHIST];

	/* state flags */
  //	int16_t bgHangover;	/* counter; number of frames after last speech frame */

//} Bgn_scdState;


typedef struct {
	int32_t hangCount;	/* counter; */
	/* history vector of past synthesis speech energy */
	int32_t cbGainHistory[L_CBGAINHIST];
	int16_t hangVar;		/* counter; */

} Cb_gain_averageState;


typedef struct {
	int32_t lsp_meanSave[M];	// Averaged LSPs saved for efficiency

} lsp_avgState;


typedef struct {
	int32_t past_r_q[M];	/* Past quantized prediction error, Q15 */
	int32_t past_lsf_q[M];	/* Past dequantized lsfs, Q15 */

} D_plsfState;
typedef struct {
	int32_t pbuf[5];
	int32_t past_gain_pit;
	int32_t prev_gp;

} ec_gain_pitchState;
typedef struct {
	int32_t gbuf[5];
	int32_t past_gain_code;
	int32_t prev_gc;

} ec_gain_codeState;
typedef struct {
	/*
	 * normal MA predictor memory, Q10
	 * (contains 20*log10(quaErr))
	 */
	int32_t past_qua_en[4];

	/*
	 * MA predictor memory for MR122 mode, Q10
	 * (contains log2(quaErr))
	 */
	int32_t past_qua_en_MR122[4];

} gc_predState;
typedef struct {
	int32_t gainMem[PHDGAINMEMSIZE];
	int32_t prevCbGain;
	int32_t prevState;
	int16_t lockFull;
	int16_t onset;

} ph_dispState;


 /*
typedef struct {
	enum DTXStateType dtxGlobalState;	// contains previous state //

	int32_t log_en;
	int32_t old_log_en;
	int32_t pn_seed_rx;
	int32_t lsp[M];
	int32_t lsp_old[M];
	int32_t lsf_hist[M * DTX_HIST_SIZE];
	int32_t lsf_hist_mean[M * DTX_HIST_SIZE];
	int32_t log_en_hist[DTX_HIST_SIZE];
	int32_t true_sid_period_inv;
	int16_t since_last_sid;
	int16_t lsf_hist_ptr;
	int16_t log_pg_mean;
	int16_t log_en_hist_ptr;
	int16_t log_en_adjust;
	int16_t dtxHangoverCount;
	int16_t decAnaElapsedCount;
	int16_t sid_frame;
	int16_t valid_data;
	int16_t dtxHangoverAdded;

	// updated in main decoder //
	int16_t data_updated;	// marker to know if CNI data is ever renewed

} dtx_decState;

 */

typedef struct {
	int32_t past_gain;

} agcState;
typedef struct {
	/* Excitation vector */
	int32_t old_exc[L_SUBFR + PIT_MAX + L_INTERPOL]; //40+143+11
	int32_t *exc;
	int32_t lsp_old[M];

	/* Filter's memory */
	int32_t mem_syn[M];

	/* pitch sharpening */
	int32_t sharp;
	int32_t old_T0;

	/* Variable holding received ltpLag, used in background noise and BFI */
	int32_t T0_lagBuff;

	/* Variables for the source characteristic detector (SCD) */
	int32_t inBackgroundNoise;
	int32_t voicedHangover;
	int32_t ltpGainHistory[9];

	/* Memories for bad frame handling */
	int32_t excEnergyHist[9];
	int16_t prev_bf;
	int16_t prev_pdf;
	int16_t state;
	int16_t nodataSeed;

	//Bgn_scdState *background_state;
	Cb_gain_averageState *Cb_gain_averState;
	lsp_avgState *lsp_avg_st;
	D_plsfState *lsfState;
	ec_gain_pitchState *ec_gain_p_st;
	ec_gain_codeState *ec_gain_c_st;
	gc_predState *pred_state;
	ph_dispState *ph_disp_st;
	//dtx_decState *dtxDecoderState;
} Decoder_amrState;
typedef struct {
	int32_t res2[L_SUBFR];
	int32_t mem_syn_pst[M];
	int32_t synth_buf[M + L_FRAME];
	int32_t preemph_state_mem_pre;
	agcState *agc_state;
} Post_FilterState;
typedef struct {
	int32_t y2_hi;
	int32_t y2_lo;
	int32_t y1_hi;
	int32_t y1_lo;
	int32_t x0;
	int32_t x1;

} Post_ProcessState;
typedef struct {
	Decoder_amrState *decoder_amrState;
	Post_FilterState *post_state;
	Post_ProcessState *postHP_state;
} Speech_Decode_FrameState;

/*
 * CodAmrReset
 *
 *
 * Parameters:
 *    state             B: state structure
 *    mode              I: AMR mode
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *    void
 */
static void Decoder_amr_reset(Decoder_amrState * state, enum Mode mode)
{
	int32_t i;

	/* Cb_gain_average_reset */
	memset(state->Cb_gain_averState->cbGainHistory, 0, L_CBGAINHIST << 2);
	state->Cb_gain_averState->hangVar = 0;
	state->Cb_gain_averState->hangCount = 0;

	/* Initialize static pointer */
	state->exc = state->old_exc + PIT_MAX + L_INTERPOL;

	/* Static vectors to zero */
	memset(state->old_exc, 0, (PIT_MAX + L_INTERPOL) << 2);

	if (mode != MRDTX)
		memset(state->mem_syn, 0, M << 2);

	/* initialize pitch sharpening */
	state->sharp = SHARPMIN;
	state->old_T0 = 40;

	/* Initialize state->lsp_old [] */
	if (mode != MRDTX) {
		state->lsp_old[0] = 30000;
		state->lsp_old[1] = 26000;
		state->lsp_old[2] = 21000;
		state->lsp_old[3] = 15000;
		state->lsp_old[4] = 8000;
		state->lsp_old[5] = 0;
		state->lsp_old[6] = -8000;
		state->lsp_old[7] = -15000;
		state->lsp_old[8] = -21000;
		state->lsp_old[9] = -26000;
	}

	/* Initialize memories of bad frame handling */
	state->prev_bf = 0;
	state->prev_pdf = 0;
	state->state = 0;
	state->T0_lagBuff = 40;
	state->inBackgroundNoise = 0;
	state->voicedHangover = 0;

	if (mode != MRDTX)
		memset(state->excEnergyHist, 0, 9 << 2);
	memset(state->ltpGainHistory, 0, 9 << 2);


	if (mode != MRDTX) {
		state->lsp_avg_st->lsp_meanSave[0] = 1384;
		state->lsp_avg_st->lsp_meanSave[1] = 2077;
		state->lsp_avg_st->lsp_meanSave[2] = 3420;
		state->lsp_avg_st->lsp_meanSave[3] = 5108;
		state->lsp_avg_st->lsp_meanSave[4] = 6742;
		state->lsp_avg_st->lsp_meanSave[5] = 8122;
		state->lsp_avg_st->lsp_meanSave[6] = 9863;
		state->lsp_avg_st->lsp_meanSave[7] = 11092;
		state->lsp_avg_st->lsp_meanSave[8] = 12714;
		state->lsp_avg_st->lsp_meanSave[9] = 13701;
	}
        

        
	memset(state->lsfState->past_r_q, 0, M << 2);

	/* Past dequantized lsfs */
	state->lsfState->past_lsf_q[0] = 1384;
	state->lsfState->past_lsf_q[1] = 2077;
	state->lsfState->past_lsf_q[2] = 3420;
	state->lsfState->past_lsf_q[3] = 5108;
	state->lsfState->past_lsf_q[4] = 6742;
	state->lsfState->past_lsf_q[5] = 8122;
	state->lsfState->past_lsf_q[6] = 9863;
	state->lsfState->past_lsf_q[7] = 11092;
	state->lsfState->past_lsf_q[8] = 12714;
	state->lsfState->past_lsf_q[9] = 13701;

	for (i = 0; i < 5; i++)
		state->ec_gain_p_st->pbuf[i] = 1640;
	state->ec_gain_p_st->past_gain_pit = 0;
	state->ec_gain_p_st->prev_gp = 16384;

	for (i = 0; i < 5; i++)
		state->ec_gain_c_st->gbuf[i] = 1;
	state->ec_gain_c_st->past_gain_code = 0;
	state->ec_gain_c_st->prev_gc = 1;

	if (mode != MRDTX) {
		for (i = 0; i < NPRED; i++) {
			state->pred_state->past_qua_en[i] = MIN_ENERGY;
			state->pred_state->past_qua_en_MR122[i] =
			    MIN_ENERGY_MR122;
		}
	}
	state->nodataSeed = 21845;

	/* Static vectors to zero */
	//memset(state->background_state->frameEnergyHist, 0, L_ENERGYHIST << 2);

	/* Initialize hangover handling */
	//state->background_state->bgHangover = 0;

	/* phDispReset */
	memset(state->ph_disp_st->gainMem, 0, PHDGAINMEMSIZE << 2);
	state->ph_disp_st->prevState = 0;
	state->ph_disp_st->prevCbGain = 0;
	state->ph_disp_st->lockFull = 0;
	state->ph_disp_st->onset = 0;	/* assume no onset in start */


        /*
	if (mode != MRDTX) {
		state->dtxDecoderState->since_last_sid = 0;
		state->dtxDecoderState->true_sid_period_inv = 8192;
		state->dtxDecoderState->log_en = 3500;
		state->dtxDecoderState->old_log_en = 3500;

		// low level noise for better performance in  DTX handover cases
		state->dtxDecoderState->pn_seed_rx = PN_INITIAL_SEED;

		// Initialize state->lsp []
		state->dtxDecoderState->lsp[0] = 30000;
		state->dtxDecoderState->lsp[1] = 26000;
		state->dtxDecoderState->lsp[2] = 21000;
		state->dtxDecoderState->lsp[3] = 15000;
		state->dtxDecoderState->lsp[4] = 8000;
		state->dtxDecoderState->lsp[5] = 0;
		state->dtxDecoderState->lsp[6] = -8000;
		state->dtxDecoderState->lsp[7] = -15000;
		state->dtxDecoderState->lsp[8] = -21000;
		state->dtxDecoderState->lsp[9] = -26000;

		// Initialize state->lsp_old []
		state->dtxDecoderState->lsp_old[0] = 30000;
		state->dtxDecoderState->lsp_old[1] = 26000;
		state->dtxDecoderState->lsp_old[2] = 21000;
		state->dtxDecoderState->lsp_old[3] = 15000;
		state->dtxDecoderState->lsp_old[4] = 8000;
		state->dtxDecoderState->lsp_old[5] = 0;
		state->dtxDecoderState->lsp_old[6] = -8000;
		state->dtxDecoderState->lsp_old[7] = -15000;
		state->dtxDecoderState->lsp_old[8] = -21000;
		state->dtxDecoderState->lsp_old[9] = -26000;
		state->dtxDecoderState->lsf_hist_ptr = 0;
		state->dtxDecoderState->log_pg_mean = 0;
		state->dtxDecoderState->log_en_hist_ptr = 0;

		// initialize decoder lsf history
		state->dtxDecoderState->lsf_hist[0] = 1384;
		state->dtxDecoderState->lsf_hist[1] = 2077;
		state->dtxDecoderState->lsf_hist[2] = 3420;
		state->dtxDecoderState->lsf_hist[3] = 5108;
		state->dtxDecoderState->lsf_hist[4] = 6742;
		state->dtxDecoderState->lsf_hist[5] = 8122;
		state->dtxDecoderState->lsf_hist[6] = 9863;
		state->dtxDecoderState->lsf_hist[7] = 11092;
		state->dtxDecoderState->lsf_hist[8] = 12714;
		state->dtxDecoderState->lsf_hist[9] = 13701;

		for (i = 1; i < DTX_HIST_SIZE; i++) {
			memcpy(&state->dtxDecoderState->lsf_hist[M * i],
			       &state->dtxDecoderState->lsf_hist[0], M << 2);
		}
		memset(state->dtxDecoderState->lsf_hist_mean, 0,
			M * DTX_HIST_SIZE << 2);

		// initialize decoder log frame energy
		for (i = 0; i < DTX_HIST_SIZE; i++) {
			state->dtxDecoderState->log_en_hist[i] =
			    state->dtxDecoderState->log_en;
		}
		state->dtxDecoderState->log_en_adjust = 0;
		state->dtxDecoderState->dtxHangoverCount = DTX_HANG_CONST;
		state->dtxDecoderState->decAnaElapsedCount = 31;
		state->dtxDecoderState->sid_frame = 0;
		state->dtxDecoderState->valid_data = 0;
		state->dtxDecoderState->dtxHangoverAdded = 0;
		state->dtxDecoderState->dtxGlobalState = DTX;
		state->dtxDecoderState->data_updated = 0;
	}

        */
	return;
}


/*
 * Lsf_lsp
 *
 *
 * Parameters:
 *    lsf               I: vector of LSFs
 *    lsp               O: vector of LSPs
 *
 * Function:
 *    Transformation lsf to lsp, order M
 *
 * Returns:
 *    void
 */
static void Lsf_lsp(int32_t lsf[], int32_t lsp[])
{
	int32_t i, ind, offset, tmp;

	for (i = 0; i < M; i++) {
		/* ind = b8-b15 of lsf[i] */
		ind = lsf[i] >> 8;

		/* offset = b0-b7  of lsf[i] */
		offset = lsf[i] & 0x00ff;

		/* lsp[i] = table[ind]+ ((table[ind+1]-table[ind])*offset) / 256 */
		tmp = ((cos_table[ind + 1] - cos_table[ind]) * offset) << 1;
		lsp[i] = cos_table[ind] + (tmp >> 9);
	}
	return;
}

/*
 * D_plsf_3
 *
 *
 * Parameters:
 *    st->past_lsf_q    I: Past dequantized LFSs
 *    st->past_r_q      B: past quantized residual
 *    mode              I: AMR mode
 *    bfi               B: bad frame indicator
 *    indice            I: quantization indices of 3 submatrices, Q0
 *    lsp1_q            O: quantized 1st LSP vector
 *
 * Function:
 *    Decodes the LSP parameters using the received quantization indices.
 *    1st order MA prediction and split by 3 vector quantization (split-VQ)
 *
 * Returns:
 *    void
 */
static void D_plsf_3(D_plsfState * st, enum Mode mode, int16_t bfi, int16_t *
		     indice, int32_t * lsp1_q)
{
	int32_t lsf1_r[M], lsf1_q[M];
	int32_t i, index, temp;
	const short *p_cb1, *p_cb2, *p_cb3, *p_dico;

        //short test[16]={0,};

	/* if bad frame */
	if (bfi != 0) {
		/* use the past LSFs slightly shifted towards their mean */
		for (i = 0; i < M; i++) {
			/* lsfi_q[i] = ALPHA*past_lsf_q[i] + ONE_ALPHA*meanLsf[i]; */
			lsf1_q[i] =
			    ((st->past_lsf_q[i] * ALPHA) >> 15) +
			    ((mean_lsf_3[i]
			      * ONE_ALPHA) >> 15);
		}

		/* estimate past quantized residual to be used in next frame */
		if (mode != MRDTX) {
			for (i = 0; i < M; i++) {
				/* temp  = meanLsf[i] +  pastR2_q[i] * pred_fac; */
				temp =
				    mean_lsf_3[i] +
				    ((st->past_r_q[i] * pred_fac[i]) >> 15);
				st->past_r_q[i] = lsf1_q[i] - temp;
			}
		} else {
			for (i = 0; i < M; i++) {
				/* temp  = meanLsf[i] +  pastR2_q[i]; */
				temp = mean_lsf_3[i] + st->past_r_q[i];
				st->past_r_q[i] = lsf1_q[i] - temp;
			}
		}
	}

	/* if good LSFs received */
	else {
		if ((mode == MR475) | (mode == MR515)) {
			/* MR475, MR515 */
			p_cb1 = dico1_lsf_3;
			p_cb2 = dico2_lsf_3;
			p_cb3 = mr515_3_lsf;

                        //memcpy(test, dico1_lsf_3, 32);
                        //i=i;
		}
                /*
                else if (mode == MR795) {
			// MR795
			p_cb1 = mr795_1_lsf;
			p_cb2 = dico2_lsf_3;
			p_cb3 = dico3_lsf_3;
		} else {
			//MR59, MR67, MR74, MR102, MRDTX
			p_cb1 = dico1_lsf_3;
			p_cb2 = dico2_lsf_3;
			p_cb3 = dico3_lsf_3;
		}
                */
		/* decode prediction residuals from 3 received indices */
		index = *indice++;
		p_dico = &p_cb1[index + index + index];
		index = *indice++;
		lsf1_r[0] = *p_dico++;
		lsf1_r[1] = *p_dico++;
		lsf1_r[2] = *p_dico++;

	       //	if ((mode == MR475) | (mode == MR515)) {
			/* MR475, MR515 only using every second entry */
			index = index << 1;
		//}
		p_dico = &p_cb2[index + index + index];
		index = *indice++;
		lsf1_r[3] = *p_dico++;
		lsf1_r[4] = *p_dico++;
		lsf1_r[5] = *p_dico;
		p_dico = &p_cb3[index << 2];
		lsf1_r[6] = *p_dico++;
		lsf1_r[7] = *p_dico++;
		lsf1_r[8] = *p_dico++;
		lsf1_r[9] = *p_dico++;

		/* Compute quantized LSFs and update the past quantized residual */
		if (mode != MRDTX) {
			for (i = 0; i < M; i++) {
				lsf1_q[i] =
				    lsf1_r[i] + (mean_lsf_3[i] +
						 ((st->past_r_q[i] *
						   pred_fac[i]) >> 15));
			}
			memcpy(st->past_r_q, lsf1_r, M << 2);
		} else {
			for (i = 0; i < M; i++) {
				lsf1_q[i] =
				    lsf1_r[i] + (mean_lsf_3[i] +
						 st->past_r_q[i]);
			}
			memcpy(st->past_r_q, lsf1_r, M << 2);
		}
	}

	/* verification that LSFs has minimum distance of LSF_GAP Hz */
	temp = LSF_GAP;

	for (i = 0; i < M; i++) {
		if (lsf1_q[i] < temp) {
			lsf1_q[i] = temp;
		}
		temp = lsf1_q[i] + LSF_GAP;
	}
	memcpy(st->past_lsf_q, lsf1_q, M << 2);

	/*  convert LSFs to the cosine domain */
	Lsf_lsp(lsf1_q, lsp1_q);
	return;
}





/* VC5.0 Global optimization does not work with this function */
#if _MSC_VER == 1100
#pragma optimize( "g", off )
#endif
/*
 * Get_lsp_pol
 *
 *
 * Parameters:
 *    lsp               I: line spectral frequencies
 *    f                 O: polynomial F1(z) or F2(z)
 *
 * Function:
 *    Find the polynomial F1(z) or F2(z) from the LSPs.
 *
 *    F1(z) = product ( 1 - 2 lsp[i] z^-1 + z^-2 )
 *             i=0,2,4,6,8
 *    F2(z) = product   ( 1 - 2 lsp[i] z^-1 + z^-2 )
 *             i=1,3,5,7,9
 *
 *    where lsp[] is the LSP vector in the cosine domain.
 *
 *    The expansion is performed using the following recursion:
 *
 *    f[0] = 1
 *    b = -2.0 * lsp[0]
 *    f[1] = b
 *    for i=2 to 5 do
 *       b = -2.0 * lsp[2*i-2];
 *       f[i] = b*f[i-1] + 2.0*f[i-2];
 *       for j=i-1 down to 2 do
 *          f[j] = f[j] + b*f[j-1] + f[j-2];
 *       f[1] = f[1] + b;
 *
 * Returns:
 *    void
 */
static void Get_lsp_pol(int32_t * lsp, int32_t * f)
{
	volatile int32_t f0, f1, f2, f3, f4, f5;
	int32_t l1, l2, l3, l4;

	/* f[0] = 1.0; */
	f0 = 16777216L;

	/* f1 = *lsp * -1024; */
	f1 = -lsp[0] << 10;
	l1 = lsp[2];
	l2 = lsp[4];
	l3 = lsp[6];
	l4 = lsp[8];
	f2 = f0 << 1;
	f2 -= (((f1 >> 16) * l1) + (((f1 & 0xFFFE) * l1) >> 16)) << 2;
	f1 -= l1 << 10;
	f3 = f1 << 1;
	f3 -= (((f2 >> 16) * l2) + (((f2 & 0xFFFE) * l2) >> 16)) << 2;
	f2 += f0;
	f2 -= (((f1 >> 16) * l2) + (((f1 & 0xFFFE) * l2) >> 16)) << 2;
	f1 -= l2 << 10;
	f4 = f2 << 1;
	f4 -= (((f3 >> 16) * l3) + (((f3 & 0xFFFE) * l3) >> 16)) << 2;
	f3 += f1;
	f3 -= (((f2 >> 16) * l3) + (((f2 & 0xFFFE) * l3) >> 16)) << 2;
	f2 += f0;
	f2 -= (((f1 >> 16) * l3) + (((f1 & 0xFFFE) * l3) >> 16)) << 2;
	f1 -= l3 << 10;
	f5 = f3 << 1;
	f5 -= (((f4 >> 16) * l4) + (((f4 & 0xFFFE) * l4) >> 16)) << 2;
	f4 += f2;
	f4 -= (((f3 >> 16) * l4) + (((f3 & 0xFFFE) * l4) >> 16)) << 2;
	f3 += f1;
	f3 -= (((f2 >> 16) * l4) + (((f2 & 0xFFFE) * l4) >> 16)) << 2;
	f2 += f0;
	f2 -= (((f1 >> 16) * l4) + (((f1 & 0xFFFE) * l4) >> 16)) << 2;
	f1 -= l4 << 10;
	f[0] = f0;
	f[1] = f1;
	f[2] = f2;
	f[3] = f3;
	f[4] = f4;
	f[5] = f5;
	return;
}

#if _MSC_VER == 1100
#pragma optimize( "", on )
#endif

/*
 * Lsp_Az
 *
 *
 * Parameters:
 *    lsp                 I: Line spectral frequencies
 *    a                   O: Predictor coefficients
 *
 * Function:
 *    Converts from the line spectral pairs (LSP) to LP coefficients,
 *    for a 10th order filter.
 *
 *    Find the coefficients of F1(z) and F2(z)
 *    Multiply F1(z) by 1+z^{-1} and F2(z) by 1-z^{-1}
 *    A(z) = ( F1(z) + F2(z) ) / 2
 *
 * Returns:
 *    void
 */
static void Lsp_Az(int32_t lsp[], int32_t a[])
{
	int32_t f1[6], f2[6];
	int32_t T0, i, j;

	Get_lsp_pol(&lsp[0], f1);
	Get_lsp_pol(&lsp[1], f2);

	for (i = 5; i > 0; i--) {
		f1[i] += f1[i - 1];
		f2[i] -= f2[i - 1];
	}
	a[0] = 4096;

	for (i = 1, j = 10; i <= 5; i++, j--) {
		T0 = f1[i] + f2[i];
		a[i] = (int16_t) (T0 >> 13);	/* emulate fixed point bug */
		if ((T0 & 4096) != 0) {
			a[i]++;
		}
		T0 = f1[i] - f2[i];
		a[j] = (int16_t) (T0 >> 13);	/* emulate fixed point bug */

		if ((T0 & 4096) != 0) {
			a[j]++;
		}
	}
	return;
}


/*
 * Log2_norm
 *
 *
 * Parameters:
 *    x                 I: input value
 *    exp               I: exponent
 *    exponent          O: Integer part of Log2. (range: 0<=val<=30)
 *    fraction          O: Fractional part of Log2. (range: 0<=val<1)
 *
 * Function:
 *    Computes log2
 *
 *    Computes log2(L_x, exp),  where   L_x is positive and
 *    normalized, and exp is the normalisation exponent
 *    If L_x is negative or zero, the result is 0.
 *
 *    The function Log2(L_x) is approximated by a table and linear
 *    interpolation. The following steps are used to compute Log2(L_x)
 *
 *    exponent = 30-normExponent
 *    i = bit25-b31 of L_x;  32<=i<=63  (because of normalization).
 *    a = bit10-b24
 *    i -=32
 *    fraction = table[i]<<16 - (table[i] - table[i+1]) * a * 2
 *
 * Returns:
 *    void
 */
static void Log2_norm(int32_t x, int32_t exp, int32_t * exponent, int32_t *
		      fraction)
{
	int32_t y, i, a;

	if (x <= 0) {
		*exponent = 0;
		*fraction = 0;
		return;
	}

	/* Extract b25-b31 */
	i = x >> 25;
	i = i - 32;

	/* Extract b10-b24 of fraction */
	a = x >> 9;
	a = a & 0xFFFE;		/* 2a */

	/* fraction */
	y = (log2_table[i] << 16) - a * (log2_table[i] - log2_table[i + 1]);
	*fraction = y >> 16;
	*exponent = 30 - exp;
	return;
}

/*
 * Log2
 *
 *
 * Parameters:
 *    x                 I: input value
 *    exponent          O: Integer part of Log2. (range: 0<=val<=30)
 *    fraction          O: Fractional part of Log2. (range: 0<=val<1)
 *
 * Function:
 *    Computes log2(L_x)
 *    If x is negative or zero, the result is 0.
 *
 * Returns:
 *    void
 */
static void Log2(int32_t x, int32_t * exponent, int32_t * fraction)
{
	int tmp, exp = 0;

	if (x != 0) {
		tmp = x;
		while (!((tmp & 0x80000000) ^ ((tmp & 0x40000000) << 1))) {
			exp++;
			tmp = tmp << 1;
		}
	}
	Log2_norm(x << exp, exp, exponent, fraction);
}

/*
 * Pow2
 *
 *
 * Parameters:
 *    exponent          I: Integer part. (range: 0<=val<=30)
 *    fraction          O: Fractional part. (range: 0.0<=val<1.0)
 *
 * Function:
 *    pow(2.0, exponent.fraction)
 *
 *    The function Pow2(L_x) is approximated by a table and linear interpolation.
 *
 *    i = bit10-b15 of fraction, 0 <= i <= 31
 *    a = biT0-b9   of fraction
 *    x = table[i]<<16 - (table[i] - table[i+1]) * a * 2
 *    x = L_x >> (30-exponent) (with rounding)
 *
 * Returns:
 *    result (range: 0<=val<=0x7fffffff)
 */
static int32_t Pow2(int32_t exponent, int32_t fraction)
{
	int32_t i, a, tmp, x, exp;

	/* Extract b10-b16 of fraction */
	i = fraction >> 10;

	/* Extract b0-b9 of fraction */
	a = (fraction << 5) & 0x7fff;

	/* table[i] << 16 */
	x = pow2_table[i] << 16;

	/* table[i] - table[i+1] */
	tmp = pow2_table[i] - pow2_table[i + 1];

	/* L_x -= tmp*a*2 */
	x -= (tmp * a) << 1;

	if (exponent >= -1) {
		exp = (30 - exponent);

		/* Rounding */
		if ((x & ((int32_t) 1 << (exp - 1))) != 0) {
			x = (x >> exp) + 1;
		} else
			x = x >> exp;
	} else
		x = 0;
	return (x);
}





/*
 * Syn_filt
 *
 *
 * Parameters:
 *    a                 I: prediction coefficients [M+1]
 *    x                 I: input signal
 *    y                 O: output signal
 *    lg                I: size of filtering
 *    mem               B: memory associated with this filtering
 *    update            I: 0=no update, 1=update of memory.
 *
 * Function:
 *    Perform synthesis filtering through 1/A(z).
 *
 * Returns:
 *    void
 */
static int32_t Syn_filt(int32_t a[], int32_t x[], int32_t y[], int32_t lg,
		       int32_t mem[]
		       , int32_t update)
{
	int32_t tmp[50];		/* malloc is slow */
	int32_t s, a0, overflow = 0;
	int32_t *yy, *yy_limit;

	/* Copy mem[] to yy[] */
	memcpy(tmp, mem, 40);
	yy = tmp + M;
	yy_limit = yy + lg;
	a0 = a[0];

	/* Do the filtering. */
	while (yy < yy_limit) {

		s = *x++ * a0;
		s -= yy[-1] * a[1];
		s -= yy[-2] * a[2];
		s -= yy[-3] * a[3];
		s -= yy[-4] * a[4];
		s -= yy[-5] * a[5];
		s -= yy[-6] * a[6];
		s -= yy[-7] * a[7];
		s -= yy[-8] * a[8];
		s -= yy[-9] * a[9];
		s -= yy[-10] * a[10];
		if (labs(s) < 0x7ffffff)
			*yy = (s + 0x800L) >> 12;
		else if (s > 0) {
			*yy = 32767;
			overflow = 1;
		} else {
			*yy = -32768;
			overflow = 1;
		}
		yy++;
	}
	memcpy(y, &tmp[M], lg << 2);

	/* Update of memory if update==1 */
	if (update) {
		memcpy(mem, &y[lg - M], 40);
	}
	return overflow;
}


/*
 * lsp_avg
 *
 *
 * Parameters:
 *    st->lsp_meanSave  B: LSP averages
 *    lsp               I: LSPs
 *
 * Function:
 *    Calculate the LSP averages
 *
 * Returns:
 *    void
 */
static void lsp_avg(lsp_avgState * st, int32_t * lsp)
{
	int32_t i, tmp;

	for (i = 0; i < M; i++) {
		/* mean = 0.84*mean */
		tmp = (st->lsp_meanSave[i] << 16);
		tmp -= (EXPCONST * st->lsp_meanSave[i]) << 1;

		/* Add 0.16 of newest LSPs to mean */
		tmp += (EXPCONST * lsp[i]) << 1;

		/* Save means */
		tmp += 0x00008000L;
		st->lsp_meanSave[i] = tmp >> 16;
	}
	return;
}

/*
 * Syn_filt_overflow
 *
 *
 * Parameters:
 *    a                 I: prediction coefficients [M+1]
 *    x                 I: input signal
 *    y                 O: output signal
 *    lg                I: size of filtering
 *    mem               B: memory associated with this filtering
 *    update            I: 0=no update, 1=update of memory.
 *
 * Function:
 *    Perform synthesis filtering through 1/A(z).
 *    Saturate after every multiplication.
 * Returns:
 *    void
 */
static void Syn_filt_overflow(int32_t a[], int32_t x[], int32_t y[], int32_t lg,
			      int32_t mem[]
			      , int32_t update)
{
	int32_t tmp[50];		/* malloc is slow */
	int32_t i, j, s, a0;
	int32_t *yy;

	/* Copy mem[] to yy[] */
	memcpy(tmp, mem, 40);
	yy = tmp + M;
	a0 = a[0];

	/* Do the filtering. */
	for (i = 0; i < lg; i++) {
		s = x[i] * a0;

		for (j = 1; j <= M; j++) {
			s -= a[j] * yy[-j];
			if (s > 1073741823) {
				s = 1073741823;
			} else if (s < -1073741824) {
				s = -1073741824;
			}
		}

		if (labs(s) < 0x7FFE800)
			*yy = (s + 0x800L) >> 12;
		else if (s > 0) {
			*yy = 32767;
		} else {
			*yy = -32768;
		}
		yy++;
	}
	memcpy(y, &tmp[M], lg << 2);

	/* Update of memory if update==1 */
	if (update) {
		memcpy(mem, &y[lg - M], 40);
	}
	return;
}



/*
 * Int_lpc_1to3
 *
 *
 * Parameters:
 *    lsp_old           I: LSP vector at the 4th subframe of past frame    [M]
 *    lsp_new           I: LSP vector at the 4th subframe of present frame [M]
 *    Az                O: interpolated LP parameters in all subframes
 *                                                                   [AZ_SIZE]
 *
 * Function:
 *    Interpolates the LSPs and converts to LPC parameters to get a different
 *    LP filter in each subframe.
 *
 *    The 20 ms speech frame is divided into 4 subframes.
 *    The LSPs are quantized and transmitted at the 4th
 *    subframes (once per frame) and interpolated at the
 *    1st, 2nd and 3rd subframe.
 *
 * Returns:
 *    void
 */
static void Int_lpc_1to3(int32_t lsp_old[], int32_t lsp_new[], int32_t Az[])
{
	int32_t lsp[M];
	int32_t i;

	for (i = 0; i < 10; i++) {
		lsp[i] = (lsp_new[i] >> 2) + (lsp_old[i] - (lsp_old[i] >> 2));
	}

	/* Subframe 1 */
	Lsp_Az(lsp, Az);
	Az += MP1;

	for (i = 0; i < 10; i++) {
		lsp[i] = (lsp_old[i] >> 1) + (lsp_new[i] >> 1);
	}

	/* Subframe 2 */
	Lsp_Az(lsp, Az);
	Az += MP1;

	for (i = 0; i < 10; i++) {
		lsp[i] = (lsp_old[i] >> 2) + (lsp_new[i] - (lsp_new[i] >> 2));
	}

	/* Subframe 3 */
	Lsp_Az(lsp, Az);
	Az += MP1;

	/* Subframe 4 */
	Lsp_Az(lsp_new, Az);
	return;
}

/*
 * D_plsf_5
 *
 *
 * Parameters:
 *    st->past_lsf_q I: Past dequantized LFSs
 *    st->past_r_q      B: past quantized residual
 *    bfi               B: bad frame indicator
 *    indice            I: quantization indices of 3 submatrices, Q0
 *    lsp1_q            O: quantized 1st LSP vector
 *    lsp2_q            O: quantized 2nd LSP vector
 *
 * Function:
 *    Decodes the 2 sets of LSP parameters in a frame
 *    using the received quantization indices.
 *
 * Returns:
 *    void
 */



/*
 * Dec_lag3
 *
 *
 * Parameters:
 *    index             I: received pitch index
 *    t0_min            I: minimum of search range
 *    t0_max            I: maximum of search range
 *    i_subfr           I: subframe flag
 *    T0_prev           I: integer pitch delay of last subframe used
 *                         in 2nd and 4th subframes
 *    T0                O: integer part of pitch lag
 *    T0_frac           O : fractional part of pitch lag
 *    flag4             I : flag for encoding with 4 bits
 * Function:
 *    Decoding of fractional pitch lag with 1/3 resolution.
 *    Extract the integer and fraction parts of the pitch lag from
 *    the received adaptive codebook index.
 *
 *    The fractional lag in 1st and 3rd subframes is encoded with 8 bits
 *    while that in 2nd and 4th subframes is relatively encoded with 4, 5
 *    and 6 bits depending on the mode.
 *
 * Returns:
 *    void
 */
static void Dec_lag3(int32_t index, int32_t t0_min, int32_t t0_max, int32_t i_subfr,
		     int32_t T0_prev, int32_t * T0, int32_t * T0_frac,
		     int32_t flag4)
{
	int32_t i, tmp_lag;

	/* if 1st or 3rd subframe */
	if (i_subfr == 0) {
		if (index < 197) {
			*T0 = (((index + 2) * 10923) >> 15) + 19;
			i = *T0 + *T0 + *T0;
			*T0_frac = (index - i) + 58;
		} else {
			*T0 = index - 112;
			*T0_frac = 0;
		}
	}

	/* 2nd or 4th subframe */
	else {
		if (flag4 == 0) {
			/* 'normal' decoding: either with 5 or 6 bit resolution */
			i = (((index + 2) * 10923) >> 15) - 1;
			*T0 = i + t0_min;
			i = i + i + i;
			*T0_frac = (index - 2) - i;
		} else {
			/* decoding with 4 bit resolution */
			tmp_lag = T0_prev;

			if ((tmp_lag - t0_min) > 5)
				tmp_lag = t0_min + 5;

			if ((t0_max - tmp_lag) > 4)
				tmp_lag = t0_max - 4;

			if (index < 4) {
				i = (tmp_lag - 5);
				*T0 = i + index;
				*T0_frac = 0;
			} else {
				if (index < 12) {
					i = (((index - 5) * 10923) >> 15) - 1;
					*T0 = i + tmp_lag;
					i = i + i + i;
					*T0_frac = (index - 9) - i;
				} else {
					i = (index - 12) + tmp_lag;
					*T0 = i + 1;
					*T0_frac = 0;
				}
			}
		}		/* end if (decoding with 4 bit resolution) */
	}
	return;
}

/*
 * Pred_lt_3or6_40
 *
 *
 * Parameters:
 *    exc               B: excitation buffer
 *    T0                I: integer pitch lag
 *    frac              I: fraction of lag
 *    flag3             I: if set, upsampling rate = 3 (6 otherwise)
 *
 * Function:
 *    Compute the result of long term prediction with fractional
 *    interpolation of resolution 1/3 or 1/6. (Interpolated past excitation).
 *
 *    Once the fractional pitch lag is determined,
 *    the adaptive codebook vector v(n) is computed by interpolating
 *    the past excitation signal u(n) at the given integer delay k
 *    and phase (fraction)  :
 *
 *          9                       9
 *    v(n) = SUM[ u(n-k-i) * b60(t+i*6) ] + SUM[ u(n-k+1+i) * b60(6-t+i*6) ],
 *          i=0                       i=0
 *    n = 0, ...,39, t = 0, ...,5.
 *
 *    The interpolation filter b60 is based on a Hamming windowed sin(x)/x
 *    function truncated at � 59 and padded with zeros at � 60 (b60(60)=0)).
 *    The filter has a cut-off frequency (-3 dB) at 3 600 Hz in
 *    the over-sampled domain.
 *
 * Returns:
 *    void
 */
static void Pred_lt_3or6_40(int32_t exc[], int32_t T0, int32_t frac, int32_t flag3)
{
	int32_t s, i;
	int32_t *x0, *x1, *x2;
	const short *c1, *c2;

	x0 = &exc[-T0];
	frac = -frac;

	if (flag3 != 0) {
		frac <<= 1;	/* inter_3l[k] = inter6[2*k] -> k' = 2*k */
	}

	if (frac < 0) {
		frac += 6;
		x0--;
	}
	c1 = &inter6[frac];
	c2 = &inter6[6 - frac];

	for (i = 0; i < 40; i++) {
		x1 = x0++;
		x2 = x0;
		s = x1[0] * c1[0];
		s += x1[-1] * c1[6];
		s += x1[-2] * c1[12];
		s += x1[-3] * c1[18];
		s += x1[-4] * c1[24];
		s += x1[-5] * c1[30];
		s += x1[-6] * c1[36];
		s += x1[-7] * c1[42];
		s += x1[-8] * c1[48];
		s += x1[-9] * c1[54];
		s += x2[0] * c2[0];
		s += x2[1] * c2[6];
		s += x2[2] * c2[12];
		s += x2[3] * c2[18];
		s += x2[4] * c2[24];
		s += x2[5] * c2[30];
		s += x2[6] * c2[36];
		s += x2[7] * c2[42];
		s += x2[8] * c2[48];
		s += x2[9] * c2[54];
		exc[i] = (s + 0x4000) >> 15;

	}
}

/*
 * decode_2i40_9bits
 *
 *
 * Parameters:
 *    subNr             I: subframe number
 *    sign              I: signs of 2 pulses
 *    index             I: Positions of the 2 pulses
 *    cod               O: algebraic (fixed) codebook excitation
 *
 * Function:
 *    Algebraic codebook decoder
 *
 * Returns:
 *    void
 */
static void decode_2i40_9bits(int32_t subNr, int32_t sign, int32_t index, int32_t
			      cod[])
{
	int32_t pos[2];
	int32_t i, j, k;

	/* Decode the positions */
	/* table bit  is the MSB */
	j = (index & 64) >> 6;
	i = index & 7;

	/* pos0 =i*5+startPos[j*8+subNr*2] */
	i = (i + (i << 2));
	k = startPos[(j << 3) + (subNr << 1)];
	pos[0] = i + k;
	index = index >> 3;
	i = index & 7;

	/* pos1 =i*5+startPos[j*8+subNr*2+1] */
	i = (i + (i << 2));
	k = startPos[((j << 3) + (subNr << 1)) + 1];
	pos[1] = (int16_t) (i + k);

	/* decode the signs  and build the codeword */
	memset(cod, 0, L_SUBFR << 2);

	for (j = 0; j < 2; j++) {
		i = sign & 1;
		sign = sign >> 1;

		if (i != 0) {
			cod[pos[j]] = 8191;	/* +1.0 */
		} else {
			cod[pos[j]] = -8192;	/* -1.0 */
		}
	}
	return;
}



/*
 * gmed_n
 *
 *
 * Parameters:
 *    ind               I: values
 *    n                 I: The number of gains (odd)
 *
 * Function:
 *    Calculates N-point median.
 *
 * Returns:
 *    index of the median value
 */
static int32_t gmed_n(int32_t ind[], int32_t n)
{
	int32_t tmp[NMAX], tmp2[NMAX];
	int32_t max, medianIndex, i, j, ix = 0;

	for (i = 0; i < n; i++) {
		tmp2[i] = ind[i];
	}

	for (i = 0; i < n; i++) {
		max = -32767;

		for (j = 0; j < n; j++) {
			if (tmp2[j] >= max) {
				max = tmp2[j];
				ix = j;
			}
		}
		tmp2[ix] = -32768;
		tmp[i] = ix;
	}
	medianIndex = tmp[(n >> 1)];
	return (ind[medianIndex]);
}


/*
 * ec_gain_pitch_update
 *
 *
 * Parameters:
 *    st->prev_gp       B: previous pitch gain
 *    st->past_gain_pit O: past gain
 *    st->pbuf          B: past gain buffer
 *    bfi               I: bad frame indicator
 *    prev_bf           I: previous frame was bad
 *    gain_pitch        B: pitch gain
 *
 * Function:
 *    Update the pitch gain concealment state
 *    Limit gain_pitch if the previous frame was bad
 *
 * Returns:
 *    gain
 */
static void ec_gain_pitch_update(ec_gain_pitchState * st, int32_t bfi,
				 int32_t prev_bf, int32_t * gain_pitch)
{
	if (bfi == 0) {
		if (prev_bf != 0) {
			if (*gain_pitch > st->prev_gp) {
				*gain_pitch = st->prev_gp;
			}
		}
		st->prev_gp = *gain_pitch;
	}
	st->past_gain_pit = *gain_pitch;

	/* if (st->past_gain_pit > 1.0) */
	if (st->past_gain_pit > 16384) {
		st->past_gain_pit = 16384;
	}
	st->pbuf[0] = st->pbuf[1];
	st->pbuf[1] = st->pbuf[2];
	st->pbuf[2] = st->pbuf[3];
	st->pbuf[3] = st->pbuf[4];
	st->pbuf[4] = st->past_gain_pit;
}

/*
 * gc_pred (366)
 *
 *
 * Parameters:
 *    st->past_qua_en         I: MA predictor
 *    st->past_qua_en_MR122   I: MA predictor MR122
 *    mode                    I: AMR mode
 *    code                    I: innovative codebook vector
 *    exp_gcode0              O: predicted gain factor (exponent)
 *    frac_gcode0             O: predicted gain factor (fraction)
 *    exp_en                  I: innovation energy (MR795) (exponent)
 *    frac_en                 I: innovation energy (MR795) (fraction)
 *
 * Function:
 *    MA prediction of the innovation energy
 *
 *    Mean removed innovation energy (dB) in subframe n
 *                          N-1
 *    E(n) = 10*log(gc*gc * SUM[(code(i) * code(i)]/N) - EMean
 *                          i=0
 *    N=40
 *
 *    Mean innovation energy (dB)
 *                   N-1
 *    Ei(n) = 10*log(SUM[(code(i) * code(i)]/N)
 *                   i=0
 *
 *    Predicted energy
 *              4
 *    Ep(n) = SUM[b(i) * R(n-i)]
 *            i=1
 *    b = [0.68 0.58 0.34 0.19]
 *    R(k) is quantified prediction error at subframe k
 *
 *    E_Mean = 36 dB (MR122)
 *
 *    Predicted gain gc is found by
 *
 *    gc = POW[10, 0.05 * (Ep(n) + EMean - Ei)]
 *
 * Returns:
 *    void
 */
static void gc_pred(gc_predState * st, enum Mode mode, int32_t * code, int32_t *
		    exp_gcode0, int32_t * frac_gcode0, int32_t * exp_en,
		    int32_t * frac_en)
{
	int32_t exp, frac, ener_code = 0, i = 0;

	/* energy of code:
	 * ener_code = sum(code[i]^2)
	 */
	while (i < L_SUBFR) {
		ener_code += code[i] * code[i];
		i++;
	}

	if ((0x3fffffff <= ener_code) | (ener_code < 0))
		ener_code = MAX_32;
	else
		ener_code <<= 1;
       /*
	if (mode == MR122) {
		int32_t ener;

		// ener_code = ener_code / lcode; lcode = 40; 1/40 = 26214 Q20
		ener_code = ((ener_code + 0x00008000L) >> 16) * 52428;

		// Q9  * Q20 -> Q30
		// energy of code:
		// * ener_code(Q17) = 10 * Log10(energy) / constant
		// *                = 1/2 * Log2(energy)
		// * constant = 20*Log10(2)

		// ener_code = 1/2 * Log2(ener_code); Note: Log2=log2+30
		Log2(ener_code, &exp, &frac);
		ener_code = ((exp - 30) << 16) + (frac << 1);

		// Q16 for log(), ->Q17 for 1/2 log()

	       //	 * predicted energy:
	       //	 * ener(Q24) = (Emean + sum{pred[i]*pastEn[i]})/constant
	       //	 *           = MEAN_ENER + sum(pred[i]*past_qua_en[i])
	       //	 * constant = 20*Log10(2)

		ener = 0;
		i = 0;

		while (i < 4) {
			ener += st->past_qua_en_MR122[i] * pred_MR122[i];
			i++;
		}
		ener <<= 1;
		ener += MEAN_ENER_MR122;

		//
		// * predicted codebook gain

		// * gc0 = Pow10( (ener*constant - ener_code*constant) / 20 )
	       //	 *     = Pow2(ener-ener_code)
	       //	 *     = Pow2(int(d)+frac(d))

		ener = (ener - ener_code) >> 1;	//Q16 
		*exp_gcode0 = ener >> 16;
		*frac_gcode0 = (ener >> 1) - (*exp_gcode0 << 15);
	}

	//all modes except 12.2 
	else

        */
        {


		int32_t tmp, gcode0;
		int exp_code;

		/*
		 * Compute: meansEner - 10log10(ener_code/ LSufr)
		 */
		exp_code = 0;
		if (ener_code != 0) {
			while (!(ener_code & 0x40000000)) {
				exp_code++;
				ener_code = ener_code << 1;
			}
		}

		/* Log2 = log2 + 27 */
		Log2_norm(ener_code, exp_code, &exp, &frac);

		/* fact = 10/log2(10) = 3.01 = 24660 Q13 */
		/* Q0.Q15 * Q13 -> Q14 */
		tmp = (exp * (-49320)) + (((frac * (-24660)) >> 15) << 1);

		/*
		 * tmp = meansEner - 10log10(ener_code/L_SUBFR)
		 *       = meansEner - 10log10(ener_code) + 10log10(L_SUBFR)
		 *       = K - fact * Log2(ener_code)
		 *     = K - fact * log2(ener_code) - fact*27
		 *
		 *   ==> K = meansEner + fact*27 + 10log10(L_SUBFR)
		 *
		 *   meansEner =       33    =  540672    Q14  (MR475, MR515, MR59)
		 *   meansEner =       28.75 =  471040    Q14  (MR67)
		 *   meansEner =       30    =  491520    Q14  (MR74)
		 *   meansEner =       36    =  589824    Q14  (MR795)
		 *   meansEner =       33    =  540672    Q14  (MR102)
		 *   10log10(L_SUBFR) = 16.02 =  262481.51 Q14
		 *   fact * 27                = 1331640    Q14
		 *   -----------------------------------------
		 *   (MR475, MR515, MR59)   K = 2134793.51 Q14 ~= 16678 * 64 * 2
		 *   (MR67)                 K = 2065161.51 Q14 ~= 32268 * 32 * 2
		 *   (MR74)                 K = 2085641.51 Q14 ~= 32588 * 32 * 2
		 *   (MR795)                K = 2183945.51 Q14 ~= 17062 * 64 * 2
		 *   (MR102)                K = 2134793.51 Q14 ~= 16678 * 64 * 2
		 */

             /*
		if (mode == MR102) {
			// mean = 33 dB
			tmp += 2134784;	// Q14 
		} else if (mode == MR795) {
			// mean = 36 dB
			tmp += 2183936;	// Q14


		       //	 * ener_code  = <xn xn> * 2^27*2^exp_code
		       //	 * frac_en    = ener_code / 2^16
		       //	 *            = <xn xn> * 2^11*2^exp_code
			// * <xn xn>    = <xn xn>*2^11*2^exp * 2^exp_en
		       //	 *           := frac_en            * 2^exp_en
		       //	 *
		       //	 * ==> exp_en = -11-exp_code;

			*frac_en = ener_code >> 16;
			*exp_en = -11 - exp_code;
		} else if (mode == MR74) {
			// mean = 30 dB
			tmp += 2085632;	// Q14
		} else if (mode == MR67) {
			// mean = 28.75 dB
			tmp += 2065152;	//Q14
		} else

              */
                {	/* MR59, MR515, MR475 */

			/* mean = 33 dB */
			tmp += 2134784;	/* Q14 */
		}

		/*
		 * Compute gcode0
		 * = Sum(i=0,3) pred[i]*past_qua_en[i] - ener_code + meanEner
		 */
		tmp = tmp << 9;	/* Q23 */

		/* Q13 * Q10 -> Q23 */
		i = 0;

		while (i < 4) {
			tmp += pred[i] * st->past_qua_en[i];
			i++;
		}
		gcode0 = tmp >> 15;	/* Q8  */

		/*
		 * gcode0 = pow(10.0, gcode0/20)
		 *        = pow(2, 3.3219*gcode0/20)
		 *        = pow(2, 0.166*gcode0)
		 */
		/* 5439 Q15 = 0.165985                                        */
		/* (correct: 1/(20*log10(2)) 0.166096 = 5443 Q15)             */
		/* For IS641 bitexactness */
	       //	if (mode == MR74) {
	       //		/*/ Q8 * Q15 -> Q24
		//	tmp = gcode0 * 10878;
	       //	} else

                {
			/* Q8 * Q15 -> Q24 */
			tmp = gcode0 * 10886;
		}
		tmp = tmp >> 9;	/* -> Q15 */

		/* -> Q0.Q15 */
		*exp_gcode0 = tmp >> 15;
		*frac_gcode0 = tmp - (*exp_gcode0 * 32768);
	}
}

/*
 * gc_pred_update
 *
 *
 * Parameters:
 *    st->past_qua_en         B: MA predictor
 *    st->past_qua_en_MR122   B: MA predictor MR122
 *    qua_ener_MR122          I: quantized energy for update (log2(quaErr))
 *    qua_ener                I: quantized energy for update (20*log10(quaErr))
 *
 * Function:
 *    Update MA predictor with last quantized energy
 *
 * Returns:
 *    void
 */
static void gc_pred_update(gc_predState * st, int32_t qua_ener_MR122,
			   int32_t qua_ener)
{
	int32_t i;

	for (i = 3; i > 0; i--) {
		st->past_qua_en[i] = st->past_qua_en[i - 1];
		st->past_qua_en_MR122[i] = st->past_qua_en_MR122[i - 1];
	}
	st->past_qua_en_MR122[0] = qua_ener_MR122;	/* log2 (quaErr), Q10 */
	st->past_qua_en[0] = qua_ener;	/* 20*log10(quaErr), Q10 */
}

/*
 * Dec_gain
 *
 *
 * Parameters:
 *    pred_state->past_qua_en       B: MA predictor
 *    pred_state->past_qua_en_MR122 B: MA predictor MR122
 *    mode                          I: AMR mode
 *    index                         I: index of quantization
 *    code                          I: Innovative vector
 *    evenSubfr                     I: Flag for even subframes
 *    gain_pit                      O: Pitch gain
 *    gain_cod                      O: Code gain
 *
 * Function:
 *    Decode the pitch and codebook gains
 *
 * Returns:
 *    void
 */
static void Dec_gain(gc_predState * pred_state, enum Mode mode, int32_t index,
		     int32_t code[], int32_t evenSubfr, int32_t * gain_pit,
		     int32_t * gain_cod)
{
	int32_t frac, gcode0, exp, qua_ener, qua_ener_MR122, g_code, tmp;
	const short *p;

	/* Read the quantized gains (table depends on mode) */
	index = index << 2;

        /*
	if ((mode == MR102) || (mode == MR74) || (mode == MR67)) {
		p = &table_gain_highrates[index];
		*gain_pit = *p++;
		g_code = *p++;
		qua_ener_MR122 = *p++;
		qua_ener = *p;
	} else
        */


        {
	       //	if (mode == MR475) {
			index = index + ((1 - evenSubfr) << 1);
			p = &table_gain_MR475[index];
			*gain_pit = *p++;
			g_code = *p++;

			/*
			 * calculate predictor update values (not stored in 4.75
			 * quantizer table to save space):
			 *   qua_ener       = log2(g)
			 *   qua_ener_MR122 = 20*log10(g)
			 */
			/* Log2(x Q12) = log2(x) + 12 */
			Log2(g_code, &exp, &frac);
			exp = exp - 12;
			tmp = frac >> 5;

			if ((frac & ((int16_t) 1 << 4)) != 0) {
				tmp++;
			}
			qua_ener_MR122 = tmp + (exp << 10);

			/* 24660 Q12 ~= 6.0206 = 20*log10(2) */
			tmp = exp * 49320;
			tmp += (((frac * 24660) >> 15) << 1);

			/* Q12 * Q0 = Q13 -> Q10 */
			qua_ener = ((tmp << 13) + 0x00008000L) >> 16;
              /*

		}
                else {
			p = &table_gain_lowrates[index];
			*gain_pit = *p++;
			g_code = *p++;
			qua_ener_MR122 = *p++;
			qua_ener = *p;
		}
                */
	}

	/*
	 * predict codebook gain
	 * gc0 = Pow2(int(d)+frac(d))
	 *     = 2^exp + 2^frac
	 * gcode0 (Q14) = 2^14*2^frac = gc0 * 2^(14-exp)
	 */
	gc_pred(pred_state, mode, code, &exp, &frac, NULL, NULL);
	gcode0 = Pow2(14, frac);

	/*
	 * read quantized gains, update table of past quantized energies
	 * st->past_qua_en(Q10) = 20 * Log10(gFac) / constant
	 *                      = Log2(gFac)
	 *                      = qua_ener
	 * constant = 20*Log10(2)
	 */
	if (exp < 11) {
		*gain_cod = (g_code * gcode0) >> (25 - exp);
	} else {
		tmp = ((g_code * gcode0) << (exp - 9));

		if ((tmp >> (exp - 9)) != (g_code * gcode0)) {
			*gain_cod = 0x7FFF;
		} else {
			*gain_cod = tmp >> 16;
		}
	}

	/* update table of past quantized energies */
	gc_pred_update(pred_state, qua_ener_MR122, qua_ener);
	return;
}



/*
 * ec_gain_code_update
 *
 *
 * Parameters:
 *    st->gbuf             B: last five gains
 *    st->past_gain_code   O: past gain
 *    st->prev_gc          B  previous gain
 *    bfi                  I: bad indicator
 *    prev_bf              I: previous frame bad indicator
 *    gain_code            O: decoded innovation gain
 *
 * Function:
 *    Update the codebook gain concealment state
 *
 * Returns:
 *    void
 */
static void ec_gain_code_update(ec_gain_codeState * st, int16_t bfi,
				int16_t prev_bf, int32_t * gain_code)
{
	/* limit gain_code by previous good gain if previous frame was bad */
	if (bfi == 0) {
		if (prev_bf != 0) {
			if (*gain_code > st->prev_gc) {
				*gain_code = st->prev_gc;
			}
		}
		st->prev_gc = *gain_code;
	}

	/* update EC states: previous gain, gain buffer */
	st->past_gain_code = *gain_code;
	st->gbuf[0] = st->gbuf[1];
	st->gbuf[1] = st->gbuf[2];
	st->gbuf[2] = st->gbuf[3];
	st->gbuf[3] = st->gbuf[4];
	st->gbuf[4] = *gain_code;
	return;
}


 
/*
 * Int_lsf
 *
 *
 * Parameters:
 *    lsf_old           I: LSF vector at the 4th subframe of past frame
 *    lsf_new           I: LSF vector at the 4th subframe of present frame
 *    i_subfr           I: current subframe
 *    lsf_out           O: interpolated LSF parameters for current subframe
 *
 * Function:
 *    Interpolates the LSFs for selected subframe
 *
 *    The LSFs are interpolated at the 1st, 2nd and 3rd
 *    ubframe and only forwarded at the 4th subframe.
 *
 *    sf1:  3/4 F0 + 1/4 F1
 *    sf2:  1/2 F0 + 1/2 F1
 *    sf3:  1/4 F0 + 3/4 F1
 *    sf4:  F1
 *
 * Returns:
 *    void
 */
static void Int_lsf(int32_t lsf_old[], int32_t lsf_new[], int i_subfr, int32_t
		    lsf_out[])
{
	int32_t i;

	switch (i_subfr) {
	case 0:
		for (i = 0; i < 10; i++) {
			lsf_out[i] =
			    lsf_old[i] - (lsf_old[i] >> 2) + (lsf_new[i] >> 2);
		}
		break;

	case 40:
		for (i = 0; i < 10; i++) {
			lsf_out[i] = (lsf_old[i] >> 1) + (lsf_new[i] >> 1);
		}
		break;

	case 80:
		for (i = 0; i < 10; i++) {
			lsf_out[i] = (lsf_old[i] >> 2) - (lsf_new[i] >> 2) +
			    lsf_new[i];
		}
		break;

	case 120:
		memcpy(lsf_out, lsf_new, M << 2);
		break;
	}
}

/*
 * Cb_gain_average
 *
 *
 * Parameters:
 *    st->cbGainHistory B: codebook gain history
 *    st->hangCount     B: hangover counter
 *    mode              I: AMR mode
 *    gain_code         I: codebook gain
 *    lsp               I: The LSP for the current frame
 *    lspAver           I: The average of LSP for 8 frames
 *    bfi               I: bad frame indication
 *    prev_bf           I: previous bad frame indication
 *    pdfi              I: potential degraded bad frame indication
 *    prev_pdf          I: previous potential degraded bad frame indication
 *    inBackgroundNoise I: background noise decision
 *    voicedHangover    I: number of frames after last voiced frame
 *
 * Function:
 *    The mixed codebook gain, used to make codebook gain more smooth in background
 *
 *
 * Returns:
 *    void
 */
static int32_t Cb_gain_average(Cb_gain_averageState * st, enum Mode mode, int32_t
			      gain_code, int32_t lsp[], int32_t lspAver[],
			      int16_t bfi, int16_t prev_bf, int16_t pdfi,
			      int16_t prev_pdf, int32_t inBackgroundNoise,
			      int32_t voicedHangover)
{
	int32_t tmp[M];
	int32_t i, cbGainMix, tmp_diff, bgMix, cbGainMean, sum, diff, tmp1, tmp2;
	int shift1, shift2, shift;

	/* set correct cbGainMix for MR74, MR795, MR122 */
	cbGainMix = gain_code;

	/*
	 * Store list of CB gain needed in the CB gain averaging                                           *
	 */
	st->cbGainHistory[0] = st->cbGainHistory[1];
	st->cbGainHistory[1] = st->cbGainHistory[2];
	st->cbGainHistory[2] = st->cbGainHistory[3];
	st->cbGainHistory[3] = st->cbGainHistory[4];
	st->cbGainHistory[4] = st->cbGainHistory[5];
	st->cbGainHistory[5] = st->cbGainHistory[6];
	st->cbGainHistory[6] = gain_code;

	/* compute lsp difference */
	for (i = 0; i < M; i++) {
		tmp1 = labs(lspAver[i] - lsp[i]);
		shift1 = 0;
		if (tmp1 != 0) {
			while (!(tmp1 & 0x2000)) {
				shift1++;
				tmp1 = tmp1 << 1;
			}
		}
		tmp2 = lspAver[i];
		shift2 = 0;
		if (tmp2 != 0) {
			while (!(tmp2 & 0x4000)) {
				shift2++;
				tmp2 = tmp2 << 1;
			}
		} else
			tmp2 = 0x4000;
		tmp[i] = (tmp1 << 15) / tmp2;
		shift = 2 + shift1 - shift2;

		if (shift >= 0) {
			tmp[i] = tmp[i] >> shift;
		} else {
			tmp[i] = tmp[i] << -(shift);
		}
	}
	diff =
	    *tmp + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] +
	    tmp[7] + tmp[8] + tmp[9];

	/* saturate */
	if (diff > 32767) {
		diff = 32767;
	}

	/* Compute hangover */
	st->hangVar += 1;

	if (diff <= 5325) {
		st->hangVar = 0;
	}

	if (st->hangVar > 10) {
		/* Speech period, reset hangover variable */
		st->hangCount = 0;
	}

	/* Compute mix constant (bgMix) */
	bgMix = 8192;

	/* MR475, MR515, MR59, MR67, MR102 */
	//if ((mode <= MR67) | (mode == MR102)) {
		/* disable mix if too short time since */
		if ((st->hangCount >= 40) & (diff <= 5325)) {	/* 0.65 in Q13 */
			/* if errors and presumed noise make smoothing probability stronger */
			if (((((pdfi != 0) & (prev_pdf != 0)) | (bfi !=
								 0) | (prev_bf
								       !=
								       0)) &
			     ((voicedHangover > 1)) & (inBackgroundNoise !=
						       0) & (mode < MR67))) {
				/* bgMix = min(0.25, max(0.0, diff-0.55)) / 0.25; */
				tmp_diff = diff - 4506;	/* 0.55 in Q13 */

				/* max(0.0, diff-0.55) */
				tmp1 = 0;

				if (tmp_diff > 0) {
					tmp1 = tmp_diff;
				}

				/* min(0.25, tmp1) */
				if (2048 >= tmp1) {
					bgMix = tmp1 << 2;
				}
			} else {
				/* bgMix = min(0.25, max(0.0, diff-0.40)) / 0.25; */
				tmp_diff = diff - 3277;	/* 0.4 in Q13 */

				/* max(0.0, diff-0.40) */
				tmp1 = 0;

				if (tmp_diff > 0) {
					tmp1 = tmp_diff;
				}

				/* min(0.25, tmp1) */
				if (2048 >= tmp1) {
					bgMix = tmp1 << 2;
				}
			}
		}

		/*
		 * Smoothen the cb gain trajectory
		 * smoothing depends on mix constant bgMix
		 */
		sum =
		    st->cbGainHistory[2] + st->cbGainHistory[3] +
		    st->cbGainHistory[4] + st->cbGainHistory[5] +
		    st->cbGainHistory[6];

		if (sum > 163822) {
			cbGainMean = 32767;
		} else {
			cbGainMean = (3277 * sum + 0x00002000L) >> 14;	/* Q1 */
		}

		/* more smoothing in error and bg noise (NB no DFI used  here) */
		if (((bfi != 0) | (prev_bf != 0)) & (inBackgroundNoise !=
						     0) & (mode < MR67)) {
			sum =
			    9362 * (st->cbGainHistory[0] +
				    st->cbGainHistory[1] +
				    st->cbGainHistory[2] +
				    st->cbGainHistory[3] +
				    st->cbGainHistory[4] +
				    st->cbGainHistory[5] +
				    st->cbGainHistory[6]);
			cbGainMean = (sum + 0x00008000L) >> 16;	/* Q1 */
		}

		/* cbGainMix = bgMix*cbGainMix + (1-bgMix)*cbGainMean; */
		sum = bgMix * cbGainMix;	/* sum in Q14 */
		sum += cbGainMean << 13;
		sum -= bgMix * cbGainMean;
		cbGainMix = (sum + 0x00001000L) >> 13;

		/* Q1 */
       //	}
	st->hangCount += 1;
	if (st->hangCount & 0x80000000)
		st->hangCount = 40;
	return cbGainMix;
}

/*
 * ph_disp
 *
 *
 * Parameters:
 *    state->gainMem    B: LTP gain memory
 *    state->prevCbGain B: Codebook gain memory
 *    mode              I: AMR mode
 *    x                 B: LTP excitation signal -> total excitation signal
 *    cbGain            I: Codebook gain
 *    ltpGain           I: LTP gain
 *    inno              B: Innovation vector
 *    pitch_fac         I: pitch factor used to scale the LTP excitation
 *    tmp_shift         I: shift factor applied to sum of scaled LTP ex & innov.
 *                         before rounding
 *
 * Function:
 *    Adaptive phase dispersion; forming of total excitation
 *
 *
 * Returns:
 *    void
 */
static void ph_disp(ph_dispState * state, enum Mode mode, int32_t x[],
		    int32_t cbGain, int32_t ltpGain, int32_t inno[],
		    int32_t pitch_fac, int32_t tmp_shift)
{
	int32_t inno_sav[L_SUBFR], ps_poss[L_SUBFR];
	int32_t i, i1, impNr, temp1, temp2, j, nze, nPulse, ppos;
	const short *ph_imp;	/* Pointer to phase dispersion filter */

	/* Update LTP gain memory */
	state->gainMem[4] = state->gainMem[3];
	state->gainMem[3] = state->gainMem[2];
	state->gainMem[2] = state->gainMem[1];
	state->gainMem[1] = state->gainMem[0];
	state->gainMem[0] = ltpGain;

	/* basic adaption of phase dispersion */
	/* no dispersion */
	impNr = 2;

	/* if (ltpGain < 0.9) */
	if (ltpGain < PHDTHR2LTP) {
		/* maximum dispersion */
		impNr = 0;

		/* if (ltpGain > 0.6 */
		if (ltpGain > PHDTHR1LTP) {
			/* medium dispersion */
			impNr = 1;
		}
	}

	/* onset indicator */
	/* onset = (cbGain  > onFact * cbGainMem[0]) */
	temp1 = ((state->prevCbGain * ONFACTPLUS1) + 0x1000) >> 13;

	if (cbGain > temp1) {
		state->onset = ONLENGTH;
	} else {
		if (state->onset > 0) {
			state->onset--;
		}
	}

	/*
	 * if not onset, check ltpGain buffer and use max phase dispersion if
	 * half or more of the ltpGain-parameters say so
	 */
	if (state->onset == 0) {
		/* Check LTP gain memory and set filter accordingly */
		i1 = 0;

		for (i = 0; i < PHDGAINMEMSIZE; i++) {
			if (state->gainMem[i] < PHDTHR1LTP) {
				i1++;
			}
		}

		if (i1 > 2) {
			impNr = 0;
		}
	}

	/* Restrict decrease in phase dispersion to one step if not onset */
	if ((impNr > (state->prevState + 1)) & (state->onset == 0)) {
		impNr--;
	}

	/* if onset, use one step less phase dispersion */
	if ((impNr < 2) & (state->onset > 0)) {
		impNr++;
	}

	/* disable for very low levels */
	if (cbGain < 10) {
		impNr = 2;
	}

	if (state->lockFull == 1) {
		impNr = 0;
	}

	/* update static memory */
	state->prevState = impNr;
	state->prevCbGain = cbGain;

	/*
	 * do phase dispersion for all modes but 12.2 and 7.4;
	 * don't modify the innovation if impNr >=2 (= no phase disp)
	 */
	if ((mode != MR122) & (mode != MR102) & (mode != MR74) & (impNr < 2)
	    ) {
		/*
		 * track pulse positions, save innovation,
		 * and initialize new innovation
		 */
		nze = 0;

		for (i = 0; i < L_SUBFR; i++) {
			if (inno[i] != 0) {
				ps_poss[nze] = i;
				nze++;
			}
		}
		memcpy(inno_sav, inno, L_SUBFR << 2);
		memset(inno, 0, L_SUBFR << 2);

		/* Choose filter corresponding to codec mode and dispersion criterium */
		ph_imp = ph_imp_mid;

		if (impNr == 0) {
			ph_imp = ph_imp_low;
		}

	       //	if (mode == MR795) {
	       //		ph_imp = ph_imp_mid_MR795;
               //
	       //		if (impNr == 0) {
	       //			ph_imp = ph_imp_low_MR795;
	       //		}
	       //	}

		/* Do phase dispersion of innovation */
		for (nPulse = 0; nPulse < nze; nPulse++) {
			ppos = ps_poss[nPulse];

			/* circular convolution with impulse response */
			j = 0;

			for (i = ppos; i < L_SUBFR; i++) {
				/* inno[i1] += inno_sav[ppos] * ph_imp[i1-ppos] */
				temp1 = (inno_sav[ppos] * ph_imp[j++]) >> 15;
				inno[i] = inno[i] + temp1;
			}

			for (i = 0; i < ppos; i++) {
				/* inno[i] += inno_sav[ppos] * ph_imp[L_SUBFR-ppos+i] */
				temp1 = (inno_sav[ppos] * ph_imp[j++]) >> 15;
				inno[i] = inno[i] + temp1;
			}
		}
	}

	/*
	 * compute total excitation for synthesis part of decoder
	 * (using modified innovation if phase dispersion is active)
	 */
	for (i = 0; i < L_SUBFR; i++) {
		/* x[i] = gain_pit*x[i] + cbGain*code[i]; */
		temp1 = x[i] * pitch_fac + inno[i] * cbGain;
		temp2 = temp1 << tmp_shift;
		x[i] = (temp2 + 0x4000) >> 15;
		if (labs(x[i]) > 32767) {
			if ((temp1 ^ temp2) & 0x80000000) {
				x[i] = (temp1 & 0x80000000) ? -32768 : 32767;
			} else {
				x[i] = (temp2 & 0x80000000) ? -32768 : 32767;
			}
		}
	}
	return;
}

/*
 * sqrt_l_exp
 *
 *
 * Parameters:
 *    x                 I: input value
 *    exp               O: right shift to be applied to result
 *
 * Function:
 *    Sqrt with exponent value.
 *
 *    y = sqrt(x)
 *    x = f * 2^-e,   0.5 <= f < 1   (normalization)
 *    y = sqrt(f) * 2^(-e/2)
 *
 *    a) e = 2k   --> y = sqrt(f)   * 2^-k
 *       (k = e div 2, 0.707 <= sqrt(f) < 1)
 *    b) e = 2k+1 --> y = sqrt(f/2) * 2^-k
 *       (k = e div 2, 0.5 <= sqrt(f/2) < 0.707)
 *
 *
 * Returns:
 *    y                 output value
 */
static int32_t sqrt_l_exp(int32_t x, int32_t * exp)
{
	int32_t y, a, i, tmp;
	int e;

	if (x <= (int32_t) 0) {
		*exp = 0;
		return (int32_t) 0;
	}
	e = 0;
	if (x != 0) {
		tmp = x;
		while (!(tmp & 0x40000000)) {
			e++;
			tmp = tmp << 1;
		}
	}
	e = e & 0xFFFE;
	x = (x << e);
	*exp = (int16_t) e;
	x = (x >> 9);
	i = (int16_t) (x >> 16);
	x = (x >> 1);
	a = x & (int16_t) 0x7fff;
	i = (i - 16);
	y = (sqrt_table[i] << 16);
	tmp = (sqrt_table[i] - sqrt_table[i + 1]);
	y -= (tmp * a) << 1;
	return (y);
}

/*
 * Ex_ctrl
 *
 *
 * Parameters:
 *    excitation        B: Current subframe excitation
 *    excEnergy         I: Exc. Energy, sqrt(totEx*totEx)
 *    exEnergyHist      I: History of subframe energies
 *    voicedHangover    I: number of frames after last voiced frame
 *    prevBFI           I: Set i previous bad frame indicators
 *    carefulFlag       I: Restrict dymamic in scaling
 *
 * Function:
 *    Charaterice synthesis speech and detect background noise
 *
 * Returns:
 *    background noise decision; 0 = no bgn, 1 = bgn
 */
static int16_t Ex_ctrl(int32_t excitation[], int32_t excEnergy, int32_t
		      exEnergyHist[], int32_t voicedHangover, int16_t prevBFI,
		      int16_t carefulFlag)
{
	int32_t i, testEnergy, scaleFactor, avgEnergy, prevEnergy, T0;
	int exp;

	/* get target level */
	avgEnergy = gmed_n(exEnergyHist, 9);
	prevEnergy = (exEnergyHist[7] + exEnergyHist[8]) >> 1;

	if (exEnergyHist[8] < prevEnergy) {
		prevEnergy = exEnergyHist[8];
	}

	/* upscaling to avoid too rapid energy rises  for some cases */
	if ((excEnergy < avgEnergy) & (excEnergy > 5)) {
		/* testEnergy = 4*prevEnergy; */
		testEnergy = prevEnergy << 2;

		if ((voicedHangover < 7) || prevBFI != 0) {
			/* testEnergy = 3*prevEnergy */
			testEnergy = testEnergy - prevEnergy;
		}

		if (avgEnergy > testEnergy) {
			avgEnergy = testEnergy;
		}

		/* scaleFactor=avgEnergy/excEnergy in Q0 */
		exp = 0;
		if (excEnergy != 0) {
			while (!(excEnergy & 0x4000)) {
				exp++;
				excEnergy = excEnergy << 1;
			}
		} else
			excEnergy = 0x4000;
		excEnergy = 536838144 / excEnergy;
		T0 = (avgEnergy * excEnergy) << 1;
		T0 = (T0 >> (20 - exp));

		if (T0 > 32767) {
			/* saturate  */
			T0 = 32767;
		}
		scaleFactor = T0;

		/* test if scaleFactor > 3.0 */
		if ((carefulFlag != 0) & (scaleFactor > 3072)) {
			scaleFactor = 3072;
		}

		/* scale the excitation by scaleFactor */
		for (i = 0; i < L_SUBFR; i++) {
			T0 = (scaleFactor * excitation[i]) << 1;
			T0 = (T0 >> 11);
			excitation[i] = T0;
		}
	}
	return 0;
}

/*
 * Inv_sqrt
 *
 *
 * Parameters:
 *    x                 I: input value
 *
 * Function:
 *    1/sqrt(x)
 *
 * Returns:
 *    y                 1/sqrt(x)
 */
static int32_t Inv_sqrt(int32_t x)
{
	int i, a, tmp, exp;
	int32_t y;

	if (x <= (int32_t) 0)
		return ((int32_t) 0x3fffffffL);
	exp = 0;
	while (!(x & 0x40000000)) {
		exp++;
		x = x << 1;
	}

	/* x is normalized */
	exp = (30 - exp);

	/* If exponent even -> shift right */
	if ((exp & 1) == 0) {
		x = (x >> 1);
	}
	exp = (exp >> 1);
	exp = (exp + 1);
	x = (x >> 9);

	/* Extract b25-b31 */
	i = (int16_t) (x >> 16);

	/* Extract b10-b24 */
	x = (x >> 1);
	a = x & (int16_t) 0x7fff;
	i = (i - 16);

	/* table[i] << 16 */
	y = inv_sqrt_table[i] << 16;

	/* table[i] - table[i+1]) */
	tmp = (inv_sqrt_table[i] - inv_sqrt_table[i + 1]);

	/* y -= tmp*a*2 */
	y -= (tmp * a) << 1;

	/* denormalization */
	y = (y >> exp);
	return (y);
}

/*
 * energy_old
 *
 *
 * Parameters:
 *    in                I: input value
 *
 * Function:
 *    Energy of signal
 *
 * Returns:
 *    Energy
 */
static int32_t energy_old(int32_t in[])
{
	int32_t temp, i, sum = 0;

	for (i = 0; i < L_SUBFR; i += 8) {
		temp = in[i] >> 2;
		sum += temp * temp;
		temp = in[i + 1] >> 2;
		sum += temp * temp;
		temp = in[i + 2] >> 2;
		sum += temp * temp;
		temp = in[i + 3] >> 2;
		sum += temp * temp;
		temp = in[i + 4] >> 2;
		sum += temp * temp;
		temp = in[i + 5] >> 2;
		sum += temp * temp;
		temp = in[i + 6] >> 2;
		sum += temp * temp;
		temp = in[i + 7] >> 2;
		sum += temp * temp;
	}

	if (sum & 0xC0000000) {
		return 0x7FFFFFFF;
	}
	return (sum << 1);
}

/*
 * energy_new
 *
 *
 * Parameters:
 *    in                I: input value
 *
 * Function:
 *    Energy of signal
 *
 * Returns:
 *    Energy
 */
static int32_t energy_new(int32_t in[])
{
	int32_t i, s = 0, overflow = 0;

	s += in[0] * in[0];
	for (i = 1; i < L_SUBFR; i += 3) {
		s += in[i] * in[i];
		s += in[i + 1] * in[i + 1];
		s += in[i + 2] * in[i + 2];

		if (s & 0xC0000000) {
			overflow = 1;
			break;
		}
	}

	/* check for overflow */
	if (overflow) {
		s = energy_old(in);
	} else {
		s = (s >> 3);
	}
	return s;
}

/*
 * agc2
 *
 *
 * Parameters:
 *    sig_in            I: Post_Filter input signal
 *    sig_out           B: Post_Filter output signal
 *
 * Function:
 *    Scales the excitation on a subframe basis
 *
 * Returns:
 *    Energy
 */
static void agc2(int32_t * sig_in, int32_t * sig_out)
{
	int32_t s;
	int i, exp;
	int16_t gain_in, gain_out, g0;

	/* calculate gain_out with exponent */
	s = energy_new(sig_out);

	if (s == 0) {
		return;
	}
	exp = 0;
	while (!(s & 0x20000000)) {
		exp++;
		s = s << 1;
	}

	gain_out = (int16_t) ((s + 0x00008000L) >> 16);

	/* calculate gain_in with exponent */
	s = energy_new(sig_in);

	if (s == 0) {
		g0 = 0;
	} else {
		i = 0;
		while (!(s & 0x40000000)) {
			i++;
			s = s << 1;
		}

		if (s < 0x7fff7fff)
			gain_in = (int16_t) ((s + 0x00008000L) >> 16);
		else
			gain_in = 32767;
		exp = (exp - i);

		/*
		 * g0 = sqrt(gain_in/gain_out);
		 */
		/* s = gain_out / gain_in */
		s = (gain_out << 15) / gain_in;
		s = (s << 7);

		if (exp > 0)
			s = (s >> exp);
		else
			s = (s << (-exp));
		s = Inv_sqrt(s);
		g0 = (int16_t) (((s << 9) + 0x00008000L) >> 16);
	}

	/* sig_out(n) = gain(n) * sig_out(n) */
	for (i = 0; i < L_SUBFR; i++) {
		sig_out[i] = (sig_out[i] * g0) >> 12;
	}
	return;
}


/*
 * Decoder_amr
 *
 *
 * Parameters:
 *    st                B: State variables
 *    mode              I: AMR mode
 *    parm              I: vector of synthesis parameters
 *    frame_type        I: received frame type
 *    synth             O: synthesis speech
 *    A_t               O: decoded LP filter in 4 subframes
 *
 * Function:
 *    Speech decoder routine
 *
 * Returns:
 *    void
 */
static void Decoder_amr(Decoder_amrState * st, enum Mode mode, int16_t parm[],
			enum RXFrameType frame_type, int32_t synth[],
			int32_t A_t[])
{
	/* LSPs */
	int32_t lsp_new[M];
	//int32_t lsp_mid[M];

	/* LSFs */
	int32_t prev_lsf[M];
	int32_t lsf_i[M];

	/* Algebraic codevector */
	int32_t code[L_SUBFR];

	/* excitation */
	int32_t excp[L_SUBFR];
	int32_t exc_enhanced[L_SUBFR];

	/* Scalars */
	int32_t i, i_subfr, overflow, T0_frac, index, temp, temp2, subfrNr,
	    excEnergy;
	int32_t gain_code, gain_code_mix, pit_sharp, pit_flag, pitch_fac, t0_min,
	    t0_max;
	int32_t gain_pit = 0, evenSubfr = 0, T0 = 0, index_mr475 = 0;
	int32_t *Az;		/* Pointer on A_t */
	int16_t flag4, carefulFlag;
	int16_t delta_frc_low, delta_frc_range, tmp_shift;
	int16_t bfi = 0, pdfi = 0;
	/* bad frame indication flag, potential degraded bad frame flag */

       //	enum DTXStateType newDTXState;	/* SPEECH , DTX, DTX_MUTE */

	/* find the new  DTX state  SPEECH OR DTX */
	//newDTXState = rx_dtx_handler(st->dtxDecoderState, frame_type);

	/* DTX actions */
       /*
	if (newDTXState != SPEECH) {
		Decoder_amr_reset(st, MRDTX);
		dtx_dec(st->dtxDecoderState, st->mem_syn, st->lsfState,
			st->pred_state, st->Cb_gain_averState, newDTXState,
			mode, parm, synth, A_t);

		// update average lsp
		Lsf_lsp(st->lsfState->past_lsf_q, st->lsp_old);
		lsp_avg(st->lsp_avg_st, st->lsfState->past_lsf_q);
		goto theEnd;
	}
        */
	/* SPEECH action state machine  */

       /*
	if (table_speech_bad[frame_type]) {
		bfi = 1;

		if (frame_type != RX_SPEECH_BAD) {
			Build_CN_param(&st->nodataSeed, mode, parm);
		}
	} else if (frame_type == RX_SPEECH_DEGRADED) {
		pdfi = 1;
	}

	if (bfi != 0) {
		st->state += 1;
	} else if (st->state == 6) {
		st->state = 5;
	} else {
		st->state = 0;
	}

	if (st->state > 6) {
		st->state = 6;
	}
        */
	/*
	 * If this frame is the first speech frame after CNI period,
	 * set the BFH state machine to an appropriate state depending
	 * on whether there was DTX muting before start of speech or not
	 * If there was DTX muting, the first speech frame is muted.
	 * If there was no DTX muting, the first speech frame is not
	 * muted. The BFH state machine starts from state 5, however, to
	 * keep the audible noise resulting from a SID frame which is
	 * erroneously interpreted as a good speech frame as small as
	 * possible (the decoder output in this case is quickly muted)
	 */

         /*
	if (st->dtxDecoderState->dtxGlobalState == DTX) {
		st->state = 5;
		st->prev_bf = 0;
	} else if (st->dtxDecoderState->dtxGlobalState == DTX_MUTE) {
		st->state = 5;
		st->prev_bf = 1;
	}
         */
	/* save old LSFs for CB gain smoothing */
	memcpy(prev_lsf, st->lsfState->past_lsf_q, M << 2);

	/*
	 * decode LSF parameters and generate interpolated lpc coefficients
	 * for the 4 subframes
	 */
	//if (mode != MR122) {
		D_plsf_3(st->lsfState, mode, bfi, parm, lsp_new);

		/* Advance synthesis parameters pointer */
		parm += 3;
		Int_lpc_1to3(st->lsp_old, lsp_new, A_t);
	//}

       // else {
       //		D_plsf_5(st->lsfState, bfi, parm, lsp_mid, lsp_new);

		/* Advance synthesis parameters pointer */
	//	parm += 5;
	 //	Int_lpc_1and3(st->lsp_old, lsp_mid, lsp_new, A_t);
       //	}

	/* update the LSPs for the next frame */
	memcpy(st->lsp_old, lsp_new, M << 2);

	/*
	 * Loop for every subframe in the analysis frame
	 *
	 * The subframe size is L_SUBFR and the loop is repeated
	 * L_FRAME/L_SUBFR times                                                                 *
	 *  - decode the pitch delay
	 *  - decode algebraic code
	 *  - decode pitch and codebook gains
	 *  - find the excitation and compute synthesis speech
	 */
	/* pointer to interpolated LPC parameters */
	Az = A_t;
	evenSubfr = 0;
	subfrNr = -1;

	for (i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR) {
		subfrNr += 1;
		evenSubfr = 1 - evenSubfr;

		/* flag for first and 3th subframe */
		pit_flag = i_subfr;

		if (i_subfr == L_FRAME_BY2) {
			//if ((mode != MR475) & (mode != MR515)) {
			 //	pit_flag = 0;
			//}
		}

		/* pitch index */
		index = *parm++;

		/*
		 * decode pitch lag and find adaptive codebook vector.
		 */
	      //	if (mode != MR122) {
			/*
			 * flag4 indicates encoding with 4 bit resolution;
			 * this is needed for mode MR475, MR515, MR59 and MR67
			 */
			flag4 = 0;

			if ((mode == MR475) || (mode == MR515) || (mode == MR59)
			    || (mode == MR67)) {
				flag4 = 1;
			}

			/*
			 * get ranges for the t0_min and t0_max
			 * only needed in delta decoding
			 */
			delta_frc_low = 5;
			delta_frc_range = 9;

			//if (mode == MR795) {
			//	delta_frc_low = 10;
			//	delta_frc_range = 19;
			//}
			t0_min = st->old_T0 - delta_frc_low;

			if (t0_min < PIT_MIN) {
				t0_min = PIT_MIN;
			}
			t0_max = t0_min + delta_frc_range;

			if (t0_max > PIT_MAX) {
				t0_max = PIT_MAX;
				t0_min = t0_max - delta_frc_range;
			}
			Dec_lag3(index, t0_min, t0_max, pit_flag, st->old_T0,
				 &T0, &T0_frac, flag4);
			st->T0_lagBuff = T0;

			if (bfi != 0) {
				if (st->old_T0 < PIT_MAX) {
					/* Graceful pitch degradation */
					st->old_T0 += 1;
				}
				T0 = st->old_T0;
				T0_frac = 0;

				if ((st->inBackgroundNoise !=
				     0) & (st->voicedHangover > 4) & ((mode ==
								       MR475)
								      || (mode
									  ==
									  MR515)
								      || (mode
									  ==
									  MR59)))
				{
					T0 = st->T0_lagBuff;
				}
			}
			Pred_lt_3or6_40(st->exc, T0, T0_frac, 1);
	       //	}

               /*
                else {
			Dec_lag6(index, PIT_MIN_MR122, PIT_MAX, pit_flag, &T0,
				 &T0_frac);

			if ((bfi != 0) || ((pit_flag != 0) & (index > 60))) {
				st->T0_lagBuff = T0;
				T0 = st->old_T0;
				T0_frac = 0;
			}
			Pred_lt_3or6_40(st->exc, T0, T0_frac, 0);
		}
               */


		/*
		 * (MR122 only: Decode pitch gain.)
		 * Decode innovative codebook.
		 * set pitch sharpening factor
		 */
		/* MR475, MR515 */
		//if ((mode == MR475) || (mode == MR515)) {
			/* index of position */
			index = *parm++;

			/* signs */
			i = *parm++;
			decode_2i40_9bits(subfrNr, i, index, code);
			pit_sharp = st->sharp << 1;
	       //	}
                /*
		// MR59
		else if (mode == MR59) {
			// index of position //
			index = *parm++;

			// signs
			i = *parm++;
			decode_2i40_11bits(i, index, code);
			pit_sharp = st->sharp << 1;
		}

		// MR67
		else if (mode == MR67) {
			// index of position
			index = *parm++;

			// signs
			i = *parm++;
			decode_3i40_14bits(i, index, code);
			pit_sharp = st->sharp << 1;
		}

		// MR74, MR795
		else if (mode <= MR795) {
			// index of position
			index = *parm++;

			// signs
			i = *parm++;
			decode_4i40_17bits(i, index, code);
			pit_sharp = st->sharp << 1;
		}

		// MR102
		else if (mode == MR102) {
			decode_8i40_31bits(parm, code);
			parm += 7;
			pit_sharp = st->sharp << 1;
		}

		///MR122
		else {
			index = *parm++;

			if (bfi != 0) {
				ec_gain_pitch(st->ec_gain_p_st, st->state,
					      &gain_pit);
			} else {
				gain_pit = d_gain_pitch(mode, index);
			}
			ec_gain_pitch_update(st->ec_gain_p_st, bfi, st->prev_bf,
					     &gain_pit);
			decode_10i40_35bits(parm, code);
			parm += 10;


		      //	 * pit_sharp = gain_pit;
		      //	 * if (pit_sharp > 1.0) pit_sharp = 1.0;
		       //
			pit_sharp = gain_pit;

			if (pit_sharp > 16383)
				pit_sharp = 32767;
			else
				pit_sharp *= 2;
		}

                */

		/*
		 * Add the pitch contribution to code[].
		 */
		for (i = T0; i < L_SUBFR; i++) {
			temp = (code[i - T0] * pit_sharp) >> 15;
			code[i] = code[i] + temp;
		}

		/*
		 * Decode codebook gain (MR122) or both pitch
		 * gain and codebook gain (all others)
		 * Update pitch sharpening "sharp" with quantized gain_pit
		 */
	       //	if (mode == MR475) {

			/* read and decode pitch and code gain */
			if (evenSubfr != 0) {
				/* index of gain(s) */
				index_mr475 = *parm++;
			}

		       //	if (bfi == 0) {
				Dec_gain(st->pred_state, mode, index_mr475,
					 code, evenSubfr, &gain_pit,
					 &gain_code);
		       //	}
                       // else {
		       //		ec_gain_pitch(st->ec_gain_p_st, st->state,
		       //			      &gain_pit);
		       //		ec_gain_code(st->ec_gain_c_st, st->pred_state,
		       //			     st->state, &gain_code);
		       //	}

			ec_gain_pitch_update(st->ec_gain_p_st, bfi, st->prev_bf,
					     &gain_pit);
			ec_gain_code_update(st->ec_gain_c_st, bfi, st->prev_bf,
					    &gain_code);
			pit_sharp = gain_pit;

			if (pit_sharp > SHARPMAX) {
				pit_sharp = SHARPMAX;
			}


	       //	}


                /*
                else if ((mode <= MR74) || (mode == MR102)) {
			// read and decode pitch and code gain
			// index of gain(s)
			index = *parm++;

			if (bfi == 0) {
				Dec_gain(st->pred_state, mode, index, code,
					 evenSubfr, &gain_pit, &gain_code);
			} else {
				ec_gain_pitch(st->ec_gain_p_st, st->state,
					      &gain_pit);
				ec_gain_code(st->ec_gain_c_st, st->pred_state,
					     st->state, &gain_code);
			}
			ec_gain_pitch_update(st->ec_gain_p_st, bfi, st->prev_bf,
					     &gain_pit);
			ec_gain_code_update(st->ec_gain_c_st, bfi, st->prev_bf,
					    &gain_code);
			pit_sharp = gain_pit;

			if (pit_sharp > SHARPMAX) {
				pit_sharp = SHARPMAX;
			}

			if (mode == MR102) {
				if (st->old_T0 > (L_SUBFR + 5)) {
					pit_sharp = pit_sharp >> 2;
				}
			}
		} else {
			// read and decode pitch gain
			// index of gain(s)
			index = *parm++;

			if (mode == MR795) {
				// decode pitch gain
				if (bfi != 0) {
					ec_gain_pitch(st->ec_gain_p_st,
						      st->state, &gain_pit);
				} else {
					gain_pit = d_gain_pitch(mode, index);
				}
				ec_gain_pitch_update(st->ec_gain_p_st, bfi,
						     st->prev_bf, &gain_pit);

				// read and decode code gain
				index = *parm++;

				if (bfi == 0) {
					d_gain_code(st->pred_state, mode, index,
						    code, &gain_code);
				} else {
					ec_gain_code(st->ec_gain_c_st,
						     st->pred_state, st->state,
						     &gain_code);
				}
				ec_gain_code_update(st->ec_gain_c_st, bfi,
						    st->prev_bf, &gain_code);
				pit_sharp = gain_pit;

				if (pit_sharp > SHARPMAX) {
					pit_sharp = SHARPMAX;
				}
			} else {	//MR122

				if (bfi == 0) {
					d_gain_code(st->pred_state, mode, index,
						    code, &gain_code);
				} else {
					ec_gain_code(st->ec_gain_c_st,
						     st->pred_state, st->state,
						     &gain_code);
				}
				ec_gain_code_update(st->ec_gain_c_st, bfi,
						    st->prev_bf, &gain_code);
				pit_sharp = gain_pit;
			}

                       


		}

                */

		/*
		 * store pitch sharpening for next subframe
		 * (for modes which use the previous pitch gain for
		 *  pitch sharpening in the search phase)
		 * do not update sharpening in even subframes for MR475
		 */
	       //	if ((mode != MR475) || evenSubfr == 0) {
                if (evenSubfr == 0) {
			st->sharp = gain_pit;

			if (st->sharp > SHARPMAX) {
				st->sharp = SHARPMAX;
			}
		}

		if (pit_sharp > 16383)
			pit_sharp = 32767;
		else
			pit_sharp *= 2;

		if (pit_sharp > 16384) {
			for (i = 0; i < L_SUBFR; i++) {
				temp = (st->exc[i] * pit_sharp) >> 15;
				temp2 = (temp * gain_pit) << 1;

				if (mode == MR122) {
					temp2 = (temp2 >> 1);
				}
				excp[i] = (temp2 + 0x00008000L) >> 16;
			}
		}

		/*
		 * Store list of LTP gains needed in the source
		 * characteristic detector (SCD)
		 */
		if (bfi == 0) {
			for (i = 0; i < 8; i++) {
				st->ltpGainHistory[i] =
				    st->ltpGainHistory[i + 1];
			}
			st->ltpGainHistory[8] = gain_pit;
		}

		/*
		 * Limit gain_pit if in background noise and BFI
		 * for MR475, MR515, MR59
		 */
		if ((st->prev_bf != 0
		     || bfi != 0) & (st->inBackgroundNoise != 0) & ((mode ==
								     MR475)
								    || (mode ==
									MR515)
								    || (mode ==
									MR59)))
		{
			/* if (gain_pit > 0.75) in Q14 */
			if (gain_pit > 12288)
				/* gain_pit = (gain_pit-0.75)/2.0 + 0.75; */
				gain_pit = ((gain_pit - 12288) >> 1) + 12288;

			/* if (gain_pit > 0.90) in Q14 */
			if (gain_pit > 14745) {
				gain_pit = 14745;
			}
		}

		/*
		 * Calculate CB mixed gain
		 */
		Int_lsf(prev_lsf, st->lsfState->past_lsf_q, i_subfr, lsf_i);
		gain_code_mix =
		    Cb_gain_average(st->Cb_gain_averState, mode, gain_code,
				    lsf_i, st->lsp_avg_st->lsp_meanSave, bfi,
				    st->prev_bf, pdfi, st->prev_pdf,
				    st->inBackgroundNoise, st->voicedHangover);

		/* make sure that MR74, MR795, MR122 have original codeGain */
		/* MR74, MR795, MR122 */
		//if ((mode > MR67) & (mode != MR102)) {
		//	gain_code_mix = gain_code;
	       //	}

		/*
		 * Find the total excitation.
		 * Find synthesis speech corresponding to st->exc[].
		 */
		/* MR475, MR515, MR59, MR67, MR74, MR795, MR102 */
		//if (mode <= MR102) {
			pitch_fac = gain_pit;
			tmp_shift = 1;
		//}

		/* MR122 */
		//else {
		//	pitch_fac = gain_pit >> 1;
		//	tmp_shift = 2;
	       //	}

		/*
		 * copy unscaled LTP excitation to exc_enhanced (used in phase
		 * dispersion below) and compute total excitation for LTP feedback
		 */
		memcpy(exc_enhanced, st->exc, L_SUBFR << 2);

		for (i = 0; i < L_SUBFR; i++) {
			/* st->exc[i] = gain_pit*st->exc[i] + gain_code*code[i]; */
			temp = (st->exc[i] * pitch_fac) + (code[i] * gain_code);
			temp2 = (temp << tmp_shift);
			if (((temp2 >> 1) ^ temp2) & 0x40000000) {
				if ((temp ^ temp2) & 0x80000000) {
					temp2 =
					    (temp & 0x80000000) ? (-1073741824L)
					    : 1073725439;
				} else {
					temp2 =
					    (temp2 & 0x80000000)
					    ? (-1073741824L) : 1073725439;
				}
			}
			st->exc[i] = (temp2 + 0x00004000L) >> 15;
		}
		/*
		 * Adaptive phase dispersion
		 */

		/* free phase dispersion adaption */
		st->ph_disp_st->lockFull = 0;

		if (((mode == MR475) || (mode == MR515)
		     || (mode == MR59)) & (st->voicedHangover >
					   3) & (st->inBackgroundNoise !=
						 0) & (bfi != 0)) {
			/*
			 * Always Use full Phase Disp.
			 * if error in bg noise
			 */
			st->ph_disp_st->lockFull = 1;
		}

		/*
		 * apply phase dispersion to innovation (if enabled) and
		 * compute total excitation for synthesis part
		 */
		ph_disp(st->ph_disp_st, mode, exc_enhanced, gain_code_mix,
			gain_pit, code, pitch_fac, tmp_shift);

		/*
		 * The Excitation control module are active during BFI.
		 * Conceal drops in signal energy if in bg noise.
		 */
		temp2 = 0;

		for (i = 0; i < L_SUBFR; i++) {
			temp2 += (exc_enhanced[i] * exc_enhanced[i]);
		}

		if (temp2 > 0x3FFFFFFF) {
			excEnergy = 11584;
		} else {
			temp2 = sqrt_l_exp(temp2, &temp);
			temp2 = (temp2 >> ((temp >> 1) + 15));
			excEnergy = temp2 >> 2;
		}

		if (((mode == MR475) || (mode == MR515)
		     || (mode == MR59)) & (st->voicedHangover >
					   5) & (st->inBackgroundNoise !=
						 0) & (st->state <
						       4) & (((pdfi != 0) &
							      (st->prev_pdf !=
							       0)) || bfi != 0
							     || st->prev_bf !=
							     0)) {
			carefulFlag = 0;

			if ((pdfi != 0) & (bfi == 0)) {
				carefulFlag = 1;
			}
			Ex_ctrl(exc_enhanced, excEnergy, st->excEnergyHist,
				st->voicedHangover, st->prev_bf, carefulFlag);
		}

		if ((st->inBackgroundNoise != 0) & (bfi != 0 || st->prev_bf != 0) & (st->state < 4)) {;	/* do nothing! */
		} else {
			/* Update energy history for all modes */
			for (i = 0; i < 8; i++) {
				st->excEnergyHist[i] = st->excEnergyHist[i + 1];
			}
			st->excEnergyHist[8] = excEnergy;
		}

		/*
		 * Excitation control module end.
		 */
		if (pit_sharp > 16384) {
			for (i = 0; i < L_SUBFR; i++) {
				excp[i] = excp[i] + exc_enhanced[i];
				if (labs(excp[i]) > 32767)
					excp[i] =
					    (excp[i] & 0x80000000) ? -32768 :
					    32767;
			}
			agc2(exc_enhanced, excp);
			overflow =
			    Syn_filt(Az, excp, &synth[i_subfr], L_SUBFR,
				     st->mem_syn, 0);
		} else {
			overflow =
			    Syn_filt(Az, exc_enhanced, &synth[i_subfr], L_SUBFR,
				     st->mem_syn, 0);
		}

		if (overflow) {
			for (i = 0; i < PIT_MAX + L_INTERPOL + L_SUBFR; i++) {
				st->old_exc[i] = st->old_exc[i] >> 2;
			}

			for (i = 0; i < L_SUBFR; i++) {
				exc_enhanced[i] = exc_enhanced[i] >> 2;
			}
			Syn_filt_overflow(Az, exc_enhanced, &synth[i_subfr],
					  L_SUBFR, st->mem_syn, 1);
		} else {
			memcpy(st->mem_syn, &synth[i_subfr + 30], 40);
		}

		/*
		 * Update signal for next frame.
		 * -> shift to the left by L_SUBFR  st->exc[]
		 */
		memcpy(&st->old_exc[0], &st->old_exc[L_SUBFR],
		       (PIT_MAX + L_INTERPOL) << 2);

		/* interpolated LPC parameters for next subframe */
		Az += MP1;

		/* store T0 for next subframe */
		st->old_T0 = T0;
	}

	/*
	 * Call the Source Characteristic Detector which updates
	 * st->inBackgroundNoise and st->voicedHangover.
	 */

      /*
	st->inBackgroundNoise =
	    Bgn_scd(st->background_state, &(st->ltpGainHistory[0]), &(synth[0]),
		    &(st->voicedHangover));
	dtx_dec_activity_update(st->dtxDecoderState, st->lsfState->past_lsf_q,
				synth);

	// store bfi for next subframe
	st->prev_bf = bfi;
	st->prev_pdf = pdfi;

	//
      //	 * Calculate the LSF averages on the eight
       //	 * previous frames
       */
      	lsp_avg(st->lsp_avg_st, st->lsfState->past_lsf_q);



// theEnd:
	//st->dtxDecoderState->dtxGlobalState = newDTXState;
	return;
}

/*
 * Residu40
 *
 *
 * Parameters:
 *    a                 I: prediction coefficients
 *    x                 I: speech signal
 *    y                 O: residual signal
 *
 * Function:
 *    The LP residual is computed by filtering the input
 *    speech through the LP inverse filter a(z)
 *
 * Returns:
 *    void
 */
static void Residu40(int32_t a[], int32_t x[], int32_t y[])
{
	int32_t s, i, j;

	for (i = 0; i < 40; i++) {
		s = a[0] * x[i] + a[1] * x[i - 1] + a[2] * x[i - 2] +
		    a[3] * x[i - 3];
		s += a[4] * x[i - 4] + a[5] * x[i - 5] + a[6] * x[i - 6] +
		    a[7] * x[i - 7];
		s += a[8] * x[i - 8] + a[9] * x[i - 9] + a[10] * x[i - 10];
		y[i] = (s + 0x800) >> 12;
		if (abs(y[i]) > 32767) {
			/* go to safe mode */
			for (i = 0; i < 40; i++) {
				s = a[0] * x[i];
				for (j = 1; j <= 10; j++) {
					s += a[j] * x[i - j];
					if (s > 1073741823) {
						s = 1073741823;
					} else if (s < -1073741824) {
						s = -1073741824;
					}
				}
				y[i] = (s + 0x800) >> 12;
				if (abs(y[i]) > 32767)
					y[i] =
					    (y[i] & 0x80000000) ? -32768 :
					    32767;
			}
			return;
		}

	}
	return;
}

/*
 * agc
 *
 *
 * Parameters:
 *    st->past_gain     B: gain memory
 *    sig_in            I: Post_Filter input signal
 *    sig_out           B: Post_Filter output signal
 *    agc_fac           I: AGC factor
 *
 * Function:
 *    Scales the Post_Filter output on a subframe basis
 *
 * Returns:
 *    void
 */
static void agc(agcState * st, int32_t * sig_in, int32_t * sig_out,
		int16_t agc_fac)
{
	int32_t s, gain_in, gain_out, g0, gain;
	int exp, i;

	/* calculate gain_out with exponent */
	s = energy_new(sig_out);

	if (s == 0) {
		st->past_gain = 0;
		return;
	}
	exp = 0;
	i = s;
	while (!(i & 0x40000000)) {
		exp++;
		i = i << 1;
	}
	exp -= 1;
	if (exp & 0x80000000) {
		s >>= 1;
	} else {
		s <<= exp;
	}
	gain_out = (s + 0x00008000L) >> 16;

	/* calculate gain_in with exponent */
	s = energy_new(sig_in);

	if (s == 0) {
		g0 = 0;
	} else {
		i = 0;
		while (!(s & 0x40000000)) {
			i++;
			s = s << 1;
		}
		s = s + 0x00008000L;

		if (s >= 0)
			gain_in = s >> 16;
		else
			gain_in = 32767;
		exp = (exp - i);

		/*
		 * g0 = (1-agc_fac) * sqrt(gain_in/gain_out);
		 */
		/* s = gain_out / gain_in */
		s = (gain_out << 15) / gain_in;
		exp = 7 - exp;

		if (exp > 0) {
			if (exp > 31) {
				if (s) {
					s = 2147483647;
				}
			} else {
				s = s << exp;
			}
		} else
			s = (s >> (-exp));
		if (s < 0)
			s = 2147483647;
		s = Inv_sqrt(s);
		i = ((s << 9) + 0x00008000L) >> 16;
		if (i & 0xFFFF8000)
			i = 32767;

		/* g0 = i * (1-agc_fac) */
		g0 = (i * (32767 - agc_fac)) >> 15;
	}

	/*
	 * compute gain[n] = agc_fac * gain[n-1] + (1-agc_fac) * sqrt(gain_in/gain_out)
	 * sig_out[n] = gain[n] * sig_out[n]
	 */
	gain = st->past_gain;

	for (i = 0; i < L_SUBFR; i++) {
		gain = (gain * agc_fac) >> 15;
		gain = gain + g0;
		sig_out[i] = (sig_out[i] * gain) >> 12;
		if (labs(sig_out[i]) > 32767)
			sig_out[i] = (sig_out[i] & 0x8000000) ? -32768 : 32767;
	}
	st->past_gain = gain;
	return;
}

/*
 * Post_Filter
 *
 *
 * Parameters:
 *    st                B: post filter states
 *    mode              I: AMR mode
 *    syn               B: synthesis speech
 *    Az_4              I: interpolated LPC parameters in all subfr.
 *
 * Function:
 *    Post_Filtering of synthesis speech.
 *
 *    inverse filtering of syn[] through A(z/0.7) to get res2[]
 *    tilt compensation filtering; 1 - MU*k*z^-1
 *    synthesis filtering through 1/A(z/0.75)
 *    adaptive gain control
 *
 * Returns:
 *    void
 */
static void Post_Filter(Post_FilterState * st, enum Mode mode, int32_t * syn,
			int32_t * Az_4)
{
	int32_t h[22], Ap3[MP1], Ap4[MP1];	/* bandwidth expanded LP parameters */
	int32_t tmp, i_subfr, i, temp1, temp2, overflow = 0;
	int32_t *Az, *p1, *p2, *syn_work = &st->synth_buf[M];
	const short *pgamma3 = &gamma3[0];
	const short *pgamma4 = &gamma4_gamma3_MR122[0];

	/*
	 * Post filtering
	 */
	memcpy(syn_work, syn, L_FRAME << 2);
	Az = Az_4;

	//if ((mode == MR122) || (mode == MR102)) {
	//	pgamma3 = &gamma4_gamma3_MR122[0];
	 //	pgamma4 = &gamma4_MR122[0];
       //	}

	for (i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR) {
		/* Find weighted filter coefficients Ap3[] and Ap[4] */
		Ap3[0] = Az[0];
		Ap4[0] = Az[0];

		for (i = 1; i <= 10; i++) {
			Ap3[i] = (Az[i] * pgamma3[i - 1] + 0x4000) >> 15;
			Ap4[i] = (Az[i] * pgamma4[i - 1] + 0x4000) >> 15;
		}

		/* filtering of synthesis speech by A(z/0.7) to find res2[] */
		Residu40(Ap3, &syn_work[i_subfr], st->res2);

		/* tilt compensation filter */
		/* impulse response of A(z/0.7)/A(z/0.75) */
		memcpy(h, Ap3, MP1 << 2);
		memset(&h[M + 1], 0, (22 - M - 1) << 2);
		Syn_filt(Ap4, h, h, 22, &h[M + 1], 0);

		/* 1st correlation of h[] */
		tmp = 16777216 + h[1] * h[1];

		for (i = 2; i < 22; i++) {
			tmp += h[i] * h[i];
			if (tmp > 0x3FFF8000)
				break;
		}
		temp1 = tmp >> 15;
		if (temp1 & 0xFFFF8000)
			temp1 = 32767;

		tmp = h[0] * h[1];

		for (i = 1; i < 21; i++) {
			tmp += h[i] * h[i + 1];
			if (abs(tmp) > 1073741823)
				tmp = 1073741823;
		}
		temp2 = tmp >> 15;

		if (temp2 <= 0) {
			temp2 = 0;
		} else {
			tmp = temp2 * 26214;
			temp2 = (tmp & 0xffff8000) / temp1;
		}

		/* preemphasis */
		p1 = st->res2 + 39;
		p2 = p1 - 1;
		tmp = *p1;

		do {
			*p1 = *p1 - ((temp2 * *p2--) >> 15);
			if (abs(*p1) > 32767) {
				*p1 = (*p1 & 0x80000000) ? -32768 : 32767;
			}
			p1--;
			*p1 = *p1 - ((temp2 * *p2--) >> 15);
			if (abs(*p1) > 32767) {
				*p1 = (*p1 & 0x80000000) ? -32768 : 32767;
			}
			p1--;
			*p1 = *p1 - ((temp2 * *p2--) >> 15);
			if (abs(*p1) > 32767) {
				*p1 = (*p1 & 0x80000000) ? -32768 : 32767;
			}
			p1--;
		} while (p1 > st->res2);
		*p1 = *p1 - ((temp2 * st->preemph_state_mem_pre) >> 15);
		if (abs(*p1) > 32767) {
			*p1 = (*p1 & 0x80000000) ? -32768 : 32767;
		}
		st->preemph_state_mem_pre = tmp;

		/* filtering through  1/A(z/0.75) */
		overflow =
		    Syn_filt(Ap4, st->res2, &syn[i_subfr], L_SUBFR,
			     st->mem_syn_pst, 0);
		if (overflow) {
			Syn_filt_overflow(Ap4, st->res2, &syn[i_subfr], L_SUBFR,
					  st->mem_syn_pst, 1);
		} else {
			memcpy(st->mem_syn_pst, &syn[i_subfr + 30], 40);
		}

		/* scale output to input */
		agc(st->agc_state, &syn_work[i_subfr], &syn[i_subfr], AGC_FAC);
		Az += MP1;
	}

	/* update syn_work[] buffer */
	memcpy(&syn_work[-M], &syn_work[L_FRAME - M], M << 2);
	return;
}

/*
 * Post_Process
 *
 *
 * Parameters:
 *    st                B: post filter states
 *    signal            B: signal
 *
 * Function:
 *    Postprocessing of input speech.
 *
 *    2nd order high pass filtering with cut off frequency at 60 Hz.
 *    Multiplication of output by two.
 *
 *
 * Returns:
 *    void
 */
static void Post_Process(Post_ProcessState * st, int32_t signal[])
{
	int32_t x2, tmp, i = 0;
	int32_t mask = 0x40000000;

	do {
		x2 = st->x1;
		st->x1 = st->x0;
		st->x0 = signal[i];

		/*
		 * y[i] = b[0]*x[i]*2 + b[1]*x[i-1]*2 + b140[2]*x[i-2]/2
		 *                    + a[1]*y[i-1] + a[2] * y[i-2];
		 */
		tmp =
		    (st->y1_hi * 15836) +
		    (((st->y1_lo * 15836) & (int32_t) 0xffff8000) >> 15);
		tmp +=
		    (st->y2_hi * -7667) +
		    (((st->y2_lo * (-7667)) & (int32_t) 0xffff8000) >> 15);
		tmp += st->x0 * 7699;
		tmp += st->x1 * -15398;
		if (((tmp >> 1) ^ tmp) & mask)
			tmp = (tmp & 0x80000000) ? -1073741824 : 1073741823;

		tmp += x2 * 7699;
		if (((tmp >> 1) ^ tmp) & mask)
			tmp = (tmp & 0x80000000) ? -1073741824 : 1073741823;

		tmp = tmp << 1;
		if (((tmp >> 1) ^ tmp) & mask)
			tmp = (tmp & 0x80000000) ? -1073741824 : 1073741823;

		tmp = tmp << 1;
		if (((tmp >> 1) ^ tmp) & mask)
			tmp = (tmp & 0x80000000) ? -1073741824 : 1073741823;

		if (labs(tmp) < 536862720) {
			signal[i++] = (tmp + 0x00002000L) >> 14;
		} else if (tmp > 0) {
			signal[i++] = 32767;
		} else {
			signal[i++] = -32768;
		}
		st->y2_hi = st->y1_hi;
		st->y2_lo = st->y1_lo;
		st->y1_hi = tmp >> 15;
		st->y1_lo = ((tmp << 1) - (st->y1_hi << 16)) >> 1;
	} while (i < 160);
	return;
}

/*
 * Speech_Decode_Frame
 *
 *
 * Parameters:
 *    st                B: decoder memory
 *    mode              I: AMR mode
 *    parm              I: speech parameters
 *    frame_type        I: Frame type
 *    synth             O: synthesis speech

 * Function:
 *    Decode one frame
 *
 * Returns:
 *    void
 */
void Speech_Decode_Frame(void *st, enum Mode mode, int16_t * parm, enum
			 RXFrameType frame_type, int16_t * synth)
{
	int32_t Az_dec[AZ_SIZE];	/* Decoded Az for post-filter in 4 subframes */
	int32_t synth_speech[L_FRAME];
	int32_t i;

	/* Synthesis */
	Decoder_amr(((Speech_Decode_FrameState *) st)->decoder_amrState, mode,
		    parm, frame_type, synth_speech, Az_dec);
	Post_Filter(((Speech_Decode_FrameState *) st)->post_state, mode,
		    synth_speech, Az_dec);

	/* post HP filter, and 15->16 bits */
	Post_Process(((Speech_Decode_FrameState *) st)->postHP_state,
		     synth_speech);

	for (i = 0; i < L_FRAME; i++) {
#ifndef NO13BIT
		/* Truncate to 13 bits */
		synth[i] = (int16_t) (synth_speech[i] & 0xfff8);
#else
		synth[i] = (int16_t) (synth_speech[i]);
#endif
	}

	return;
}

/*
 * Decoder_amr_exit
 *
 *
 * Parameters:
 *    state                I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
static void Decoder_amr_exit(Decoder_amrState ** state)
{
	if (state == NULL || *state == NULL)
		return;
	if ((*state)->lsfState)
		free((*state)->lsfState);
	if ((*state)->ec_gain_p_st)
		free((*state)->ec_gain_p_st);
	if ((*state)->ec_gain_c_st)
		free((*state)->ec_gain_c_st);
	if ((*state)->pred_state)
		free((*state)->pred_state);
	//if ((*state)->background_state)
	 //	free((*state)->background_state);
	if ((*state)->ph_disp_st)
		free((*state)->ph_disp_st);
	if ((*state)->Cb_gain_averState)
		free((*state)->Cb_gain_averState);
	if ((*state)->lsp_avg_st)
		free((*state)->lsp_avg_st);
	//if ((*state)->dtxDecoderState)
	//	free((*state)->dtxDecoderState);

	/* deallocate memory */
	free(*state);
	*state = NULL;
	return;
}

/*
 * Post_Filter_exit
 *
 *
 * Parameters:
 *    state                I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
static void Post_Filter_exit(Post_FilterState ** state)
{
	if (state == NULL || *state == NULL)
		return;
	free((*state)->agc_state);

	/* deallocate memory */
	free(*state);
	*state = NULL;
	return;
}

/*
 * Post_Process_reset
 *
 *
 * Parameters:
 *    state             B: state structure
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *    -1 failure
 */
static int Post_Process_reset(Post_ProcessState * state)
{
	if ((Post_ProcessState *) state == NULL) {
		//fprintf(stderr, "Post_Process_reset: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	state->y2_hi = 0;
	state->y2_lo = 0;
	state->y1_hi = 0;
	state->y1_lo = 0;
	state->x0 = 0;
	state->x1 = 0;
	return 0;
}

/*
 * Post_Process_exit
 *
 *
 * Parameters:
 *    state                I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
static void Post_Process_exit(Post_ProcessState ** state)
{
	if (state == NULL || *state == NULL)
		return;

	/* deallocate memory */
	free(*state);
	*state = NULL;
	return;
}

/*
 * Decoder_amr_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    success = 0
 */
static int Decoder_amr_init(Decoder_amrState ** state)
{
	Decoder_amrState *s;

	if ((Decoder_amrState *) state == NULL) {
		mtrap(__LINE__);
	}
	*state = NULL;

	/* allocate memory */
	if ((s = calloc(1, sizeof(Decoder_amrState))) == NULL) {
		mtrap(__LINE__);
	}
	//dec_mem+=sizeof(Decoder_amrState);
        memadd(sizeof(Decoder_amrState), 2);
	/* DPlsf_init */
	/* allocate memory */
	if ((s->lsfState = calloc(1, sizeof(D_plsfState))) == NULL) {
		mtrap(__LINE__);
	}
	//dec_mem+=sizeof(D_plsfState);
        memadd(sizeof(D_plsfState), 3);
	/* ecGainPitchInit */
	/* allocate memory */
	if ((s->ec_gain_p_st =
	     calloc(1, sizeof(ec_gain_pitchState))) ==
	    NULL) {
		mtrap(__LINE__);
		goto lfree;
	}
  //dec_mem+=sizeof(ec_gain_pitchState);
  memadd(sizeof(ec_gain_pitchState), 4);

	/* ecGainCodeInit */
	/* allocate memory */
	if ((s->ec_gain_c_st =
	     calloc(1, sizeof(ec_gain_codeState))) == NULL) {
		mtrap(__LINE__);
		goto lfree;
	}
 // dec_mem+=sizeof(ec_gain_codeState);
  memadd(sizeof(ec_gain_codeState), 5);


	/* gcPredInit */
	/* allocate memory */
	if ((s->pred_state = calloc(1, sizeof(gc_predState)))
	    == NULL) {
		mtrap(__LINE__);
	}
  //dec_mem+=sizeof(gc_predState);
  memadd(sizeof(gc_predState), 6);

	/* Cb_gain_averageInit */
	/* allocate memory */
	if ((s->Cb_gain_averState =
	     calloc(1, sizeof(Cb_gain_averageState))) ==
	    NULL) {
		mtrap(__LINE__);
	}
	//dec_mem+=sizeof(Cb_gain_averageState);
         memadd(sizeof(Cb_gain_averageState), 7);

	memset(s->Cb_gain_averState->cbGainHistory, 0, L_CBGAINHIST << 2);

	/* Initialize hangover handling */
	s->Cb_gain_averState->hangVar = 0;
	s->Cb_gain_averState->hangCount = 0;

	/* lsp_avgInit */
	/* allocate memory */


	if ((s->lsp_avg_st = calloc(1, sizeof(lsp_avgState)))
	    == NULL) {
		mtrap(__LINE__);
	}
  //dec_mem+=sizeof(lsp_avgState);
  memadd(sizeof(lsp_avgState), 8);

	/* Bgn_scdInit */
	/* allocate memory */
       //	if ((s->background_state =
       //	     calloc(1, sizeof(Bgn_scdState))) == NULL) {
	//	mtrap(__LINE__);
	//}
  //dec_mem+=sizeof(Bgn_scdState);
	/* phDispInit */
	/* allocate memory */
	if ((s->ph_disp_st = calloc(1, sizeof(ph_dispState)))
	    == NULL) {
		mtrap(__LINE__);
	}
  //dec_mem+=sizeof(ph_dispState);
  memadd(sizeof(ph_dispState), 9);

	/* dtxDecInit */
	/* allocate memory */
	//if ((s->dtxDecoderState = calloc(1, sizeof(dtx_decState))) == NULL) {
	//	mtrap(__LINE__);
	//}
	//dec_mem+=sizeof(dtx_decState);
	Decoder_amr_reset(s, (enum Mode)0);
	*state = s;
	return 0;
lfree:
	if (s->ph_disp_st)
		free(s->ph_disp_st);
	//if (s->background_state)
	//	free(s->background_state);
	if (s->lsp_avg_st)
		free(s->lsp_avg_st);
	if (s->Cb_gain_averState)
		free(s->Cb_gain_averState);
	if (s->pred_state)
		free(s->pred_state);
	if (s->ec_gain_c_st)
		free(s->ec_gain_c_st);
	if (s->ec_gain_p_st)
		free(s->ec_gain_p_st);
	if (s->lsfState)
		free(s->lsfState);
	if (s)
		free(s);
	return -1;
}

/*
 * Post_Filter_reset
 *
 *
 * Parameters:
 *    state             B: state structure
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *    -1 failure
 */
static int Post_Filter_reset(Post_FilterState * state)
{
	if ((Post_FilterState *) state == NULL) {
		//fprintf(stderr, "Post_Filter_reset: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	state->preemph_state_mem_pre = 0;
	state->agc_state->past_gain = 4096;
	memset(state->mem_syn_pst, 0, M << 2);
	memset(state->res2, 0, L_SUBFR << 2);
	memset(state->synth_buf, 0, (L_FRAME + M) << 2);
	return 0;
}

/*
 * Post_Filter_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    success = 0
 */
static int Post_Filter_init(Post_FilterState ** state)
{
	Post_FilterState *s;

	if ((Post_FilterState *) state == NULL) {
		//fprintf(stderr, "F057:invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	*state = NULL;

	/* allocate memory */
	if ((s = calloc(1, sizeof(Post_FilterState))) ==
	    NULL) {
		//fprintf(stderr, "F057:can not malloc filter structure\n");
		//return -1;
				mtrap(__LINE__);
	}
			dec_mem+=sizeof(Post_FilterState);
	s->agc_state = NULL;

	/* allocate memory */
	if ((s->agc_state = calloc(1, sizeof(agcState))) == NULL) {
		//fprintf(stderr, "agcInit: can not malloc state structure\n");
		mtrap(__LINE__);
		//free(s);
		//return -1;
	}
	dec_mem+=sizeof(agcState);
	Post_Filter_reset(s);
	*state = s;
	return 0;
}

/*
 * Post_Process_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    success = 0
 */
static int Post_Process_init(Post_ProcessState ** state)
{
	Post_ProcessState *s;

	if ((Post_ProcessState *) state == NULL) {
		//fprintf(stderr, "Post_Process_init: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	*state = NULL;

	/* allocate memory */
	if ((s = calloc(1, sizeof(Post_ProcessState))) ==
	    NULL) {
		//fprintf(stderr,
		//	"Post_Process_init: can not malloc state structure\n");
		//return -1;
				mtrap(__LINE__);
	}
			dec_mem+=sizeof(Post_ProcessState);
	Post_Process_reset(s);
	*state = s;
	return 0;
}

/*
 * Speech_Decode_Frame_exit
 *
 *
 * Parameters:
 *    state                I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
void Speech_Decode_Frame_exit(void *_st)
{
	Speech_Decode_FrameState *st = _st;
	
	if (_st == NULL)
		return;

	

	Decoder_amr_exit(&st->decoder_amrState);
	Post_Filter_exit(&st->post_state);
	Post_Process_exit(&st->postHP_state);

	/* deallocate memory */
	free(st);
	return;
}

/*
 * Speech_Decode_Frame_reset
 *
 *
 * Parameters:
 *    state             B: state structure
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *    -1 = failure
 */
int Speech_Decode_Frame_reset(void **st)
{
	Speech_Decode_FrameState *state;

	if (st == NULL || *st == NULL)
		return (-1);
	state = (Speech_Decode_FrameState *) st;
	Decoder_amr_reset(state->decoder_amrState, (enum Mode)0);
	Post_Filter_reset(state->post_state);
	Post_Process_reset(state->postHP_state);
	return 0;
}

/*
 * Speech_Decode_Frame_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    success = 0
 */
void *Speech_Decode_Frame_init()
{
	Speech_Decode_FrameState *s;

	/* allocate memory */
	if ((s = calloc(1, sizeof(Speech_Decode_FrameState))) == NULL) {
		//fprintf(stderr,
		//	"Speech_Decode_Frame_init: can not malloc state "
		//	"structure\n");
		//return NULL;
		mtrap(__LINE__);
	}
	dec_mem+=sizeof(Speech_Decode_FrameState);
	s->decoder_amrState = NULL;
	s->post_state = NULL;
	s->postHP_state = NULL;

	if (Decoder_amr_init(&s->decoder_amrState)
	    || Post_Filter_init(&s->post_state)
	    || Post_Process_init(&s->postHP_state)) {
		Speech_Decode_Frame_exit(s);
		return NULL;
	}
	return s;
}
