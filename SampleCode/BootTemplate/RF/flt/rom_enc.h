/* vim: set tabstop=4:softtabstop=4:shiftwidth=4:noexpandtab */

/*
 * ===================================================================
 *  TS 26.104
 *  REL-5 V5.4.0 2004-03
 *  REL-6 V6.1.0 2004-03
 *  3GPP AMR Floating-point Speech Codec
 * ===================================================================
 *
 */

/*
 * rom_enc.h
 *
 *
 * Project:
 *    AMR Floating-Point Codec
 *
 * Contains:
 *    This file contains all the tables needed by AMR encoder functions.
 *
 */

#pragma once

#ifndef _ROM_ENC_H_
#define _ROM_ENC_H_

#include <stdint.h>
#include "sp_enc.h"

#define M               10	/* Order of LP filter */
#define MP1             (M+1)	/* Order of LP filter + 1 */
#define L_WINDOW        240	/* Window size in LP analysis */
#define L_NEXT          40	/* Overhead in LP analysis */
#define LTPG_MEM_SIZE   5	/* number of stored past LTP coding gains + 1 */
#define N_FRAME         7	/* old pitch gains in average calculation */
#define DTX_HIST_SIZE   8	/* DTX history size */
#define L_TOTAL         320	/* Total size of speech buffer. */
#define L_FRAME         160	/* Frame size */
#define L_FRAME_BY2     80	/* Frame size divided by 2 */
#define L_SUBFR         40	/* Subframe size */
#define L_CODE          40	/* codevector length */
#define PIT_MAX         143	/* Maximum pitch lag */
#define PIT_MIN         20	/* Minimum pitch lag */
#define PIT_MIN_MR122   18	/* Minimum pitch lag (MR122 mode) */
#define L_INTERPOL      (10+1)	/* Length of filter for interpolation */
#define NPRED           4	/* number of prediction taps */
#define SHARPMIN        0	/* Minimum value of pitch sharpening */
#define MAX_PRM_SIZE    57	/* max. num. of params */
#define L_INTER_SRCH    4	/* Length of filter for CL LTP search interpolation */
#define GP_CLIP         0.95F	/* Pitch gain clipping */
#define UP_SAMP_MAX     6
#define NB_TRACK        5	/* number of tracks */
#define NB_TRACK_MR102  4	/* number of tracks mode mr102 */
#define STEP            5	/* codebook step size */
#define STEP_MR102      4	/* codebook step size mode mr102 */
#define NC              M/2	/* Order of LP filter divided by 2 */

#define SCALE_LSP_FREQ  (float)(4000.0/3.141592654)
#define SCALE_FREQ_LSP  (float)(3.141592654/4000.0)
#define SLOPE1_WGHT_LSF (float)((3.347-1.8)/(450.0-0.0))
#define SLOPE2_WGHT_LSF (float)((1.8-1.0)/(1500.0-450.0))

#define FRAME_LEN 160		/* Length (samples) of the input frame          */
#define COMPLEN 9		/* Number of sub-bands used by VAD              */
#define INV_COMPLEN 3641	/* 1.0/COMPLEN*2^15                             */
#define LOOKAHEAD 40		/* length of the lookahead used by speech coder */

#define UNITY 512		/* Scaling used with SNR calculation            */
#define UNIRSHFT 6		/* = log2(MAX_16/UNITY)                         */

#define TONE_THR 0.65F		/* Threshold for tone detection   */

/* Constants for background spectrum update */
#define ALPHA_UP1   (float)(1.0 - 0.95)	/* Normal update, upwards:   */
#define ALPHA_DOWN1 (float)(1.0 - 0.936)	/* Normal update, downwards  */
#define ALPHA_UP2   (float)(1.0 - 0.985)	/* Forced update, upwards    */
#define ALPHA_DOWN2 (float)(1.0 - 0.943)	/* Forced update, downwards  */
#define ALPHA3      (float)(1.0 - 0.95)	/* Update downwards          */
#define ALPHA4      (float)(1.0 - 0.9)	/* For stationary estimation */
#define ALPHA5      (float)(1.0 - 0.5)	/* For stationary estimation */

/* Constants for VAD threshold */
#define VAD_THR_HIGH 1260	/* Highest threshold                 */
#define VAD_THR_LOW  720	/* Lowest threshold                  */
#define VAD_P1 0		/* Noise level for highest threshold */
#define VAD_P2 6300		/* Noise level for lowest threshold  */
#define VAD_SLOPE (float)(VAD_THR_LOW-VAD_THR_HIGH)/(float)(VAD_P2-VAD_P1)

/* Parameters for background spectrum recovery function */
#define STAT_COUNT 20		/* threshold of stationary detection counter         */
#define STAT_COUNT_BY_2 10	/* threshold of stationary detection counter         */
#define CAD_MIN_STAT_COUNT 5	/* threshold of stationary detection counter         */

#define STAT_THR_LEVEL 184	/* Threshold level for stationarity detection        */
#define STAT_THR 1000		/* Threshold for stationarity detection              */

/* Limits for background noise estimate */
#define NOISE_MIN 40		/* minimum */
#define NOISE_MAX 16000		/* maximum */
#define NOISE_INIT 150		/* initial */

/* Constants for VAD hangover addition */
#define HANG_NOISE_THR 100
#define BURST_LEN_HIGH_NOISE 4
#define HANG_LEN_HIGH_NOISE 7
#define BURST_LEN_LOW_NOISE 5
#define HANG_LEN_LOW_NOISE 4

/* Thresholds for signal power */
#define VAD_POW_LOW (short)15000/2	/* If input power is lower,                    */
				       /*     VAD is set to 0                         */
#define POW_PITCH_THR (short)343040/2	/* If input power is lower, pitch              */
				       /*     detection is ignored                    */

#define POW_COMPLEX_THR (short)15000/2	/* If input power is lower, complex            */
				       /* flags  value for previous frame  is un-set  */
/*
 * VAD Constants
 */

/* Constants for the filter bank */
#define LEVEL_SHIFT 0		/* scaling                                  */
#define COEFF3   (float)13363/32768	/* coefficient for the 3rd order filter     */
#define COEFF5_1 (float)21955/32768	/* 1st coefficient the for 5th order filter */
#define COEFF5_2 (float)6390/32768	/* 2nd coefficient the for 5th order filter */

/* Constants for pitch detection */
#define LTHRESH 4
#define NTHRESH 4

/* Constants for complex signal VAD  */
#define CVAD_THRESH_ADAPT_HIGH  0.6F	/* threshold for adapt stopping high */
#define CVAD_THRESH_ADAPT_LOW  0.5F	/* threshold for adapt stopping low */
#define CVAD_THRESH_IN_NOISE  0.65F * 32768.0F	/* threshold going into speech on
						   a short term basis */
#define CVAD_THRESH_HANG  0.70F	/* threshold */
#define CVAD_HANG_LIMIT  (int16_t)(100)	/* 2 second estimation time */
#define CVAD_HANG_LENGTH  (int16_t)(250)	/* 5 second hangover */

#define CVAD_LOWPOW_RESET 0.40F	/* init in low power segment */
#define CVAD_MIN_CORR 0.40F	/* lowest adaptation value */

#define CVAD_BURST 20		/* speech burst length for speech reset */
#define CVAD_ADAPT_SLOW 1.0F - 0.98F	/* threshold for slow adaption */
#define CVAD_ADAPT_FAST 1.0F - 0.92F	/* threshold for fast adaption */
#define CVAD_ADAPT_REALLY_FAST 1.0F - 0.80F	/* threshold for really fast adaption */


//==========================================

#define NB_QUA_PITCH 16
#define NB_QUA_CODE 32
#define PAST_RQ_INIT_SIZE 8
#define DICO1_SIZE_3  256
#define DICO2_SIZE_3  512
#define DICO3_SIZE_3  512
#define MR515_3_SIZE  128
#define MR795_1_SIZE  512
#define DICO1_SIZE_5 128
#define DICO2_SIZE_5 256
#define DICO3_SIZE_5 256
#define DICO4_SIZE_5 256
#define DICO5_SIZE_5 64
#define MR475_VQ_SIZE 256
#define VQ_SIZE_HIGHRATES 128
#define VQ_SIZE_LOWRATES 64
#define DTX_VQ_SIZE 47

/*
extern int8_t f_trackTable[4 * 5];
extern const float f_gamma1[M];
extern const float f_gamma1_12k2[M];
extern const float f_gamma2[M];
extern const float f_b60[UP_SAMP_MAX * (L_INTERPOL - 1) + 1];
extern const short f_inter6[61];
extern const int16_t f_startPos1[2];
extern const int16_t f_startPos2[4];
extern const int16_t f_startPos[2 * 4 * 2];
extern const float f_qua_gain_pitch[NB_QUA_PITCH];
extern const float f_qua_gain_pitch_MR122[NB_QUA_PITCH];
extern const float f_gain_factor[NB_QUA_CODE];
extern const int8_t f_gray[8];
extern const float f_grid[61];
extern const float f_b24[UP_SAMP_MAX * L_INTER_SRCH + 1];
extern const float f_lag_wind[M];
extern const float f_lsp_init_data[M];
extern const float f_past_rq_init[80];
extern const float f_mean_lsf_3[10];
extern const float f_mean_lsf_5[10];
extern const float f_pred_fac[10];
extern const float f_dico1_lsf_3[768];
extern const float f_dico2_lsf_3[1536];
extern const float f_dico3_lsf_3[2048];
extern const float f_mr515_3_lsf[512];
extern const float f_mr795_1_lsf[1536];
extern const float f_dico1_lsf_5[DICO1_SIZE_5 * 4];
extern const float f_dico2_lsf_5[DICO2_SIZE_5 * 4];
extern const float f_dico3_lsf_5[DICO3_SIZE_5 * 4];
extern const float f_dico4_lsf_5[DICO4_SIZE_5 * 4];
extern const float f_dico4_lsf_5[DICO4_SIZE_5 * 4];
extern const float f_dico5_lsf_5[DICO5_SIZE_5 * 4];
extern const float f_table_gain_MR475[MR475_VQ_SIZE * 4];
extern const float f_table_highrates[VQ_SIZE_HIGHRATES * 2];
extern const float f_table_lowrates[VQ_SIZE_LOWRATES * 2];
extern const short f_qua_gain_code_MR122[NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES+(MR475_VQ_SIZE * 2)+DTX_VQ_SIZE + 1];
extern const short f_qua_gain_code[NB_QUA_CODE+VQ_SIZE_HIGHRATES+VQ_SIZE_LOWRATES + (MR475_VQ_SIZE * 2)+DTX_VQ_SIZE + 3];
extern const float f_window_200_40[240];
extern const float f_window_232_8[240];
extern const float f_window_160_80[240];
extern const float f_corrweight[251];
//extern const struct mode_dep_parm[8];
extern const short f_log2_table[33];
extern const short f_pow2_table[33];
*/


extern int8_t f_trackTable[20]; 
extern const float f_gamma1[10];
extern const float f_gamma2[10];
extern const float f_b60[61];
extern const short f_inter6[61];
extern const int16_t f_startPos[16];
extern const float f_gain_factor[32];
extern const float f_grid[61];
extern const float f_b24[25];
extern const float f_lag_wind[10];
extern const float f_lsp_init_data[10];
extern const float f_mean_lsf_3[10];
extern const float f_pred_fac[10];
extern const float f_dico1_lsf_3[768];
extern const float f_dico2_lsf_3[1536];
extern const float f_mr515_3_lsf[512];
extern const float f_table_gain_MR475[1024];
extern const short f_qua_gain_code[786];
extern const float f_window_200_40[240];
extern const short f_log2_table[33];
extern const short f_pow2_table[33];

//total 19098 bytes

#endif /* _ROM_ENC_H_ */

