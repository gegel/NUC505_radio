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
 * rom_dec.h
 *
 *
 * Project:
 *    AMR Floating-Point Codec
 *
 * Contains:
 *    This file contains all the tables needed by AMR decoder functions.
 *
 */

#pragma once

#ifndef _ROM_DEC_H_
#define _ROM_DEC_H_

/*
 * include files
 */
#include <stdint.h>
#include "interf_rom.h"

/*
 * definition of constants
 */
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
#define NPRED 4			/* number of prediction taps */
#define SHARPMIN  0		/* Minimum value of pitch sharpening */
#define MAX_PRM_SIZE    57	/* max. num. of params */
#define L_INTER_SRCH    4	/* Length of filter for CL LTP search interpolation */
#define GP_CLIP         0.95F	/* Pitch gain clipping */
#define UP_SAMP_MAX     6
#define NB_TRACK        5	/* number of tracks */
#define NB_TRACK_MR102  4	/* number of tracks mode mr102 */
#define STEP            5	/* codebook step size */
#define STEP_MR102      4	/* codebook step size mode mr102 */
#define NC              M/2	/* Order of LP filter divided by 2 */

/* vad */
#define COMPLEN               9	/* Number of sub-bands used by VAD */
#define L_ENERGYHIST          60
#define L_CBGAINHIST          7
#define PHDGAINMEMSIZE        5
#define MIN_ENERGY            -14336	/* 14 Q10 */
#define MIN_ENERGY_MR122      -2381	/* 14 / (20*log10(2)) Q10 */
#define PN_INITIAL_SEED       0x70816958L	/* Pseudo noise generator seed value  */
#define MIN_16                (int16_t)-32768
#define MAX_16                (int16_t)0x7fff
#define MAX_32                (short)0x7fffffffL
#define EXPCONST              5243	/* 0.16 in Q15 */
#define DTX_MAX_EMPTY_THRESH  50
#define DTX_ELAPSED_FRAMES_THRESH (24 + 7 -1)
#define LSF_GAP               205	/* Minimum distance between LSF after quantization; 50 Hz = 205 */
#define LSP_PRED_FAC_MR122    21299	/* MR122 LSP prediction factor (0.65 Q15) */
#define POS_CODE              8191
#define NEG_CODE              8191
#define NMAX                  9	/* largest N used in median calculation */
#define MEAN_ENER_MR122       783741L	/* 36/(20*log10(2)) (Q17) */
#define SHARPMAX              13017	/* Maximum value of pitch sharpening */
#define FRAMEENERGYLIMIT      17578	/* 150 */
#define LOWERNOISELIMIT       20	/*  5 */
#define UPPERNOISELIMIT       1953	/* 50 */
#define AZ_SIZE               (4*M+4)	/* Size of array of LP filters in 4 subfr.s */
#define AGC_FAC               29491	/* Factor for automatic gain control 0.9 */
#define PHDGAINMEMSIZE        5
#define PHDTHR1LTP            9830	/* 0.6 in Q14 */
#define PHDTHR2LTP            14746	/* 0.9 in Q14 */
#define ONFACTPLUS1           16384	/* 2.0 in Q13 */
#define ONLENGTH              2
#define DTX_HANG_CONST        7	/* yields eight frames of SP HANGOVER */

/* number of parameters */
#define PRMNO_MR475 17
#define PRMNO_MR515 19
#define PRMNO_MR59  19
#define PRMNO_MR67  19
#define PRMNO_MR74  19
#define PRMNO_MR795 23
#define PRMNO_MR102 39
#define PRMNO_MR122 57
#define PRMNO_MRDTX 5


//===========================

#define DICO1_SIZE_5  128
#define DICO2_SIZE_5  256
#define DICO3_SIZE_5  256
#define DICO4_SIZE_5  256
#define DICO5_SIZE_5  64

#define MR475_VQ_SIZE 256

#define VQ_SIZE_HIGHRATES 128

#define VQ_SIZE_LOWRATES 64
#define NB_QUA_PITCH 16
#define NB_QUA_CODE 32
#define PAST_RQ_INIT_SIZE 8

#define ALPHA     29491
#define ONE_ALPHA 3277

#define ALPHA_122     31128
#define ONE_ALPHA_122 1639

#define DICO1_SIZE_3  256
#define DICO2_SIZE_3  512
#define DICO3_SIZE_3  512

#define MR515_3_SIZE  128

#define MR795_1_SIZE  512


typedef struct
{
 short tab1[8];
 short tab2[16];
}dec_r;

 void init_dec_rom(void);

extern const short cdown[7];
extern const short pdown[7];
extern const short pred[4];
extern const short gamma4_gamma3_MR122[10];
extern const short gamma3[10];
extern const short sqrt_table[49];
extern const short inv_sqrt_table[49];
extern const short log2_table[33];
extern const short pow2_table[33];
extern const short cos_table[65];
extern const short ph_imp_low[40];
extern const short ph_imp_mid[40];
extern const short mean_lsf_3[10];
extern const short pred_fac[10];
extern const short dico1_lsf_3[768];
extern const short dico2_lsf_3[1536];
extern const short mr515_3_lsf[512];
extern const short table_gain_MR475[1024];
extern const short inter6[61];
extern const short window_200_40[240];
extern const int8_t startPos[16];

//total 9048

/*
extern const short* cdown;
extern const short* pdown;
extern const short* pred;
extern const short* gamma4_gamma3_MR122;
extern const short* gamma3;
extern const short* sqrt_table;
extern const short* inv_sqrt_table;
extern const short* log2_table;
extern const short* pow2_table;
extern const short* cos_table;
extern const short* ph_imp_low;
extern const short* ph_imp_mid;
extern const short* mean_lsf_3;
extern const short* pred_fac;
extern const short* dico1_lsf_3;
extern const short* dico2_lsf_3;
extern const short* mr515_3_lsf;
extern const short* table_gain_MR475;
extern const short* inter6;
extern const short* window_200_40;
extern const int8_t* startPos;
 */

#endif /* _ROM_DEC_H_ */

