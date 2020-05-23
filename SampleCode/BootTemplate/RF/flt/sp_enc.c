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
 * sp_enc.c
 *
 *
 * Project:
 *    AMR Floating-Point Codec
 *
 * Contains:
 *    This module contains all the functions needed encoding 160
 *    16-bit speech samples to AMR encoder parameters.
 *
 */
 


/////////////#include "arm_math.h"
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "platform.h"
//#include <strings.h>
//#include <memory.h>

#include <math.h>
#include <float.h>

#include "sp_enc.h"
#include "rom_enc.h"
//#include <ophtools.h>



 struct {
	int8_t max_frac_lag;	/* lag up to which fractional lags are used */
	int8_t flag3;		/* enable 1/3 instead of 1/6 fract. resolution */
	int8_t first_frac;	/* first fractional to check */
	int8_t last_frac;	/* last fractional to check */
	int8_t delta_int_low;	/* integer lag below TO to start search from */
	int8_t delta_int_range;	/* integer range around T0 */
	int8_t delta_frc_low;	/* fractional below T0 */
	int8_t delta_frc_range;	/* fractional range around T0 */
	int8_t pit_min;		/* minimum pitch */

} mode_dep_parm[1] =


        {
	       {84, 1, -2, 2, 5, 10, 5, 9, PIT_MIN} //MR475
        };



volatile int16_t ssss=0;
int enc_mem=0;

void mtrap(int line);
void memadd(int len, unsigned char id);

;

/*
 * Definition of structures used in encoding process
 */
    typedef struct {
	float y2;
	float y1;
	float x0;
	float x1;

} Pre_ProcessState;

#ifdef VAD2
 
/* Defines for VAD2 */
#define	FRM_LEN1		80
#define	DELAY0			24
#define	FFT_LEN1		128

#define	UPDATE_CNT_THLD1	50

#define	INIT_FRAMES		4

#define	CNE_SM_FAC1		0.1F
#define	CEE_SM_FAC1		0.55F

#define	HYSTER_CNT_THLD1	6	/* forced update constants... */
#define	HIGH_ALPHA1		0.9F
#define	LOW_ALPHA1		0.7F
#define	ALPHA_RANGE1		(HIGH_ALPHA1-LOW_ALPHA1)

#define NORM_ENRG		(4.0F)	/* account for div by 2 by the HPF */
#define	MIN_CHAN_ENRG		(0.0625F / NORM_ENRG)
#define	INE			(16.0F / NORM_ENRG)
#define	NOISE_FLOOR		(1.0F / NORM_ENRG)

#define	PRE_EMP_FAC1		(-0.8F)

#define	NUM_CHAN		16
#define	LO_CHAN			0
#define	HI_CHAN			15
#define	UPDATE_THLD		35

#define	SINE_START_CHAN		2
#define	P2A_THRESH		10.0F
#define	DEV_THLD1		28.0F

/* Defines for the FFT function */
#define	SIZE			128
#define	SIZE_BY_TWO		64
#define	NUM_STAGE		6

#define	PI			3.141592653589793F

#define	TRUE			1
#define	FALSE			0

/* Macros */
#define	min(a,b)		((a)<(b)?(a):(b))
#define	max(a,b)		((a)>(b)?(a):(b))
#define	square(a)		((a)*(a))

/* structures */
typedef struct {
	float pre_emp_mem;
	int16_t update_cnt;
	int16_t hyster_cnt;
	int16_t last_update_cnt;
	float ch_enrg_long_db[NUM_CHAN];
	int32_t Lframe_cnt;
	float ch_enrg[NUM_CHAN];
	float ch_noise[NUM_CHAN];
	float tsnr;
	int16_t hangover;
	int16_t burstcount;
	int16_t fupdate_flag;
	float negSNRvar;
	float negSNRbias;
	float R0;
	float Rmax;
	int16_t LTP_flag;
} vadState;
#else
typedef struct {
	float bckr_est[COMPLEN];	/* background noise estimate */
	float ave_level[COMPLEN];

	/* averaged input components for stationary estimation */
	float old_level[COMPLEN];	/* input levels of the previous frame */
	float sub_level[COMPLEN];

	/* input levels calculated at the end of a frame (lookahead) */
	float a_data5[3][2];	/* memory for the filter bank */
	float a_data3[5];	/* memory for the filter bank */
	float best_corr_hp;	/* FIP filtered value */

	/* counts length of a speech burst incl HO addition */
	float corr_hp_fast;	/* filtered value */
	int32_t vadreg;		/* flags for intermediate VAD decisions */
	int32_t pitch;		/* flags for pitch detection */
	int32_t oldlag_count, oldlag;	/* variables for pitch detection */
	int32_t complex_high;	/* flags for complex detection */
	int32_t complex_low;	/* flags for complex detection */
	int32_t complex_warning;	/* complex background warning */
	int32_t tone;		/* flags for tone detection */
	int16_t burst_count;	/* counts length of a speech burst */
	int16_t hang_count;	/* hangover counter */
	int16_t stat_count;	/* stationary counter */
	int16_t complex_hang_count;	/* complex hangover counter, used by VAD */
	int16_t complex_hang_timer;	/* hangover initiator, used by CAD */
	int16_t speech_vad_decision;	/* final decision */
	int16_t sp_burst_count;

} vadState;
#endif
#define DTX_HIST_SIZE 8
#define DTX_ELAPSED_FRAMES_THRESH (24 + 7 -1)
#define DTX_HANG_CONST 7	/* yields eight frames of SP HANGOVER */


/*
typedef struct {
	float lsp_hist[M * DTX_HIST_SIZE];
	float log_en_hist[DTX_HIST_SIZE];
	int32_t init_lsf_vq_index;
	int16_t hist_ptr;
	int16_t log_en_index;
	int16_t lsp_index[3];

	// DTX handler stuff 
	int16_t dtxHangoverCount;
	int16_t decAnaElapsedCount;

} dtx_encState;

 */

typedef struct {
	/* gain history */
	float gp[N_FRAME];

	/* counters */
	int16_t count;

} tonStabState;
typedef struct {
	int32_t past_qua_en[4];

	/* normal MA predictor memory, (contains 20*log10(qua_err)) */
} gc_predState;

typedef struct {
	float prev_alpha;	/* previous adaptor output, */
	float prev_gc;	/* previous code gain, */
	float ltpg_mem[LTPG_MEM_SIZE];	/* LTP coding gain history, */
	int16_t onset;		/* onset state, */

	/* (ltpg_mem[0] not used for history) */
} gain_adaptState;
typedef struct {

	float sf0_target_en;
	float sf0_coeff[5];
	int32_t sf0_gcode0_exp;
	int32_t sf0_gcode0_fra;
	int16_t *gain_idx_ptr;

	gc_predState *gc_predSt;
	gc_predState *gc_predUncSt;
	gain_adaptState *adaptSt;
} gainQuantState;
typedef struct {
	int32_t T0_prev_subframe;	/* integer pitch lag of previous sub-frame */

} Pitch_frState;
typedef struct {
	Pitch_frState *pitchSt;
} clLtpState;
typedef struct {
	float ada_w;
	int32_t old_T0_med;
	int16_t wght_flg;

} pitchOLWghtState;
typedef struct {
	float past_rq[M];	/* Past quantized prediction error */

} Q_plsfState;
typedef struct {
	/* Past LSPs */
	float lsp_old[M];
	float lsp_old_q[M];

	/* Quantization state */
	Q_plsfState *qSt;
} lspState;
typedef struct {
	float old_A[M + 1];	/* Last A(z) for case of unstable filter */

} LevinsonState;
typedef struct {
	LevinsonState *LevinsonSt;
} lpcState;
typedef struct {
	/* Speech vector */
	float old_speech[L_TOTAL];
	float *speech, *p_window, *p_window_12k2;
	float *new_speech;	/* Global variable */

	/* Weight speech vector */
	float old_wsp[L_FRAME + PIT_MAX];
	float *wsp;

	/* OL LTP states */
	int32_t old_lags[5];
	float ol_gain_flg[2];

	/* Excitation vector */
	float old_exc[L_FRAME + PIT_MAX + L_INTERPOL];
	float *exc;

	/* Zero vector */
	float ai_zero[L_SUBFR + MP1];
	float *zero;

	/* Impulse response vector */
	float *h1;
	float hvec[L_SUBFR * 2];

	/* Substates */
	lpcState *lpcSt;
	lspState *lspSt;
	clLtpState *clLtpSt;
	gainQuantState *gainQuantSt;
	pitchOLWghtState *pitchOLWghtSt;
	tonStabState *tonStabSt;
	vadState *vadSt;

	int32_t dtx;

	//dtx_encState *dtxEncSt;

	/* Filter's memory */
	float mem_syn[M], mem_w0[M], mem_w[M];
	float mem_err[M + L_SUBFR], *error;
	float sharp;

} cod_amrState;
typedef struct {
	cod_amrState *cod_amr_state;
	Pre_ProcessState *pre_state;

	int32_t dtx;

} Speech_Encode_FrameState;

/*
 * Dotproduct40
 *
 *
 * Parameters:
 *    x                 I: First input
 *    y                 I: Second input
 * Function:
 *    Computes dot product size 40
 *
 * Returns:
 *    acc                dot product
 */


static float Dotproduct40(float * x, float * y)
{
	float acc;

	acc = x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
	acc += x[4] * y[4] + x[5] * y[5] + x[6] * y[6] + x[7] * y[7];
	acc += x[8] * y[8] + x[9] * y[9] + x[10] * y[10] + x[11] * y[11];
	acc += x[12] * y[12] + x[13] * y[13] + x[14] * y[14] + x[15] * y[15];
	acc += x[16] * y[16] + x[17] * y[17] + x[18] * y[18] + x[19] * y[19];
	acc += x[20] * y[20] + x[21] * y[21] + x[22] * y[22] + x[23] * y[23];
	acc += x[24] * y[24] + x[25] * y[25] + x[26] * y[26] + x[27] * y[27];
	acc += x[28] * y[28] + x[29] * y[29] + x[30] * y[30] + x[31] * y[31];
	acc += x[32] * y[32] + x[33] * y[33] + x[34] * y[34] + x[35] * y[35];
	acc += x[36] * y[36] + x[37] * y[37] + x[38] * y[38] + x[39] * y[39];
	return (acc);
}






/*
 * Autocorr
 *
 *
 * Parameters:
 *    x                 I: Input signal
 *    r                 O: Autocorrelations
 *    wind              I: Window for LPC analysis
 * Function:
 *    Calculate autocorrelation with window, LPC order = M
 *
 * Returns:
 *    void
 */
static void Autocorr(float x[], float r[], const float wind[])
{
	int32_t i, j;		/* Counters */
	float y[L_WINDOW + M + 1];	/* Windowed signal */
	float sum;		/* temp */

	/*
	 * Windowing of signal
	 */
	for (i = 0; i < L_WINDOW; i++) {
		y[i] = x[i] * wind[i];
	}

	/*
	 * Zero remaining memory
	 */
	memset(&y[L_WINDOW], 0, 44);

	/*
	 * Autocorrelation
	 */
	for (i = 0; i <= M; i++) {
		sum = 0;

		for (j = 0; j < L_WINDOW; j += 40) {
			sum += Dotproduct40(&y[j], &y[j + i]);
		}
		r[i] = (float) sum;
	}
}

/*
 * Levinson
 *
 *
 * Parameters:
 *    old_A             I: Vector of old LP coefficients [M+1]
 *    r                 I: Vector of autocorrelations    [M+1]
 *    a                 O: LP coefficients               [M+1]
 *    rc                O: Reflection coefficients       [4]
 * Function:
 *    Levinson-Durbin algorithm
 *
 * Returns:
 *    void
 *
 */
static void Levinson(float * old_A, float * r, float * A, float * rc)
{
	float sum, at, err;
	int32_t l, j, i;
	float rct[M];		/* temporary reflection coefficients  0,...,m-1 */

	rct[0] = (-r[1]) / r[0];
	A[0] = 1.0F;
	A[1] = rct[0];
	err = r[0] + r[1] * rct[0];

	if (err <= 0.0F)
		err = 0.01F;

	for (i = 2; i <= M; i++) {
		sum = 0.0F;

		for (j = 0; j < i; j++)
			sum += r[i - j] * A[j];
		rct[i - 1] = (-sum) / (err);

		for (j = 1; j <= (i / 2); j++) {
			l = i - j;
			at = A[j] + rct[i - 1] * A[l];
			A[l] += rct[i - 1] * A[j];
			A[j] = at;
		}
		A[i] = rct[i - 1];
		err += rct[i - 1] * sum;

		if (err <= 0.0F)
			err = 0.01F;
	}
	memcpy(rc, rct, 4 * sizeof(float));
	memcpy(old_A, A, MP1 * sizeof(float));
}

/*
 * lpc
 *
 *
 * Parameters:
 *    old_A             O: Vector of old LP coefficients [M+1]
 *    x                 I: Input signal
 *    x_12k2            I: Input signal 12.2k
 *    a                 O: predictor coefficients
 *    mode              I: AMR mode
 * Function:
 *    LP analysis
 *
 *    In 12.2 kbit/s mode linear prediction (LP) analysis is performed
 *    twice per speech frame using the auto-correlation approach with
 *    30 ms asymmetric windows. No lookahead is used in
 *    the auto-correlation computation.
 *
 *    In other modes analysis is performed once per speech frame
 *    using the auto-correlation approach with 30 ms asymmetric windows.
 *    A lookahead of 40 samples (5 ms) is used in the auto-correlation computation.
 *
 *    The auto-correlations of windowed speech are converted to the LP
 *    coefficients using the Levinson-Durbin algorithm.
 *    Then the LP coefficients are transformed to the Line Spectral Pair
 *    (LSP) domain  for quantization and interpolation purposes.
 *    The interpolated quantified and unquantized filter coefficients
 *    are converted back to the LP filter coefficients
 *    (to construct the synthesis and weighting filters at each subframe).
 *
 * Returns:
 *    void
 *
 */
static void lpc(float * old_A, float x[], float x_12k2[], float a[],
		enum Mode mode)
{
	int32_t i;
	float r[MP1];
	float rc[4];
       /*
	if (mode == MR122) {
		Autocorr(x_12k2, r, f_window_160_80);


		 // Lag windowing

		for (i = 1; i <= M; i++) {
			r[i] = r[i] * f_lag_wind[i - 1];
		}
		r[0] *= 1.0001F;

		if (r[0] < 1.0F)
			r[0] = 1.0F;


		 // Levinson Durbin

		Levinson(old_A, r, &a[MP1], rc);


		 // Autocorrelations

		Autocorr(x_12k2, r, f_window_232_8);


		 // Lag windowing

		for (i = 1; i <= M; i++) {
			r[i] = r[i] * f_lag_wind[i - 1];
		}
		r[0] *= 1.0001F;

		if (r[0] < 1.0F)
			r[0] = 1.0F;


		 // Levinson Durbin

		Levinson(old_A, r, &a[MP1 * 3], rc);
	} else


         */
        {
		/*
		 * Autocorrelations
		 */
		Autocorr(x, r, f_window_200_40);

		/*
		 * a 60 Hz bandwidth expansion is used by lag windowing
		 * the auto-correlations. Further, auto-correlation[0] is
		 * multiplied by the white noise correction factor 1.0001
		 * which is equivalent to adding a noise floor at -40 dB.
		 */
		for (i = 1; i <= M; i++) {
			r[i] = r[i] * f_lag_wind[i - 1];
		}
		r[0] *= 1.0001F;

		if (r[0] < 1.0F)
			r[0] = 1.0F;

		/*
		 * Levinson Durbin
		 */
		Levinson(old_A, r, &a[MP1 * 3], rc);
	}
}

/*
 * Chebps
 *
 *
 * Parameters:
 *    x                 I: Cosine value for the freqency given
 *    f                 I: angular freqency
 * Function:
 *    Evaluates the Chebyshev polynomial series
 *
 * Returns:
 *    result of polynomial evaluation
 */
static float Chebps(float x, float f[])
{
	float b0, b1, b2, x2;
	int32_t i;

	x2 = 2.0F * x;
	b2 = 1.0F;
	b1 = x2 + f[1];

	for (i = 2; i < 5; i++) {
		b0 = x2 * b1 - b2 + f[i];
		b2 = b1;
		b1 = b0;
	}
	return (x * b1 - b2 + f[i]);
}

/*
 * Az_lsp
 *
 *
 * Parameters:
 *    a                 I: Predictor coefficients              [MP1]
 *    lsp               O: Line spectral pairs                 [M]
 *    old_lsp           I: Old lsp, in case not found 10 roots [M]
 *
 * Function:
 *    LP to LSP conversion
 *
 *    The LP filter coefficients A[], are converted to the line spectral pair
 *    (LSP) representation for quantization and interpolation purposes.
 *
 * Returns:
 *    void
 */
static void Az_lsp(float a[], float lsp[], float old_lsp[])
{
	int32_t i, j, nf, ip;
	float xlow, ylow, xhigh, yhigh, xmid, ymid, xint;
	float y;
	float *coef;
	float f1[6], f2[6];

	/*
	 *  find the sum and diff. pol. F1(z) and F2(z)
	 */
	f1[0] = 1.0F;
	f2[0] = 1.0F;

	for (i = 0; i < (NC); i++) {
		f1[i + 1] = a[i + 1] + a[M - i] - f1[i];
		f2[i + 1] = a[i + 1] - a[M - i] + f2[i];
	}
	f1[NC] *= 0.5F;
	f2[NC] *= 0.5F;

	/*
	 * find the LSPs using the Chebychev pol. evaluation
	 */
	nf = 0;			/* number of found frequencies */
	ip = 0;			/* indicator for f1 or f2 */
	coef = f1;
	xlow = f_grid[0];
	ylow = Chebps(xlow, coef);
	j = 0;

	while ((nf < M) && (j < 60)) {
		j++;
		xhigh = xlow;
		yhigh = ylow;
		xlow = f_grid[j];
		ylow = Chebps(xlow, coef);

		if (ylow * yhigh <= 0) {
			/* divide 4 times the interval */
			for (i = 0; i < 4; i++) {
				xmid = (xlow + xhigh) * 0.5F;
				ymid = Chebps(xmid, coef);

				if (ylow * ymid <= 0.0F) {
					yhigh = ymid;
					xhigh = xmid;
				} else {
					ylow = ymid;
					xlow = xmid;
				}
			}

			/*
			 * Linear interpolation
			 * xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow)
			 */
			y = yhigh - ylow;

			if (y == 0) {
				xint = xlow;
			} else {
				y = (xhigh - xlow) / (yhigh - ylow);
				xint = xlow - ylow * y;
			}
			lsp[nf] = xint;
			xlow = xint;
			nf++;

			if (ip == 0) {
				ip = 1;
				coef = f2;
			} else {
				ip = 0;
				coef = f1;
			}
			ylow = Chebps(xlow, coef);
		}
	}

	/* Check if M roots found */
	if (nf < M) {
		memcpy(lsp, old_lsp, M << 2);
	}
	return;
}

/*
 * Get_lsp_pol
 *
 *
 * Parameters:
 *    lsp                 I: line spectral frequencies
 *    f                   O: polynomial F1(z) or F2(z)
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
static void Get_lsp_pol(float * lsp, float * f)
{
	int32_t i, j;
	float T0;

	f[0] = 1.0F;
	f[1] = -2.0F * lsp[0];

	for (i = 2; i <= 5; i++) {
		T0 = -2.0F * lsp[2 * i - 2];
		f[i] = (float) (T0 * f[i - 1] + 2.0F * f[i - 2]);

		for (j = i - 1; j >= 2; j--) {
			f[j] = f[j] + T0 * f[j - 1] + f[j - 2];
		}
		f[1] = f[1] + T0;
	}
	return;
}

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
 * Returns:
 *    void
 */
static void Lsp_Az(float lsp[], float a[])
{
	float f1[6], f2[6];
	int32_t i, j;

	Get_lsp_pol(&lsp[0], f1);
	Get_lsp_pol(&lsp[1], f2);

	for (i = 5; i > 0; i--) {
		f1[i] += f1[i - 1];
		f2[i] -= f2[i - 1];
	}
	a[0] = 1;

	for (i = 1, j = 10; i <= 5; i++, j--) {
		a[i] = (float) ((f1[i] + f2[i]) * 0.5F);
		a[j] = (float) ((f1[i] - f2[i]) * 0.5F);
	}
	return;
}


/*
 * Lsp_lsf
 *
 *
 * Parameters:
 *    lsp               I: LSP vector
 *    lsf               O: LSF vector
 *
 * Function:
 *    Transformation lsp to lsf, LPC order M
 *
 * Returns:
 *    void
 */
static void Lsp_lsf(float lsp[], float lsf[])
{
	int32_t i;

	for (i = 0; i < M; i++) {
		lsf[i] = (float) (acosf(lsp[i]) * SCALE_LSP_FREQ);
	}
	return;
}

/*
 * Lsf_wt
 *
 *
 * Parameters:
 *    lsf               I: LSF vector
 *    wf                O: square of weighting factors
 *
 * Function:
 *    Compute LSF weighting factors
 *
 * Returns:
 *    void
 */
static void Lsf_wt(float * lsf, float * wf)
{
	float temp;
	int32_t i;

	wf[0] = lsf[1];

	for (i = 1; i < 9; i++) {
		wf[i] = lsf[i + 1] - lsf[i - 1];
	}
	wf[9] = 4000.0F - lsf[8];

	for (i = 0; i < 10; i++) {
		if (wf[i] < 450.0F) {
			temp = 3.347F - SLOPE1_WGHT_LSF * wf[i];
		} else {
			temp = 1.8F - SLOPE2_WGHT_LSF * (wf[i] - 450.0F);
		}
		wf[i] = temp * temp;
	}
	return;
}

/*
 * Reorder_lsf
 *
 *
 * Parameters:
 *    lsf               B: vector of LSFs
 *    min_dist          I: minimum required distance
 *
 * Function:
 *    Make sure that the LSFs are properly ordered and to keep a certain minimum
 *    distance between adjacent LSFs. LPC order = M.
 *
 * Returns:
 *    void
 */
static void Reorder_lsf(float * lsf, float min_dist)
{
	float lsf_min;
	int32_t i;

	lsf_min = min_dist;

	for (i = 0; i < M; i++) {
		if (lsf[i] < lsf_min) {
			lsf[i] = lsf_min;
		}
		lsf_min = lsf[i] + min_dist;
	}
}

/*
 * Lsf_lsp
 *
 *
 * Parameters:
 *    lsf               I: vector of LSFs
 *    lsp	            O: vector of LSPs
 *
 * Function:
 *    Transformation lsf to lsp, order M
 *
 * Returns:
 *    void
 */
static void Lsf_lsp(float lsf[], float lsp[])
{
	int32_t i;

	for (i = 0; i < M; i++) {
		lsp[i] = (float) cosf(SCALE_FREQ_LSP * lsf[i]);
	}
	return;
}

/*
 * Vq_subvec3
 *
 *
 * Parameters:
 *    lsf_r1            I: 1st LSF residual vector
 *    dico              I: quantization codebook
 *    wf1               I: 1st LSF weighting factors
 *    dico_size         I: size of quantization codebook
 *    use_half          I: use every second entry in codebook
 *
 * Function:
 *    Quantization of a 3 dimensional subvector
 *
 * Returns:
 *    index             quantization index
 */
static int16_t Vq_subvec3(float * lsf_r1, const float * dico, float * wf1,
			 int16_t dico_size, int32_t use_half)
{
	float dist, dist_min;
	float temp;
	const float *p_dico;
	int32_t i, index = 0;

	dist_min = FLT_MAX;
	p_dico = dico;

	if (use_half == 0) {
		for (i = 0; i < dico_size; i++) {
			temp = lsf_r1[0] - *p_dico++;
			temp *= wf1[0];
			dist = temp * temp;
			temp = lsf_r1[1] - *p_dico++;
			temp *= wf1[1];
			dist += temp * temp;
			temp = lsf_r1[2] - *p_dico++;
			temp *= wf1[2];
			dist += temp * temp;

			if (dist < dist_min) {
				dist_min = dist;
				index = i;
			}
		}
		p_dico = &dico[(3 * index)];
	} else {
		for (i = 0; i < dico_size; i++) {
			temp = lsf_r1[0] - *p_dico++;
			temp *= wf1[0];
			dist = temp * temp;
			temp = lsf_r1[1] - *p_dico++;
			temp *= wf1[1];
			dist += temp * temp;
			temp = lsf_r1[2] - *p_dico++;
			temp *= wf1[2];
			dist += temp * temp;

			if (dist < dist_min) {
				dist_min = dist;
				index = i;
			}
			p_dico = p_dico + 3;
		}
		p_dico = &dico[6 * index];
	}

	/* Reading the selected vector */
	lsf_r1[0] = *p_dico++;
	lsf_r1[1] = *p_dico++;
	lsf_r1[2] = *p_dico++;
	return (int16_t) index;
}

/*
 * Vq_subvec4
 *
 *
 * Parameters:
 *    lsf_r1            I: 1st LSF residual vector
 *    dico              I: quantization codebook
 *    wf1               I: 1st LSF weighting factors
 *    dico_size         I: size of quantization codebook
 *
 * Function:
 *    Quantization of a 4 dimensional subvector
 *
 * Returns:
 *    index             quantization index
 */
static int16_t Vq_subvec4(float * lsf_r1, const float * dico, float * wf1,
			 int16_t dico_size)
{
	float dist, dist_min;
	float temp;
	const float *p_dico;
	int32_t i, index = 0;

	dist_min = FLT_MAX;
	p_dico = dico;

	for (i = 0; i < dico_size; i++) {
		temp = lsf_r1[0] - *p_dico++;
		temp *= wf1[0];
		dist = temp * temp;
		temp = lsf_r1[1] - *p_dico++;
		temp *= wf1[1];
		dist += temp * temp;
		temp = lsf_r1[2] - *p_dico++;
		temp *= wf1[2];
		dist += temp * temp;
		temp = lsf_r1[3] - *p_dico++;
		temp *= wf1[3];
		dist += temp * temp;

		if (dist < dist_min) {
			dist_min = dist;
			index = i;
		}
	}

	/* Reading the selected vector */
	p_dico = &dico[index << 2];
	lsf_r1[0] = *p_dico++;
	lsf_r1[1] = *p_dico++;
	lsf_r1[2] = *p_dico++;
	lsf_r1[3] = *p_dico++;
	return (int16_t) index;
}

/*
 * Q_plsf_3
 *
 *
 * Parameters:
 *    mode              I: AMR mode
 *    past_rq           B: past quantized residual
 *    lsp1              I: 1st LSP vector
 *    lsp1_q            O: quantized 1st LSP vector
 *    indice            I: quantization indices of 5 matrices and
 *                         one sign for 3rd
 *    pred_init_i       O: init index for MA prediction in DTX mode
 *
 * Function:
 *    Quantization of LSF parameters with 1st order MA prediction and
 *    split by 3 vector quantization (split-VQ)
 *
 * Returns:
 *    void
 */
static void Q_plsf_3(enum Mode mode, float * past_rq, float * lsp1,
		     float * lsp1_q, int16_t * indice, int32_t * pred_init_i)
{
	float lsf1[M], wf1[M], lsf_p[M], lsf_r1[M];
	float lsf1_q[M];
//	float pred_init_err;
//	float min_pred_init_err;
//	float temp_r1[M];
//	float temp_p[M];
	int32_t i;

	/* convert LSFs to normalize frequency domain */
	Lsp_lsf(lsp1, lsf1);

	/* compute LSF weighting factors */
	Lsf_wt(lsf1, wf1);

	/* Compute predicted LSF and prediction error */
	if (mode != MRDTX)
        {
		for (i = 0; i < M; i++) {
			lsf_p[i] = f_mean_lsf_3[i] + past_rq[i] * f_pred_fac[i];
			lsf_r1[i] = lsf1[i] - lsf_p[i];
		}
	}

        /*
        else {
		//
		 // DTX mode, search the init vector that yields
		// lowest prediction resuidual energy
		 //
		*pred_init_i = 0;
		min_pred_init_err = FLT_MAX;

		for (j = 0; j < PAST_RQ_INIT_SIZE; j++) {
			pred_init_err = 0;

			for (i = 0; i < M; i++) {
				temp_p[i] =
				    f_mean_lsf_3[i] + f_past_rq_init[j * M + i];
				temp_r1[i] = lsf1[i] - temp_p[i];
				pred_init_err += temp_r1[i] * temp_r1[i];
			}	// next i

			if (pred_init_err < min_pred_init_err) {
				min_pred_init_err = pred_init_err;
				memcpy(lsf_r1, temp_r1, M << 2);
				memcpy(lsf_p, temp_p, M << 2);
				memcpy(past_rq, &f_past_rq_init[j * M], M << 2);
				*pred_init_i = j;
			}
		}



	}

         */

	/* Split-VQ of prediction error */
	/* MR475, MR515 */
	if ((mode == MR475)) {
		indice[0] =
		    Vq_subvec3(&lsf_r1[0], f_dico1_lsf_3, &wf1[0], DICO1_SIZE_3,
			       0);
		indice[1] =
		    Vq_subvec3(&lsf_r1[3], f_dico2_lsf_3, &wf1[3],
			       DICO2_SIZE_3 / 2, 1);
		indice[2] =
		    Vq_subvec4(&lsf_r1[6], f_mr515_3_lsf, &wf1[6], MR515_3_SIZE);
	}

        /*
	// MR795
	else if (mode == MR795) {
		indice[0] =
		    Vq_subvec3(&lsf_r1[0], f_mr795_1_lsf, &wf1[0], MR795_1_SIZE,
			       0);
		indice[1] =
		    Vq_subvec3(&lsf_r1[3], f_dico2_lsf_3, &wf1[3], DICO2_SIZE_3,
			       0);
		indice[2] =
		    Vq_subvec4(&lsf_r1[6], f_dico3_lsf_3, &wf1[6], DICO3_SIZE_3);
	}

	// MR59, MR67, MR74, MR102 , MRDTX
	else {
		indice[0] =
		    Vq_subvec3(&lsf_r1[0], f_dico1_lsf_3, &wf1[0], DICO1_SIZE_3,
			       0);
		indice[1] =
		    Vq_subvec3(&lsf_r1[3], f_dico2_lsf_3, &wf1[3], DICO2_SIZE_3,
			       0);
		indice[2] =
		    Vq_subvec4(&lsf_r1[6], f_dico3_lsf_3, &wf1[6], DICO3_SIZE_3);
	}

         */


	/* Compute quantized LSFs and update the past quantized residual */
	for (i = 0; i < M; i++) {
		lsf1_q[i] = lsf_r1[i] + lsf_p[i];
		past_rq[i] = lsf_r1[i];
	}

	/* verification that LSFs has mimimum distance of LSF_GAP 50 Hz */
	Reorder_lsf(lsf1_q, 50.0F);

	/*  convert LSFs to the cosine domain */
	Lsf_lsp(lsf1_q, lsp1_q);
}


/*
 * Int_lpc_1to3_2
 *
 *
 * Parameters:
 *    lsp_old           I: LSP vector at the 4th subfr. of past frame      [M]
 *    lsp_new           I: LSP vector at the 4th subframe of present frame [M]
 *    az                O: interpolated LP parameters in subframes 1, 2 and 3
 *                                                                   [AZ_SIZE]
 *
 * Function:
 *    Interpolation of the LPC parameters.
 *
 * Returns:
 *    void
 */
static void Int_lpc_1to3_2(float lsp_old[], float lsp_new[], float az[])
{
	float lsp[M];
	int32_t i;

	for (i = 0; i < M; i += 2) {
		lsp[i] = lsp_new[i] * 0.25F + lsp_old[i] * 0.75F;
		lsp[i + 1] = lsp_new[i + 1] * 0.25F + lsp_old[i + 1] * 0.75F;
	}

	/* Subframe 1 */
	Lsp_Az(lsp, az);
	az += MP1;

	for (i = 0; i < M; i += 2) {
		lsp[i] = (lsp_old[i] + lsp_new[i]) * 0.5F;
		lsp[i + 1] = (lsp_old[i + 1] + lsp_new[i + 1]) * 0.5F;
	}

	/* Subframe 2 */
	Lsp_Az(lsp, az);
	az += MP1;

	for (i = 0; i < M; i += 2) {
		lsp[i] = lsp_old[i] * 0.25F + lsp_new[i] * 0.75F;
		lsp[i + 1] = lsp_old[i + 1] * 0.25F + lsp_new[i + 1] * 0.75F;
	}

	/* Subframe 3 */
	Lsp_Az(lsp, az);
	return;
}

/*
 * Int_lpc_1to3
 *
 *
 * Parameters:
 *    lsp_old           I: LSP vector at the 4th subfr. of past frame      [M]
 *    lsp_new           I: LSP vector at the 4th subframe of present frame [M]
 *    az                O: interpolated LP parameters in all subframes
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
static void Int_lpc_1to3(float lsp_old[], float lsp_new[], float az[])
{
	float lsp[M];
	int32_t i;

	for (i = 0; i < M; i++) {
		lsp[i] = lsp_new[i] * 0.25F + lsp_old[i] * 0.75F;
	}

	/* Subframe 1 */
	Lsp_Az(lsp, az);
	az += MP1;

	for (i = 0; i < M; i++) {
		lsp[i] = (lsp_old[i] + lsp_new[i]) * 0.5F;
	}

	/* Subframe 2 */
	Lsp_Az(lsp, az);
	az += MP1;

	for (i = 0; i < M; i++) {
		lsp[i] = lsp_old[i] * 0.25F + lsp_new[i] * 0.75F;
	}

	/* Subframe 3 */
	Lsp_Az(lsp, az);
	az += MP1;

	/* Subframe 4 */
	Lsp_Az(lsp_new, az);
	return;
}

/*
 * lsp
 *
 *
 * Parameters:
 *    req_mode          I: requested mode
 *    used_mode         I: used mode
 *    lsp_old           B: old LSP vector
 *    lsp_old_q         B: old quantized LSP vector
 *    past_rq           B: past quantized residual
 *    az                B: interpolated LP parameters
 *    azQ               O: quantization interpol. LP parameters
 *    lsp_new           O: new lsp vector
 *    anap              O: analysis parameters
 *
 * Function:
 *    From A(z) to lsp. LSP quantization and interpolation
 *
 * Returns:
 *    void
 */
static void lsp(enum Mode req_mode, enum Mode used_mode, float * lsp_old,
		float * lsp_old_q, float * past_rq, float az[],
		float azQ[], float lsp_new[], int16_t ** anap)
{
	float lsp_new_q[M];	/* LSPs at 4th subframe */
//	float lsp_mid[M], lsp_mid_q[M];	/* LSPs at 2nd subframe */
	int32_t pred_init_i;	/* init index for MA prediction in DTX mode */


        /*
	if (req_mode == MR122)
        {
		Az_lsp(&az[MP1], lsp_mid, lsp_old);
		Az_lsp(&az[MP1 * 3], lsp_new, lsp_mid);

		//
		 // Find interpolated LPC parameters in all subframes
		 // (both quantized and unquantized).
		 // The interpolated parameters are in array A_t[] of size (M+1)*4
		 // and the quantized interpolated parameters are in array Aq_t[]

		Int_lpc_1and3_2(lsp_old, lsp_mid, lsp_new, az);

		if (used_mode != MRDTX) {
		// LSP quantization (lsp_mid[] and lsp_new[] jointly quantized)
			Q_plsf_5(past_rq, lsp_mid, lsp_new, lsp_mid_q,
				 lsp_new_q, *anap);
			Int_lpc_1and3(lsp_old_q, lsp_mid_q, lsp_new_q, azQ);

			// Advance analysis parameters pointer
			(*anap) += 5;
		}
	} else
        */


        {
		/* From A(z) to lsp */
		Az_lsp(&az[MP1 * 3], lsp_new, lsp_old);

		/*
		 * Find interpolated LPC parameters in all subframes
		 * (both quantized and unquantized).
		 * The interpolated parameters are in array A_t[] of size (M+1)*4
		 * and the quantized interpolated parameters are in array Aq_t[]
		 */
		Int_lpc_1to3_2(lsp_old, lsp_new, az);

		/* LSP quantization */
		if (used_mode != MRDTX) {
			Q_plsf_3(req_mode, past_rq, lsp_new, lsp_new_q, *anap,
				 &pred_init_i);
			Int_lpc_1to3(lsp_old_q, lsp_new_q, azQ);

			/* Advance analysis parameters pointer */
			(*anap) += 3;
		}
	}

	/* update the LSPs for the next frame */
	memcpy(lsp_old, lsp_new, M << 2);
	memcpy(lsp_old_q, lsp_new_q, M << 2);
}

/*
 * check_lsp
 *
 *
 * Parameters:
 *    count          B: counter for resonance
 *    lsp            B: LSP vector
 *
 * Function:
 *    Check the LSP's to detect resonances
 *
 *    Resonances in the LPC filter are monitored to detect possible problem
 *    areas where divergence between the adaptive codebook memories in
 *    the encoder and the decoder could cause unstable filters in areas
 *    with highly correlated continuos signals. Typically, this divergence
 *    is due to channel errors.
 *    The monitoring of resonance signals is performed using unquantized LSPs
 *    q(i), i = 1,...,10. The algorithm utilises the fact that LSPs are
 *    closely located at a peak in the spectrum. First, two distances,
 *    dist 1 and dist 2 ,are calculated in two different regions,
 *    defined as
 *
 *    dist1 = min[q(i) - q(i + 1)],  i = 4,...,8
 *    dist2 = min[q(i) - q(i + 1)],  i = 2,3
 *
 *    Either of these two minimum distance conditions must be fulfilled
 *    to classify the frame as a resonance frame and increase the resonance
 *    counter.
 *
 *    if(dist1 < TH1) || if (dist2 < TH2)
 *       counter++
 *    else
 *       counter = 0
 *
 *    TH1 = 0.046
 *    TH2 = 0.018, q(2) > 0.98
 *    TH2 = 0.024, 0.93 < q(2) <= 0.98
 *    TH2 = 0.018, otherwise
 *
 *    12 consecutive resonance frames are needed to indicate possible
 *    problem conditions, otherwise the LSP_flag is cleared.
 *
 * Returns:
 *    resonance flag
 */
static int16_t check_lsp(int16_t * count, float * lsp)
{
	float dist, dist_min1, dist_min2, dist_th;
	int32_t i;

	/*
	 * Check for a resonance:
	 * Find minimum distance between lsp[i] and lsp[i+1]
	 */
	dist_min1 = FLT_MAX;

	for (i = 3; i < 8; i++) {
		dist = lsp[i] - lsp[i + 1];

		if (dist < dist_min1) {
			dist_min1 = dist;
		}
	}
	dist_min2 = FLT_MAX;

	for (i = 1; i < 3; i++) {
		dist = lsp[i] - lsp[i + 1];

		if (dist < dist_min2) {
			dist_min2 = dist;
		}
	}

	if (lsp[1] > 0.98F) {
		dist_th = 0.018F;
	} else if (lsp[1] > 0.93F) {
		dist_th = 0.024F;
	} else {
		dist_th = 0.034F;
	}

	if ((dist_min1 < 0.046F) || (dist_min2 < dist_th)) {
		*count += 1;
	} else {
		*count = 0;
	}

	/* Need 12 consecutive frames to set the flag */
	if (*count >= 12) {
		*count = 12;
		return 1;
	} else {
		return 0;
	}
}

/*
 * Weight_Ai
 *
 *
 * Parameters:
 *    a                 I: LPC coefficients                    [M+1]
 *    fac               I: Spectral expansion factors.         [M+1]
 *    a_exp             O: Spectral expanded LPC coefficients  [M+1]
 *
 * Function:
 *    Spectral expansion of LP coefficients
 *
 * Returns:
 *    void
 */
static void Weight_Ai(float a[], const float fac[], float a_exp[])
{
	int32_t i;

	a_exp[0] = a[0];

	for (i = 1; i <= M; i++) {
		a_exp[i] = a[i] * fac[i - 1];
	}
	return;
}

/*
 * Residu
 *
 *
 * Parameters:
 *    a                 I: prediction coefficients
 *    x                 I: speech signal
 *    y                 O: residual signal
 *
 * Function:
 *    Computes the LTP residual signal.
 *
 * Returns:
 *    void
 */
static void Residu(float a[], float x[], float y[])
{
	float s;
	int32_t i;

	for (i = 0; i < L_SUBFR; i += 4) {
		s = x[i] * a[0];
		s += x[i - 1] * a[1];
		s += x[i - 2] * a[2];
		s += x[i - 3] * a[3];
		s += x[i - 4] * a[4];
		s += x[i - 5] * a[5];
		s += x[i - 6] * a[6];
		s += x[i - 7] * a[7];
		s += x[i - 8] * a[8];
		s += x[i - 9] * a[9];
		s += x[i - 10] * a[10];
		y[i] = s;
		s = x[i + 1] * a[0];
		s += x[i] * a[1];
		s += x[i - 1] * a[2];
		s += x[i - 2] * a[3];
		s += x[i - 3] * a[4];
		s += x[i - 4] * a[5];
		s += x[i - 5] * a[6];
		s += x[i - 6] * a[7];
		s += x[i - 7] * a[8];
		s += x[i - 8] * a[9];
		s += x[i - 9] * a[10];
		y[i + 1] = s;
		s = x[i + 2] * a[0];
		s += x[i + 1] * a[1];
		s += x[i] * a[2];
		s += x[i - 1] * a[3];
		s += x[i - 2] * a[4];
		s += x[i - 3] * a[5];
		s += x[i - 4] * a[6];
		s += x[i - 5] * a[7];
		s += x[i - 6] * a[8];
		s += x[i - 7] * a[9];
		s += x[i - 8] * a[10];
		y[i + 2] = s;
		s = x[i + 3] * a[0];
		s += x[i + 2] * a[1];
		s += x[i + 1] * a[2];
		s += x[i] * a[3];
		s += x[i - 1] * a[4];
		s += x[i - 2] * a[5];
		s += x[i - 3] * a[6];
		s += x[i - 4] * a[7];
		s += x[i - 5] * a[8];
		s += x[i - 6] * a[9];
		s += x[i - 7] * a[10];
		y[i + 3] = s;
	}
	return;
}

/*
 * Syn_filt
 *
 *
 * Parameters:
 *    a                 I: prediction coefficients [M+1]
 *    x                 I: input signal
 *    y                 O: output signal
 *    mem               B: memory associated with this filtering
 *    update            I: 0=no update, 1=update of memory.
 *
 * Function:
 *    Perform synthesis filtering through 1/A(z).
 *
 * Returns:
 *    void
 */
static void Syn_filt(float a[], float x[], float y[], float mem[],
		     int16_t update)
{
	float tmp[50];
	float sum;
	float *yy;
	int32_t i;

	/* Copy mem[] to yy[] */
	yy = tmp;

	for (i = 0; i < M; i++) {
		*yy++ = mem[i];
	}

	/* Do the filtering. */
	for (i = 0; i < L_SUBFR; i = i + 4) {
		sum = x[i] * a[0];
		sum -= a[1] * yy[-1];
		sum -= a[2] * yy[-2];
		sum -= a[3] * yy[-3];
		sum -= a[4] * yy[-4];
		sum -= a[5] * yy[-5];
		sum -= a[6] * yy[-6];
		sum -= a[7] * yy[-7];
		sum -= a[8] * yy[-8];
		sum -= a[9] * yy[-9];
		sum -= a[10] * yy[-10];
		*yy++ = sum;
		y[i] = (float) yy[-1];
		sum = x[i + 1] * a[0];
		sum -= a[1] * yy[-1];
		sum -= a[2] * yy[-2];
		sum -= a[3] * yy[-3];
		sum -= a[4] * yy[-4];
		sum -= a[5] * yy[-5];
		sum -= a[6] * yy[-6];
		sum -= a[7] * yy[-7];
		sum -= a[8] * yy[-8];
		sum -= a[9] * yy[-9];
		sum -= a[10] * yy[-10];
		*yy++ = sum;
		y[i + 1] = (float) yy[-1];
		sum = x[i + 2] * a[0];
		sum -= a[1] * yy[-1];
		sum -= a[2] * yy[-2];
		sum -= a[3] * yy[-3];
		sum -= a[4] * yy[-4];
		sum -= a[5] * yy[-5];
		sum -= a[6] * yy[-6];
		sum -= a[7] * yy[-7];
		sum -= a[8] * yy[-8];
		sum -= a[9] * yy[-9];
		sum -= a[10] * yy[-10];
		*yy++ = sum;
		y[i + 2] = (float) yy[-1];
		sum = x[i + 3] * a[0];
		sum -= a[1] * yy[-1];
		sum -= a[2] * yy[-2];
		sum -= a[3] * yy[-3];
		sum -= a[4] * yy[-4];
		sum -= a[5] * yy[-5];
		sum -= a[6] * yy[-6];
		sum -= a[7] * yy[-7];
		sum -= a[8] * yy[-8];
		sum -= a[9] * yy[-9];
		sum -= a[10] * yy[-10];
		*yy++ = sum;
		y[i + 3] = (float) yy[-1];
	}

	/* Update of memory if update==1 */
	if (update != 0) {
		for (i = 0; i < M; i++) {
			mem[i] = y[30 + i];
		}
	}
	return;
}

/*
 * pre_big
 *
 *
 * Parameters:
 *    mode              I: AMR mode
 *    gamma1            I: spectral exp. factor 1
 *    gamma1_12k2       I: spectral exp. factor 1 for modes above MR795
 *    gamma2            I: spectral exp. factor 2
 *    A_t               I: A(z) unquantized, for 4 subframes
 *    frame_offset      I: frameoffset, 1st or second big_sbf
 *    speech            I: speech
 *    mem_w             B: synthesis filter memory state
 *    wsp               O: weighted speech
 *
 * Function:
 *    Big subframe (2 subframes) preprocessing
 *
 *    Open-loop pitch analysis is performed in order to simplify the pitch
 *    analysis and confine the closed-loop pitch search to a small number of
 *    lags around the open-loop estimated lags.
 *    Open-loop pitch estimation is based on the weighted speech signal Sw(n)
 *    which is obtained by filtering the input speech signal through
 *    the weighting filter
 *
 *    W(z) = A(z/g1) / A(z/g2)
 *
 *    That is, in a subframe of size L, the weighted speech is given by:
 *
 *                    10                           10
 *    Sw(n) = S(n) + SUM[a(i) * g1(i) * S(n-i)] - SUM[a(i) * g2(i) * Sw(n-i)],
 *                   i=1                          i=1
 *    n = 0, ..., L-1
 *
 * Returns:
 *    void
 */
static int32_t pre_big(enum Mode mode, const float gamma1[], const float gamma2[], float A_t[],
		      int16_t frame_offset, float speech[], float mem_w[],
		      float wsp[])
{
	float Ap1[MP1], Ap2[MP1];
	int32_t offset, i;

	/* A(z) with spectral expansion */
	const float *g1;

       //	g1 = gamma1_12k2;

	//if (mode <= MR795) {
		g1 = gamma1;
       //	}
	offset = 0;

	if (frame_offset > 0) {
		offset = MP1 << 1;
	}

	/* process two subframes (which form the "big" subframe) */
	for (i = 0; i < 2; i++) {
		/* a(i) * g1(i) */
		Weight_Ai(&A_t[offset], g1, Ap1);

		/* a(i) * g2(i) */
		Weight_Ai(&A_t[offset], gamma2, Ap2);

		/*
		 *       10
		 *  S(n) + SUM[a(i) * g1(i) * S(n-i)]
		 *       i=1
		 */
		Residu(Ap1, &speech[frame_offset], &wsp[frame_offset]);

		/*
		 *          10                            10
		 *  S(n) + SUM[a(i) * g1(i) * S(n-i)]    SUM[a(i) * g2(i) * Sn(n-i)]
		 *         i=1                           i=1
		 */
		Syn_filt(Ap2, &wsp[frame_offset], &wsp[frame_offset], mem_w, 1);
		offset += MP1;
		frame_offset += L_SUBFR;
	}
	return 0;
}

/*
 * comp_corr
 *
 *
 * Parameters:
 *    sig               I: signal
 *    L_frame           I: length of frame to compute pitch
 *    lag_max           I: maximum lag
 *    lag_min           I: minimum lag
 *    corr              O: correlation of selected lag
 *
 * Function:
 *    Calculate all correlations in a given delay range.
 *
 * Returns:
 *    void
 */
static void comp_corr(float sig[], int32_t L_frame, int32_t lag_max, int32_t
		      lag_min, float corr[])
{
	int32_t i, j;
	float *p, *p1;
	float T0;

	for (i = lag_max; i >= lag_min; i--) {
		p = sig;
		p1 = &sig[-i];
		T0 = 0.0F;

		for (j = 0; j < L_frame; j = j + 40, p += 40, p1 += 40) {
			T0 +=
			    p[0] * p1[0] + p[1] * p1[1] + p[2] * p1[2] +
			    p[3] * p1[3];
			T0 +=
			    p[4] * p1[4] + p[5] * p1[5] + p[6] * p1[6] +
			    p[7] * p1[7];
			T0 +=
			    p[8] * p1[8] + p[9] * p1[9] + p[10] * p1[10] +
			    p[11] * p1[11];
			T0 +=
			    p[12] * p1[12] + p[13] * p1[13] + p[14] * p1[14] +
			    p[15] * p1[15];
			T0 +=
			    p[16] * p1[16] + p[17] * p1[17] + p[18] * p1[18] +
			    p[19] * p1[19];
			T0 +=
			    p[20] * p1[20] + p[21] * p1[21] + p[22] * p1[22] +
			    p[23] * p1[23];
			T0 +=
			    p[24] * p1[24] + p[25] * p1[25] + p[26] * p1[26] +
			    p[27] * p1[27];
			T0 +=
			    p[28] * p1[28] + p[29] * p1[29] + p[30] * p1[30] +
			    p[31] * p1[31];
			T0 +=
			    p[32] * p1[32] + p[33] * p1[33] + p[34] * p1[34] +
			    p[35] * p1[35];
			T0 +=
			    p[36] * p1[36] + p[37] * p1[37] + p[38] * p1[38] +
			    p[39] * p1[39];
		}
		corr[-i] = T0;
	}
	return;
}

/*
 * vad_tone_detection
 *
 *
 * Parameters:
 *    st->tone          B: flags indicating presence of a tone
 *    T0                I: autocorrelation maxima
 *    t1                I: energy
 *
 * Function:
 *    Set tone flag if pitch gain is high.
 *    This is used to detect signaling tones and other signals
 *    with high pitch gain.
 *
 * Returns:
 *    void
 */

static void vad_tone_detection(vadState * st, float T0, float t1)
{
	if ((t1 > 0) && (T0 > t1 * TONE_THR)) {
		st->tone = st->tone | 0x00004000;
	}
}


/*
 * Lag_max
 *
 *
 * Parameters:
 *    vadSt          B: vad structure
 *    corr           I: correlation vector
 *    sig            I: signal
 *    L_frame        I: length of frame to compute pitch
 *    lag_max        I: maximum lag
 *    lag_min        I: minimum lag
 *    cor_max        O: maximum correlation
 *    dtx            I: dtx on/off
 *
 * Function:
 *    Compute the open loop pitch lag.
 *
 * Returns:
 *    p_max             lag found
 */

static int16_t Lag_max(vadState * vadSt, float corr[], float sig[], int16_t
		      L_frame, int32_t lag_max, int32_t lag_min,
		      float * cor_max, int32_t dtx)

{
	float max, T0;
	float *p;
	int32_t i, j, p_max;

	max = -FLT_MAX;
	p_max = lag_max;

	for (i = lag_max, j = (PIT_MAX - lag_max - 1); i >= lag_min; i--, j--) {
		if (corr[-i] >= max) {
			max = corr[-i];
			p_max = i;
		}
	}

	/* compute energy for normalization */
	T0 = 0.0F;
	p = &sig[-p_max];

	for (i = 0; i < L_frame; i++, p++) {
		T0 += *p * *p;
	}

	if (dtx) {

		/* check tone */
		vad_tone_detection(vadSt, max, T0);

	}

	if (T0 > 0.0F)
		T0 = 1.0F / (float) sqrtf(T0);
	else
		T0 = 0.0F;

	/* max = max/sqrtf(energy) */
	max *= T0;
	*cor_max = max;
	return ((int16_t) p_max);
}

/*
 * hp_max
 *
 *
 * Parameters:
 *    corr           I: correlation vector
 *    sig            I: signal
 *    L_frame        I: length of frame to compute pitch
 *    lag_max        I: maximum lag
 *    lag_min        I: minimum lag
 *    cor_hp_max     O: max high-pass filtered correlation
 *
 * Function:
 *    Find the maximum correlation of scal_sig[] in a given delay range.
 *
 *    The correlation is given by
 *       cor[t] = <scal_sig[n],scal_sig[n-t]>,  t=lag_min,...,lag_max
 *    The functions outputs the maximum correlation after normalization
 *    and the corresponding lag.
 *
 * Returns:
 *    void
 */

static void hp_max(float corr[], float sig[], int32_t L_frame, int32_t
		   lag_max, int32_t lag_min, float * cor_hp_max)
{
	float T0, t1, max;
	float *p, *p1;
	int32_t i;

	max = -FLT_MAX;

	for (i = lag_max - 1; i > lag_min; i--) {
		/* high-pass filtering */
		T0 = ((corr[-i] * 2) - corr[-i - 1]) - corr[-i + 1];
		T0 = (float) fabsf(T0);

		if (T0 >= max) {
			max = T0;
		}
	}

	/* compute energy */
	p = sig;
	p1 = &sig[0];
	T0 = 0;

	for (i = 0; i < L_frame; i++, p++, p1++) {
		T0 += *p * *p1;
	}
	p = sig;
	p1 = &sig[-1];
	t1 = 0;

	for (i = 0; i < L_frame; i++, p++, p1++) {
		t1 += *p * *p1;
	}

	/* high-pass filtering */
	T0 = T0 - t1;
	T0 = (float) fabsf(T0);

	/* max/T0 */
	if (T0 != 0) {
		*cor_hp_max = max / T0;
	} else {
		*cor_hp_max = 0;
	}
}


/*
 * vad_tone_detection_update
 *
 *
 * Parameters:
 *    st->tone          B: flags indicating presence of a tone
 *    one_lag_per_frame I: 1 open-loop lag is calculated per each frame
 *
 * Function:
 *    Update the tone flag register.
 *
 * Returns:
 *    void
 */

static void vad_tone_detection_update(vadState * st, int16_t one_lag_per_frame)
{
	/* Shift tone flags right by one bit */
	st->tone = st->tone >> 1;

	/*
	 * If open-loop lag is calculated only once in each frame,
	 * do extra update and assume that the other tone flag
	 * of the frame is one.
	 */
	if (one_lag_per_frame != 0) {
		st->tone = st->tone >> 1;
		st->tone = st->tone | 0x00002000;
	}
}


/*
 * Pitch_ol
 *
 *
 * Parameters:
 *    mode           I: AMR mode
 *    vadSt          B: VAD state struct
 *    signal         I: signal used to compute the open loop pitch
 *                                                 [[-pit_max]:[-1]]
 *    pit_min        I: minimum pitch lag
 *    pit_max        I: maximum pitch lag
 *    L_frame        I: length of frame to compute pitch
 *    dtx            I: DTX flag
 *    idx            I: frame index
 *
 * Function:
 *    Compute the open loop pitch lag.
 *
 *    Open-loop pitch analysis is performed twice per frame (each 10 ms)
 *    to find two estimates of the pitch lag in each frame.
 *    Open-loop pitch analysis is performed as follows.
 *    In the first step, 3 maxima of the correlation:
 *
 *          79
 *    O(k) = SUM Sw(n)*Sw(n-k)
 *          n=0
 *
 *    are found in the three ranges:
 *       pit_min     ...      2*pit_min-1
 *       2*pit_min   ...      4*pit_min-1
 *       4*pit_min   ...      pit_max
 *
 *    The retained maxima O(t(i)), i = 1, 2, 3, are normalized by dividing by
 *
 *    SQRT[SUM[POW(Sw(n-t(i)), 2]], i = 1, 2, 3,
 *         n
 *
 *    respectively.
 *    The normalized maxima and corresponding delays are denoted by
 *    (M(i), t(i)), i = 1, 2, 3. The winner, Top, among the three normalized
 *    correlations is selected by favouring the delays with the values
 *    in the lower range. This is performed by weighting the normalized
 *    correlations corresponding to the longer delays. The best
 *    open-loop delay Top is determined as follows:
 *
 *    Top = t(1)
 *    M(Top) = M(1)
 *    if M(2) > 0.85 * M(Top)
 *       M(Top) = M(2)
 *       Top = t(2)
 *    end
 *    if M(3) > 0.85 * M(Top)
 *       M(Top) = M(3)
 *       Top = t(3)
 *    end
 *
 * Returns:
 *    void
 */
static int32_t Pitch_ol(enum Mode mode, vadState * vadSt, float signal[],
		       int32_t pit_min, int32_t pit_max, int16_t L_frame,
		       int32_t dtx, int16_t idx)
{

	float corr[PIT_MAX + 1];
	float max1, max2, max3, p_max1, p_max2, p_max3;
	float *corr_ptr;
	int32_t i, j;

	float corr_hp_max;



       //	if (dtx) {
		/* update tone detection */
		if ((mode == MR475) || (mode == MR515)) {
			vad_tone_detection_update(vadSt, 1);
		} else {
			vad_tone_detection_update(vadSt, 0);
		}
       //	}


	corr_ptr = &corr[pit_max];

	/*        79             */
	/* O(k) = SUM Sw(n)*Sw(n-k)   */
	/*        n=0               */
	comp_corr(signal, L_frame, pit_max, pit_min, corr_ptr);


	/* Find a maximum for each section.  */
	/* Maxima 1  */
	j = pit_min << 2;
	p_max1 =
	    Lag_max(vadSt, corr_ptr, signal, L_frame, pit_max, j, &max1, dtx);

	/* Maxima 2 */
	i = j - 1;
	j = pit_min << 1;
	p_max2 = Lag_max(vadSt, corr_ptr, signal, L_frame, i, j, &max2, dtx);

	/* Maxima 3 */
	i = j - 1;
	p_max3 =
	    Lag_max(vadSt, corr_ptr, signal, L_frame, i, pit_min, &max3, dtx);

       //	if (dtx) {
		if (idx == 1) {
			/* calculate max high-passed filtered correlation of all lags */
			hp_max(corr_ptr, signal, L_frame, pit_max, pit_min,
			       &corr_hp_max);

			/* update complex background detector */
			vadSt->best_corr_hp = corr_hp_max * 0.5F;
		}
       //	}


	/* The best open-loop delay */
	if ((max1 * 0.85F) < max2) {
		max1 = max2;
		p_max1 = p_max2;

	}

	if ((max1 * 0.85F) < max3) {
		p_max1 = p_max3;

	}

	return (int32_t) p_max1;
}




/*
 * ol_ltp
 *
 *
 * Parameters:
 *    mode              I: AMR mode
 *    vadSt             B: VAD state struct
 *    wsp               I: signal used to compute the OL pitch
 *    T_op              O: open loop pitch lag
 *    ol_gain_flg       I: OL gain flag
 *    old_T0_med        O: old Cl lags median
 *    wght_flg          I: weighting function flag
 *    ada_w             B:
 *    old_lags          I: history with old stored Cl lags
 *    ol_gain_flg       I: OL gain flag
 *    dtx               I: DTX flag
 *    idx               I: frame index
 *
 * Function:
 *    Compute the open loop pitch lag.
 *
 *    Open-loop pitch analysis is performed in order to simplify
 *    the pitch analysis and confine the closed-loop pitch search to
 *    a small number of lags around the open-loop estimated lags.
 *    Open-loop pitch estimation is based on the weighted speech signal Sw(n)
 *    which is obtained by filtering the input speech signal through
 *    the weighting filter W(z) = A(z/g1) / A(z/g2). That is,
 *    in a subframe of size L, the weighted speech is given by:
 *
 *                10
 *    Sw(n) = S(n) + SUM[ a(i) * g1(i) * S(n-i) ]
 *                i=1
 *                   10
 *                - SUM[ a(i) * g2(i) * Sw(n-i) ], n = 0, ..., L-1
 *                  i=1
 *
 * Returns:
 *    void
 */
static void ol_ltp(enum Mode mode, vadState * vadSt, float wsp[],
		   int32_t * T_op, float ol_gain_flg[], int32_t * old_T0_med,
		   int16_t * wght_flg, float * ada_w, int32_t * old_lags,
		   int32_t dtx, int16_t idx)
{
	//if (mode != MR102) {
		ol_gain_flg[0] = 0;
		ol_gain_flg[1] = 0;
       //	}

       //	if ((mode == MR475) || (mode == MR515)) {
		*T_op =
		    Pitch_ol(mode, vadSt, wsp, PIT_MIN, PIT_MAX, L_FRAME, dtx,
			     idx);
       //}

       /*
        else {
		if (mode <= MR795) {
			*T_op =
			    Pitch_ol(mode, vadSt, wsp, PIT_MIN, PIT_MAX,
				     L_FRAME_BY2, dtx, idx);
		} else if (mode == MR102) {
			*T_op =
			    Pitch_ol_wgh(old_T0_med, wght_flg, ada_w, vadSt,
					 wsp, old_lags, ol_gain_flg, idx, dtx);
		} else {
			*T_op =
			    Pitch_ol(mode, vadSt, wsp, PIT_MIN_MR122, PIT_MAX,
				     L_FRAME_BY2, dtx, idx);
		}

	}

        */
}

/*
 * subframePreProc
 *
 *
 * Parameters:
 *    mode           I: AMR mode
 *    gamma1         I: spectral exp. factor 1
 *    gamma1_12k2    I: spectral exp. factor 1 for EFR
 *    gamma2         I: spectral exp. factor 2
 *    A              I: A(z) unquantized for the 4 subframes
 *    Aq             I: A(z)   quantized for the 4 subframes
 *    speech         I: speech segment
 *    mem_err        I: pointer to error signal
 *    mem_w0         I: memory of weighting filter
 *    zero           I: pointer to zero vector
 *    ai_zero        O: history of weighted synth. filter
 *    exc            O: long term prediction residual
 *    h1             O: impulse response
 *    xn             O: target vector for pitch search
 *    res2           O: long term prediction residual
 *    error          O: error of LPC synthesis filter
 *
 * Function:
 *    Subframe preprocessing
 *
 *    Impulse response computation:
 *       The impulse response, h(n), of the weighted synthesis filter
 *
 *       H(z) * W(z) = A(z/g1) / ( A'(z) * A(z/g2) )
 *
 *       is computed each subframe. This impulse response is needed for
 *       the search of adaptive and fixed codebooks. The impulse response h(n)
 *       is computed by filtering the vector of coefficients of
 *       the filter A(z/g1) extended by zeros through the two filters
 *       1/A'(z) and 1/A(z/g2).
 *
 *    Target signal computation:
 *       The target signal for adaptive codebook search is usually computed
 *       by subtracting the zero input response of
 *       the weighted synthesis filter H(z) * W(z) from the weighted
 *       speech signal Sw(n). This is performed on a subframe basis.
 *       An equivalent procedure for computing the target signal is
 *       the filtering of the LP residual signal res(n) through
 *       the combination of synthesis filter 1/A'(z) and
 *       the weighting filter A(z/g1)/A(z/g2). After determining
 *       the excitation for the subframe, the initial states of
 *       these filters are updated by filtering the difference between
 *       the LP residual and excitation.
 *
 *       The residual signal res(n) which is needed for finding
 *       the target vector is also used in the adaptive codebook search
 *       to extend the past excitation buffer. This simplifies
 *       the adaptive codebook search procedure for delays less than
 *       the subframe size of 40. The LP residual is given by:
 *
 *                        10
 *       res(n) = S(n) + SUM[A'(i)* S(n-i)
 *                       i=1
 *
 * Returns:
 *    void
 */
static void subframePreProc(enum Mode mode, const float gamma1[],
                            const float gamma2[],
			    float * A, float * Aq, float * speech,
			    float * mem_err, float * mem_w0, float * zero,
			    float ai_zero[], float * exc, float h1[],
			    float xn[], float res2[], float error[])
{
	float Ap1[MP1];	/* weighted LPC coefficients */
	float Ap2[MP1];	/* weighted LPC coefficients */
	const float *g1;

	/* mode specific pointer to gamma1 values */
	g1 = gamma1;

	//if ((mode == MR122) || (mode == MR102)) {
	 //	g1 = gamma1_12k2;
	//}

	/* Find the weighted LPC coefficients for the weighting filter. */
	Weight_Ai(A, g1, Ap1);
	Weight_Ai(A, gamma2, Ap2);

	/*
	 * Compute impulse response, h1[],
	 * of weighted synthesis filter A(z/g1)/A(z/g2)
	 */
	memcpy(ai_zero, Ap1, MP1 << 2);
	Syn_filt(Aq, ai_zero, h1, zero, 0);
	Syn_filt(Ap2, h1, h1, zero, 0);

	/*
	 * Find the target vector for pitch search:
	 */
	/* LP residual */
	Residu(Aq, speech, res2);
	memcpy(exc, res2, L_SUBFR << 2);

	/* Synthesis filter */
	Syn_filt(Aq, exc, error, mem_err, 0);
	Residu(Ap1, error, xn);

	/* target signal xn[] */
	Syn_filt(Ap2, xn, xn, mem_w0, 0);
}

/*
 * getRange
 *
 *
 * Parameters:
 *    T0                I: integer pitch
 *    delta_low         I: search start offset
 *    delta_range       I: search range
 *    pitmin            I: minimum pitch
 *    pitmax            I: maximum pitch
 *    T0_min            I: search range minimum
 *    T0_max            I: search range maximum
 *
 * Function:
 *    Sets range around open-loop pitch or integer pitch of last subframe
 *
 *    Takes integer pitch T0 and calculates a range around it with
 *    T0_min = T0-delta_low and T0_max = (T0-delta_low) + delta_range
 *    T0_min and T0_max are bounded by pitmin and pitmax
 *
 * Returns:
 *    void
 */
static void getRange(int32_t T0, int16_t delta_low, int16_t delta_range,
		     int16_t pitmin, int16_t pitmax, int32_t * T0_min,
		     int32_t * T0_max)
{
	*T0_min = T0 - delta_low;

	if (*T0_min < pitmin) {
		*T0_min = pitmin;
	}
	*T0_max = *T0_min + delta_range;

	if (*T0_max > pitmax) {
		*T0_max = pitmax;
		*T0_min = *T0_max - delta_range;
	}
}

/*
 * Norm_Corr
 *
 *
 * Parameters:
 *    exc         I: excitation buffer                      [L_SUBFR]
 *    xn          I: target vector                          [L_SUBFR]
 *    h           I: impulse response of synthesis and weighting filters
 *                                                          [L_SUBFR]
 *    t_min       I: interval to compute normalized correlation
 *    t_max       I: interval to compute normalized correlation
 *    corr_norm   O: Normalized correlation                 [wT_min-wT_max]
 *
 * Function:
 *    Normalized correlation
 *
 *    The closed-loop pitch search is performed by minimizing
 *    the mean-square weighted error between the original and
 *    synthesized speech. This is achieved by maximizing the term:
 *
 *            39                           39
 *    R(k) = SUM[ X(n) * Yk(n)) ] / SQRT[ SUM[ Yk(n) * Yk(n)] ]
 *           n=0                          n=0
 *
 *    where X(n) is the target signal and Yk(n) is the past filtered
 *    excitation at delay k (past excitation convolved with h(n) ).
 *    The search range is limited around the open-loop pitch.
 *
 *    The convolution Yk(n) is computed for the first delay t_min in
 *    the searched range, and for the other delays in the search range
 *    k = t_min + 1, ..., t_max, it is updated using the recursive relation:
 *
 *    Yk(n) = Yk-1(n-1) + u(-k) * h(n),
 *
 *    where u(n), n = -( 143 + 11 ), ..., 39, is the excitation buffer.
 *    Note that in search stage, the samples u(n), n = 0, ..., 39,
 *    are not known, and they are needed for pitch delays less than 40.
 *    To simplify the search, the LP residual is copied to u(n) in order
 *    to make the relation in above equation valid for all delays.
 *
 * Returns:
 *    void
 */
static void Norm_Corr(float exc[], float xn[], float h[], int32_t t_min,
		      int32_t t_max, float corr_norm[])
{
	float exc_temp[L_SUBFR];
	float *p_exc;
	float corr, norm;
	float sum;
	int32_t i, j, k;

	k = -t_min;
	p_exc = &exc[-t_min];

	/* compute the filtered excitation for the first delay t_min */
	/* convolution Yk(n) */
	for (j = 0; j < L_SUBFR; j++) {
		sum = 0;

		for (i = 0; i <= j; i++) {
			sum += p_exc[i] * h[j - i];
		}
		exc_temp[j] = sum;
	}

	/* loop for every possible period */
	for (i = t_min; i <= t_max; i++) {
		/*        39                     */
		/* SQRT[ SUM[ Yk(n) * Yk(n)] ]   */
		/*       n=0                     */
		norm = (float) Dotproduct40(exc_temp, exc_temp);

		if (norm == 0)
			norm = 1.0;
		else
			norm = (float) (1.0F / (sqrtf(norm)));

		/*        39                  */
		/* SQRT[ SUM[ X(n) * Yk(n)] ] */
		/*       n=0                  */
		corr = (float) Dotproduct40(xn, exc_temp);

		/* R(k) */
		corr_norm[i] = corr * norm;

		/* modify the filtered excitation exc_tmp[] for the next iteration */
		if (i != t_max) {
			k--;

			for (j = L_SUBFR - 1; j > 0; j--) {
				/* Yk(n) = Yk-1(n-1) + u(-k) * h(n) */
				exc_temp[j] = exc_temp[j - 1] + exc[k] * h[j];
			}
			exc_temp[0] = exc[k];
		}
	}
}

/*
 * Interpol_3or6
 *
 *
 * Parameters:
 *    x                 I: input vector
 *    frac              I: fraction  (-2..2 for 3*, -3..3 for 6*)
 *    flag3             I: if set, upsampling rate = 3 (6 otherwise)
 *
 * Function:
 *    Interpolating the normalized correlation with 1/3 or 1/6 resolution.
 *
 *    The interpolation is performed using an FIR filter b24
 *    based on a Hamming windowed sin(x)/x function truncated at ±23
 *    and padded with zeros at ±24 (b24(24) = 0). The filter has its
 *    cut-off frequency (-3 dB) at 3 600 Hz in the over-sampled domain.
 *    The interpolated values of R(k) for the fractions -3/6 to 3/6
 *    are obtained using the interpolation formula:
 *
 *              3                            3
 *    R(k)t = SUM[ R(k-i) * b24(t+i*6) ] + SUM [ R(k+1+i) * b24(6-t+i*6) ],
 *            i=0                          i=0
 *    t = 0, ..., 5,
 *
 *    where t = 0, ..., 5, corresponds to the fractions
 *    0, 1/6, 2/6, 3/6, -2/6, and -1/6, respectively. Note that it is
 *    necessary to compute the correlation terms using a range t_min - 4,
 *    t_max + 4, to allow for the proper interpolation.
 *
 * Returns:
 *    s                 interpolated value
 */
static float Interpol_3or6(float * x, int32_t frac, int16_t flag3)
{
	float s;
	float *x1, *x2;
	const float *c1, *c2;
	int32_t i, k;

	if (flag3 != 0) {
		/* inter_3[k] = b60[2*k] -> k' = 2*k */
		frac <<= 1;
	}

	if (frac < 0) {
		frac += UP_SAMP_MAX;
		x--;
	}
	x1 = &x[0];
	x2 = &x[1];
	c1 = &f_b24[frac];
	c2 = &f_b24[UP_SAMP_MAX - frac];
	s = 0;

	for (i = 0, k = 0; i < L_INTER_SRCH; i++, k += UP_SAMP_MAX) {
		/* R(k-i) * b24(t+i*6) */
		s += x1[-i] * c1[k];

		/* R(k+1+i) * b24(6-t+i*6) */
		s += x2[i] * c2[k];
	}
	return s;
}

/*
 * searchFrac
 *
 *
 * Parameters:
 *    lag               B: integer pitch
 *    frac              B: start point of search - fractional pitch
 *    last_frac         I: endpoint of search
 *    corr              I: normalized correlation
 *    flag3             I: if set, upsampling rate = 3 (6 otherwise)
 *
 * Function:
 *    Find fractional pitch
 *
 *    The function interpolates the normalized correlation at the
 *    fractional positions around lag T0. The position at which the
 *    interpolation function reaches its maximum is the fractional pitch.
 *    Starting point of the search is frac, end point is last_frac.
 *    frac is overwritten with the fractional pitch.
 *
 * Returns:
 *    void
 */
static void searchFrac(int32_t * lag, int32_t * frac, int16_t last_frac, float
		       corr[], int16_t flag3)
{
	float max, corr_int;
	int32_t i;

	/*
	 * Test the fractions around T0 and choose the one which maximizes
	 * the interpolated normalized correlation.
	 */
	max = Interpol_3or6(&corr[*lag], *frac, flag3);

	for (i = *frac + 1; i <= last_frac; i++) {
		corr_int = Interpol_3or6(&corr[*lag], i, flag3);

		if (corr_int > max) {
			max = corr_int;
			*frac = i;
		}
	}

	if (flag3 == 0) {
		/* Limit the fraction value in the interval [-2,-1,0,1,2,3] */
		if (*frac == -3) {
			*frac = 3;
			*lag -= 1;
		}
	} else {
		/* limit the fraction value between -1 and 1 */
		if (*frac == -2) {
			*frac = 1;
			*lag -= 1;
		}

		if (*frac == 2) {
			*frac = -1;
			*lag += 1;
		}
	}
}

/*
 * Enc_lag3
 *
 *
 * Parameters:
 *    T0             I: Pitch delay
 *    T0_frac        I: Fractional pitch delay
 *    T0_prev        I: Integer pitch delay of last subframe
 *    T0_min         I: minimum of search range
 *    T0_max         I: maximum of search range
 *    delta_flag     I: Flag for 1st (or 3rd) subframe
 *    flag4          I: Flag for encoding with 4 bits
 *
 * Function:
 *    Encoding of fractional pitch lag with 1/3 resolution.
 *
 * Returns:
 *    index             index of encoding
 */
static int32_t Enc_lag3(int32_t T0, int32_t T0_frac, int32_t T0_prev, int32_t T0_min,
		       int32_t T0_max, int16_t delta_flag, int16_t flag4)
{
	int32_t index, i, tmp_ind, uplag, tmp_lag;

	/* if 1st or 3rd subframe */
	if (delta_flag == 0) {
		/* encode pitch delay (with fraction) */
		if (T0 <= 85) {
			index = T0 * 3 - 58 + T0_frac;
		} else {
			index = T0 + 112;
		}
	}

	/* if second or fourth subframe */
	else {
		if (flag4 == 0) {
			/* 'normal' encoding: either with 5 or 6 bit resolution */
			index = 3 * (T0 - T0_min) + 2 + T0_frac;
		} else {
			/* encoding with 4 bit resolution */
			tmp_lag = T0_prev;

			if ((tmp_lag - T0_min) > 5)
				tmp_lag = T0_min + 5;

			if ((T0_max - tmp_lag) > 4)
				tmp_lag = T0_max - 4;
			uplag = T0 + T0 + T0 + T0_frac;
			i = tmp_lag - 2;
			tmp_ind = i + i + i;

			if (tmp_ind >= uplag) {
				index = (T0 - tmp_lag) + 5;
			} else {
				i = tmp_lag + 1;
				i = i + i + i;

				if (i > uplag) {
					index = (uplag - tmp_ind) + 3;
				} else {
					index = (T0 - tmp_lag) + 11;
				}
			}
		}		/* end if (encoding with 4 bit resolution) */
	}			/* end if (second of fourth subframe) */
	return index;
}

/*
 * Pitch_fr
 *
 *
 * Parameters:
 *    T0_prev_subframe  B: integer pitch lag of previous sub-frame
 *    mode              I: codec mode
 *    T_op              I: open-loop pitch estimations for
 *                         the 2 big subframes [2]
 *    exc               I: excitation buffer
 *    xn                I: target vector
 *    h                 I: impulse response of synthesis
 *                         and weighting filters
 *    i_subfr           I: subframe number
 *    pit_frac          O: pitch period (fractional)
 *    resu3             O: subsample resolution 1/3 (=1) or 1/6 (=0)
 *    ana_index         O: index of encoding
 *
 * Function:
 *    Closed-loop pitch search
 *
 *    In the first and third subframes, a fractional pitch delay is used
 *    with resolutions: 1/6 in the range [17 3/6, 94 3/6] and integers only
 *    in the range [95, 143]. For the second and fourth subframes,
 *    a pitch resolution of 1/6 is always used in
 *    the range [T1 - 5 3/6, T1 + 4 /3/6], where T1 is nearest integer to
 *    the fractional pitch lag of the previous (1st or 3rd) subframe,
 *    bounded by 18...143.
 *
 *    Closed-loop pitch analysis is performed around
 *    the open-loop pitch estimates on a subframe basis.
 *    In the first (and third) subframe the range Top±3,
 *    bounded by 18...143, is searched. For the other subframes,
 *    closed-loop pitch analysis is performed around the integer pitch
 *    selected in the previous subframe, as described above.
 *    The pitch delay is encoded with 9 bits in the first and
 *    third subframes and the relative delay of the other subframes
 *    is encoded with 6 bits.
 *
 *    The closed-loop pitch search is performed by minimizing
 *    the mean-square weighted error between the original and
 *    synthesized speech. This is achieved by maximizing the term:
 *
 *            39                           39
 *    R(k) = SUM[ X(n) * Yk(n)) ] / SQRT[ SUM[ Yk(n) * Yk(n)] ]
 *           n=0                          n=0
 *
 *    where X(n) is the target signal and Yk(n) is the past filtered
 *    excitation at delay k (past excitation convolved with h(n) ).
 *
 *    Once the optimum integer pitch delay is determined, the fractions
 *    from -3/6 to 3/6 with a step of 1/6 around that integer are tested.
 *    The fractional pitch search is performed by interpolating
 *    the normalized correlation R(k) and searching for its maximum.
 *    The interpolation is performed using an FIR filter b24
 *    based on a Hamming windowed sin(x)/x function truncated at ±23
 *    and padded with zeros at ±24 (b24(24) = 0). The filter has its
 *    cut-off frequency (-3 dB) at 3 600 Hz in the over-sampled domain.
 *    The interpolated values of R(k) for the fractions -3/6 to 3/6
 *    are obtained using the interpolation formula:
 *
 *              3                            3
 *    R(k)t = SUM[ R(k-i) * b24(t+i*6) ] + SUM [ R(k+1+i) * b24(6-t+i*6) ],
 *            i=0                          i=0
 *    t = 0, ..., 5,
 *
 *    where t = 0, ..., 5, corresponds to the fractions
 *    0, 1/6, 2/6, 3/6, -2/6, and -1/6, respectively. Note that it is
 *    necessary to compute the correlation terms using a range t_min -4,
 *    t_max + 4, to allow for the proper interpolation.
 *
 * Returns:
 *    lag             closed-loop pitch lag
 */
static int32_t Pitch_fr(int32_t * T0_prev_subframe, enum Mode mode, int32_t T_op[],
		       float exc[], float xn[], float h[], int16_t i_subfr,
		       int32_t * pit_frac, int16_t * resu3, int32_t * ana_index)
{
	float corr_v[40];
	float max;
	float *corr;
	int32_t i, t_min, t_max, T0_min, T0_max;
	int32_t lag, frac, tmp_lag;
	int16_t max_frac_lag, flag3, flag4, last_frac;
	int16_t delta_int_low, delta_int_range, delta_frc_low, delta_frc_range;
	int16_t pit_min;
	int16_t frame_offset;
	int16_t delta_search;

	/* set mode specific variables */
	max_frac_lag = mode_dep_parm[mode].max_frac_lag;
	flag3 = mode_dep_parm[mode].flag3;
	frac = mode_dep_parm[mode].first_frac;
	last_frac = mode_dep_parm[mode].last_frac;
	delta_int_low = mode_dep_parm[mode].delta_int_low;
	delta_int_range = mode_dep_parm[mode].delta_int_range;
	delta_frc_low = mode_dep_parm[mode].delta_frc_low;
	delta_frc_range = mode_dep_parm[mode].delta_frc_range;
	pit_min = mode_dep_parm[mode].pit_min;

	/* decide upon full or differential search */
	delta_search = 1;

	if ((i_subfr == 0) || (i_subfr == L_FRAME_BY2)) {
		/* Subframe 1 and 3 */
		if (((mode != MR475) && (mode != MR515)) || (i_subfr !=
							     L_FRAME_BY2)) {
			/*
			 * set T0_min, T0_max for full search
			 * this is *not* done for mode MR475, MR515 in subframe 3
			 */
			delta_search = 0;	/* no differential search */

			/*
			 * calculate index into T_op which contains the open-loop
			 * pitch estimations for the 2 big subframes
			 */
			frame_offset = 1;

			if (i_subfr == 0)
				frame_offset = 0;

			/*
			 * get T_op from the corresponding half frame and
			 * set T0_min, T0_max
			 */
			getRange(T_op[frame_offset], delta_int_low,
				 delta_int_range, pit_min, PIT_MAX, &T0_min,
				 &T0_max);
		} else {
			/* mode MR475, MR515 and 3. Subframe: delta search as well */
			getRange(*T0_prev_subframe, delta_frc_low,
				 delta_frc_range, pit_min, PIT_MAX, &T0_min,
				 &T0_max);
		}
	} else {
		/*
		 * for Subframe 2 and 4
		 * get range around T0 of previous subframe for delta search
		 */
		getRange(*T0_prev_subframe, delta_frc_low, delta_frc_range,
			 pit_min, PIT_MAX, &T0_min, &T0_max);
	}

	/* Find interval to compute normalized correlation */
	t_min = T0_min - L_INTER_SRCH;
	t_max = T0_max + L_INTER_SRCH;
	corr = &corr_v[-t_min];

	/* Compute normalized correlation between target and filtered excitation */
	Norm_Corr(exc, xn, h, t_min, t_max, corr);

	/* Find integer pitch */
	max = corr[T0_min];
	lag = T0_min;

	for (i = T0_min + 1; i <= T0_max; i++) {
		if (corr[i] >= max) {
			max = corr[i];
			lag = i;
		}
	}

	/* Find fractional pitch   */
	if ((delta_search == 0) && (lag > max_frac_lag)) {
		/*
		 * full search and integer pitch greater than max_frac_lag
		 * fractional search is not needed, set fractional to zero
		 */
		frac = 0;
	} else {
		/*
		 * if differential search AND mode MR475 OR MR515 OR MR59 OR MR67
		 * then search fractional with 4 bits resolution
		 */
		if ((delta_search != 0)
		    && ((mode == MR475) || (mode == MR515) || (mode == MR59)
			|| (mode == MR67))) {
			/*
			 * modify frac or last_frac according to position of last
			 * integer pitch: either search around integer pitch,
			 * or only on left or right side
			 */
			tmp_lag = *T0_prev_subframe;

			if ((tmp_lag - T0_min) > 5)
				tmp_lag = T0_min + 5;

			if ((T0_max - tmp_lag) > 4)
				tmp_lag = T0_max - 4;

			if ((lag == tmp_lag) || (lag == (tmp_lag - 1))) {
				/* normal search in fractions around T0 */
				searchFrac(&lag, &frac, last_frac, corr, flag3);
			} else if (lag == (tmp_lag - 2)) {
				/* limit search around T0 to the right side */
				frac = 0;
				searchFrac(&lag, &frac, last_frac, corr, flag3);
			} else if (lag == (tmp_lag + 1)) {
				/* limit search around T0 to the left side */
				last_frac = 0;
				searchFrac(&lag, &frac, last_frac, corr, flag3);
			} else {
				/* no fractional search */
				frac = 0;
			}
		} else
			/* test the fractions around T0 */
			searchFrac(&lag, &frac, last_frac, corr, flag3);
	}

	/*
	 *  encode pitch
	 */
	if (flag3 != 0) {
		/*
		 * flag4 indicates encoding with 4 bit resolution;
		 * this is needed for mode MR475, MR515 and MR59
		 */
		flag4 = 0;

		//if ((mode == MR475) || (mode == MR515) || (mode == MR59)
		 //   || (mode == MR67)) {
			flag4 = 1;
	       //	}

		/* encode with 1/3 subsample resolution */
		*ana_index =
		    Enc_lag3(lag, frac, *T0_prev_subframe, T0_min, T0_max,
			     delta_search, flag4);
	}

        //else {
		/* encode with 1/6 subsample resolution */
	 //	*ana_index = Enc_lag6(lag, frac, T0_min, delta_search);
	//}

	/*
	 *  update state variables
	 */
	*T0_prev_subframe = lag;

	/*
	 * update output variables
	 */
	*resu3 = flag3;
	*pit_frac = frac;
	return (lag);
}

/*
 * Pred_lt_3or6
 *
 *
 * Parameters:
 *    exc      B: excitation buffer
 *    T0       I: integer pitch lag
 *    frac     I: fraction of lag
 *    flag3    I: if set, upsampling rate = 3 (6 otherwise)
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
 *            9                              9
 *    v(n) = SUM[ u(n-k-i) * b60(t+i*6) ] + SUM[ u(n-k+1+i) * b60(6-t+i*6) ],
 *           i=0                            i=0
 *    n = 0, ...,39, t = 0, ...,5.
 *
 *    The interpolation filter b60 is based on a Hamming windowed sin(x)/x
 *    function truncated at ± 59 and padded with zeros at ± 60 (b60(60)=0)).
 *    The filter has a cut-off frequency (-3 dB) at 3 600 Hz in
 *    the over-sampled domain.
 *
 * Returns:
 *    void
 */
static void Pred_lt_3or6(float exc[], int32_t T0, int32_t frac, int16_t flag3)
{
	float s;
	float *x0, *x1, *x2;
	const float *c1, *c2;
	int32_t j;

	x0 = &exc[-T0];
	frac = -frac;

	if (flag3 != 0) {
		/* inter_3l[k] = b60[2*k] -> k' = 2*k */
		frac <<= 1;
	}

	if (frac < 0) {
		frac += UP_SAMP_MAX;
		x0--;
	}
	c1 = &f_b60[frac];
	c2 = &f_b60[UP_SAMP_MAX - frac];

	for (j = 0; j < L_SUBFR; j += 4) {
		x1 = x0++;
		x2 = x0;
		s = x1[0] * c1[0] + x2[0] * c2[0];
		s += x1[-1] * c1[6] + x2[1] * c2[6];
		s += x1[-2] * c1[12] + x2[2] * c2[12];
		s += x1[-3] * c1[18] + x2[3] * c2[18];
		s += x1[-4] * c1[24] + x2[4] * c2[24];
		s += x1[-5] * c1[30] + x2[5] * c2[30];
		s += x1[-6] * c1[36] + x2[6] * c2[36];
		s += x1[-7] * c1[42] + x2[7] * c2[42];
		s += x1[-8] * c1[48] + x2[8] * c2[48];
		s += x1[-9] * c1[54] + x2[9] * c2[54];
		exc[j] = (float) floorf(s + 0.5F);
		x1 = x0++;
		x2 = x0;
		s = x1[0] * c1[0] + x2[0] * c2[0];
		s += x1[-1] * c1[6] + x2[1] * c2[6];
		s += x1[-2] * c1[12] + x2[2] * c2[12];
		s += x1[-3] * c1[18] + x2[3] * c2[18];
		s += x1[-4] * c1[24] + x2[4] * c2[24];
		s += x1[-5] * c1[30] + x2[5] * c2[30];
		s += x1[-6] * c1[36] + x2[6] * c2[36];
		s += x1[-7] * c1[42] + x2[7] * c2[42];
		s += x1[-8] * c1[48] + x2[8] * c2[48];
		s += x1[-9] * c1[54] + x2[9] * c2[54];
		exc[j + 1] = (float) floorf(s + 0.5F);
		x1 = x0++;
		x2 = x0;
		s = x1[0] * c1[0] + x2[0] * c2[0];
		s += x1[-1] * c1[6] + x2[1] * c2[6];
		s += x1[-2] * c1[12] + x2[2] * c2[12];
		s += x1[-3] * c1[18] + x2[3] * c2[18];
		s += x1[-4] * c1[24] + x2[4] * c2[24];
		s += x1[-5] * c1[30] + x2[5] * c2[30];
		s += x1[-6] * c1[36] + x2[6] * c2[36];
		s += x1[-7] * c1[42] + x2[7] * c2[42];
		s += x1[-8] * c1[48] + x2[8] * c2[48];
		s += x1[-9] * c1[54] + x2[9] * c2[54];
		exc[j + 2] = (float) floorf(s + 0.5F);
		x1 = x0++;
		x2 = x0;
		s = x1[0] * c1[0] + x2[0] * c2[0];
		s += x1[-1] * c1[6] + x2[1] * c2[6];
		s += x1[-2] * c1[12] + x2[2] * c2[12];
		s += x1[-3] * c1[18] + x2[3] * c2[18];
		s += x1[-4] * c1[24] + x2[4] * c2[24];
		s += x1[-5] * c1[30] + x2[5] * c2[30];
		s += x1[-6] * c1[36] + x2[6] * c2[36];
		s += x1[-7] * c1[42] + x2[7] * c2[42];
		s += x1[-8] * c1[48] + x2[8] * c2[48];
		s += x1[-9] * c1[54] + x2[9] * c2[54];
		exc[j + 3] = (float) floorf(s + 0.5F);
	}
	return;
}

static void Pred_lt_3or6_fixed(int32_t exc[], int32_t T0, int32_t frac,
			       int32_t flag3)
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
	c1 = &f_inter6[frac];
	c2 = &f_inter6[6 - frac];

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
 * G_pitch
 *
 *
 * Parameters:
 *    xn       I: Pitch target
 *    y1       I: Filtered adaptive codebook
 *    gCoeff   O: Correlations need for gain quantization
 *
 * Function:
 *    Compute the pitch (adaptive codebook) gain.
 *
 *    The adaptive codebook gain is given by
 *
 *       g = <x[], y[]> / <y[], y[]>
 *
 *    where x[] is the target vector, y[] is the filtered adaptive
 *    codevector, and <> denotes dot product.
 *
 * Returns:
 *    gain              gain saturated to 1.2
 */
static float G_pitch(float xn[], float y1[], float gCoeff[])
{
	float gain, sum;

	/* Compute scalar product <y1[],y1[]> */
	sum = (float) Dotproduct40(y1, y1);

	/* Avoid case of all zeros */
	sum += 0.01F;
	gCoeff[0] = sum;

	/* Compute scalar product <xn[],y1[]> */
	sum = (float) Dotproduct40(xn, y1);
	gCoeff[1] = sum;

	/* compute gain = xy/yy */
	gain = (float) (gCoeff[1] / gCoeff[0]);

	/* if(gain >1.2) gain = 1.2 */
	if (gain < 0.0F)
		gain = 0.0F;

	if (gain > 1.2F)
		gain = 1.2F;
	return (gain);
}

/*
 * check_gp_clipping
 *
 *
 * Parameters:
 *    gp          I: old pitch gains
 *    g_pitch     I: pitch gain
 *
 * Function:
 *    Verify that the sum of the last (N_FRAME+1) pitch gains is under
 *    a certain threshold.
 *
 * Returns:
 *    True or false
 */
static int16_t check_gp_clipping(float * gp, float g_pitch)
{
	float sum;
	int32_t i;

	sum = g_pitch;

	for (i = 0; i < N_FRAME; i++) {
		sum += gp[i];
	}

	if (sum > 7.6F /*FGP_CLIP */ ) {
		return 1;
	} else {
		return 0;
	}
}

/*
 * cl_ltp
 *
 *
 * Parameters:
 *    T0_prev_subframe  B: Integer pitch lag of previous sub-frame
 *    gp                I: Gain history
 *    mode              I: Coder mode
 *    frame_offset      I: Offset to subframe
 *    T_op              I: Open loop pitch lags
 *    h1                I: Impulse response vector
 *    exc               B: Excitation vector
 *    res2              B: Long term prediction residual
 *    xn                I: Target vector for pitch search
 *    lsp_flag          I: LSP resonance flag
 *    xn2               O: Target vector for codebook search
 *    y1                O: Filtered adaptive excitation
 *    T0                O: Pitch delay (integer part)
 *    T0_frac           O: Pitch delay (fractional part)
 *    gain_pit          O: Pitch gain
 *    gCoeff[]          O: Correlations between xn, y1, & y2
 *    anap              O: Analysis parameters
 *    gp_limit          O: pitch gain limit
 *
 * Function:
 *    Closed-loop ltp search
 *
 *    Adaptive codebook search is performed on a subframe basis.
 *    It consists of performing closed-loop pitch search, and then computing
 *    the adaptive codevector by interpolating the past excitation at
 *    the selected fractional pitch lag.
 *    The adaptive codebook parameters (or pitch parameters) are
 *    the delay and gain of the pitch filter. In the adaptive codebook approach
 *    for implementing the pitch filter, the excitation is repeated for delays
 *    less than the subframe length. In the search stage, the excitation is
 *    extended by the LP residual to simplify the closed-loop search.
 *
 * Returns:
 *    void
 */
static void cl_ltp(int32_t * T0_prev_subframe, float * gp, enum Mode mode,
		   int16_t frame_offset, int32_t T_op[], float * h1,
		   float * exc, float res2[], float xn[], int16_t lsp_flag,
		   float xn2[], float y1[], int32_t * T0, int32_t * T0_frac,
		   float * gain_pit, float gCoeff[], int16_t ** anap,
		   float * gp_limit)
{
	float s;
	int32_t i, n;
	int16_t gpc_flag, resu3;	/* flag for upsample resolution */

	int32_t exc_tmp[314];
	int32_t *exc_tmp_p;

	exc_tmp_p = exc_tmp + PIT_MAX + L_INTERPOL;

	/* Closed-loop fractional pitch search */
	*T0 = Pitch_fr(T0_prev_subframe, mode, T_op, exc, xn, h1, frame_offset,
		       T0_frac, &resu3, &i);
	*(*anap)++ = (int16_t) i;

	/*
	 * Compute the adaptive codebook vector
	 * using fixed point. This is required
	 * to maintain encoder/decoder excitation
	 * syncronisation
	 */
	for (i = -(PIT_MAX + L_INTERPOL); i < 40; i++)
		exc_tmp_p[i] = (int32_t) exc[i];

	Pred_lt_3or6_fixed(exc_tmp_p, *T0, *T0_frac, resu3);

	for (i = -(PIT_MAX + L_INTERPOL); i < 40; i++)
		exc[i] = (float) exc_tmp_p[i];

	/*
	 *   Convolve to get filtered adaptive codebook vector
	 *  y[n] = sum_{i=0}^{n} x[i] h[n-i], n=0,...,L-1
	 */
	for (n = 0; n < L_SUBFR; n++) {
		s = 0;

		for (i = 0; i <= n; i++) {
			s += exc[i] * h1[n - i];
		}
		y1[n] = s;
	}

	/* The adaptive codebook gain */
	*gain_pit = G_pitch(xn, y1, gCoeff);

	/* check if the pitch gain should be limit due to resonance in LPC filter */
	gpc_flag = 0;
	*gp_limit = 2.0F;

	if ((lsp_flag != 0) && (*gain_pit > 0.95F)) {
		gpc_flag = check_gp_clipping(gp, *gain_pit);
	}

	/*
	 * special for the MR475, MR515 mode; limit the gain to 0.85 to
	 * cope with bit errors in the decoder in a better way.
	 */
	//if ((mode == MR475) || (mode == MR515)) {


		if (*gain_pit > 0.85F) {
			*gain_pit = 0.85F;
		}

		if (gpc_flag != 0)
			*gp_limit = GP_CLIP;
       //	}

        /*
        else {
		if (gpc_flag != 0) {
			*gp_limit = GP_CLIP;
			*gain_pit = GP_CLIP;
		}




		//
		 // 12k2 gain_pit is quantized here and not in gainQuant.
		 //
		if (mode == MR122) {

			*(*anap)++ =
			    q_gain_pitch(MR122, *gp_limit, gain_pit, NULL,
					 NULL);
		}
	}

        */

	/*
	 * Update target vector for codebook search
	 * Find LTP residual
	 */
	for (i = 0; i < L_SUBFR; i++) {
		xn2[i] = xn[i] - y1[i] * *gain_pit;
		res2[i] = res2[i] - exc[i] * *gain_pit;
	}
}

/*
 * DotProduct
 *
 *
 * Parameters:
 *    x                 I: first input
 *    y                 I: second input
 *    len               I: length of product
 *
 * Function:
 *    Computes dot product
 *
 * Returns:
 *    acc               dot product
 */
static float DotProduct(float * x, float * y, int32_t len)
{
	int32_t i;
	float acc;

	acc = 0.0F;

	for (i = 0; i < len; i++)
		acc += x[i] * y[i];
	return (acc);
}

/*
 * cor_h_x
 *
 *
 * Parameters:
 *    h                 I: impulse response of weighted synthesis filter
 *    x                 I: target
 *    dn                O: correlation between target and impulse response
 *
 * Function:
 *    Computes correlation between target signal and impulse response.
 *
 * Returns:
 *    void
 */
static void cor_h_x(float h[], float x[], float dn[])
{
	int32_t i;

	dn[0] = (float) Dotproduct40(h, x);

	for (i = 1; i < L_CODE; i++)
		dn[i] = (float) DotProduct(h, &x[i], L_CODE - i);
}

/*
 * set_sign
 *
 *
 * Parameters:
 *    dn                B: correlation between target and h[]
 *    sign              O: sign of dn[]
 *    dn2               O: maximum of correlation in each track
 *    n                 I: # of maximum correlations in dn2[]
 *
 * Function:
 *    Builds sign[] vector.
 *
 * Returns:
 *    void
 */
static void set_sign(float dn[], float sign[], float dn2[], int16_t n)
{
	float val, min;
	int32_t i, j, k, pos = 0;

	/* set sign according to dn[] */
	for (i = 0; i < L_CODE; i++) {
		val = dn[i];

		if (val >= 0) {
			sign[i] = 1.0F;
		} else {
			sign[i] = -1.0F;
			val = -val;
		}

		/* modify dn[] according to the fixed sign */
		dn[i] = val;
		dn2[i] = val;
	}

	/* keep 8-n maximum positions/8 of each track and store it in dn2[] */
	for (i = 0; i < NB_TRACK; i++) {
		for (k = 0; k < (8 - n); k++) {
			min = FLT_MAX;

			for (j = i; j < L_CODE; j += STEP) {
				if (dn2[j] >= 0) {
					val = dn2[j] - min;

					if (val < 0) {
						min = dn2[j];
						pos = j;
					}
				}
			}
			dn2[pos] = -1.0F;
		}
	}
	return;
}

/*
 * cor_h
 *
 *
 * Parameters:
 *    h                I: h[]
 *    sign             I: sign information
 *    rr               O: correlations
 *
 * Function:
 *    Computes correlations of h[] needed for the codebook search,
 *    and includes the sign information into the correlations.
 *
 * Returns:
 *    void
 */
static void cor_h(float h[], float sign[], float rr[][L_CODE])
{
	float sum;
	float *prr, *ph, *ph_max;
	float *rrj, *rri, *signi, *signj;
	int32_t ii, total_loops, four_loops;

	sum = 0.0F;

	/* Compute diagonal matrix of autocorrelation of h */
	rr[0][0] = (float) Dotproduct40(h, h);
	prr = &rr[39][39];
	ph = &h[0];
	ph_max = ph + 39;

	/*
	 * speed optimization of code:
	 * for (k=0; k<m; k++)
	 * {
	 * sum += h[k]*h[k];
	 * rr[i][i] = sum;
	 * i--;
	 * }
	 */
	do {
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
		sum += *ph * *ph;
		ph++;
		*prr = sum;
		prr -= 41;
	} while (ph < ph_max);

	/*
	 * Compute upper & bottom symmetric matrix of autocorrelation of h.
	 */
	/* speed optimization of code:
	 * for (ii=1; ii<L_CODE; ii++)
	 * {
	 * j = m;
	 * i = j - ii;
	 * sum = (float)0.0;
	 * for ( k = 0; k < (L_CODE-ii); k++ )
	 * {
	 * sum += h[k]*h[k+ii];
	 * rr[i][j] = rr[j][i] = (sum*sign[i]*sign[j]);
	 * i--; j--;
	 * }
	 * }
	 */
	ii = 1;

	for (total_loops = 9; total_loops >= 0; total_loops--) {
		rrj = rri = &rr[39][39];
		rrj -= ii;
		rri = (rri - 40 * ii);
		signi = signj = &sign[39];
		signi -= ii;
		sum = 0.0F;
		ph = &h[0];

		for (four_loops = 0; four_loops < total_loops; four_loops++) {
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
		}
		sum += *ph * *(ph + ii);
		ph++;
		*rri = *rrj = sum * *signi * *signj;
		rri -= 41;
		rrj -= 41;
		signi--;
		signj--;
		sum += *ph * *(ph + ii);
		ph++;
		*rri = *rrj = sum * *signi * *signj;
		rri -= 41;
		rrj -= 41;
		signi--;
		signj--;
		sum += *ph * *(ph + ii);
		*rri = *rrj = sum * *signi * *signj;
		ii++;
		rrj = rri = &rr[39][39];
		rrj -= ii;
		rri = (rri - 40 * ii);
		signi = signj = &sign[39];
		signi -= ii;
		sum = 0.0F;
		ph = &h[0];

		for (four_loops = 0; four_loops < total_loops; four_loops++) {
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
		}
		sum += *ph * *(ph + ii);
		ph++;
		*rri = *rrj = sum * *signi * *signj;
		rri -= 41;
		rrj -= 41;
		signi--;
		signj--;
		sum += *ph * *(ph + ii);
		*rri = *rrj = sum * *signi * *signj;
		ii++;
		rrj = rri = &rr[39][39];
		rrj -= ii;
		rri = (rri - 40 * ii);
		signi = signj = &sign[39];
		signi -= ii;
		sum = 0.0F;
		ph = &h[0];

		for (four_loops = 0; four_loops < total_loops; four_loops++) {
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
		}
		sum += *ph * *(ph + ii);
		*rri = *rrj = sum * *signi * *signj;
		ii++;
		rrj = rri = &rr[39][39];
		rrj -= ii;
		rri = (rri - 40 * ii);
		signi = signj = &sign[39];
		signi -= ii;
		sum = 0.0F;
		ph = &h[0];

		for (four_loops = 0; four_loops < total_loops; four_loops++) {
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * (*signi) * (*signj);
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
			sum += *ph * *(ph + ii);
			ph++;
			*rri = *rrj = sum * *signi * *signj;
			rri -= 41;
			rrj -= 41;
			signi--;
			signj--;
		}
		ii++;
	}
	return;
}

/*
 * search_2i40_9bits
 *
 *
 * Parameters:
 *    subNr             I: subframe number
 *    dn                I: correlation between target and h[]
 *    rr                I: matrix of autocorrelation
 *    codvec            O: algebraic codebook vector
 *
 * Function:
 *    Search the best codevector; determine positions of the 2 pulses
 *    in the 40-sample frame.
 *
 *    First subframe:
 *    first    i0 :  0, 5, 10, 15, 20, 25, 30, 35.
 *             i1 :  2, 7, 12, 17, 22, 27, 32, 37.
 *    second   i0 :  1, 6, 11, 16, 21, 26, 31, 36.
 *             i1 :  3, 8, 13, 18, 23, 28, 33, 38.
 *
 *    Second subframe:
 *    first    i0 :  0, 5, 10, 15, 20, 25, 30, 35.
 *             i1 :  3, 8, 13, 18, 23, 28, 33, 38.
 *    second   i0 :  2, 7, 12, 17, 22, 27, 32, 37.
 *             i1 :  4, 9, 14, 19, 24, 29, 34, 39.
 *
 *    Third subframe:
 *    first    i0 :  0, 5, 10, 15, 20, 25, 30, 35.
 *             i1 :  2, 7, 12, 17, 22, 27, 32, 37.
 *    second   i0 :  1, 6, 11, 16, 21, 26, 31, 36.
 *             i1 :  4, 9, 14, 19, 24, 29, 34, 39.
 *
 *    Fourth subframe:
 *    first    i0 :  0, 5, 10, 15, 20, 25, 30, 35.
 *             i1 :  3, 8, 13, 18, 23, 28, 33, 38.
 *    second   i0 :  1, 6, 11, 16, 21, 26, 31, 36.
 *             i1 :  4, 9, 14, 19, 24, 29, 34, 39.
 *
 * Returns:
 *    void
 */
static void search_2i40_9bits(int16_t subNr, float dn[], float rr[][L_CODE],
			      int32_t codvec[])
{
	float ps0, ps1, psk, alp, alp0, alp1, alpk, sq, sq1;
	int32_t i0, i1, ix, i;
	int16_t ipos[2];
	int16_t track1;

	psk = -1;
	alpk = 1;

	for (i = 0; i < 2; i++) {
		codvec[i] = i;
	}

	/* main loop: try 2x4  tracks        */
	for (track1 = 0; track1 < 2; track1++) {
		ipos[0] = f_startPos[(subNr << 1) + (track1 << 3)];
		ipos[1] = f_startPos[(subNr << 1) + 1 + (track1 << 3)];

		/* i0 loop: try 8 positions       */
		for (i0 = ipos[0]; i0 < L_CODE; i0 += STEP) {
			ps0 = dn[i0];
			alp0 = rr[i0][i0];

			/* i1 loop: 8 positions */
			sq = -1;
			alp = 1;
			ix = ipos[1];

			for (i1 = ipos[1]; i1 < L_CODE; i1 += STEP) {
				ps1 = ps0 + dn[i1];
				alp1 = alp0 + rr[i1][i1] + 2.0F * rr[i0][i1];
				sq1 = ps1 * ps1;

				if ((alp * sq1) > (sq * alp1)) {
					sq = sq1;
					alp = alp1;
					ix = i1;
				}
			}

			/* memorise codevector if this one is better than the last one */
			if ((alpk * sq) > (psk * alp)) {
				psk = sq;
				alpk = alp;
				codvec[0] = i0;
				codvec[1] = ix;
			}
		}
	}
	return;
}

/*
 * build_code_2i40_9bits
 *
 *
 * Parameters:
 *    subNr             I: subframe number
 *    codvec            I: position of pulses
 *    dn_sign           I: sign of pulses
 *    cod               O: algebraic codebook vector
 *    h                 I: impulse response of weighted synthesis filter
 *    y                 I: filtered innovative code
 *    anap              O: analysis parameters
 *
 * Function:
 *    Builds the codeword, the filtered codeword and index of the
 *    codevector, based on the signs and positions of 2 pulses.
 *
 * Returns:
 *    void
 */
static void build_code_2i40_9bits(int16_t subNr, int32_t codvec[], float
				  dn_sign[], float cod[], float h[],
				  float y[], int16_t * anap)
{
	float s;
	float *p0, *p1;
	int32_t _sign[2];
	int32_t i, j, k, track, index, indx = 0, rsign = 0;
	int8_t first, *pt;

	pt = &f_trackTable[subNr + (subNr << 2)];
	memset(cod, 0, 160);

	for (k = 0; k < 2; k++) {
		/* read pulse position */
		i = codvec[k];

		/* read sign */
		j = (int32_t) dn_sign[i];

		/* index = pos/5 */
		index = i / 5;

		/* track = pos%5 */
		track = i % 5;
		first = pt[track];

		if (first == 0) {
			if (k == 0) {
				/*  position of 1st pulse   */
				track = 0;
			} else {
				track = 1;

				/*  position of 2nd pulse   */
				index <<= 3;
			}
		} else {
			if (k == 0) {
				track = 0;

				/*  position of 1st pulse, subset 2 */
				index += 64;
			} else {
				track = 1;
				index <<= 3;
			}
		}

		if (j > 0) {
			cod[i] = 0.9998779296875F;
			_sign[k] = 1;

			/*     sign information */
			rsign = rsign + (1 << track);
		} else {
			cod[i] = -1;
			_sign[k] = -1;
		}
		indx = indx + index;
	}
	p0 = h - codvec[0];
	p1 = h - codvec[1];

	for (i = 0; i < L_CODE; i++) {
		s = *p0++ * _sign[0];
		s += *p1++ * _sign[1];
		y[i] = s;
	}
	anap[0] = (int16_t) indx;
	anap[1] = (int16_t) rsign;
}

/*
 * code_2i40_9bits
 *
 *
 * Parameters:
 *    subNr             I: subframe number
 *    x                 I: target vector
 *    h                 I: impulse response of weighted synthesis filter
 *    T0                I: Pitch lag
 *    pitch_sharp       I: Last quantized pitch gain
 *    code              O: innovative codebook
 *    y                 O: filtered fixed codebook excitation
 *    anap              O: analysis parameters
 *
 * Function:
 *    Searches a 9 bit algebraic codebook containing 2 pulses
 *    in a frame of 40 samples.
 *
 *    The code length is 40, containing 2 nonzero pulses: i0...i1.
 *    All pulses can have two possible amplitudes: +1 or -1.
 *    Pulse i0 can have 8 possible positions, pulse i1 can have
 *    8 positions. Also coded is which track pair should be used,
 *    i.e. first or second pair. Where each pair contains 2 tracks.
 *
 * Returns:
 *    void
 */
static void code_2i40_9bits(int16_t subNr, float x[], float h[], int32_t T0,
			    float pitch_sharp, float code[], float y[],
			    int16_t * anap)
{
	float rr[L_CODE][L_CODE];
	float dn[L_CODE], dn_sign[L_CODE], dn2[L_CODE];
	int32_t codvec[2];
	int32_t i;

	if ((T0 < L_CODE) && (pitch_sharp != 0.0F))
		for (i = T0; i < L_CODE; i++) {
			h[i] += h[i - T0] * pitch_sharp;
		}
	cor_h_x(h, x, dn);
	set_sign(dn, dn_sign, dn2, 8);
	cor_h(h, dn_sign, rr);
	search_2i40_9bits(subNr, dn, rr, codvec);
	build_code_2i40_9bits(subNr, codvec, dn_sign, code, h, y, anap);

	/*
	 * Compute innovation vector gain.
	 * Include fixed-gain pitch contribution into code[].
	 */
	if ((T0 < L_CODE) && (pitch_sharp != 0.0F))
		for (i = T0; i < L_CODE; i++) {
			code[i] += code[i - T0] * pitch_sharp;
		}
}

/*
 * cbsearch
 *
 *
 * Parameters:
 *    mode              I: AMR mode
 *    subnr             I: Subframe
 *    x                 I: Target vector
 *    h                 B: Impulse response of weighted synthesis filter
 *    T0                I: Pitch lag
 *    pitch_sharp       I: Last quantized pitch gain
 *    gain_pit          I: Algebraic codebook gain
 *    code              O: Innovative codebook
 *    y                 O: Filtered fixed codebook excitation
 *    res2              I: residual after long term prediction
 *    anap              O: Signs and positions of the pulses
 *
 * Function:
 *    Innovative codebook search (find index and gain)
 *
 * Returns:
 *    void
 */
static void cbsearch(enum Mode mode, int16_t subnr, float x[],
		     float h[], int32_t T0, float pitch_sharp,
		     float gain_pit, float code[], float y[],
		     float * res2, int16_t ** anap)
{
	//switch (mode) {
	//case MR475:
	//case MR515:
		code_2i40_9bits(subnr, x, h, T0, pitch_sharp, code, y, *anap);
		(*anap) += 2;
	 //	break;
       /*
	case MR59:
		code_2i40_11bits(x, h, T0, pitch_sharp, code, y, *anap);
		(*anap) += 2;
		break;
	case MR67:
		code_3i40_14bits(x, h, T0, pitch_sharp, code, y, *anap);
		(*anap) += 2;
		break;
	case MR74:
	case MR795:
		code_4i40_17bits(x, h, T0, pitch_sharp, code, y, *anap);
		(*anap) += 2;
		break;
	case MR102:
		code_8i40_31bits(x, res2, h, T0, pitch_sharp, code, y, *anap);
		*anap += 7;
		break;
	default:
		code_10i40_35bits(x, res2, h, T0, gain_pit, code, y, *anap);
		*anap += 10;
	}
       */
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
	y = (f_log2_table[i] << 16) - a * (f_log2_table[i] - f_log2_table[i + 1]);
	*fraction = y >> 16;
	*exponent = 30 - exp;
	return;
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
	x = f_pow2_table[i] << 16;

	/* table[i] - table[i+1] */
	tmp = f_pow2_table[i] - f_pow2_table[i + 1];

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
 * gc_pred
 *
 *
 * Parameters:
 *    past_qua_en       I: MA predictor
 *    mode              I: AMR mode
 *    code              I: innovative codebook vector
 *    gcode0            O: predicted gain factor
 *    en                I: innovation energy (only calculated for MR795)
 *
 * Function:
 *    MA prediction of the innovation energy
 *
 *    Mean removed innovation energy (dB) in subframe n
 *                          N-1
 *    E(n) = 10*log(gc*gc * SUM[(code(i) * code(i)]/N) - E_mean
 *                          i=0
 *    N=40
 *
 *    Mean innovation energy (dB)
 *                   N-1
 *    Ei(n) = 10*log(SUM[(code(i) * code(i)]/N)
 *                   i=0
 *
 *    Predicted energy
 *             4
 *    Ep(n) = SUM[b(i) * R(n-i)]
 *            i=1
 *    b = [0.68 0.58 0.34 0.19]
 *    R(k) is quantified prediction error at subframe k
 *
 *    E_Mean = 36 dB (MR122)
 *
 *    Predicted gain gc is found by
 *
 *    gc = POW[10, 0.05 * (Ep(n) + E_mean - Ei)]
 *
 * Returns:
 *    void
 */
static void gc_pred(int32_t * past_qua_en, enum Mode mode, float * code,
		    int32_t * gcode0_exp, int32_t * gcode0_fra, float * en)
{
	float ener_code;
	int32_t exp, frac, ener,  tmp;
	int exp_code;

	/* energy of code */
	ener_code = Dotproduct40(code, code);
        /*
	if (mode == MR122) {

		ener = (int32_t) (ener_code * 33554432);
		//ener_code = ener_code / lcode; lcode = 40; 1/40 = 26214 Q20
		ener = ((ener + 0x00008000L) >> 16) * 52428;

		Log2(ener, &exp, &frac);
		ener = ((exp - 30) << 16) + (frac << 1);

		ener_tmp = 44 * f_qua_gain_code_MR122[past_qua_en[0]];
		ener_tmp += 37 * f_qua_gain_code_MR122[past_qua_en[1]];
		ener_tmp += 22 * f_qua_gain_code_MR122[past_qua_en[2]];
		ener_tmp += 12 * f_qua_gain_code_MR122[past_qua_en[3]];

		ener_tmp = ener_tmp << 1;
		ener_tmp += 783741L;

		//
		// * predicted codebook gain
		// * gc0 = Pow10( (ener*constant - ener_code*constant) / 20 )
		// *     = Pow2(ener-ener_code)
		// *     = Pow2(int(d)+frac(d))
		//
		ener = (ener_tmp - ener) >> 1;	// Q16
		*gcode0_exp = ener >> 16;
		*gcode0_fra = (ener >> 1) - (*gcode0_exp << 15);
	} else

        */

        {
		ener = (int32_t) (ener_code * 134217728);
		if (ener < 0)
			ener = 0x7fffffff;

		frexpf((float) ener, &exp_code);
		exp_code = 31 - exp_code;
		ener <<= exp_code;

		Log2_norm(ener, exp_code, &exp, &frac);

		tmp = (exp * (-49320)) + (((frac * (-24660)) >> 15) << 1);

	       //	if (mode == MR102) {
			/* mean = 33 dB */
		 //	tmp += 2134784;	/* Q14 */
	//} else if (mode == MR795) {
			/* mean = 36 dB */
	       //		tmp += 2183936;	/* Q14 */

		//	*en = (float) ener_code;
	       //	} else if (mode == MR74) {
			/* mean = 30 dB */
	       //		tmp += 2085632;	/* Q14 */
	       //	} else if (mode == MR67) {
			/* mean = 28.75 dB */
		//	tmp += 2065152;	/* Q14 */
	       //	} else
                {	/* MR59, MR515, MR475 */

			/* mean = 33 dB */
			tmp += 2134784;	/* Q14 */
		}

		tmp = tmp << 9;

		tmp += 5571 * f_qua_gain_code[past_qua_en[0]];
		tmp += 4751 * f_qua_gain_code[past_qua_en[1]];
		tmp += 2785 * f_qua_gain_code[past_qua_en[2]];
		tmp += 1556 * f_qua_gain_code[past_qua_en[3]];

		tmp = tmp >> 15;	/* Q8  */

		/*
		 * gcode0 = pow(10.0, gcode0/20)
		 *        = pow(2, 3.3219*gcode0/20)
		 *        = pow(2, 0.166*gcode0)
		 */
		/* 5439 Q15 = 0.165985                                        */
		/* (correct: 1/(20*log10(2)) 0.166096 = 5443 Q15)             */
		/* For IS641 bitexactness */
		//if (mode == MR74) {
			/* Q8 * Q15 -> Q24 */
		 //	tmp = tmp * 10878;
		//} else
               // {
			/* Q8 * Q15 -> Q24 */
			tmp = tmp * 10886;
	       //	}
		tmp = tmp >> 9;	/* -> Q15 */

		*gcode0_exp = tmp >> 15;
		*gcode0_fra = tmp - (*gcode0_exp * 32768);
	}
}

/*
 * calc_filt_energies
 *
 *
 * Parameters:
 *    mode              I: AMR mode
 *    xn                I: LTP target vector
 *    xn2               I: CB target vector
 *    y1                I: Adaptive codebook
 *    y2                I: Filtered innovative vector
 *    gCoeff            I: Correlations <xn y1> <y1 y1>
 *    coeff             O: energy coefficients
 *    cod_gain          O: optimum codebook gain
 *
 * Function:
 *    Calculation of several energy coefficients for filtered excitation signals
 *
 *    Compute coefficients need for the quantization and the optimum
 *    codebook gain gcu (for MR475 only).
 *
 *       coeff[0] =    y1 y1
 *       coeff[1] = -2 xn y1
 *       coeff[2] =    y2 y2
 *       coeff[3] = -2 xn y2
 *       coeff[4] =  2 y1 y2
 *
 *
 *       gcu = <xn2, y2> / <y2, y2> (0 if <xn2, y2> <= 0)
 *
 *    Product <y1 y1> and <xn y1> have been computed in G_pitch() and
 *    are in vector gCoeff[].
 *
 * Returns:
 *    void
 */
static void calc_filt_energies(enum Mode mode, float xn[], float xn2[],
			       float y1[], float y2[], float gCoeff[],
			       float coeff[], float * cod_gain)
{
	float sum, ener_init = 0.01F;

       //	if ((mode == MR795) || (mode == MR475))
		ener_init = 0;
	coeff[0] = gCoeff[0];
	coeff[1] = -2.0F * gCoeff[1];

	/* Compute scalar product <y2[],y2[]> */
	sum = (float) Dotproduct40(y2, y2);
	sum += ener_init;
	coeff[2] = sum;

	/* Compute scalar product -2*<xn[],y2[]> */
	sum = (float) Dotproduct40(xn, y2);
	sum += ener_init;
	coeff[3] = -2.0F * sum;

	/* Compute scalar product 2*<y1[],y2[]> */
	sum = (float) Dotproduct40(y1, y2);
	sum += ener_init;
	coeff[4] = 2.0F * sum;

	//if ((mode == MR475) || (mode == MR795)) {
		/* Compute scalar product <xn2[],y2[]> */
		sum = (float) Dotproduct40(xn2, y2);

		if (sum <= 0) {
			*cod_gain = 0;
		} else {
			/*
			 * gcu = <xn2, y2> / <y2, y2>
			 */
			*cod_gain = sum / coeff[2];
		}
       //	}
}

/*
 * MR475_update_unq_pred
 *
 *
 * Parameters:
 *    past_qua_en       I: MA predictor memory, 20*log10(qua_err)
 *    gcode0            I: predicted CB gain
 *    cod_gain          I: optimum codebook gain
 *
 * Function:
 *    Use optimum codebook gain and update "unquantized"
 *    gain predictor with the (bounded) prediction error
 *
 *    Prediction error is given by:
 *
 *       R(n) = E(n) - E_pred(n) = 20 * log(cf),
 *
 *    where correction factor cf between the gain Gc and
 *    the estimated pne Gc' is given by:
 *
 *       cf = Gc/Gc'.
 *
 * Returns:
 *    void
 */
static void MR475_update_unq_pred(int32_t * past_qua_en, float gcode0, float
				  cod_gain)
{
	float qua_ener, pred_err_fact;
	int32_t i, index, energy, max, s;

	if (cod_gain <= 0) {
		/*MIN_QUA_ENER */
		qua_ener = -32.0F;
	} else {
		if (gcode0 != 0) {
			pred_err_fact = cod_gain / gcode0;
		} else {
			pred_err_fact = 10.0F;
		}

		if (pred_err_fact < 0.0251189F) {
			/*MIN_QUA_ENER */
			qua_ener = -32.0F;
		} else if (pred_err_fact > 7.8125F) {
			/*MAX_QUA_ENER */
			qua_ener = 17.8558F;
		} else {
			qua_ener = (float) (20.0F * log10f(pred_err_fact));
		}
	}
	energy = (int32_t) (qua_ener * 1024 + 0.5F);
	max = abs(energy - f_qua_gain_code[0]);
	index = 0;
	/* find match from table */
	for (i = 1;
	     i <
	     NB_QUA_CODE + VQ_SIZE_HIGHRATES + VQ_SIZE_LOWRATES +
	     MR475_VQ_SIZE * 2 + 3; i++) {
		s = abs(energy - f_qua_gain_code[i]);
		if (s < max) {
			max = s;
			index = i;
			if (s == 0) {
				break;
			}
		}
	}
	/* update MA predictor memory */
	for (i = 3; i > 0; i--) {
		past_qua_en[i] = past_qua_en[i - 1];
	}
	past_qua_en[0] = index;
}

/*
 * MR475_gain_quant
 *
 *
 * Parameters:
 *    past_qua_en          B: MA predictor memory, 20*log10(qua_err)
 *    sf0_gcode0_exp(fra)  I: predicted CB gain from subframe 0 (or 2)
 *    sf0_coeff            I: energy coeff. from subframe 0 (or 2)
 *    sf0_target_en        I: target energy from subframe 0 (or 2)
 *    sf1_code_nosharp     I: innovative codebook vector (L_SUBFR)
 *                            (without pitch sharpening)
 *                            from subframe 1 (or 3)
 *    sf1_gcode0_exp(fra)  I: predicted CB gain from subframe 1 (or 3)
 *    sf1_coeff            I: energy coeff. subframe 1 (or 3)
 *    sf1_target_en        I: target energy from subframe 1 (or 3)
 *    gp_limit             I: pitch gain limit
 *    sf0_gain_pit         O: Pitch gain subframe 0 (or 2)
 *    sf0_gain_cod         O: Code gain subframe 0 (or 2)
 *    sf1_gain_pit         O: Pitch gain subframe 1 (or 3)
 *    sf1_gain_cod         O: Code gain subframe 1 (or 3)
 *
 * Function:
 *    Quantization of pitch and codebook gains for two subframes
 *    (using predicted codebook gain)
 *
 * Returns:
 *    index             index of quantization
 */
static int16_t MR475_gain_quant(int32_t * past_qua_en, int32_t sf0_gcode0_exp,
			       int32_t sf0_gcode0_fra, float sf0_coeff[],
			       float sf0_target_en,
			       float sf1_code_nosharp[],
			       int32_t sf1_gcode0_exp, int32_t sf1_gcode0_fra,
			       float sf1_coeff[], float sf1_target_en,
			       float gp_limit, float * sf0_gain_pit,
			       float * sf0_gain_cod, float * sf1_gain_pit,
			       float * sf1_gain_cod)
{
	float temp, temp2, g_pitch, g2_pitch, g_code, g2_code, g_pit_cod,
	    dist_min, sf0_gcode0, sf1_gcode0;
	const float *p;
	int32_t i, tmp, g_code_tmp, gcode0, index = 0;

	sf0_gcode0 = (float) Pow2(sf0_gcode0_exp, sf0_gcode0_fra);
	sf1_gcode0 = (float) Pow2(sf1_gcode0_exp, sf1_gcode0_fra);

	if ((sf0_target_en * 2.0F) < sf1_target_en) {
		sf0_coeff[0] *= 2.0F;
		sf0_coeff[1] *= 2.0F;
		sf0_coeff[2] *= 2.0F;
		sf0_coeff[3] *= 2.0F;
		sf0_coeff[4] *= 2.0F;
	} else if (sf0_target_en > (sf1_target_en * 4.0F)) {
		sf1_coeff[0] *= 2.0F;
		sf1_coeff[1] *= 2.0F;
		sf1_coeff[2] *= 2.0F;
		sf1_coeff[3] *= 2.0F;
		sf1_coeff[4] *= 2.0F;
	}

	/*
	 * Codebook search:
	 * For each pair (g_pitch, g_fac) in the table calculate the
	 * terms t[0..4] and sum them up; the result is the mean squared
	 * error for the quantized gains from the table. The index for the
	 * minimum MSE is stored and finally used to retrieve the quantized
	 * gains
	 */
	dist_min = FLT_MAX;
	p = &f_table_gain_MR475[0];

	for (i = 0; i < MR475_VQ_SIZE; i++) {
		/* subframe 0 (and 2) calculations */
		g_pitch = *p++;
		g_code = *p++;
		g_code *= sf0_gcode0;
		g2_pitch = g_pitch * g_pitch;
		g2_code = g_code * g_code;
		g_pit_cod = g_code * g_pitch;
		temp = sf0_coeff[0] * g2_pitch;
		temp += sf0_coeff[1] * g_pitch;
		temp += sf0_coeff[2] * g2_code;
		temp += sf0_coeff[3] * g_code;
		temp += sf0_coeff[4] * g_pit_cod;
		temp2 = g_pitch - gp_limit;

		/* subframe 1 (and 3) calculations */
		g_pitch = *p++;
		g_code = *p++;

		if (temp2 <= 0 && (g_pitch <= gp_limit)) {
			g_code *= sf1_gcode0;
			g2_pitch = g_pitch * g_pitch;
			g2_code = g_code * g_code;
			g_pit_cod = g_code * g_pitch;
			temp += sf1_coeff[0] * g2_pitch;
			temp += sf1_coeff[1] * g_pitch;
			temp += sf1_coeff[2] * g2_code;
			temp += sf1_coeff[3] * g_code;
			temp += sf1_coeff[4] * g_pit_cod;

			/*
			 * store table index if MSE for this index is lower
			 * than the minimum MSE seen so far
			 */
			if (temp < dist_min) {
				dist_min = temp;
				index = i;
			}
		}
	}

	/*
	 *  read quantized gains and update MA predictor memories
	 *
	 * for subframe 0, the pre-calculated gcode0 is the same
	 * as one calculated from the "real" predictor using quantized gains
	 */
	tmp = index << 2;
	p = &f_table_gain_MR475[tmp];
	*sf0_gain_pit = *p++;
	g_code_tmp = (int32_t) (*p++ * 4096 + 0.5F);

	gcode0 = Pow2(14, sf0_gcode0_fra);
	if (sf0_gcode0_exp < 11) {
		*sf0_gain_cod =
		    (float) ((g_code_tmp * gcode0) >> (25 - sf0_gcode0_exp));
	} else {
		i = ((g_code_tmp * gcode0) << (sf0_gcode0_exp - 9));

		if ((i >> (sf0_gcode0_exp - 9)) != (g_code_tmp * gcode0)) {
			*sf0_gain_cod = (float) 0x7FFF;
		} else {
			*sf0_gain_cod = (float) (i >> 16);
		}
	}

	*sf0_gain_cod *= 0.5F;

	for (i = 3; i > 0; i--) {
		past_qua_en[i] = past_qua_en[i - 1];
	}
	past_qua_en[0] =
	    NB_QUA_CODE + VQ_SIZE_HIGHRATES + VQ_SIZE_LOWRATES + (index << 1);

	/*
	 * calculate new predicted gain for subframe 1 (this time using
	 * the real, quantized gains)
	 */
	gc_pred(past_qua_en, MR475, sf1_code_nosharp, &sf1_gcode0_exp,
		&sf1_gcode0_fra, &sf0_gcode0);

	tmp += 2;
	p = &f_table_gain_MR475[tmp];
	*sf1_gain_pit = *p++;
	g_code_tmp = (int32_t) (*p++ * 4096 + 0.5F);

	gcode0 = Pow2(14, sf1_gcode0_fra);
	if (sf1_gcode0_exp < 11) {
		*sf1_gain_cod =
		    (float) ((g_code_tmp * gcode0) >> (25 - sf1_gcode0_exp));
	} else {
		i = ((g_code_tmp * gcode0) << (sf1_gcode0_exp - 9));

		if ((i >> (sf1_gcode0_exp - 9)) != (g_code_tmp * gcode0)) {
			*sf1_gain_cod = (float) 0x7FFF;
		} else {
			*sf1_gain_cod = (float) (i >> 16);
		}
	}

	*sf1_gain_cod *= 0.5F;

	for (i = 3; i > 0; i--) {
		past_qua_en[i] = past_qua_en[i - 1];
	}
	past_qua_en[0] =
	    NB_QUA_CODE + VQ_SIZE_HIGHRATES + VQ_SIZE_LOWRATES + (index << 1) +
	    1;

	return (int16_t) index;
}




/*
 * gainQuant
 *
 *
 * Parameters:
 *    gcode0            I: predicted CB gain
 *    coeff             I: energy coefficients (5)
 *    gp_limit          I: pitch gain limit
 *    gain_pit          O: Pitch gain
 *    gain_cod          O: Code gain
 *    qua_ener          O: quantized energy error,
 *    mode              I: AMR mode
 *    even_subframe     I: even subframe indicator flag
 *    past_qua_en       B: past quantized energies [4]
 *    past_qua_en_unq   B: past energies [4]
 *    sf0_coeff         B: energy coefficients subframe 0 (or 2)
 *    sf0_target_en     B: target energy from subframe 0 (or 2)
 *    sf0_gcode0        B: predicted gain factor subframe 0 (or 2)
 *    gain_idx_ptr      B: gain index pointer
 *    sf0_gain_pit      B: Pitch gain subframe 0 (or 2)
 *    sf0_gain_cod      B: Code gain subframe 0 (or 2)
 *    res               I: LP residual
 *    exc               I: LTP excitation (unfiltered)
 *    code              I: innovative codebook vector
 *    xn                I: Target vector
 *    xn2               I: CB target vector
 *    y1                I: Adaptive codebook
 *    y2                I: Filtered innovative vector
 *    gCoeff            I: Correlations <xn y1> <y1 y1>
 *    gp_limit          I: pitch gain limit
 *    gain_pit          O: Pitch gain
 *    gain_cod          O: Code gain
 *    prev_gc           B: B: previous code gain
 *    onset             B: onset indicator
 *    ltpg_mem          B: stored past LTP coding gains
 *    prev_alpha        B: previous gain adaptation factor
 *    anap              B: Index of quantization
 *
 * Function:
 *    Quantization of gains
 *
 * Returns:
 *    index             index of quantization
 */
static void gainQuant(enum Mode mode, int32_t even_subframe, int32_t *
		      past_qua_en, int32_t * past_qua_en_unq,
		      float * sf0_coeff, float * sf0_target_en,
		      int32_t * sf0_gcode0_exp, int32_t * sf0_gcode0_fra,
		      int16_t ** gain_idx_ptr, float * sf0_gain_pit,
		      float * sf0_gain_cod, float * res, float * exc,
		      float code[], float xn[], float xn2[], float y1[],
		      float y2[], float gCoeff[], float gp_limit,
		      float * gain_pit, float * gain_cod, float * prev_gc,
		      int16_t * onset, float * ltpg_mem, float * prev_alpha,
		      int16_t ** anap)
{
	float coeff[5];
	float gcode0, cod_gain, en = 0;
	int32_t exp, frac; //, qua_ener_index;

       //	if (mode == MR475) {
		if (even_subframe != 0) {
			/*
			 * save position in output parameter stream and current
			 * state of codebook gain predictor
			 */
			*gain_idx_ptr = (*anap)++;
			past_qua_en_unq[0] = past_qua_en[0];
			past_qua_en_unq[1] = past_qua_en[1];
			past_qua_en_unq[2] = past_qua_en[2];
			past_qua_en_unq[3] = past_qua_en[3];

			/*
			 * predict codebook gain (using "unquantized" predictor)
			 * (note that code[] is unsharpened in MR475)
			 */
			gc_pred(past_qua_en, mode, code, sf0_gcode0_exp,
				sf0_gcode0_fra, &en);
			gcode0 =
			    (float) Pow2(*sf0_gcode0_exp, *sf0_gcode0_fra);

			/*
			 * calculate energy coefficients for quantization
			 * and store them in state structure (will be used
			 * in next subframe when real quantizer is run)
			 */
			calc_filt_energies(mode, xn, xn2, y1, y2, gCoeff,
					   sf0_coeff, &cod_gain);

			/* store optimum codebook gain */
			*gain_cod = cod_gain;
			*sf0_target_en = (float) Dotproduct40(xn, xn);

			/*
			 * calculate optimum codebook gain and update
			 * "unquantized" predictor
			 */
			MR475_update_unq_pred(past_qua_en_unq, gcode0,
					      cod_gain);

			/* the real quantizer is not run here... */
		} else {
			/*
			 * predict codebook gain (using "unquantized" predictor)
			 * (note that code[] is unsharpened in MR475)
			 */
			gc_pred(past_qua_en_unq, mode, code, &exp, &frac, &en);

			/* calculate energy coefficients for quantization */
			calc_filt_energies(mode, xn, xn2, y1, y2, gCoeff, coeff,
					   &cod_gain);
			en = (float) Dotproduct40(xn, xn);

			/* run real (4-dim) quantizer and update real gain predictor */
			**gain_idx_ptr =
			    MR475_gain_quant(past_qua_en, *sf0_gcode0_exp,
					     *sf0_gcode0_fra, sf0_coeff,
					     *sf0_target_en, code, exp, frac,
					     coeff, en, gp_limit, sf0_gain_pit,
					     sf0_gain_cod, gain_pit, gain_cod);
		}
       //	}


        /*
        else {
		//
		 // predict codebook gain and quantize
		 //  (also compute normalized CB innovation energy for MR795)

		gc_pred(past_qua_en, mode, code, &exp, &frac, &en);

		if (mode == MR122) {
			//
		       //	 * Compute the innovative codebook gain.
		      //	 * The innovative codebook gain is given by
		      //	 * g = <xn2[], y2[]> / <y2[], y2[]>
		      //	 * where xn2[] is the target vector,
		       //	 * y2[] is the filtered innovative
		       //	 * codevector


			gcode0 = (float) Pow2(exp, frac);
			// saturation at decoder
			if (gcode0 > 2047.9375F)
				gcode0 = 2047.9375F;

			*gain_cod =
			    (float) (Dotproduct40(xn2, y2) /
				       (Dotproduct40(y2, y2) + 0.01F));

			if (*gain_cod < 0)
				*gain_cod = 0.0F;
			*(*anap)++ =
			    q_gain_code(gcode0, gain_cod, &qua_ener_index);
		} else {
			// calculate energy coefficients for quantization
			calc_filt_energies(mode, xn, xn2, y1, y2, gCoeff, coeff,
					   &cod_gain);

			if (mode == MR795) {
				MR795_gain_quant(prev_gc, onset, ltpg_mem,
						 prev_alpha, res, exc, code,
						 coeff, en, exp, frac, cod_gain,
						 gp_limit, gain_pit, gain_cod,
						 &qua_ener_index, anap);
			} else {

				*(*anap)++ =
				    Qua_gain(mode, exp, frac, coeff, gp_limit,
					     gain_pit, gain_cod,
					     &qua_ener_index);
			}
		}

		//
		 // update table of past quantized energies
		 //
		for (i = 3; i > 0; i--) {
			past_qua_en[i] = past_qua_en[i - 1];
		}
		past_qua_en[0] = qua_ener_index;

	}

        */
}

/*
 * subframePostProc
 *
 *
 * Parameters:
 *    speech            I: Speech segment
 *    i_subfr           I: Subframe number
 *    gain_pit          I: Pitch gain
 *    gain_code         I: Decoded innovation gain
 *    a_q               I: A(z) quantized for the 4 subframes
 *    synth             I: Local synthesis
 *    xn                I: Target vector for pitch search
 *    code              I: Fixed codebook exitation
 *    y1                I: Filtered adaptive exitation
 *    y2                I: Filtered fixed codebook excitation
 *    mem_syn           B: memory of synthesis filter
 *    mem_err           O: pointer to error signal
 *    mem_w0            O: memory of weighting filter
 *    exc               O: long term prediction residual
 *    sharp             O: pitch sharpening value
 *
 * Function:
 *    Subframe post processing
 *
 *    Memory update (all modes)
 *    An update of the states of the synthesis and weighting filters is needed
 *   in order to compute the target signal in the next subframe.
 *   After the two gains are quantified, the excitation signal, u(n),
 *   in the present subframe is found by:
 *
 *   u(n) = Gp_q * v(n) + Gc_q * c(n), n = 0, ..., 39,
 *
 *   where Gp_q and Gc_q are the quantified adaptive and fixed codebook gains,
 *   respectively, v(n) the adaptive codebook vector
 *   (interpolated past excitation), and c(n) is the fixed codebook vector
 *   (algebraic code including pitch sharpening). The states of the filters
 *   can be updated by filtering the signal res_lp(n) - u(n)
 *   (difference between residual and excitation) through the filters
 *   1 / A_q(z) and A(z/g1) / A(z/g2) for the 40-sample subframe and saving
 *   the states of the filters. This would require 3 filterings.
 *   A simpler approach which requires only one filtering is as follows.
 *   The local synthesized speech, S_s(n), is computed by filtering
 *   the excitation signal through 1 / A_q(z). The output of the filter
 *   due to the input res_lp(n) - u(n) is equivalent to e(n) = S(n) - S_s(n).
 *   So the states of the synthesis filter 1 / A_q(z) are given by
 *   e(n), n = 30, ..., 39. Updating the states of the filter A(z/g1) / A(z/g2)
 *   can be done by filtering the error signal e(n) through this filter
 *   to find the perceptually weighted error ew(n). However, the signal ew(n)
 *   can be equivalently found by:
 *
 *   ew(n) = x(n) - Gp_q * y(n) - Gc_q(n) * z(n)
 *
 *   Since the signals x(n), y(n) and z(n) are available, the states of
 *   the weighting filter are updated by computing ew(n) for n = 30, ..., 39.
 *   This saves two filterings.
 *
 * Returns:
 *    void
 */
static void subframePostProc(float * speech, int16_t i_subfr, float gain_pit,
			     float gain_code, float * a_q, float synth[],
			     float xn[], float code[], float y1[],
			     float y2[], float * mem_syn, float * mem_err,
			     float * mem_w0, float * exc, float * sharp)
{
	int32_t i, j;

	/*
	 * Update pitch sharpening "sharp" with quantized gain_pit
	 */
	*sharp = gain_pit;
	if (*sharp > 0.794556F) {
		*sharp = 0.794556F;
	}

	/* Find the total excitation */
	for (i = 0; i < L_SUBFR; i += 4) {
		exc[i + i_subfr] =
		    (float)
		    floorf((gain_pit * exc[i + i_subfr] + gain_code * code[i]) +
			  0.5F);
		exc[i + i_subfr + 1] =
		    (float)
		    floorf((gain_pit * exc[i + i_subfr + 1] +
			   gain_code * code[i + 1]) + 0.5F);
		exc[i + i_subfr + 2] =
		    (float)
		    floorf((gain_pit * exc[i + i_subfr + 2] +
			   gain_code * code[i + 2]) + 0.5F);
		exc[i + i_subfr + 3] =
		    (float)
		    floorf((gain_pit * exc[i + i_subfr + 3] +
			   gain_code * code[i + 3]) + 0.5F);
	}

	/* The local synthesis speech */
	Syn_filt(a_q, &exc[i_subfr], &synth[i_subfr], mem_syn, 1);

	for (i = L_SUBFR - M, j = 0; i < L_SUBFR; i++, j++) {
		/* e(n) = S(n) - S_s(n) */
		mem_err[j] = speech[i_subfr + i] - synth[i_subfr + i];

		/* ew(n) = x(n) - Gp_q * y(n) - Gc_q(n) * z(n) */
		mem_w0[j] = xn[i] - y1[i] * gain_pit - y2[i] * gain_code;
	}
}

/*
 * Convolve
 *
 *
 * Parameters:
 *    x                 I: First input
 *    h                 I: second input
 *    y                 O: output
 *
 * Function:
 *    Convolution
 *
 * Returns:
 *    void
 */
static void Convolve(float x[], float h[], float y[])
{
	int32_t i, n;
	float s;

	for (n = 0; n < L_SUBFR; n++) {
		s = 0.0F;

		for (i = 0; i <= n; i++) {
			s += x[i] * h[n - i];
		}
		y[n] = s;
	}
	return;
}

/*
 * complex_estimate_adapt
 *
 *
 * Parameters:
 *    st->best_corr_hp  I: complex background detector
 *    st->corr_hp_fast  B: long term complex signal estimate
 *    low_power         I: very low level flag of the input frame
 *
 * Function:
 *    Update/adapt of complex signal estimate
 *
 * Returns:
 *    void
 */

static void complex_estimate_adapt(vadState * st, int16_t low_power)
{
	float alpha;

	/* adapt speed on own state */
	/* decrease */
	if (st->best_corr_hp < st->corr_hp_fast) {
		/* low state */
		if (st->corr_hp_fast < CVAD_THRESH_ADAPT_HIGH) {
			alpha = CVAD_ADAPT_FAST;
		}

		/* high state */
		else {
			alpha = CVAD_ADAPT_REALLY_FAST;
		}
	}

	/* increase */
	else {
		if (st->corr_hp_fast < CVAD_THRESH_ADAPT_HIGH) {
			alpha = CVAD_ADAPT_FAST;
		} else {
			alpha = CVAD_ADAPT_SLOW;
		}
	}
	st->corr_hp_fast =
	    st->corr_hp_fast - alpha * st->corr_hp_fast +
	    alpha * st->best_corr_hp;

	if (st->corr_hp_fast < CVAD_MIN_CORR) {
		st->corr_hp_fast = CVAD_MIN_CORR;
	}

	if (low_power != 0) {
		st->corr_hp_fast = CVAD_MIN_CORR;
	}
}


/*
 * complex_vad
 *
 *
 * Parameters:
 *    st->complex_high        B: 1 if (corr_hp_fast > CVAD_THRESH_ADAPT_HIGH)
 *    st->complex_low         B: 1 if (corr_hp_fast > CVAD_THRESH_ADAPT_LOW)
 *    low_power               I: flag power of the input frame
 *    st->best_corr_hp        I: complex background detector
 *    st->corr_hp_fast        B: long term complex signal estimate
 *    st->complex_hang_timer  B: complex hang timer
 *
 *
 * Function:
 *    Complex background decision
 *
 * Returns:
 *    void
 */

static int32_t complex_vad(vadState * st, int16_t low_power)
{
	st->complex_high = st->complex_high >> 1;
	st->complex_low = st->complex_low >> 1;

	if (low_power == 0) {
		if (st->corr_hp_fast > CVAD_THRESH_ADAPT_HIGH) {
			st->complex_high = st->complex_high | 0x00004000;
		}

		if (st->corr_hp_fast > CVAD_THRESH_ADAPT_LOW) {
			st->complex_low = st->complex_low | 0x00004000;
		}
	}

	if (st->corr_hp_fast > CVAD_THRESH_HANG) {
		st->complex_hang_timer += 1;
	} else {
		st->complex_hang_timer = 0;
	}
	return (int16_t) (((st->complex_high & 0x00007f80) == 0x00007f80)
			 || ((st->complex_low & 0x00007fff) == 0x00007fff));
}


/*
 * complex_vad
 *
 *
 * Parameters:
 *    st->complex_warning  I: flags for complex detection
 *    st->stat_count       B: stationary counter
 *    st->pitch            I: flags for pitch detection
 *    st->tone             I: flags indicating presence of a tone
 *    st->vadreg           I: intermediate VAD flags
 *    level                I: sub-band levels of the input frame
 *    st->ave_level        B: Average amplitude estimate
 *
 * Function:
 *    Control update of the background noise estimate
 *
 * Returns:
 *    void
 */

static void update_cntrl(vadState * st, float level[])
{
	float stat_rat, num, denom;
	float alpha;
	int32_t i;

	/*
	 * handle highband complex signal input  separately
	 * if ther has been highband correlation for some time
	 * make sure that the VAD update speed is low for a while
	 */
	if (st->complex_warning != 0) {
		if (st->stat_count < CAD_MIN_STAT_COUNT) {
			st->stat_count = CAD_MIN_STAT_COUNT;
		}
	}

	/*
	 * NB stat_count is allowed to be decreased by one below again
	 * deadlock in speech is not possible unless the signal is very
	 * complex and need a high rate
	 * if fullband pitch or tone have been detected for a while, initialize stat_count
	 */
	if (((st->pitch & 0x6000) == 0x6000) || ((st->tone & 0x00007c00) ==
						 0x7c00)) {
		st->stat_count = STAT_COUNT;
	} else {
		/* if 8 last vad-decisions have been "0", reinitialize stat_count */
		if ((st->vadreg & 0x7f80) == 0) {
			st->stat_count = STAT_COUNT;
		} else {
			stat_rat = 0;

			for (i = 0; i < COMPLEN; i++) {
				if (level[i] > st->ave_level[i]) {
					num = level[i];
					denom = st->ave_level[i];
				} else {
					num = st->ave_level[i];
					denom = level[i];
				}

				/* Limit nimimum value of num and denom to STAT_THR_LEVEL */
				if (num < STAT_THR_LEVEL) {
					num = STAT_THR_LEVEL;
				}

				if (denom < STAT_THR_LEVEL) {
					denom = STAT_THR_LEVEL;
				}
				stat_rat += num / denom * 64;
			}

			/* compare stat_rat with a threshold and update stat_count */
			if (stat_rat > STAT_THR) {
				st->stat_count = STAT_COUNT;
			} else {
				if ((st->vadreg & 0x4000) != 0) {
					if (st->stat_count != 0) {
						st->stat_count -= 1;
					}
				}
			}
		}
	}

	/* Update average amplitude estimate for stationarity estimation */
	alpha = ALPHA4;

	if (st->stat_count == STAT_COUNT) {
		alpha = 1.0F;
	} else if ((st->vadreg & 0x4000) == 0) {
		alpha = ALPHA5;
	}

	for (i = 0; i < COMPLEN; i++) {
		st->ave_level[i] += alpha * (level[i] - st->ave_level[i]);
	}
}


/*
 * noise_estimate_update
 *
 *
 * Parameters:
 *    st                      B: State struct
 *    level                   I: sub-band levels of the input frame
 *    st->vadreg              I: intermediate VAD flags
 *    st->pitch               I: flags for pitch detection
 *    st->complex_hang_count  I: signal is too complex for VAD
 *    st->stat_count          B: stationary counter
 *    st->old_level           B: signal levels of the previous frame
 *    st->bckr_est            B: noise estimate
 *
 * Function:
 *    Update of background noise estimate
 *
 * Returns:
 *    void
 */

static void noise_estimate_update(vadState * st, float level[])
{
	float alpha_up, alpha_down, bckr_add;
	int32_t i;

	/* Control update of bckr_est[] */
	update_cntrl(st, level);

	/* Choose update speed */
	bckr_add = 2;

	if (((0x7800 & st->vadreg) == 0) && ((st->pitch & 0x7800) == 0)
	    && (st->complex_hang_count == 0)) {
		alpha_up = ALPHA_UP1;
		alpha_down = ALPHA_DOWN1;
	} else {
		if ((st->stat_count == 0) && (st->complex_hang_count == 0)) {
			alpha_up = ALPHA_UP2;
			alpha_down = ALPHA_DOWN2;
		} else {
			alpha_up = 0;
			alpha_down = ALPHA3;
			bckr_add = 0;
		}
	}

	/* Update noise estimate (bckr_est) */
	for (i = 0; i < COMPLEN; i++) {
		float temp;

		temp = st->old_level[i] - st->bckr_est[i];

		/* update downwards */
		if (temp < 0) {
			st->bckr_est[i] =
			    (-2 + (st->bckr_est[i] + (alpha_down * temp)));

			/* limit minimum value of the noise estimate to NOISE_MIN */
			if (st->bckr_est[i] < NOISE_MIN) {
				st->bckr_est[i] = NOISE_MIN;
			}
		}

		/* update upwards */
		else {
			st->bckr_est[i] =
			    (bckr_add + (st->bckr_est[i] + (alpha_up * temp)
			     ));

			/* limit maximum value of the noise estimate to NOISE_MAX */
			if (st->bckr_est[i] > NOISE_MAX) {
				st->bckr_est[i] = NOISE_MAX;
			}
		}
	}

	/* Update signal levels of the previous frame (old_level) */
	for (i = 0; i < COMPLEN; i++) {
		st->old_level[i] = level[i];
	}
}


/*
 * hangover_addition
 *
 *
 * Parameters:
 *    noise_level             I: average level of the noise estimates
 *    low_power               I: flag power of the input frame
 *    st->burst_count         O: counter for the length of speech bursts
 *    st->hang_count          O: hangover counter
 *    st->complex_hang_count  B: signal is too complex for VAD
 *    st->complex_hang_timer  B: complex hang timer
 *    st->vadreg              I: intermediate VAD flags
 *    st->corr_hp_fast        I: long term complex signal estimate
 *
 * Function:
 *    Add hangover for complex signal or after speech bursts
 *
 * Returns:
 *    VAD_flag indicating final VAD decision
 */

static int16_t hangover_addition(vadState * st, float noise_level, int16_t
				low_power)
{
	int16_t hang_len, burst_len;

	/*
	 * Calculate burst_len and hang_len
	 * burst_len: number of consecutive intermediate vad flags with "1"-decision
	 * required for hangover addition
	 * hang_len:  length of the hangover
	 */
	if (noise_level > HANG_NOISE_THR) {
		burst_len = BURST_LEN_HIGH_NOISE;
		hang_len = HANG_LEN_HIGH_NOISE;
	} else {
		burst_len = BURST_LEN_LOW_NOISE;
		hang_len = HANG_LEN_LOW_NOISE;
	}

	/*
	 * if the input power (pow_sum) is lower than a threshold, clear
	 * counters and set VAD_flag to "0"  "fast exit"
	 */
	if (low_power != 0) {
		st->burst_count = 0;
		st->hang_count = 0;
		st->complex_hang_count = 0;
		st->complex_hang_timer = 0;
		return 0;
	}

	if (st->complex_hang_timer > CVAD_HANG_LIMIT) {
		if (st->complex_hang_count < CVAD_HANG_LENGTH) {
			st->complex_hang_count = CVAD_HANG_LENGTH;
		}
	}

	/* long time very complex signal override VAD output function */
	if (st->complex_hang_count != 0) {
		st->burst_count = BURST_LEN_HIGH_NOISE;
		st->complex_hang_count -= 1;
		return 1;
	} else {
		/* let hp_corr work in from a noise_period indicated by the VAD */
		if (((st->vadreg & 0x3ff0) == 0) && (st->corr_hp_fast >
						     CVAD_THRESH_IN_NOISE)) {
			return 1;
		}
	}

	/* update the counters (hang_count, burst_count) */
	if ((st->vadreg & 0x4000) != 0) {
		st->burst_count += 1;

		if (st->burst_count >= burst_len) {
			st->hang_count = hang_len;
		}
		return 1;
	} else {
		st->burst_count = 0;

		if (st->hang_count > 0) {
			st->hang_count -= 1;
			return 1;
		}
	}
	return 0;
}


/*
 * vad_decision
 *
 *
 * Parameters:
 *    st                      B: State struct
 *    level                   I: sub-band levels of the input frame
 *    pow_sum                 I: power of the input frame
 *    st->bckr_est            I: background noise components
 *    st->vadreg              I: intermediate VAD flags
 *    st->complex_warning     O: flags for complex detection
 *    st->speech_vad_decision O: speech VAD flag
 *
 * Function:
 *    Calculates VAD_flag
 *
 * Returns:
 *    VAD_flag indicating final VAD decision
 */

static int16_t vad_decision(vadState * st, float level[COMPLEN], float
			   pow_sum)
{
	float snr_sum, temp, vad_thr, noise_level;
	int32_t i;
	int16_t low_power_flag;

	/*
	 * Calculate squared sum of the input levels (level)
	 * divided by the background noise components (bckr_est).
	 */
	snr_sum = 0;

	for (i = 0; i < COMPLEN; i++) {
		temp = level[i] / st->bckr_est[i];
		snr_sum += temp * temp;
	}
	snr_sum = snr_sum * 56.8889F;

	/* Calculate average level of estimated background noise */
	noise_level =
	    st->bckr_est[0] + st->bckr_est[1] + st->bckr_est[2] +
	    st->bckr_est[3] + st->bckr_est[4] + st->bckr_est[5] +
	    st->bckr_est[6] + st->bckr_est[7] + st->bckr_est[8];
	noise_level = noise_level * 0.111111F;

	/* Calculate VAD threshold */
	vad_thr = VAD_SLOPE * (noise_level - VAD_P1) + VAD_THR_HIGH;

	if (vad_thr < VAD_THR_LOW) {
		vad_thr = VAD_THR_LOW;
	}

	/* Shift VAD decision register */
	st->vadreg >>= 1;

	/* Make intermediate VAD decision */
	if (snr_sum > vad_thr) {
		st->vadreg = st->vadreg | 0x4000;
	}

	/*
	 * primary vad decision made
	 * check if the input power (pow_sum) is lower than a threshold"
	 */
	if (pow_sum < VAD_POW_LOW) {
		low_power_flag = 1;
	} else {
		low_power_flag = 0;
	}

	/*
	 * update complex signal estimate st->corr_hp_fast and hangover reset timer using
	 * low_power_flag and corr_hp_fast and various adaptation speeds
	 */
	complex_estimate_adapt(st, low_power_flag);

	/* check multiple thresholds of the st->corr_hp_fast value */
	st->complex_warning = complex_vad(st, low_power_flag);

	/* Update speech subband vad background noise estimates */
	noise_estimate_update(st, level);

	/*
	 *  Add speech and complex hangover and return speech VAD_flag
	 *  long term complex hangover may be added
	 */
	st->speech_vad_decision =
	    hangover_addition(st, noise_level, low_power_flag);
	return (st->speech_vad_decision);
}


/*
 * level_calculation
 *
 *
 * Parameters:
 *    data              I: signal buffer
 *    sub_level         B: level calculate at the end of the previous frame/
 *                         level of signal calculated from the last
 *                         (count2 - count1) samples
 *    count1            I: number of samples to be counted
 *    count2            I: number of samples to be counted
 *    ind_m             I: step size for the index of the data buffer
 *    ind_a             I: starting index of the data buffer
 *    scale             I: scaling for the level calculation
 *
 * Function:
 *    Calculate signal level in a sub-band.
 *
 *    Level is calculated by summing absolute values of the input data.
 *
 * Returns:
 *    signal level
 */

static float level_calculation(float data[], float * sub_level, int16_t
				 count1, int16_t count2, int16_t ind_m,
				 int16_t ind_a, int16_t scale)
{
	float level, temp1;
	int32_t i;

	temp1 = 0;

	for (i = count1; i < count2; i++) {
		temp1 += (float) fabsf(data[ind_m * i + ind_a]);
	}
	level = temp1 + *sub_level;
	*sub_level = temp1;

	for (i = 0; i < count1; i++) {
		level += (float) fabsf(data[ind_m * i + ind_a]);
	}
	return (scale * level);
}


/*
 * filter3
 *
 *
 * Parameters:
 *    in0               B: input values; output low-pass part
 *    in1               B: input values; output high-pass part
 *    data              B: updated filter memory
 *
 * Function:
 *    Third-order half-band lowpass/highpass filter pair.
 *
 * Returns:
 *   void
 */

static void filter3(float * in0, float * in1, float * data)
{
	float temp1, temp2;

	temp1 = *in1 - (COEFF3 * *data);
	temp2 = *data + (COEFF3 * temp1);
	*data = temp1;
	*in1 = (*in0 - temp2) * 0.5F;
	*in0 = (*in0 + temp2) * 0.5F;
}


/*
 * filter5
 *
 *
 * Parameters:
 *    in0               B: input values; output low-pass part
 *    in1               B: input values; output high-pass part
 *    data              B: updated filter memory
 *
 * Function:
 *    Fifth-order half-band lowpass/highpass filter pair.
 *
 * Returns:
 *   void
 */

static void filter5(float * in0, float * in1, float data[])
{
	float temp0, temp1, temp2;

	temp0 = *in0 - (COEFF5_1 * data[0]);
	temp1 = data[0] + (COEFF5_1 * temp0);
	data[0] = temp0;
	temp0 = *in1 - (COEFF5_2 * data[1]);
	temp2 = data[1] + (COEFF5_2 * temp0);
	data[1] = temp0;
	*in0 = (temp1 + temp2) * 0.5F;
	*in1 = (temp1 - temp2) * 0.5F;
}


/*
 * first_filter_stage
 *
 *
 * Parameters:
 *    in                I: input signal
 *    out               O: output values,
 *                         every other output is low-pass part and
 *                         every other output is high-pass part
 *    data              B: updated filter memory
 *
 * Function:
 *    Calculate 5th order half-band lowpass/highpass filter pair
 *
 * Returns:
 *   void
 */
#ifndef VAD2
static void first_filter_stage(float in[], float out[], float data[])
{
	float temp0, temp1, temp2, temp3;
	float data0, data1;
	int32_t i;

	data0 = data[0];
	data1 = data[1];

	for (i = 0; i < L_SUBFR; i++) {
		temp0 = (in[4 * i + 0] * 0.25F) - (COEFF5_1 * data0);
		temp1 = data0 + (COEFF5_1 * temp0);
		temp3 = (in[4 * i + 1] * 0.25F) - (COEFF5_2 * data1);
		temp2 = data1 + (COEFF5_2 * temp3);
		out[4 * i + 0] = temp1 + temp2;
		out[4 * i + 1] = temp1 - temp2;
		data0 = (in[4 * i + 2] * 0.25F) - (COEFF5_1 * temp0);
		temp1 = temp0 + (COEFF5_1 * data0);
		data1 = (in[4 * i + 3] * 0.25F) - (COEFF5_2 * temp3);
		temp2 = temp3 + (COEFF5_2 * data1);
		out[4 * i + 2] = temp1 + temp2;
		out[4 * i + 3] = temp1 - temp2;
	}
	data[0] = data0;
	data[1] = data1;
}
#endif

/*
 * filter_bank
 *
 *
 * Parameters:
 *    in                I: input frame
 *    st->a_data5       B: filter memory
 *    st->a_data3       B: filter memory
 *    st->sub_level     B: level memory
 *    level             O: signal levels at each band
 *
 * Function:
 *    Divides input signal into 9-bands and calcultes level of the signal in each band
 *
 * Returns:
 *    void
 */
#ifndef VAD2
static void filter_bank(vadState * st, float in[], float level[])
{
	int32_t i;
	float tmp_buf[FRAME_LEN];

	/* calculate the filter bank */
	first_filter_stage(in, tmp_buf, st->a_data5[0]);

	for (i = 0; i < FRAME_LEN / 4; i++) {
		filter5(&tmp_buf[4 * i], &tmp_buf[4 * i + 2], st->a_data5[1]);
		filter5(&tmp_buf[4 * i + 1], &tmp_buf[4 * i + 3],
			st->a_data5[2]);
	}

	for (i = 0; i < FRAME_LEN / 8; i++) {
		filter3(&tmp_buf[8 * i + 0], &tmp_buf[8 * i + 4],
			&st->a_data3[0]);
		filter3(&tmp_buf[8 * i + 2], &tmp_buf[8 * i + 6],
			&st->a_data3[1]);
		filter3(&tmp_buf[8 * i + 3], &tmp_buf[8 * i + 7],
			&st->a_data3[4]);
	}

	for (i = 0; i < FRAME_LEN / 16; i++) {
		filter3(&tmp_buf[16 * i + 0], &tmp_buf[16 * i + 8],
			&st->a_data3[2]);
		filter3(&tmp_buf[16 * i + 4], &tmp_buf[16 * i + 12],
			&st->a_data3[3]);
	}

	/* calculate levels in each frequency band */
	/* 3000 - 4000 Hz */
	level[8] =
	    level_calculation(tmp_buf, &st->sub_level[8], FRAME_LEN / 4 - 8,
			      FRAME_LEN / 4, 4, 1, 1);

	/* 2500 - 3000 Hz */
	level[7] =
	    level_calculation(tmp_buf, &st->sub_level[7], FRAME_LEN / 8 - 4,
			      FRAME_LEN / 8, 8, 7, 2);

	/* 2000 - 2500 Hz */
	level[6] =
	    level_calculation(tmp_buf, &st->sub_level[6], FRAME_LEN / 8 - 4,
			      FRAME_LEN / 8, 8, 3, 2);

	/* 1500 - 2000 Hz */
	level[5] =
	    level_calculation(tmp_buf, &st->sub_level[5], FRAME_LEN / 8 - 4,
			      FRAME_LEN / 8, 8, 2, 2);

	/* 1000 - 1500 Hz */
	level[4] =
	    level_calculation(tmp_buf, &st->sub_level[4], FRAME_LEN / 8 - 4,
			      FRAME_LEN / 8, 8, 6, 2);

	/* 750 - 1000 Hz */
	level[3] =
	    level_calculation(tmp_buf, &st->sub_level[3], FRAME_LEN / 16 - 2,
			      FRAME_LEN / 16, 16, 4, 2);

	/* 500 - 750 Hz */
	level[2] =
	    level_calculation(tmp_buf, &st->sub_level[2], FRAME_LEN / 16 - 2,
			      FRAME_LEN / 16, 16, 12, 2);

	/* 250 - 500 Hz */
	level[1] =
	    level_calculation(tmp_buf, &st->sub_level[1], FRAME_LEN / 16 - 2,
			      FRAME_LEN / 16, 16, 8, 2);

	/* 0 - 250 Hz */
	level[0] =
	    level_calculation(tmp_buf, &st->sub_level[0], FRAME_LEN / 16 - 2,
			      FRAME_LEN / 16, 16, 0, 2);
}
#endif

/*
 * vad
 *
 *
 * Parameters:
 *    in_buf            I: samples of the input frame
 *    st                B: State struct
 *    st->pitch         B: flags for pitch detection
 *    st->complex_low   B: complex flag
 *
 * Function:
 *    Voice Activity Detection (VAD)
 *
 * Returns:
 *    VAD Decision, 1 = speech, 0 = noise
 */

static int16_t vad(vadState * st, float in_buf[])
{
	float level[COMPLEN];
	float pow_sum;
	int32_t i;

	/* Calculate power of the input frame. */
	pow_sum = 0L;

	for (i = -40; i < 120; i += 8) {
		pow_sum += in_buf[i] * in_buf[i];
		pow_sum += in_buf[i + 1] * in_buf[i + 1];
		pow_sum += in_buf[i + 2] * in_buf[i + 2];
		pow_sum += in_buf[i + 3] * in_buf[i + 3];
		pow_sum += in_buf[i + 4] * in_buf[i + 4];
		pow_sum += in_buf[i + 5] * in_buf[i + 5];
		pow_sum += in_buf[i + 6] * in_buf[i + 6];
		pow_sum += in_buf[i + 7] * in_buf[i + 7];
	}

	/*
	 * If input power is very low, clear pitch flag of the current frame
	 */
	if (pow_sum < POW_PITCH_THR) {
		st->pitch = (int16_t) (st->pitch & 0x3fff);
	}

	/*
	 * If input power is very low, clear complex flag of the "current" frame
	 */
	if (pow_sum < POW_COMPLEX_THR) {
		st->complex_low = (int16_t) (st->complex_low & 0x3fff);
	}

	/*
	 * Run the filter bank which calculates signal levels at each band
	 */
	filter_bank(st, in_buf, level);
	return (vad_decision(st, level, pow_sum));
}


/*
 * vad_pitch_detection
 *
 *
 * Parameters:
 *    st->oldlag        B: old LTP lag
 *    T_op              I: speech encoder open loop lags
 *    st->pitch         B: flags for pitch detection
 *    st                B: State struct
 *    st->pitch         B: flags for pitch detection
 *    st->oldlag_count  B: lag count
 *
 * Function:
 *    Test if signal contains pitch or other periodic component.
 *
 * Returns:
 *    Boolean voiced / unvoiced decision in state variable
 */

static void vad_pitch_detection(vadState * st, int32_t T_op[])
{
	int32_t lagcount, i;

	lagcount = 0;

	for (i = 0; i < 2; i++) {
		if (abs(st->oldlag - T_op[i]) < LTHRESH) {
			lagcount += 1;
		}

		/* Save the current LTP lag */
		st->oldlag = T_op[i];
	}

	/*
	 * Make pitch decision.
	 * Save flag of the pitch detection to the variable pitch.
	 */
	st->pitch = st->pitch >> 1;

	if ((st->oldlag_count + lagcount) >= NTHRESH) {
		st->pitch = st->pitch | 0x4000;
	}

	/* Update oldlagcount */
	st->oldlag_count = lagcount;
}



/*
 * cod_amr
 *
 *
 * Parameters:
 *    st          B: state structure
 *    mode        I: encoder mode
 *    new_speech  I: input speech frame, size L_FRAME
 *    st          B: State struct
 *    ana         O: Analysis parameters
 *    used_mode   B: In: -1 forces VAD on, Out:used encoder mode
 *    synth       O: local synthesis, size L_FRAME
 *
 * Function:
 *    GSM adaptive multi rate speech encoder
 *
 * Returns:
 *    void
 */
static void cod_amr(cod_amrState * st, enum Mode mode, float new_speech[],
		    int16_t ana[], enum Mode *used_mode, float synth[])
{
	/* LPC coefficients */
	float A_t[(MP1) * 4];	/* A(z) unquantized for the 4 subframes */
	float Aq_t[(MP1) * 4];	/* A(z)   quantized for the 4 subframes */
	float *A, *Aq;	/* Pointer on Aq_t */
	float lsp_new[M];

	/* Other vectors */
	float xn[L_SUBFR];	/* Target vector for pitch search */
	float xn2[L_SUBFR];	/* Target vector for codebook search */
	float code[L_SUBFR];	/* Fixed codebook excitation */
	float y1[L_SUBFR];	/* Filtered adaptive excitation */
	float y2[L_SUBFR];	/* Filtered fixed codebook excitation */
	float gCoeff[3];	/* Correlations between xn, y1, & y2: */
	float res[L_SUBFR];	/* Short term (LPC) prediction residual */
	float res2[L_SUBFR];	/* Long term (LTP) prediction residual */

	/* Vector and scalars needed for the MR475 */
	float xn_sf0[L_SUBFR];	/* Target vector for pitch search */
	float y2_sf0[L_SUBFR];	/* Filtered codebook innovation */
	float code_sf0[L_SUBFR];	/* Fixed codebook excitation */
	float h1_sf0[L_SUBFR];	/* The impulse response of sf0 */
	float mem_syn_save[M];	/* Filter memory */
	float mem_w0_save[M];	/* Filter memory */
	float mem_err_save[M];	/* Filter memory */
	float sharp_save = 0;	/* Sharpening */
	float gain_pit_sf0;	/* Quantized pitch gain for sf0 */
	float gain_code_sf0;	/* Quantized codebook gain for sf0 */
	int16_t i_subfr_sf0 = 0;	/* Position in exc[] for sf0 */

	/* Scalars & Flags */
	float gain_pit, gain_code;
	float gp_limit;	/* pitch gain limit value */
	int32_t T0_sf0 = 0;	/* Integer pitch lag of sf0 */
	int32_t T0_frac_sf0 = 0;	/* Fractional pitch lag of sf0 */
	int32_t T0, T0_frac;
	int32_t T_op[2];
	int32_t evenSubfr;
	int32_t i;
	int16_t i_subfr, subfrNr;
	int16_t lsp_flag = 0;	/* indicates resonance in LPC filter */
	//int16_t compute_sid_flag;
	int16_t vad_flag;

	memcpy(st->new_speech, new_speech, L_FRAME << 2);

	//if (st->dtx) {

		/* Find VAD decision (option 1) */
		vad_flag = vad(st->vadSt, st->new_speech);

		/* force VAD on   */
	       //	if (((int16_t) (*used_mode)) < 0)	//to fucking reference coder: enum Mode not be negative!!!
	       //		vad_flag = 1;

	       if(vad_flag) *used_mode = mode;  else *used_mode = MRDTX;
               //compute_sid_flag = 0;
		/* NB! used_mode may change here */
	       //	compute_sid_flag =
	       //	    tx_dtx_handler(vad_flag, &st->dtxEncSt->decAnaElapsedCount,
		//		   &st->dtxEncSt->dtxHangoverCount, used_mode);


	//}

        //else {
	//	compute_sid_flag = 0;
	//	*used_mode = mode;
	//}

	/*
	 * Perform LPC analysis:
	 * Autocorrelation + Lag windowing.
	 * Levinson-durbin algorithm to find a[].
	 * Convert a[] to lsp[].
	 * Quantize and code the LSPs.
	 * find the interpolated LSPs and convert to a[] for all
	 * subframes (both quantized and unquantized).
	 */
	/* LP analysis */
	lpc(st->lpcSt->LevinsonSt->old_A, st->p_window, st->p_window_12k2, A_t,
	    mode);

	/*
	 * The LP filter coefficients, are converted to
	 * the line spectral pair (LSP) representation for
	 * quantization and interpolation purposes.
	 */
	lsp(mode, *used_mode, st->lspSt->lsp_old, st->lspSt->lsp_old_q,
	    st->lspSt->qSt->past_rq, A_t, Aq_t, lsp_new, &ana);

	/* Buffer lsp's and energy */
	//dtx_buffer(&st->dtxEncSt->hist_ptr, st->dtxEncSt->lsp_hist, lsp_new,
	//	   st->new_speech, st->dtxEncSt->log_en_hist);

	if (*used_mode == MRDTX) {
	      //	dtx_enc(&st->dtxEncSt->log_en_index, st->dtxEncSt->log_en_hist,
	      //		st->dtxEncSt->lsp_hist, st->dtxEncSt->lsp_index,
	      //		&st->dtxEncSt->init_lsf_vq_index, compute_sid_flag,
	      //		&st->lspSt->qSt->past_rq[0],
	      //		st->gainQuantSt->gc_predSt->past_qua_en, &ana);
		memset(st->old_exc, 0, (PIT_MAX + L_INTERPOL) << 2);
		memset(st->mem_w0, 0, M << 2);
		memset(st->mem_err, 0, M << 2);
		memset(st->zero, 0, L_SUBFR << 2);
		memset(st->hvec, 0, L_SUBFR << 2);
		memset(st->lspSt->qSt->past_rq, 0, M << 2);
		memcpy(st->lspSt->lsp_old, lsp_new, M << 2);
		memcpy(st->lspSt->lsp_old_q, lsp_new, M << 2);

		/* Reset clLtp states */
		st->clLtpSt->pitchSt->T0_prev_subframe = 0;
		st->sharp = 0;
	} else {
		/* check resonance in the filter */
		lsp_flag = check_lsp(&st->tonStabSt->count, st->lspSt->lsp_old);
	}


	for (subfrNr = 0, i_subfr = 0; subfrNr < 2; subfrNr++, i_subfr +=
	     L_FRAME_BY2) {
		/*
		 * Pre-processing on 80 samples
		 * Find the weighted input speech for the whole speech frame
		 */
		pre_big(mode, f_gamma1, f_gamma2, A_t, i_subfr,
			st->speech, st->mem_w, st->wsp);

		/* Find open loop pitch lag for two subframes */
             /*
		if ((mode != MR475) && (mode != MR515)) {
			ol_ltp(mode, st->vadSt, &st->wsp[i_subfr],
			       &T_op[subfrNr], st->ol_gain_flg,
			       &st->pitchOLWghtSt->old_T0_med,
			       &st->pitchOLWghtSt->wght_flg,
			       &st->pitchOLWghtSt->ada_w, st->old_lags, st->dtx,
			       subfrNr);
		}
              */

	}

	//if ((mode == MR475) || (mode == MR515))
        {
		/*
		 * Find open loop pitch lag for ONE FRAME ONLY
		 * search on 160 samples
		 */
		ol_ltp(mode, st->vadSt, &st->wsp[0], &T_op[0], st->ol_gain_flg,
		       &st->pitchOLWghtSt->old_T0_med,
		       &st->pitchOLWghtSt->wght_flg, &st->pitchOLWghtSt->ada_w,
		       st->old_lags, st->dtx, 1);
		T_op[1] = T_op[0];
	}


       //	if (st->dtx) {
		vad_pitch_detection(st->vadSt, T_op);
	//}

	if (*used_mode == MRDTX) {
		goto the_end;
	}

	/*
	 * Loop for every subframe in the analysis frame
	 *
	 * To find the pitch and innovation parameters. The subframe size is
	 * L_SUBFR and the loop is repeated L_FRAME/L_SUBFR times.
	 *     - find the weighted LPC coefficients
	 *     - find the LPC residual signal res[]
	 *     - compute the target signal for pitch search
	 *     - compute impulse response of weighted synthesis filter (h1[])
	 *     - find the closed-loop pitch parameters
	 *     - encode the pitch dealy
	 *     - update the impulse response h1[] by including fixed-gain pitch
	 *     - find target vector for codebook search
	 *     - codebook search
	 *     - encode codebook address
	 *     - VQ of pitch and codebook gains
	 *     - find synthesis speech
	 *     - update states of weighting filter
	 */
	/* pointer to interpolated LPC parameters */
	A = A_t;

	/* pointer to interpolated quantized LPC parameters */
	Aq = Aq_t;
	evenSubfr = 0;
	subfrNr = -1;

	for (i_subfr = 0; i_subfr < L_FRAME; i_subfr += L_SUBFR) {
		subfrNr += 1;
		evenSubfr = 1 - evenSubfr;

		if ((evenSubfr != 0) && (*used_mode == MR475)) {
			memcpy(mem_syn_save, st->mem_syn, M << 2);
			memcpy(mem_w0_save, st->mem_w0, M << 2);
			memcpy(mem_err_save, st->mem_err, M << 2);
			sharp_save = st->sharp;
		}

		/* Preprocessing of subframe */

               /*
		if (*used_mode != MR475) {
			subframePreProc(*used_mode, f_gamma1, f_gamma2,
					A, Aq, &st->speech[i_subfr],
					st->mem_err, st->mem_w0, st->zero,
					st->ai_zero, &st->exc[i_subfr], st->h1,
					xn, res, st->error);
		}

		//MR475
		else
                 */

                {
			subframePreProc(*used_mode, f_gamma1, f_gamma2,
					A, Aq, &st->speech[i_subfr],
					st->mem_err, mem_w0_save, st->zero,
					st->ai_zero, &st->exc[i_subfr], st->h1,
					xn, res, st->error);

			if (evenSubfr != 0) {
				memcpy(h1_sf0, st->h1, L_SUBFR << 2);
			}
		}

		/* copy the LP residual (res2 is modified in the CL LTP search) */
		memcpy(res2, res, L_SUBFR << 2);

		/* Closed-loop LTP search */
		cl_ltp(&st->clLtpSt->pitchSt->T0_prev_subframe,
		       st->tonStabSt->gp, *used_mode, i_subfr, T_op, st->h1,
		       &st->exc[i_subfr], res2, xn, lsp_flag, xn2, y1, &T0,
		       &T0_frac, &gain_pit, gCoeff, &ana, &gp_limit);

		/* update LTP lag history */
		if ((subfrNr == 0) && (st->ol_gain_flg[0] > 0)) {
			st->old_lags[1] = T0;
		}

		if ((subfrNr == 3) && (st->ol_gain_flg[1] > 0)) {
			st->old_lags[0] = T0;
		}

		/* Innovative codebook search (find index and gain) */
		cbsearch(*used_mode, subfrNr, xn2, st->h1, T0, st->sharp,
			 gain_pit, code, y2, res2, &ana);

		/* Quantization of gains. */
		gainQuant(*used_mode, evenSubfr,
			  st->gainQuantSt->gc_predSt->past_qua_en,
			  st->gainQuantSt->gc_predUncSt->past_qua_en,
			  st->gainQuantSt->sf0_coeff,
			  &st->gainQuantSt->sf0_target_en,
			  &st->gainQuantSt->sf0_gcode0_exp,
			  &st->gainQuantSt->sf0_gcode0_fra,
			  &st->gainQuantSt->gain_idx_ptr, &gain_pit_sf0,
			  &gain_code_sf0, res, &st->exc[i_subfr], code, xn, xn2,
			  y1, y2, gCoeff, gp_limit, &gain_pit, &gain_code,
			  &st->gainQuantSt->adaptSt->prev_gc,
			  &st->gainQuantSt->adaptSt->onset,
			  st->gainQuantSt->adaptSt->ltpg_mem,
			  &st->gainQuantSt->adaptSt->prev_alpha, &ana);

		/* update gain history */
		for (i = 0; i < N_FRAME - 1; i++) {
			st->tonStabSt->gp[i] = st->tonStabSt->gp[i + 1];
		}
		st->tonStabSt->gp[N_FRAME - 1] = gain_pit;

		/* Subframe Post Processing */
                /*
		if (*used_mode != MR475)
                {
			subframePostProc(st->speech, i_subfr, gain_pit,
					 gain_code, Aq, synth, xn, code, y1, y2,
					 st->mem_syn, st->mem_err, st->mem_w0,
					 st->exc, &st->sharp);
		} else
                 */
                {
			if (evenSubfr != 0) {
				i_subfr_sf0 = i_subfr;
				memcpy(xn_sf0, xn, L_SUBFR << 2);
				memcpy(y2_sf0, y2, L_SUBFR << 2);
				memcpy(code_sf0, code, L_SUBFR << 2);
				T0_sf0 = T0;
				T0_frac_sf0 = T0_frac;

				/* Subframe Post Porcessing */
				subframePostProc(st->speech, i_subfr, gain_pit,
						 gain_code, Aq, synth, xn, code,
						 y1, y2, mem_syn_save,
						 st->mem_err, mem_w0_save,
						 st->exc, &st->sharp);
				st->sharp = sharp_save;
			} else {
				/*
				 * update both subframes for the MR475
				 * Restore states for the MR475 mode
				 */
				memcpy(st->mem_err, mem_err_save, M << 2);

				/* re-build excitation for sf 0 */
				Pred_lt_3or6(&st->exc[i_subfr_sf0], T0_sf0,
					     T0_frac_sf0, 1);
				Convolve(&st->exc[i_subfr_sf0], h1_sf0, y1);
				Aq -= MP1;
				subframePostProc(st->speech, i_subfr_sf0,
						 gain_pit_sf0, gain_code_sf0,
						 Aq, synth, xn_sf0, code_sf0,
						 y1, y2_sf0, st->mem_syn,
						 st->mem_err, st->mem_w0,
						 st->exc, &sharp_save);

				/* overwrites sharp_save */
				Aq += MP1;

				/*
				 * re-run pre-processing to get xn right (needed by postproc)
				 * (this also reconstructs the unsharpened h1 for sf 1)
				 */
				subframePreProc(*used_mode, f_gamma1,
						f_gamma2, A, Aq,
						&st->speech[i_subfr],
						st->mem_err, st->mem_w0,
						st->zero, st->ai_zero,
						&st->exc[i_subfr], st->h1, xn,
						res, st->error);

				/* re-build excitation sf 1 (changed if lag < L_SUBFR) */
				Pred_lt_3or6(&st->exc[i_subfr], T0, T0_frac, 1);
				Convolve(&st->exc[i_subfr], st->h1, y1);
				subframePostProc(st->speech, i_subfr, gain_pit,
						 gain_code, Aq, synth, xn, code,
						 y1, y2, st->mem_syn,
						 st->mem_err, st->mem_w0,
						 st->exc, &st->sharp);
			}
		}

		/* interpolated LPC parameters for next subframe */
		A += MP1;
		Aq += MP1;
	}
 the_end:

	/* Update signal for next frame. */
	for (i = 0; i < PIT_MAX; i++) {
		st->old_wsp[i] = st->old_wsp[L_FRAME + i];
	}

	for (i = 0; i < PIT_MAX + L_INTERPOL; i++) {
		st->old_exc[i] = st->old_exc[L_FRAME + i];
	}

	for (i = 0; i < L_TOTAL - L_FRAME; i++) {
		st->old_speech[i] = st->old_speech[L_FRAME + i];
	}
}

/*
 * Pre_Process_reset
 *
 *
 * Parameters:
 *    state                O: state structure
 *
 * Function:
 *    Initializes state memory to zero
 *
 * Returns:
 *
 */
static int32_t Pre_Process_reset(Pre_ProcessState * state)
{
	if (state == (Pre_ProcessState *) NULL) {
		//fprintf(stderr, "Pre_Process_reset: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	state->y2 = 0;
	state->y1 = 0;
	state->x0 = 0;
	state->x1 = 0;
	return 0;
}

/*
 * Pre_Process_exit
 *
 *
 * Parameters:
 *    state             I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
static void Pre_Process_exit(Pre_ProcessState ** state)
{
	if (state == NULL || *state == NULL)
		return;

	/* deallocate memory */
	free(*state);
	*state = NULL;
	return;
}

/*
 * Pre_Process_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    succeed = 0
 */
static int32_t Pre_Process_init(Pre_ProcessState ** state)
{
	Pre_ProcessState *s;

	if (state == (Pre_ProcessState * *)NULL) {
		//fprintf(stderr, "Pre_Process_init: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	*state = NULL;

	/* allocate memory */
	if ((s = calloc(1, sizeof(Pre_ProcessState))) == NULL) {
		//fprintf(stderr,
		//	"Pre_Process_init: can not malloc state structure\n");
		//return -1;
		mtrap(__LINE__);
	}
	//enc_mem+=sizeof(Pre_ProcessState);
        memadd(sizeof(Pre_ProcessState), 11);

	Pre_Process_reset(s);
	*state = s;
	return 0;
}

/*
 * Pre_Process
 *
 *
 * Parameters:
 *    y2, y1, x0, x1    B: filter memory
 *    speech            I: speech vector to be processed
 *    fspeech           O: processed vector
 *    size              I: size of the vector
 *
 * Function:
 *    Pre-processing
 *
 *    Two pre-processing functions are applied prior to
 *    the encoding process: high-pass filtering and signal down-scaling.
 *    Down-scaling consists of dividing the input by a factor of 2
 *    to reduce the possibility of overflows in the fixed-point
 *    implementation. The high-pass filter serves as a precaution
 *    against undesired low frequency components. A filter with
 *    a cut off frequency of 80 Hz is used, and it is given by:
 *
 *            0.927246093 - 1.8544941z^-1 + 0.927246903z^-2
 *    H(z) = -----------------------------------------------
 *                1 - 1.906005859z^-1 + 0.911376953z^-2
 *
 *    Down-scaling and high-pass filtering are combined by dividing
 *    the coefficients at the numerator by 2.
 *
 * Returns:
 *    void
 */
static void Pre_Process(float * y2, float * y1, float * x0, float
			* x1, int16_t * speech, float * f_speech)
{
	int32_t i;
	float x2;
	float tmp;

	for (i = 0; i < 160; i++) {
		x2 = *x1;
		*x1 = *x0;
		//ssss=speech[i];
		*x0 = speech[i];   //*x0 = speech[i];
		tmp =
		    (float) (0.4636230465F * *x0 - 0.92724705F * *x1 +
			       0.4636234515F * x2 + 1.906005859F * *y1 -
			       0.911376953F * *y2);
		f_speech[i] = tmp;
		*y2 = *y1;
		*y1 = tmp;
	}

	if ((fabsf(*y1) + fabsf(*y2)) < 0.0000000001F)
		*y2 = *y1 = 0;
}

/*
 * cod_amr_reset
 *
 *
 * Parameters:
 *    s                 B: state structure
 *    dtx               I: dtx on/off
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *    void
 */
static void cod_amr_reset(cod_amrState * s, int32_t dtx)
{
	int32_t i;

	/* reset DTX */
	s->dtx = dtx;

	/* reset Pitch_frState */
	s->clLtpSt->pitchSt->T0_prev_subframe = 0;

	/* reset Q_plsfState */
	memset(s->lspSt->qSt->past_rq, 0, M * sizeof(float));
	memcpy(s->lspSt->lsp_old, f_lsp_init_data, sizeof(f_lsp_init_data));
	memcpy(s->lspSt->lsp_old_q, f_lsp_init_data, sizeof(f_lsp_init_data));

	/* reset gc_predState */
	for (i = 0; i < NPRED; i++) {
		s->gainQuantSt->gc_predSt->past_qua_en[i] =
		    NB_QUA_CODE + VQ_SIZE_HIGHRATES + VQ_SIZE_LOWRATES +
		    MR475_VQ_SIZE * 2 + DTX_VQ_SIZE;
		s->gainQuantSt->gc_predUncSt->past_qua_en[i] =
		    NB_QUA_CODE + VQ_SIZE_HIGHRATES + VQ_SIZE_LOWRATES +
		    MR475_VQ_SIZE * 2 + DTX_VQ_SIZE;
	}

	/* reset gain_adaptState */
	s->gainQuantSt->adaptSt->onset = 0;
	s->gainQuantSt->adaptSt->prev_alpha = 0.0F;
	s->gainQuantSt->adaptSt->prev_gc = 0.0F;
	memset(s->gainQuantSt->adaptSt->ltpg_mem, 0,
		LTPG_MEM_SIZE * sizeof(float));
	s->gainQuantSt->sf0_gcode0_exp = 0;
	s->gainQuantSt->sf0_gcode0_fra = 0;
	s->gainQuantSt->sf0_target_en = 0.0F;
	memset(s->gainQuantSt->sf0_coeff, 0, 5 * sizeof(float));
	s->gainQuantSt->gain_idx_ptr = NULL;

	/* reset pitchOLWghtState */
	s->pitchOLWghtSt->old_T0_med = 40;
	s->pitchOLWghtSt->ada_w = 0.0F;
	s->pitchOLWghtSt->wght_flg = 0;

	/* reset tonStabState */
	s->tonStabSt->count = 0;
	memset(s->tonStabSt->gp, 0, N_FRAME * sizeof(float));

	/* reset LevinsonState */
	s->lpcSt->LevinsonSt->old_A[0] = 1.0F;
	memset(&s->lpcSt->LevinsonSt->old_A[1], 0, M * sizeof(float));


	/* reset vadState */
	s->vadSt->oldlag_count = 0;
	s->vadSt->oldlag = 0;
	s->vadSt->pitch = 0;
	s->vadSt->tone = 0;
	s->vadSt->complex_high = 0;
	s->vadSt->complex_low = 0;
	s->vadSt->complex_hang_timer = 0;
	s->vadSt->vadreg = 0;
	s->vadSt->burst_count = 0;
	s->vadSt->hang_count = 0;
	s->vadSt->complex_hang_count = 0;

	/* initialize memory used by the filter bank */
	for (i = 0; i < 3; i++) {
		s->vadSt->a_data5[i][0] = 0;
		s->vadSt->a_data5[i][1] = 0;
	}

	for (i = 0; i < 5; i++) {
		s->vadSt->a_data3[i] = 0;
	}

	/* reset dtx_encState */
	/* initialize the rest of the memory */
	for (i = 0; i < COMPLEN; i++) {
		s->vadSt->bckr_est[i] = NOISE_INIT;
		s->vadSt->old_level[i] = NOISE_INIT;
		s->vadSt->ave_level[i] = NOISE_INIT;
		s->vadSt->sub_level[i] = 0;
	}
	s->vadSt->best_corr_hp = CVAD_LOWPOW_RESET;
	s->vadSt->speech_vad_decision = 0;
	s->vadSt->complex_warning = 0;
	s->vadSt->sp_burst_count = 0;
	s->vadSt->corr_hp_fast = CVAD_LOWPOW_RESET;

       /*
	s->dtxEncSt->hist_ptr = 0;
	s->dtxEncSt->log_en_index = 0;
	s->dtxEncSt->init_lsf_vq_index = 0;
	s->dtxEncSt->lsp_index[0] = 0;
	s->dtxEncSt->lsp_index[1] = 0;
	s->dtxEncSt->lsp_index[2] = 0;

	for (i = 0; i < DTX_HIST_SIZE; i++) {
		memcpy(&s->dtxEncSt->lsp_hist[i * M], f_lsp_init_data,
		       sizeof(float) * M);
	}
	memset(s->dtxEncSt->log_en_hist, 0, M * sizeof(float));
	s->dtxEncSt->dtxHangoverCount = DTX_HANG_CONST;
	s->dtxEncSt->decAnaElapsedCount = DTX_ELAPSED_FRAMES_THRESH;
        */
	/* init speech pointers */
	/* New speech */
	s->new_speech = s->old_speech + L_TOTAL - L_FRAME;

	/* Present frame */
	s->speech = s->new_speech - L_NEXT;
	s->p_window = s->old_speech + L_TOTAL - L_WINDOW;

	/* For LPC window                            */
	s->p_window_12k2 = s->p_window - L_NEXT;

	/* Initialize static pointers */
	s->wsp = s->old_wsp + PIT_MAX;
	s->exc = s->old_exc + PIT_MAX + L_INTERPOL;
	s->zero = s->ai_zero + MP1;
	s->error = s->mem_err + M;
	s->h1 = &s->hvec[L_SUBFR];

	/* Static vectors to zero */
	memset(s->old_speech, 0, L_TOTAL * sizeof(float));
	memset(s->old_exc, 0, (PIT_MAX + L_INTERPOL) * sizeof(float));
	memset(s->old_wsp, 0, PIT_MAX * sizeof(float));
	memset(s->mem_syn, 0, M * sizeof(float));
	memset(s->mem_w, 0, M * sizeof(float));
	memset(s->mem_w0, 0, M * sizeof(float));
	memset(s->mem_err, 0, M * sizeof(float));
	memset(s->ai_zero, 0, L_SUBFR * sizeof(float));
	memset(s->hvec, 0, L_SUBFR * sizeof(float));

	for (i = 0; i < 5; i++) {
		s->old_lags[i] = 40;
	}
	s->sharp = 0.0F;
}

/*
 * cod_amr_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *    dtx               I: dtx mode used
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    succeed = 0
 */
static int32_t cod_amr_init(cod_amrState ** state, int32_t dtx)
{
	cod_amrState *s;

	if ((s = calloc(1, sizeof(cod_amrState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
		mtrap(__LINE__);
	}
  memadd(sizeof(cod_amrState), 12);
  //enc_mem+=sizeof(cod_amrState);
	/* init clLtpState */
	if ((s->clLtpSt = calloc(1, sizeof(clLtpState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
		mtrap(__LINE__);
	}
  //enc_mem+=sizeof(clLtpState);
  memadd(sizeof(clLtpState), 13);
	/* init Pitch_frState */
	if ((s->clLtpSt->pitchSt =
	     calloc(1, sizeof(Pitch_frState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}
  //enc_mem+=sizeof(Pitch_frState);
  memadd(sizeof(Pitch_frState), 14);

	/* init lspState */
	if ((s->lspSt = calloc(1, sizeof(lspState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
		mtrap(__LINE__);
	}
  //enc_mem+=sizeof(lspState);
  memadd(sizeof(lspState), 15);

	/* init Q_plsfState */
	if ((s->lspSt->qSt = calloc(1, sizeof(Q_plsfState))) ==
	    NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				mtrap(__LINE__);
	}
  //enc_mem+=sizeof(Q_plsfState);
  memadd(sizeof(Q_plsfState), 16);
	/* init gainQuantState */
	if ((s->gainQuantSt =
	     calloc(1, sizeof(gainQuantState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}

  //enc_mem+=sizeof(gainQuantState);
  memadd(sizeof(gainQuantState), 17);

	/* init gc_predState x2 */
	if ((s->gainQuantSt->gc_predSt =
	     calloc(1, sizeof(gc_predState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}

  //enc_mem+=sizeof(gc_predState);
  memadd(sizeof(gc_predState), 18);

	if ((s->gainQuantSt->gc_predUncSt =
	     calloc(1, sizeof(gc_predState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}
  //enc_mem+=sizeof(gc_predState);
  memadd(sizeof(gc_predState), 19);

	/* init gain_adaptState */
	if ((s->gainQuantSt->adaptSt =
	     calloc(1, sizeof(gain_adaptState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}

  //enc_mem+=sizeof(gain_adaptState);
  memadd(sizeof(gain_adaptState), 20);

	/* init pitchOLWghtState */
	if ((s->pitchOLWghtSt =
	     calloc(1, sizeof(pitchOLWghtState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}
  //enc_mem+=sizeof(pitchOLWghtState);
  memadd(sizeof(pitchOLWghtState), 21);

	/* init tonStabState */
	if ((s->tonStabSt = calloc(1, sizeof(tonStabState)))
	    == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				mtrap(__LINE__);
	}
  //enc_mem+=sizeof(tonStabState);
   memadd(sizeof(tonStabState), 22);

	/* init lpcState */
	if ((s->lpcSt = calloc(1, sizeof(lpcState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
		mtrap(__LINE__);
	}
  //enc_mem+=sizeof(lpcState);
  memadd(sizeof(lpcState), 23);

	/* init LevinsonState */
	if ((s->lpcSt->LevinsonSt =
	     calloc(1, sizeof(LevinsonState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
				 mtrap(__LINE__);
	}
  //enc_mem+=sizeof(LevinsonState);
  memadd(sizeof(LevinsonState), 24);

	if ((s->vadSt = calloc(1, sizeof(vadState))) == NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;
		mtrap(__LINE__);
	}
  //enc_mem+=sizeof(vadState);
  memadd(sizeof(vadState), 25);

	/* Init dtx_encState */
	//if ((s->dtxEncSt = calloc(1, sizeof(dtx_encState))) ==
	//    NULL) {
		//fprintf(stderr, "can not malloc state structure\n");
		//goto lfree;]
		    //		mtrap(__LINE__);
	//}
	//enc_mem+=sizeof(dtx_encState);


	cod_amr_reset(s, dtx);
	*state = s;
	return 0;
/*
 lfree:
	if (s) {
		//if (s->dtxEncSt)
		 //	free(s->dtxEncSt);
		if (s->vadSt)
			free(s->vadSt);
		if (s->lpcSt) {
			if (s->lpcSt->LevinsonSt)
				free(s->lpcSt->LevinsonSt);
			free(s->lpcSt);
		}
		if (s->tonStabSt)
			free(s->tonStabSt);
		if (s->pitchOLWghtSt)
			free(s->pitchOLWghtSt);
		if (s->gainQuantSt) {
			if (s->gainQuantSt->gc_predUncSt)
				free(s->gainQuantSt->gc_predUncSt);
			if (s->gainQuantSt->gc_predSt)
				free(s->gainQuantSt->gc_predSt);
			free(s->gainQuantSt);
		}
		if (s->lspSt) {
			if (s->lspSt->qSt)
				free(s->lspSt->qSt);
			free(s->lspSt);
		}
		if (s->clLtpSt) {
			if (s->clLtpSt->pitchSt)
				free(s->clLtpSt->pitchSt);
			free(s->clLtpSt);
		}
		free(s);
	}
	
	*/
	//return -1;
}

/*
 * cod_amr_exit
 *
 *
 * Parameters:
 *    state             I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
static void cod_amr_exit(cod_amrState ** state)
{
	if (state == NULL || *state == NULL)
		return;

	/* deallocate memory */
	free((*state)->vadSt);
	free((*state)->gainQuantSt->gc_predSt);
	free((*state)->gainQuantSt->gc_predUncSt);
	free((*state)->gainQuantSt->adaptSt);
	free((*state)->clLtpSt->pitchSt);
	free((*state)->lspSt->qSt);
	free((*state)->lpcSt->LevinsonSt);
	free((*state)->lpcSt);
	free((*state)->lspSt);
	free((*state)->clLtpSt);
	free((*state)->gainQuantSt);
	free((*state)->pitchOLWghtSt);
	free((*state)->tonStabSt);
	//free((*state)->dtxEncSt);
	free(*state);
	*state = NULL;
	return;
}

/*
 * Speech_Encode_Frame_init
 *
 *
 * Parameters:
 *    state             O: state structure
 *    dtx               I: dtx mode used
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    succeed = 0
 */
void *Speech_Encode_Frame_init(int dtx)
{
	Speech_Encode_FrameState *s;

	/* allocate memory */
	if ((s =
	     calloc(1, sizeof(Speech_Encode_FrameState))) == NULL) {
		//fprintf(stderr,
		//	"Speech_Encode_Frame_init: can not malloc state "
		//	"structure\n");
		//return NULL;
				 mtrap(__LINE__);
	}
	//enc_mem+=sizeof(Speech_Encode_FrameState);
        memadd(sizeof(Speech_Encode_FrameState), 26);

	s->pre_state = NULL;
	s->cod_amr_state = NULL;
	s->dtx = dtx;

	if (Pre_Process_init(&s->pre_state) || cod_amr_init(&s->cod_amr_state,
							    dtx)) {
		Speech_Encode_Frame_exit((void **)(&s));
		return NULL;
	}
	return s;
}

/*
 * Speech_Encode_Frame_reset
 *
 *
 * Parameters:
 *    state          O: state structure
 *
 * Function:
 *    Resets state memory
 *
 * Returns:
 *
 */
int Speech_Encode_Frame_reset(void *st, int dtx)
{
	Speech_Encode_FrameState *state;
	state = (Speech_Encode_FrameState *) st;

	if ((Speech_Encode_FrameState *) state == NULL) {
		//fprintf(stderr,
		//	"Speech_Encode_Frame_reset: invalid parameter\n");
		//return -1;
		mtrap(__LINE__);
	}
	Pre_Process_reset(state->pre_state);
	cod_amr_reset(state->cod_amr_state, dtx);
	return 0;
}

/*
 * Speech_Encode_Frame_exit
 *
 *
 * Parameters:
 *    state            I: state structure
 *
 * Function:
 *    The memory used for state memory is freed
 *
 * Returns:
 *    Void
 */
void Speech_Encode_Frame_exit(void **st)
{
	if ((Speech_Encode_FrameState *) (*st) == NULL)
		return;
	Pre_Process_exit(&(((Speech_Encode_FrameState *) (*st))->pre_state));
	cod_amr_exit(&(((Speech_Encode_FrameState *) (*st))->cod_amr_state));

	/* deallocate memory */
	free(*st);
	*st = NULL;
	return;
}

/*
 * Speech_Encode_Frame
 *
 *
 * Parameters:
 *    st                B: state structure
 *    mode              I: speech coder mode
 *    new_speech        I: speech input, size L_FRAME
 *    prm               O: Analysis parameters
 *    used_mode         B: force VAD/used_mode
 * Function:
 *    Encode one frame
 *
 * Returns:
 *    Void
 */
void Speech_Encode_Frame(void *st, enum Mode mode, int16_t * new_speech, int16_t *
			 prm, enum Mode *used_mode)
{
	float syn[L_FRAME];	/* Buffer for synthesis speech */
	float speech[160];
	//int32_t i;

	Speech_Encode_FrameState *state;
	state = (Speech_Encode_FrameState *) st;

	//for (i = 0; i < 160; i++) {
	//	new_speech[i] = (int16_t) (new_speech[i] & 0xfff8);
	//}

	/* filter + downscaling */
	Pre_Process(&state->pre_state->y2, &state->pre_state->y1,
		    &state->pre_state->x0, &state->pre_state->x1, new_speech,
		    speech);

	/* Call the speech encoder */
	cod_amr(state->cod_amr_state, mode, speech, prm, used_mode, syn);

}
