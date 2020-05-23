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
 * interf_dec.c
 *
 *
 * Project:
 *     AMR Floating-Point Codec
 *
 * Contains:
 *    This module provides means to conversion from 3GPP or ETSI
 *    bitstream to AMR parameters
 */

/*
 * include files
 */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "platform.h"

//#include <memory.h>
//#include "amr_speech_importance.h"
#include "sp_dec.h"
#include "interf_rom.h"
#include "rom_dec.h"
//#include <ophtools.h>

int dif_mem=0;

void mtrap(int line);
void memadd(int len, unsigned char id);

/*
 * definition of constants
 */
#define EHF_MASK 0x0008		/* encoder homing frame pattern */
typedef

    struct {
	int reset_flag_old;	/* previous was homing frame */

	enum RXFrameType prev_ft;	/* previous frame type */
	enum Mode prev_mode;	/* previous mode */
	void *decoder_State;	/* Points decoder state */

} dec_interface_State;

/*
 * Decoder_Interface_reset
 *
 *
 * Parameters:
 *    st                O: state struct
 *
 * Function:
 *    Reset homing frame counter
 *
 * Returns:
 *    void
 */
void Decoder_Interface_reset(dec_interface_State * st)
{
	st->reset_flag_old = 1;
	st->prev_ft = RX_SPEECH_GOOD;
	st->prev_mode = MR475;	/* minimum bitrate */
}

/*
 * Decoder_Interface_init
 *
 *
 * Parameters:
 *    void
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    success           : pointer to structure
 *    failure           : NULL
 */
void *Decoder_Interface_init(void)
{
	dec_interface_State *s;

	/* allocate memory */
	if ((s = calloc(1, sizeof(dec_interface_State))) ==
	    NULL) {
		mtrap(__LINE__);
	}

	//dif_mem+=sizeof(dec_interface_State);
        memadd(sizeof(dec_interface_State), 1);

	s->decoder_State = Speech_Decode_Frame_init();

	if (s->decoder_State == NULL) {
		free(s);
		return NULL;
	}
	Decoder_Interface_reset(s);
	return s;
}

/*
 * Decoder_Interface_exit
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
void Decoder_Interface_exit(void *state)
{
	dec_interface_State *s;
	s = (dec_interface_State *) state;

	/* free memory */
	Speech_Decode_Frame_exit(s->decoder_State);
	free(s);
	s = NULL;
}



short AMR475_decode(void *st, uint8_t * serial, int16_t * synth)
{
       int32_t i, j;
       int16_t prm[PRMNO_MR475];
       int16_t bits;		//bits counter
       dec_interface_State *s;


        if(serial[0]&0x01)
        {
         memset(synth, 0, L_FRAME*sizeof(short));
         return 0;
        }


	s = (dec_interface_State *) st;	//state

        //unpack bits ro parameters
	memset(prm, 0, sizeof(prm));


		bits = 1;	//skip first bit (SID flag)
		for (i = 0; i < PRMNO_MR475; i++) {
			for (j = 0; j < bitno_MR475[i]; j++) {
				if (serial[bits >> 3] &
				    (((unsigned char) 1) << (bits & 7)))
					prm[i] |= ((unsigned short) 1 << j);
				bits++;
			}
		}

        Speech_Decode_Frame(s->decoder_State, (enum Mode)0, prm, RX_SPEECH_GOOD, synth);
        return 160;
}

