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
 * interf_enc.c
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

/*
 * include files
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "platform.h"
//#include <memory.h>
//#include "amr_speech_importance.h"
#include "sp_enc.h"
#include "interf_rom.h"
//#include <ophtools.h>
int eif_mem=0;

void mtrap(int line);
void memadd(int len, unsigned char id);

/*
 * Declare structure types
 */
/* Declaration transmitted frame types */
enum TXFrameType { TX_SPEECH_GOOD = 0,
	TX_SID_FIRST,
	TX_SID_UPDATE,
	TX_NO_DATA,
	TX_SPEECH_DEGRADED,
	TX_SPEECH_BAD,
	TX_SID_BAD,
	TX_ONSET,
	TX_N_FRAMETYPES		/* number of frame types */
};

/* Declaration of interface structure */
typedef struct {
	int16_t sid_update_counter;	/* Number of frames since last SID */
	int16_t sid_handover_debt;	/* Number of extra SID_UPD frames to schedule */
	int32_t dtx;
	enum TXFrameType prev_ft;	/* Type of the previous frame */
	void *encoderState;	/* Points encoder state structure */
} enc_interface_State;

/*
 * Sid_Sync_reset
 *
 *
 * Parameters:
 *    st                O: state structure
 *
 * Function:
 *    Initializes state memory
 *
 * Returns:
 *    void
 */
static void Sid_Sync_reset(enc_interface_State * st)
{
	st->sid_update_counter = 3;
	st->sid_handover_debt = 0;
	st->prev_ft = TX_SPEECH_GOOD;
}

/*
 * Encoder_Interface_init
 *
 *
 * Parameters:
 *    dtx               I: DTX flag
 *
 * Function:
 *    Allocates state memory and initializes state memory
 *
 * Returns:
 *    pointer to encoder interface structure
 */
void *Encoder_Interface_init(int dtx)
{
	enc_interface_State *s;

	/* allocate memory */
	if ((s = calloc(1, sizeof(enc_interface_State))) ==
	    NULL) {
		mtrap(__LINE__);
	}
        //eif_mem+=sizeof(enc_interface_State);
        memadd(sizeof(enc_interface_State), 10);


	s->encoderState = Speech_Encode_Frame_init(dtx);
	Sid_Sync_reset(s);
	s->dtx = dtx;
	return s;
}

/*
 * DecoderInterfaceExit
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
void Encoder_Interface_exit(void *state)
{
	enc_interface_State *s;
	s = (enc_interface_State *) state;

	/* free memory */
	Speech_Encode_Frame_exit(&s->encoderState);
	free(s);
	s = NULL;
}


short AMR475_encode(void *st, int16_t * speech, uint8_t * serial)
{
	int16_t prm[PRMNO_MR475]; //parameters outputed by encoder
	int i, j;
	int16_t bits;  //output bits counter
        enum Mode used_mode = (enum Mode)0; //result of encoding (VAD output)
        enc_interface_State *s = (enc_interface_State *) st;;	//pointer to state

        //encode speech frame
        Speech_Encode_Frame(s->encoderState, (enum Mode)0, speech, prm, &used_mode);
        memset(serial, 0, 12); //clear output


        //check silency
        if(used_mode)
        {
         serial[0]=1; //set silency flag
         return 0;
        }

         //convert parameters to bits
        bits = 1;	//reserve first bit for no-SID flag=0 (bit 0)
        for (i = 0; i < PRMNO_MR475; i++)	//pack speech to bits
        {
         for (j = 0; j < bitno_MR475[i]; j++)
         {
          if (prm[i] & ((unsigned short) 1 << j)) serial[bits >> 3] |=((unsigned char) 1 << (bits & 0x07));
          bits++;	//updates bits counter
         }
        }

        return 12;

}

//chack for malloc fail
void mtrap(int line) //malloc failure trap
{
  
}

void memadd(int len, unsigned char id)
{
 
}


