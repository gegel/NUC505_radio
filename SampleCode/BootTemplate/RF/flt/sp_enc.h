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
 * sp_enc.h
 *
 *
 * Project:
 *    AMR Floating-Point Codec
 *
 * Contains:
 *    Defines interface to AMR encoder
 *
 */

#pragma once

#ifndef _SP_ENC_H_
#define _SP_ENC_H_

/*
 * include files
 */
#include <stdint.h>
#include "modes.h"

/*
 * Function prototypes
 */

/*
 * initialize one instance of the speech encoder
 * Stores pointer to filter status struct in *st. This pointer has to
 * be passed to Speech_Encode_Frame in each call.
 * returns 0 on success
 */
void *Speech_Encode_Frame_init(int dtx);
/*
 * reset speech encoder (i.e. set state memory to zero)
 * returns 0 on success
 */
int Speech_Encode_Frame_reset(void *st, int dtx);

/*
 * de-initialize speech encoder (i.e. free status struct)
 * stores NULL in *st
 */
void Speech_Encode_Frame_exit(void **st);

/*
 * Encodes one speech frame
 * Returns analysis parameters
 */
void Speech_Encode_Frame(void *st, enum Mode mode, short *newSpeech,
			 short *prm, enum Mode *usedMode);

#endif /* _SP_ENC_H_ */
