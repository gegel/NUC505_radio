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
 * interf_enc.h
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

#ifndef _interf_enc_h_
#define _interf_enc_h_

/*
 * include files
 */
//#include"sp_enc.h"
#include"modes.h"
/*
 * Function prototypes
 */
/*
 * Reserve and init. memory
 */
void *Encoder_Interface_init(int dtx);

/*
 * Exit and free memory
 */
void Encoder_Interface_exit(void *state);


short AMR475_encode(void *st, short * speech, unsigned char * serial);

#endif /* _interf_enc_h_ */

