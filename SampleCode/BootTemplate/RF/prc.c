/*
Data processing procedures with RX and TX 3-level FIFOs for 12 bytes AMR data packets 
Thre are hight level for:
- modulation next packet from TX fifo to data bits will be output over SPi to radio
- demodulation bits inputted from radio over SPI to RX fifo
- encoding PCM frame recorded by audio engine to TX fifo
- decoding from RX fifo to PCM frame will be played with audio engine

Crystal clock of remote transmitter can be a little difference with local clock
Local receiver synchronized with remote transmitter by change a little receiving SPI clock
So local player must be also synchronized with remote clock in lock loop
There are two ways of this sync:
- open way: received frames mark as voice or silency by remote voice active detector,
  this fremes are queued in RX fifo, so we can add/skip silency frames in fifo near underrun/overrun
- close way: if remote speech haven't silency too long we can a little resample PCM
  before output for playing for prevent queue underrun/overrun	 
Both ways are realized in this software module

Note: all received frames marked as silency replased by exact silency, comfort noise
generatiot not uses and remove from AMR codec because code size are cryticat for our case
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "flt/interf_enc.h"
#include "flt/interf_dec.h"
#include "mdm/mdm.h"

#include "prc.h"

#define VTEST  //include voice test table with AMR 4750 coded male speech count from 1 to 10 in Russian 

#ifdef VTEST
#include "voice.h"
#endif

//----------internal procedures prototypes---------------
short prc_resample(short* src, short* dst, short r); //playing resampler

//----------global values--------------------------

//RX fifo
unsigned char prc_rbuf[FIFOCNT][PKTLENGTH]; //speech blocks fifo for receiving/playing
unsigned char prc_in_rptr=2; //receiving input pointer
unsigned char prc_out_rptr=0; //playing output pointer

//TX fifo
unsigned char prc_tbuf[FIFOCNT][PKTLENGTH]; //speech blocks fifo for grabbing/sending
unsigned char prc_in_tptr=2; //grabbing input pointer
unsigned char prc_out_tptr=0; //sending output pointer

//TX resampling buffer
short prc_pbuf[2*FRAMELENGTH]; //player's resampling buffer
short prc_pptr=0; //number of samples ready for playing

//other values
unsigned char prc_syn=0;  //demodulator's sync flag
unsigned char prc_pdif=2; //internal diff between pointers fixed at a moment
short prc_pleft=0;  //resampler's last sample
unsigned short prc_ppos=PRC_ONE; //resamplers current position

//AMR codec states
int *prc_amrenc=0; //AMR encode state
int *prc_amrdec=0; //AMR decode state

#ifdef VTEST
#define VPAUSE 100 //x20 mS
unsigned short prc_vtcnt=0;
unsigned char prc_vton=0;


//enable/disable voice test count
void prc_vt(unsigned char on)
{
	prc_vton=on;
}
#else
void prc_vt(unsigned char on)
{

}
#endif	


//------------------------------------------initializing--------------------------------------

//initialize processing
void prc_init(void)
{
 #ifdef VTEST
  prc_vtcnt=0;
  prc_vton=0;	
 #endif
	
 //clear pointers
 prc_in_rptr=2;  //input rx pointer must be diff by 2 from output pointer
 prc_out_rptr=0;
 prc_in_tptr=2; //input tx pointer must be diff by 2 from output pointer
 prc_out_tptr=0;	
 prc_pptr=0;
 prc_pdif=2; //diff must be 2
 memset(prc_rbuf, 0, sizeof(prc_rbuf)); //clear buffers
 memset(prc_tbuf, 0, sizeof(prc_tbuf));
 memset(prc_pbuf, 0, sizeof(prc_pbuf));
 mdm_init(); //initialize modem
}

//initialize amr codec
void amr_init(void)
{
 prc_amrenc = Encoder_Interface_init(1); //create AMR encoder
 prc_amrdec = Decoder_Interface_init();    //create AMR decoder		
}


//-----------------------------grabbing/sendigng for radio transmitter----------------------------------

//encode grabbed audio frame 
short prc_enc(short* pcm)
{
 short i;
 i=AMR475_encode(prc_amrenc, (short*)pcm, prc_tbuf[prc_in_tptr]);   //encode pcm to tx fifo

#ifdef VTEST
 if(prc_vton)  //check voice test active
 {
  if(prc_vtcnt<AMR_SPEECH) //replace recorded frame by speech count 
	{
		memcpy(prc_tbuf[prc_in_tptr], speech_tbl[prc_vtcnt], 12); //use table
		i=1; //set voice flag
	}
	else //or pause between conts 1-9
	{
		memset(prc_tbuf[prc_in_tptr], 0, 12); //output silency
		prc_tbuf[prc_in_tptr][0]=1; //clear voive flag in data block
		i=0; //clear voice flag
	}
	prc_vtcnt++; //move pointer to next table block
	if(prc_vtcnt>=(AMR_SPEECH+VPAUSE)) prc_vtcnt=0; //ring with table size + pause size
 }	 
#endif	
		

 prc_in_tptr++; //move TX fifo input pointer
 if(prc_in_tptr==FIFOCNT) prc_in_tptr=0; //ring tx fifo input pointer	
 return i; //returns voice flag	
}


//modulate bits for SPI transmitter
void prc_mdm(unsigned short* bits)
{
 mdm_modulate(prc_tbuf[prc_out_tptr++], (unsigned short*)bits);	//modulate bits from tx fifo
 if(prc_out_tptr==FIFOCNT) prc_out_tptr=0;	//ring tx fifo output pointer
}


//-----------------------------------receiving/playing for radio receiver--------------------------

//demodulate PCI recevied bits
unsigned short prc_dmd(unsigned short* bits)
{
	unsigned short n; //demodulator's sync flaf and rate correction value
	unsigned char s, v; //demodulator's sync and voice flags
	
	
	//demodulate bits to rx fifo, return sync
	n=(short)mdm_demodulate((unsigned short*)bits, prc_rbuf[prc_in_rptr]); 
	s=n>>15;   //sync flag (bit15, carrier detected is 1)
	v=1&prc_rbuf[prc_in_rptr][0]; //voice flag of received frame
  
	//move input pointers for received frames
	if( (!v) || (prc_pdif<=2) ) //check frame not voice or no overrun (no 3)
  {
	 if(s) prc_in_rptr++; //move rx fifo pointer (only on sync ok)
	 if(prc_in_rptr==FIFOCNT) prc_in_rptr=0; //ring rx fifo input pointer
	}
	
	//check modem's sync flag change
	if(s!=prc_syn)
	{
		if(prc_syn) //sync lost  
		{ //reset play fifo
		 memset(prc_rbuf, 0, sizeof(prc_rbuf)); 
		 memset(prc_pbuf, 0, sizeof(prc_pbuf));
		 prc_in_rptr=1; //set input rx pointer must be diff from output pointer
     prc_out_rptr=0; //clear output rx pointer
		 prc_pdif=1;	//set pointer's difference
		 prc_pleft=0;  //clear resampler's last sample
     prc_ppos=PRC_ONE; //clear resamplers current position
    }
		prc_syn=s; //set actual sync flag
	}
	
	prc_set(); //capure pointer's difference at this time point
	
	return n; //returns demodulator's sync flag and rate correction value for SPI receiver
}

//decode frame for playing with resampling
void prc_decr(short* pcm)
{
 short i, r;
 unsigned char v; //voice flag	
	
 r=prc_pdif; //get pointer's difference of rx fifo
 r-=2; //compute correction value (-1 to 1), must lock to 0
 r*=PRC_RC; //apply rate correction coefficient (not more then 64! equal (+- 0.4%))	
 
	//check modem sync flag
 if(!prc_syn)
 {
  for(i=0;i<FRAMELENGTH;i++) pcm[i]=0; //output silency on no sync
	return; // goto play_end; 
 }	 

//check we haven't samples for output frame 
 while(prc_pptr<FRAMELENGTH) 
 {	
	 //no samples yet: decode new 
	 v=AMR475_decode(prc_amrdec, prc_rbuf[prc_out_rptr], prc_pbuf+FRAMELENGTH); //decode rx fifo to buffer
 
	 //check this is voice sample or not fifo underrun
	 if( v  || (prc_pdif>=2))  
	 {
		prc_out_rptr++; //move pointer, skip for unvoiced on underrun
	  if(prc_out_rptr==FIFOCNT) prc_out_rptr=0;	//ring rx fifo output pointer
	 }
	
	 //check this is unvoice sample
	if(!v)
	{ //reset resampler
		prc_pleft=0;  //clear resampler's last sample
    prc_ppos=PRC_ONE; //clear resamplers current position
	}
	
	//resample decoded samples
	i=prc_resample(prc_pbuf+FRAMELENGTH, prc_pbuf+prc_pptr, r); //returns number of output samples only on range 159-161   
	prc_pptr+=i; //add actual number of new samples were add
	
 } //end of while

 //output frame of pcm
 memcpy(pcm, prc_pbuf, FRAMELENGTH*sizeof(short)); //output frame from resampling buffer
 
 prc_pptr-=FRAMELENGTH; //number of samples in tail
 if(prc_pptr) memcpy(prc_pbuf, prc_pbuf+FRAMELENGTH, prc_pptr*sizeof(short)); //move tail to start of buffer
 
}

//decode frame for playing without resampling
void prc_dec(short* pcm)
{
	AMR475_decode(prc_amrdec, prc_rbuf[prc_out_rptr], pcm); //decode rx fifo to output pcm
  prc_out_rptr++;
	if(prc_out_rptr==FIFOCNT) prc_out_rptr=0;	//ring rx fifo output pointer
}

//------------------------------------------rate correction for playing-----------------------------

//perform rate correction: this procedure threaded safe and must be called on some fixing time point
//for example from audio playing interupt
void prc_set(void)
{
  short r; //rate correction value
	
	r=prc_in_rptr-prc_out_rptr; //get difference between in and out pointers in receivin/playing fifo (-2 to 2)
  if(r<=0) r+=3; //match to 1 to 3   
  prc_pdif=r; //set pointer's difference value at this moment
}


//playing resampler: resample pcm from src to dst with rate correction r (+-64)
//must returns maximum 161 samples! so r cann't be less of 64!!!
short prc_resample(short* src, short* dst, short r)
{
 short i;  //samples counter
 short pcur; //currently processed sample	
 int diff;  //difference between current and last samples
 short* out = dst; //pointer to carrently output sample


 r+=PRC_ONE;  //fit signed rate correction to Q15
 for(i=0;i<FRAMELENGTH;i++) //process all samples of frame
 {
  pcur= *src; //processed sample
	diff = pcur - prc_pleft; //difference between processed and last sample
	while(prc_ppos<=PRC_ONE) //step position in loop up to next sample
	{		
	 *(out++) = prc_pleft + (diff>>31)*((abs(diff)*prc_ppos)>>14); //output sample depends difference and position
   prc_ppos+=r; //move position depend rate value
  }
  src++; //move to next sample will be processed
	prc_pleft=pcur;	//set currently processed sample as last samle
	prc_ppos-=PRC_ONE; //ring position to ONE (exactly input rate)
 }	
	
 return(short)(out-dst); //returns number of samples were outputted		
}

//get correction value for external hardware correction of playing sampling rate
signed char prc_get(void)
{
	return (signed char)prc_pdif-2; //value in range -1 to 1, must be lock to 0
	
}

//get pointers difference value at this moment
unsigned char prc_tst(void)
{
	signed char r;
	
	r=prc_in_rptr-prc_out_rptr; //pointer difference
	if(r<=0) r+=3; //in range 1-3, must be lock to 2
	return r;
}
