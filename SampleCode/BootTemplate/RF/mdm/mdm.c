/*
 Modulator and demodulator uses for pass speech over CMT2300A radio in FSK direct mode
 Frame's payload is 12 bytes (96 bits) corresponds AMR 4750bps audio codec
 Modulator append 4 sync bits after payload, gray and
 output bitstream 16*(96+4) probes (each probe is one bit, so output is two bytes
 is 0x0000 or 0xFFFF for each payload+sync bit 0 or 1). Probes will be transmitt to
 CMT2300A direct input over SPI master with 16*5000=80 000 bps bitrate.

 CMT2300A Receiver outputs data to direct output. This data captured by MC SPI master
 on 80 000 bps so we have bitstream 16 bits (probes) per one data bit
 Tis stream must be aligned to bit and frame bounaries then convert to 12 bytes
 of frame payload.  Additionaly subsequentional frame number 0-15 will be output.
 Also demodulator output sync flag 0 or 1 depends sync is compleet.

 Demodulator's input is unaligned stream of 100 bytes with probes captured by SPI
 Each processing output are aligned 12 bytes of payload, subsequentional frame number,
 sync flag and returns rate correction value in range +-6 for correct recever SPI rate
 must be same transmitter rate.
*/


#include <string.h>
#include <stdlib.h>

#include "mdm.h"

//user's settings
#define MD_ON 4 //matching counter trashold
#define MD_OFF 16 //unmatching counter trashold

//modem's parameters
#define POSINBIT 16      //number of probes in one data bit
#define BITSINFRAME 100   //number of data bits in frame
#define BYTESINFRAME 12   //number of byte in frame

//========================tables====================
//bit mask for byte
const unsigned char mask[8]={1,2,4,8,16,32,64,128};

//bit mask for word
//bit mask for byte
const unsigned short wmask[16]={1,2,4,8,16,32,64,128, 256,512,1024,2048,4096,8192,16384,32768};

//reverse bit order for 4 bits
const unsigned char rev[16]={0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
//decodes 4 probes to soft bit value
const signed char rbit[16]={-2,-1,-1,0,-1,0,0,1,-1,0,1,1,0,1,1,2};


//========================modulator=================
unsigned char tx_cnt=0; //counter of modulated frames
unsigned char tx_old=0; //for test only



//==========================demodulator=================
unsigned short md_pos_sh=0; //probes shifter
unsigned short md_posc[POSINBIT]; //probes waights

unsigned char md_bit_sh=0; //bits shifter
unsigned char md_bitc[BITSINFRAME]; //bits weights

unsigned char md_cnt=0; //counter of received frame
unsigned char md_del=0; //delta between our receive counter and frame number
unsigned char md_old=0; //number of prevoius frame

unsigned char md_on=0;  //counter of sync matches for set sync flag
unsigned char md_off=0; //counter of sync missmatch for clear sync flag
unsigned char md_syn=0; //sync flag (0 or 1)

unsigned char md_data[16]; //received byte array + statistic

volatile unsigned short md_test=0;


void mdm_init(void)
{
 //========================modulator=================
 tx_cnt=0; //counter of modulated frames
 tx_old=0; //for test only



//==========================demodulator=================
 md_pos_sh=0; //probes shifter
 memset(md_posc, 0, sizeof(md_posc));
 md_bit_sh=0; //bits shifter
 memset(md_bitc, 0, sizeof(md_bitc));
	
 md_cnt=0; //counter of received frame
 md_del=0; //delta between our receive counter and frame number
 md_old=0; //number of prevoius frame

 md_on=0;  //counter of sync matches for set sync flag
 md_off=0; //counter of sync missmatch for clear sync flag
 md_syn=0; //sync flag (0 or 1)
 memset(md_data, 0, sizeof(md_data));	
}



//-------------------------------------------------------
//Modulate payload is 96 bits (12 bytes) of data
//to 16*100 probes (100 u16 words)
//------------------------------------------------------
void mdm_modulate(unsigned char* data, unsigned short* out)
{
 unsigned short i;
 unsigned char bf[13];
 unsigned char bit;
 //unsigned short j, m; //for test only

 memcpy(bf, data, 12); //copy payload
 bf[12]=tx_cnt&0x0F; //copy counter
 if(tx_cnt&1)
 {
  bf[12]^=0x0E; //scramble counter
  for(i=0;i<12;i++) bf[i]^=0xAA; //scramble payload
 }
 else for(i=0;i<12;i++) bf[i]^=0x55; //scramble payload

 //convert bytes to probes each in one bit, 8 probes for 1 bit of payload
 for(i=0;i<BITSINFRAME;i++) //process each data bit to 8 probes
 {
   if(bf[i>>3]&mask[i&7]) bit=1; else bit=0; //get bit value
   if(bit) out[i]=0xFFFF; else out[i]=0x0000;  //convert to ideal probes

	 /*
   //=================Test ONLY!==============================
   //add noise on bit changed
   if(bit!=tx_old) //check bit changed
   {
    out[i]^=0x0002; //invert second bit after change
    if(i) out[i-1]^=0x4000; //invert second bit before change (not for first frame bit)
   }
   tx_old=bit; //save current bit for next


   //add spontaneous noise
   j=rand()%16; //0-7
   m=1;
   m<<=j;
   out[i]^=m; //invert 1 bit on byte

   j=rand()%16; //0-7
   m=1;
   m<<=j;
  // if((i<98)&&(i>17)) out[i]^=m; //invert 1 bit on byte

   j=rand()%16; //0-7
   m=1;
   m<<=j;
   out[i]^=m; //invert 1 bit on byte
 //==================end of test================================
*/

 }

 tx_cnt++;
}


 
//------------------------------------------------------
//demodulate 100 bytes (each data bit have 8 probes one bit each)
//output payload 12 bytes (96 bits) of data,
//------------------------------------------------------
unsigned short mdm_demodulate(unsigned short* in, unsigned char* data)
{
 unsigned char bit; //processed bit counter (0-99)
 unsigned char prb; //processed probes counter (0-7)

 unsigned char pos; //bit boundary lag (0-7)
 unsigned char lag; //lag of frame boundary (0-99)

 unsigned short posc[POSINBIT]; //probes weights for processed frame
 unsigned char dat; //current bit value (0 or 1)
 unsigned char cnt; //number of this frame
 unsigned char ccnt; //expected number of this frame

 unsigned char b, d, k; //unsigned
 signed char c; //signed
 short i,j; //short
 unsigned short n;

 //-----Step 1: search bit boundary-------------------------
 //search best probe position for detect data bit
 n=0;pos=0;
 for(i=0;i<POSINBIT;i++) //search on 16 possible positions
 {
  if(md_posc[i]>n) //if waight of this position the best
  {
   n=md_posc[i]; //save best weights
   pos=(unsigned char)i; //set this position the best
  }
 }

 //-----------Step 2: search frame boundary-----------------
 //search best frame boundary lag
 b=0;lag=0;
 for(i=0;i<BITSINFRAME;i++) //check all lags (on each data bit in frame)
 {
  d=md_bitc[i]>>4; //get match counter for this lag
  if(d>b) //if this lag is best
  {
   lag=(unsigned char)i; //set this lag as best
   b=d; //set best counter value for next comparing
  }
 }

//--------------------Step 3:  process all new probes in loop-------------
 //process all bits in frame
 memset(posc, 0, sizeof(posc)); //clear probes position counters for this frame
 for(bit=0;bit<BITSINFRAME;bit++) //bit loop
 {

  //------------------Substep 3.1: process probes in one bit--------------
  //process all probes
  for(prb=0;prb<POSINBIT;prb++) //pos loop
  {
   //add next probe to 16-probes shifter  (expected one data bit)
   md_pos_sh<<=1; //shift prevoius probes
   b=(!(!(in[bit]&wmask[prb]))); //get probe value 0 or 1 from input stream
   md_pos_sh|=b; //add current probe to shifter

   //detect expected bit, save only for best position, accumulate quality of each position
   c=rbit[(md_pos_sh>>4)&0x0F]; //value of probes 4-7 (in range -2 to 2)
   c+=rbit[(md_pos_sh>>8)&0x0F]; //value of probes 4-11 in range -4 to 4 (soft bit)
   if(prb==pos) dat=!(c>>7); //hard decission on best position: set current bit as sign of probes (>=0 is bit 1)
   c+=rbit[(md_pos_sh>>12)&0x0F]; //value of probes 4-15
   c+=rbit[(md_pos_sh>>0)&0x0F];  //value of probes 0-15 (in range -8 to 8)
   posc[prb]+=abs(c); //add quality to weight of current position (value in ramge 0 to 8)

  }//end op pos loop

  //---------------------Substep 3.2: save received bit------------------
  if(dat) //check received bit is 1
  {
   c=(signed char)bit-(signed char)lag-(signed char)1;
   if(c<0) c+=BITSINFRAME;  //get actual bit number (0-99) depends bit counter and lag
   md_data[c>>3]|=mask[c&7];//set received bit in bytes array
  }

  //---------------------Substep 3.3: check frame sync--------------
  //add next data bit to shifter
  md_bit_sh<<=1; //shift previout bits
  md_bit_sh|=dat; //add current bit to shifter

  //get current frame number
  b=rev[md_bit_sh&0x0F]; //try use 4 LSB as counter
  if(b&1) b^=0x0E; //descrambly counter

  //get expected frame number (incremented number of last received frame)
  d=0x0F&md_bitc[bit]; //get number of last frame for current lag
  d++; //expected number is previous number + 1
  d&=0x0F; //ring to 4 bits frame counter

  //upcount match or downcount unmatch
  k=md_bitc[bit]>>4; //get counter of match for current lag
  if(b==d) //check match of (numbers of last and current frames are subseqent)
  {
   if(k<15) k++; //increment matching counter up to 15
  }
  else if(k) k--; //not match: decrement down to 0
  md_bitc[bit]=b|((unsigned char)k<<4); //replace storage by current values


  //-----------------------Substep 3.4: check frame boundary and output received data--------
  //output frame on frame boundary
  if(lag==bit) //chech we receive last bit of frame
  {
    cnt=b; //counter of current frame
    ccnt=0x0F&(md_del+md_cnt);  //restore frame number from our counter and delta
		if(cnt&1) b=0xAA; else b=0x55; //set scrambling value
    for(i=0;i<12;i++) md_data[i]^=b; //descramble payload
    memcpy(data, md_data, BYTESINFRAME); //output frame payload 12 bytes
    memset(md_data, 0, BYTESINFRAME); //clear buffer for next
  }

 } //end of bit loop


 //---------------------Step 5: set/clear sync flag------------------
 //restore frame number
 b=0x0F&(cnt-md_cnt); //delta of frame counter for this frame
 if(b==md_old) //check delta is eqaual to last delta (so stable)
 {
  if(md_off) md_off--;  //down to 0 missmatch counnter
  md_on++; //up match counter
  if(md_on>=MD_ON) //if on trashold reached
  {
   md_on--; //fix on counter
   md_syn=1; //set syn flag
   md_del=b; //set current delta between our counter and numbers of received frames
  }
 }
 else //delta is changed: frame number is not subsequent
 {
  if(md_on) md_on--;  //down to 0 match counter
  md_off++; //up missmatch counter
  if(md_off>=MD_OFF) //if off trashold reached
  {
   md_off--;  //fix missmatch counter
   md_syn=0; //clear syn flag
  }
 }
 md_old=b; //set current delta value for next

 //--------------------Step 6: compute rate correction value for frame---------
 //add bit position to accumulator
 //search best value of bit position and get fast rate correction value
 c=0; j=0;
 for(i=0;i<POSINBIT;i++) //look for all positions for obtain bit boundary
 {
  md_posc[i]-=md_posc[i]>>5; //*32 filter
  md_posc[i]+=posc[i]; //accumulate position for averaging
  if(posc[i]>j) //check this position is best
  {
   j=posc[i]; //save position's value
   c=(signed char)i; //set best position in range 0-15
  }
 }


 //c-=8; //match position 0-15 to range -8 to 7
 //if(c<0) c++; //match in range -7 to 7 with two centering zero 
 //if(!md_syn) c=0; //skip rate correction while no sync

 //------------------------Step 7: output statistic-----------------
 //output statistic
 n=lag; n<<=8; //frame boundary lag (0-99)
 n|=(ccnt<<4); //frame number (0-15)
 n|=c; //bit position (0-15)
 if(md_syn) n|=0x8000; //sync lock flag (0/1)
 else data[0]|=1; //or set silency flag in outputted frame

 //data[12]=ccnt;  //restore frame number from our counter and delta
 //data[13]=pos; //bit boundary lag
 //data[14]=lag; //frame boundary lag
 //data[15]=md_syn; //flag of sync ok
 
 md_cnt++; //count processed frames

 return n;  //return probes rate correction value
}

