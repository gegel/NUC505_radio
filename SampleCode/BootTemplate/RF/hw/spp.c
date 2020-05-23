/*
 Procedures for working with SPI0 (RX) and SPI1(TX)
 Transmitter work in direct mode: during TX SPI1 must send bits 1/0 to input with
 bitrate of 5000Hz. SPI1 works with CLK=80KHz so each bit coded as 0x0000 or oxFFFF
 Receiver also work in direct mode so on RX output are bits 1 /0 with bitrate 5000KHz
 can be a little difference with our bitrate (cause +- 40ppm difference between crystals on RX and TX sides)
 SPI0 receive 16 probes for one data bit uses in software demodulator for CDR.
 Clock of SPI1 is 80Khz can be adjusted a little for lock with clock on transmitter side
 Rate correction value returned by demodulator only if carrier discowered.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NUC505Series.h"


#include "spp.h"

//------------SPI SEND-------------------------
unsigned int spp_tbuf[2*BLOCKLENGTH]; //u32 word buffer for spi send
volatile unsigned char spp_tptr=0; //pointer to TX word in buffer
volatile unsigned char spp_tflg=0; //flag of frame sended

//------------SPI RCVD-------------------------
unsigned int spp_rbuf[2*BLOCKLENGTH]; //u32 word buffer for spi rcvd
unsigned char spp_rptr=0; //pointer to TX/RX word in buffer
volatile unsigned char spp_rflg=0; //flag of frame ready

//-----------RCVD rate correction--------------
volatile unsigned char spp_r=0; //rate will be applied to SPI clock
volatile unsigned char spp_c=0; //counter of periods rate will be applied 
short spp_d=0;  //average rate 
unsigned char spp_n=0; //number of averaging values

//#define SPI1_TEST   //test for change rate of TX SPI1


#ifdef SPI1_TEST
volatile unsigned char spp_tspi_r=0;  //new rate for change
volatile unsigned char spp_tspi_c=0;  //counter of apply rate changing
#endif

//initialize SPI
void spp_init(void)
{
    short i;
	
	  //clear global values
	  spp_tptr=0; //pointer to TX word in buffer
    spp_tflg=0; //flag of frame sended
	  spp_rptr=0; //pointer to TX/RX word in buffer
    spp_rflg=0; //flag of frame ready
	  spp_r=0; //rate will be applied to SPI clock
    spp_c=0; //counter of periods rate will be applied 
    spp_n=0; //average rate 
    spp_d=0; //number of averaging values
	  memset(spp_tbuf, 0, sizeof(spp_tbuf)); //clear buffers
	  memset(spp_rbuf, 0, sizeof(spp_rbuf));
	
	  SYS_ResetModule(SPI0_RST); //reset receving SPI
	  SYS_ResetModule(SPI1_RST); //reset transmition SPI
	  
	  //start receiving SPI
	  SPI_Open(SPI0, SPI_MASTER, SPI_MODE_3, 32, 80000); //CLK idle low, x change on rising, rx latch in falling
    SPI_SET_LSB_FIRST(SPI0); //set LSB to MSB mode
		SPI_SET_SUSPEND_CYCLE(SPI0, 0); //no suspend beetwen words
		SPI_ENABLE_3WIRE_MODE(SPI0); //not use SS line
    SPI_ClearRxFIFO(SPI0); //clear RX FIFO
		SPI_ClearTxFIFO(SPI0); //clear RX FIFO
    SPI_SetFIFO(SPI0, 2, 6); //set FIFO levels
		NVIC_EnableIRQ(SPI0_IRQn); //enable IRQ
		SPI_EnableInt(SPI0, SPI_FIFO_RXTH_INT_MASK); //enable RX interupt
	  for(i=0;i<8;i++) SPI_WRITE_TX(SPI0, 0xFFFFFFFF); //start work 
    
		//start transmitting SPI
		SPI_Open(SPI1, SPI_MASTER, SPI_MODE_3, 32, 80000); //CLK idle low, x change on rising, rx latch in falling
		SPI_SET_LSB_FIRST(SPI1); //set LSB to MSB mode
		SPI_SET_SUSPEND_CYCLE(SPI1, 0); //no suspend beetwen words
		SPI_ENABLE_3WIRE_MODE(SPI1); //not use SS line
    SPI_ClearRxFIFO(SPI1); //clear RX FIFO
		SPI_ClearTxFIFO(SPI1); //clear RX FIFO
    SPI_SetFIFO(SPI1, 2, 6); //set FIFO levels
		NVIC_EnableIRQ(SPI1_IRQn); //enable IRQ
		SPI_EnableInt(SPI1, SPI_FIFO_TXTH_INT_MASK); //enable TX interupt
		for(i=0;i<8;i++)  SPI_WRITE_TX(SPI1, 0xFFFFFFFF);
}

//correct receiving SPI rate using value returns by demodulator
//actually our receiving SPI rate must be lock to their transmitting SPI rate
void spp_rate(unsigned short nn)
{
 if(nn&0x8000) //check demodulator sync ok 
 {
	       short j;
	       unsigned char c;
         
	       //spi rx rate correction
				 j=(short)(nn&0x0F)-8; //convert bit boundary 0-15 to range -8 to 7 
         if(j<0) j++; //set to symmetric range -7 to 7 with two centering zeroes
				 spp_d+=j; //average value
				 spp_n++;  //count averaged values
					
				 //check we have 8 averaged rate values in accumulator	
				 if(spp_n&4) 
				 {	 		
					spp_d>>=1; //get averaged value 
				  c=0; //clear timer rate
				  if(spp_d <-1) c=TMR_SPI-1; //or set timer rate for increasing bit position
          if(spp_d > 1) c=TMR_SPI+1; //or set timer rate for decreasing bit position
				  if(spp_d) spp_c=1  +  ((abs(spp_d))>>2); // spi_c=2; //1+abs(j); //set count of spi words will be rx with corrected rate
				  spp_r=c; //set rate will be applied from next spi word (in interupt)
					spp_d=0; //clear accumulator
					spp_n=0; 
				 }			 
				 
	}
	else //no sync
	{
					spp_d=0; //clear averaged values
					spp_n=0; //clear number ov values
	}
}

//poll receiving SPI have bits, returns 1 and pointer to receving bits
unsigned short* spp_rcvd(void)
{
	unsigned short* ptr=0;
	if(spp_rflg) //poll receving ready flag
	{
	 unsigned char c;
   c=spp_rflg-1; //number of double buffer's part with received bits
	 spp_rflg=0; //clear flag
   ptr=(unsigned short*)(spp_rbuf+c*BLOCKLENGTH); //set pointer to buffer with received bits
	}
	return ptr;
}

//poll sending SPI need bits, returns 1 and pointer to buffer for new bits 
unsigned short* spp_send(void)
{
	unsigned short* ptr=0;
	if(spp_tflg) //poll transmitting compleet flag
	{
	 unsigned char c;
   c=spp_tflg-1; //number of double buffer's part for new data
	 spp_tflg=0; //clear flag
   ptr=(unsigned short*)(spp_tbuf+c*BLOCKLENGTH); //set pointer to buffer for new bits for transmitt
	}
	return ptr;
}

#ifdef SPI1_TEST
//test procedure for change SPI1 TX rate
void spp_test(signed char r)
{
	spp_tspi_c=2;
  spp_tspi_r=TMR_SPI+r;
}
#endif

//SPI0 data received
void SPI0_IRQHandler(void)
  {	
	  //call PRC rate correction here!!!
		
		//apply SPI0 rate correction value
		if(spp_r) //check new rate correction value
	  {
		 SPI0->CLKDIV=spp_r;	//set new rate
	   spp_r=0;	//clearerr value (set once!)
	  }
		 
		//restore SPI0 base rate
	  if(spp_c) //check rate correction counter
	  {
	   spp_c--; //decrement counter
     if(!spp_c) SPI0->CLKDIV=TMR_SPI; //if correction interval elapser set basic rate						
	  }
		
		//keep transmitt buffer full for run master
		while(SPI_GET_TX_FIFO_FULL_FLAG(SPI0) == 0) //while transmitting buffer will be full
		{	
			SPI_WRITE_TX(SPI0, 0xFFFFFFFF);
	  }
		
		//get all received data from RX FIFO
		while(SPI_GET_RX_FIFO_EMPTY_FLAG(SPI0) == 0) //while receiving buffer will be empty
    {	
			spp_rbuf[spp_rptr] = SPI_READ_RX(SPI0); //read data from FIFO
			spp_rptr++; //move pointer
			if(spp_rptr==BLOCKLENGTH) spp_rflg=1; //check first part is ready, set flag
			else if(spp_rptr==2*BLOCKLENGTH) //check second part is ready, set flag and ring double buffer pointer
			{
				spp_rptr=0; //ring spi double buffer
				spp_rflg=2; //set flag for second part of buffer ready
			}	 
    }	
  }
	

	//SPI1 data send
	void SPI1_IRQHandler(void)
  {	
		
#ifdef SPI1_TEST		
		
		//apply SPI0 rate correction value
		if(spp_tspi_r) //check new rate correction value
	  {
		 SPI1->CLKDIV=spp_tspi_r;	//set new rate
	   spp_tspi_r=0;	//clearerr value (set once!)
	  }
		 
		//restore SPI0 base rate
	  if(spp_tspi_c) //check rate correction counter
	  {
	   spp_tspi_c--; //decrement counter
     if(!spp_tspi_c) SPI1->CLKDIV=TMR_SPI; //if correction interval elapser set basic rate						
	  }
#endif		
		//put data to FIFO up to fullfill
		while(SPI_GET_TX_FIFO_FULL_FLAG(SPI1) == 0) //while transmitting buffer will be full
    {
			SPI_WRITE_TX(SPI1, spp_tbuf[spp_tptr]); //put data
			spp_tptr++; //move pointer
			if(spp_tptr==BLOCKLENGTH) spp_tflg=1; //check first part is processed, set flag
			else if(spp_tptr==2*BLOCKLENGTH) //check second part is processed, set flag and ring double buffer pointer
			{
				spp_tptr=0; //ring transmitting double buffer
				spp_tflg=2; //set flag the second part of transmitting buffer is empy
			}	
    }
  }
	
	


