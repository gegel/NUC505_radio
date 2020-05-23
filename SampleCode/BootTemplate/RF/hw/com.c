/*
 Procedures for UART1 RX/TX. So UART1 has 64 bytes RX and TX fifos 
 there is no need to USE DMA or interupts
 Receiving event generate by timeout of received byte stream
 Receving block of bytes (1-64 byttes) put to buffer 
 and rx procedure returns pointer to this buffer
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NUC505Series.h"

#include "com.h"

//----------global variables--------------

unsigned char com_rout[COM_RBUF]; //output buffer for store received data block

//initialize UART0
void com_init(void)
{
  SYS_ResetModule(UART1_RST); //reset UART0 module
	UART_Open(UART1, COMBDR); //open UART0 with specified baudrate (interupts must be rary as possible
  UART_SetTimeoutCnt(UART1, COMTOUT);	//set commbytetimeut
  UART_EnableInt(UART1, (UART_INTEN_RDAIEN_Msk | UART_INTEN_RXTOIEN_Msk)); //enable interupts: RX fifo level reached, RX timeout
}

//send sata over UART
short com_write(unsigned char* data, short len)
{
	short i;
	
	for(i=0;i<len;i++) //put bytes to TX fifo (64 bytes total size)
	{
		if(UART_IS_TX_FULL(UART1)) break; //brak if fifo full
		UART_WRITE(UART1, data[i]); //put next byte	
	}
	
	return i; //return number of bytes put to buffer
}

//poll having received block
short com_poll(void)
{
	short i=0;
	if(UART_GET_INT_FLAG(UART1, UART_INTSTS_RXTOINT_Msk)) //chack rx timeout flag
	{
	  if(UART_IS_RX_FULL(UART1)) i=64;
		else i=(short)((UART1->FIFOSTS)>>8)&0x3F; //get rx fifo level
	}
	return i; //return number of bytes in rx fifo or 0
}

//get pointer to received data and clear flags
unsigned char* com_read(void)
{	
	short i=0;
	while(!UART_GET_RX_EMPTY(UART1)) //check rx fifo is not empty
	{
   com_rout[i++]=UART_READ(UART1); //read byte from FIFO
	}		
	return com_rout; //returns pointer to received block
}



