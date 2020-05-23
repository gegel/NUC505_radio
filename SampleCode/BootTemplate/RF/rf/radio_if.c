/*
 Radio conrol interface for Hoperf RFM300 module based on CMT2300A transiever in direct mode
 Baud rate is 5KHz, freq. deviation is 5KHz, base freq. is 868.0MHz channel step is 10KHz up.
 Work mode use direct TX/RX without CDR and without mute output gpio.
 Power save mode use 20 mS RX / 500mS sleep authomatic duty cycle with muting output gpio to 0
 RSSI source in power down mode for unmute radio output is PJD with 8 jups level.
 Tables for initialization registers for work and power down modes created with RFPDF grapic tool
 See exported text files cmt2300a_params.exp and pd_cmt2300a_params.exp for details.
 In Tx mode GPIO pin set as input, in work RX mode, sleep mode and no-inited - as output (inactive - 0)
 In work mode we can fast set TX or RX mode, new frequency channel will be applied.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gpc.h"
#include "radio.h"
#include "radio_if.h"

unsigned char rd_ch=0; //radio channel 0-255
unsigned char id=0; //device id 0-255 (not used in direct mode)

//init radio for work (0) or power_down(1)
void radio_init(unsigned char on)
{
 RF_Init(on);   //init registers from preset tables
 RF_Config(rd_ch); //config and goto sleep	
}
 
//set frequency channel (will be apply on TX/RX mode change)
unsigned char radio_ch(short ch)
{
	if(ch>=0) rd_ch=(unsigned char)ch;
	return rd_ch;
}

//set frqequence channel over UART command as string F=123<CR>
unsigned char radio_cmd(unsigned char* cmd)
{
	short i; //loop counter
	unsigned short c=0; //accumulator for atoi output
	unsigned char d; //current char
	unsigned char e=0; //compleet flag
	
	//check command header
	if((cmd[0]=='F')&&(cmd[1]=='=')) //must be F=
	{
	 	//process up to 3 digits in loop
		for(i=2;i<5;i++) //process channel digits can be 1 to 3 (value from 0 to 255)
		{
		 d=cmd[i]-'0'; //digits will be in range 0-9, all great is not digits
		 if(d>'9') break; //break on first non-digit
		 c*=10; //multiple accumulator by 10 for add new digit
  	 c+=d; //add new digit to accumulator		
		}
	  
		//set compleet flag
		if((c<256)&&(d==0xDD)) e=1; //check accu,ulated value maus be in range 0-255 end first non-digit is <CR>
	}
	
	//check compleet flag
	if(e)  //atoi compleet
	{
		rd_ch=(unsigned char)c; //set new radio channel
		e=5; //set number of chars in answer
		memcpy(cmd, "OK\r\n", e); //answer "OK"
	}
	else //atoi fail
	{
		e=8; //set number of of chars in answer
		memcpy(cmd, "ERROR\r\n", e); //answer "ERROR"
	}
	return e; //return answer length
}


//change radio mode (rx or tx) and set frequency channel (0-255)
unsigned char radio_mode(unsigned char mode)
{
 CMT2300A_GoStby(); //goto standby mode first
 if(0==CMT2300A_IsExist()) return 0; //check CMT2300 is exist
 CMT2300A_ClearInterruptFlags(); //clear interupts flags	 
 CMT2300A_SetFrequencyChannel(rd_ch);	//SET_GPIO_AFOD new radio channel
 	
 if(mode) //rx
 {
	 CMT2300A_EnableTxDin(0); //set GPIO as output
	 CMT2300A_GoRx(); //goto RX mode	
 }
 else //tx
 {
	 CMT2300A_EnableTxDin(1); //set GPIO as input
   CMT2300A_GoTx(); //goto TX mode	  
 }
 
	return 1;
}



