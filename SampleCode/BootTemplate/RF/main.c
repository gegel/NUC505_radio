/**************************************************************************//**
 * @file        main.c
 * @version     V1.00
 * $Revision:   1$
 * $Date:       15/05/20 3:00p$
 * @brief       Demonstrate simplex radio with GoperF RFM300W module based on CMT2300A transiever
 *
 *              Author: Segey Gayevsky, Gamma UA. Kiev, 2020
 * @note        The main() function cannot be debugged until C startup has completed.
 *              It is because C startup will be responsible for copying main() from ROM to RAM.
 *              So don't check the Run to main debug option or set BP before C startup is completed.
 * Copyright (C) 2014 Nuvoton Technology Corp. All rights reserved.
 ******************************************************************************/


/*
 
 Nuvoton NUC505 series have low-cost but low-speed internal SPI Flash.
 While code executed directly from this Flash overal perfomance cann't
 be more then only 1-2 MIPS. So we use bootloader places all data from 
 flash to RAM once on start. After this  all access will be to fast RAM 
 with zero weight state on maximal CPU frequence 100 MHz. This is need 
 for audio processing. AMR audio encoder require near 10 mS in this case,
 audio decoder - 2 mS. Frame duration is 20 mS so duplex can be run realtime.
 
 The only restriction is RAM size is 127K in this chip. 
 We use floating point AMR version because fixed point version too slow even with 
 using ETSI DSP for ARM. Floating point codec adopted for single precission 
 hardware of Cortex M4F. All AMR modes except 4750 bps and comfort noise generation 
 were remove from AMR codec source reducing code and tables size to useable with chip RAM.
 Stack size is 0x2000 and heap size is 0x3000 requires (codec use malloc!).
 Build result with level3 optimization for speed  is near: 
 CODE=41540 RO=31828 RW=604 ZI=24524
 
 MainInSRAM template project was used as a base: startup is placed in ROM, 
 all other code copies to RAM on boot stage. See main_on_sram.ld for details.
 Back effect is software power reset is not work so we use signal output pin
 connected to hardware reset pin for reboot chip after exiting power down.
 
 Power consumption in power down mode for NUC505 chip is near 700uA because
 RAM retention is required and there are no documented way for disabling it even 
 exiting power down only by reboot is used.
 before entering chip power down all periferal modules disabled and radio set
 to duty cycle periodic RX mode (20 mS RX to 500 mS sleep) with current is only 
 350uA (noraml RX current of CMT2300A is 8 mA).  
 
 Exiting power down avaliable by
 PTT button press (chip will be resetted immediately for entering work mode)
 or by pulse on radio output (will be detected on 8 subsequent phaze jumps of carrier).
 So there are sponaneous data packets (most transmitters works in packet mode on ISM)
 can be unexpected wakeups. So timeout from work to power down on radio and 
 user inactivity is so long (near 60 sec) we must prevent unexpected wakeups as possible.
 So after first pulse on radio out chip wake up and work on slow clock only 750 Khz
 Chip poll radio output during 1 sec (2 duty RX cycles) and count pulses on radio out
 On carrier there are manu pulses but in occasional packet there are only 1-2 pulses,
 on this case CPU back to power down for new pulse on radio out. So we have low probalility
 of false carrier detection and comparatively low  ready for RX (near 1-1.5 sec after
 transmitting started). 
 
 CDR (clock data restore) is disable in CMT2300A, so we have "raw" data on radio output.
 SPI sample this data with 80KHz rate, so we have 16 probes for one data bit (5000 pbs radio rate).
 There are software CDR in demodulation software discower bit boundary and detect soft bits (LLR)
 on best position skipping probes on sides of bit field. The result is averaged correction value for 
 lock receiving SPI rate to transmittin rate can be a little difference with +-(20+20) ppm
 of crystal resonator inaccuracy. 
 
 Other demodulator task is align data block (encoded audio frame) bits boundary.
 In packet mode this is perhaps thanks preabule and followed sync sequence. But in stream
 mode this way have big overhead. We use 4 bits of sync followed every 96 payload bits.
 So each 12 bytes compressed audio frame (20 mS duration) transmitted as 96 payload + 4 sync 
 bits with 5Khz radio data rate. Sync is frame counter from 0 to 15. There are scrambling depends
 counter;s LSB: frame scrambled with different masks and three MSB of counter inverted on LSB=1
 Demodulator search bit lag as the position of best accumulated value of incremented previous 
 counter value match current value. Besides there are mutch and unmatch counters with trasholds 
 for set/reset carrier lock flag used for mute audio output on carrier loss.
 And there are local frame counter corrected by received counters values. This counter is stable
 (not affected by noise) and can be use, for example, for example for cryptographic proprerties. 
 
*/



//-----------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NUC505Series.h"
#include "hw/gpc.h" //io hardware 
#include "rf/radio_if.h"  //radio control
#include "hw/spp.h" //spi receive/send
#include "hw/aud.h"  //audio record/play
#include "hw/com.h"  //audio record/play
#include "hw/pdc.h"  //entering power down mode
#include "prc.h"   //processing


//-------------------definitions-------------------
#define SLEEP_TOUT 3000  //x20mS timeout of inactivity for entering sleep mode 
#define PRE_TX 50 //number of silency frames transmitted on TX start
#define PRE_RX 50 //number of frames muted on RX start
#define LOG_ON //output modem statistic to UART

//------------------prototypes---------------
void SYS_Init(void);   //system initialization
void SYS_PD(void); //preparing to power down mode

//---------------global values---------------
unsigned short* uptr=0; //bits pointer
short* sptr=0; //pcm pointer

unsigned char ison=0; //flag of radio avaliable
unsigned char rch=0; //radio channel
unsigned char vol=25; //in range 0-31
unsigned char rmode=1;  //radio mode (0-TX, 1-RX)

unsigned short w; //modem state
unsigned char syn=0; //RX carrier lock flag
unsigned char vad=0; //TX voice detector flag

unsigned char cnt=0; //LED blink counter
unsigned short tmr=0; //sleep timer
 
char str[16];   //string for output statistic in test


#if defined (__GNUC__)
#define VECTOR_SIZE     32
uint32_t VectorTable[VECTOR_SIZE] __attribute__ ((aligned(256)));
#endif


//main procedure with loop
int main(void)
{
    /* Init System, IP clock and multi-function I/O */
    SYS_Init();

    /* Relocate vector table in SRAM for fast interrupt handling. */
    {		
#if defined ( __CC_ARM )
        extern uint32_t __Vectors[];
        extern uint32_t __Vectors_Size[];
        extern uint32_t Image$$ER_VECTOR2$$ZI$$Base[];

        //printf("Relocate vector table in SRAM (0x%08X) for fast interrupt handling.\n", Image$$ER_VECTOR2$$ZI$$Base);
        memcpy((void *) Image$$ER_VECTOR2$$ZI$$Base, (void *) __Vectors, (unsigned int) __Vectors_Size);
        SCB->VTOR = (uint32_t) Image$$ER_VECTOR2$$ZI$$Base;

#elif defined (__ICCARM__)
#pragma section = "VECTOR2"
        extern uint32_t __Vectors[];
        extern uint32_t __Vectors_Size[];

        //printf("Relocate vector table in SRAM (0x%08X) for fast interrupt handling.\n", __section_begin("VECTOR2"));
        memcpy((void *) __section_begin("VECTOR2"), (void *) __Vectors, (unsigned int) __Vectors_Size);
        SCB->VTOR = (uint32_t) __section_begin("VECTOR2");

#elif defined (__GNUC__)
        extern uint32_t __Vectors[];
        extern uint32_t __Vectors_Size[];
        memcpy(VectorTable, (uint32_t*)0x0, (unsigned int) __Vectors_Size);
        SCB->VTOR = (uint32_t)VectorTable;
#endif
    }

	//-------------------start initialization---------------------	
   
   radio_ch(rch); //set radio channel  
	 radio_init(1); //init radio to work mode
	 com_init(); //initialize uart
   amr_init(); //initialize codec	
	 spp_init();  //init spi 
	 aud_init();  //init audio and uart
   
	strt:	  //entering to main loop (after swithing TX/RX)
	
	 //initialize	
	 prc_init(); //reinit processing engine
	 tmr=0; //clear sleep timer
	 
	//get RX/TX mode by PTT button
	 rmode = CMT2300A_ReadBTN(); //set RX or TX mode mode by button press/release
	 if(CMT2300A_ReadTST()) //work mode
	 {
		 if(rmode) aud_vol(vol); else aud_vol(0); //mute playing on TX
		 ison=radio_mode(rmode); //apply receiving/transmitting radio mode and channel
		 
	 }
	 else //test mode
	 {
	  ison=radio_mode(0); //force TX
	  prc_vt(1); //force voice count instead speech from mike
	 }
	
	 if(ison) CMT2300A_WriteTST_1(); else CMT2300A_WriteTST_0();  //off on-board LED
	 
	 //set mike ans speacker control hw pins
	 CMT2300A_WriteSPC_0(); //disable speaker
	 if(!rmode) //transmitting
	 {
		 CMT2300A_WriteMIC_1(); //enable MIC on TX mode
   }
	 else //transmitting
	 {
		 CMT2300A_WriteMIC_0(); //receiving
	 }
	
	 
//==========================main loop===============================
	 while(rmode==CMT2300A_ReadBTN()) //loop while button state keep stable
	 //while(1)
	 {
    //--------------------Transmitting---------------------
		 
		 sptr=aud_grab(); //check audio grabbing pcm frame
     if(sptr) vad=prc_enc((short*)sptr); //encoding	grabbed pcm frame
		 
		 uptr=spp_send(); //check radio send bits
		 if(uptr) prc_mdm(uptr); //modulating new bits for next sending
		 
		//---------------------Receiving------------------------

    uptr=spp_rcvd();  //check radio receiving bits
    if(uptr) //demodulation 
		{
     unsigned char frame_num;
		 unsigned char bit_lag;			
		 unsigned char play_lag;	
			
		 w=prc_dmd(uptr); //demodulate received bits, get resulting flags	
     spp_rate(w); //correct receiving rate
		 syn=w>>15; //get sync flag
		 frame_num=(w>>4)&0xF; //get frame number
     bit_lag=w&0x0F; //get bit lag
     play_lag=prc_tst(); //get state of playing fifo

			
#ifdef LOG_ON
     //output statistic for test			
		 sprintf(str, "%X%c%X%c%d\r\n", frame_num, syn+'*', bit_lag, (char)vad+'*', play_lag);
     com_write((unsigned char*)str, 7);
#endif


			//LED indication
	   if(!rmode)  //TX
		 {                                //during TX:
			if(vad) CMT2300A_WriteLED_1(); else CMT2300A_WriteLED_0(); //set LED on speech/silency 
      tmr=0; //clear sleep timeout on sync			  
		 }                                          //during RX:
		 else if(syn) //sync OK
		 {
			if(cnt&2) CMT2300A_WriteLED_1(); else CMT2300A_WriteLED_0(); //blink fast
      CMT2300A_WriteSPC_1(); //speaker on
      tmr=0; //clear sleep timeout on sync			 
		 }
		 else //no sync
		 {
			if((cnt&0x3F)<4) CMT2300A_WriteLED_1(); else CMT2300A_WriteLED_0(); //blink slow
			CMT2300A_WriteSPC_0(); //speaker off
		 }
		 cnt++;	 
				
	   //check for sleep
     if((++tmr)>SLEEP_TOUT) //check for sleep timeout
		 {  
			 SYS_PD(); //preparing MC for power down: disable all periferal modules and interupts 
			 radio_init(0); //initialize radio to power saving mode
			 pdc_init(); //set IO pins, change system clock to internal RC and entering power down mode 
		 }
	 
		}	//end of demodulation		
		
		sptr=aud_play(); //check audio play frame		
		if(sptr) prc_decr(sptr); //decode new frame for next playing	
	
		
//-------------------Reading UART command--------------------------			
	  w=com_poll(); //poll having data block receiving over uart
		if(w)
		{
			unsigned char* cmd;
			if(w>16) w=16; //restrict length
			cmd=com_read(); //get pointer to received data
			memcpy(str, cmd, w); //copy received data to string (for test)
      w=radio_cmd((unsigned char*)str);
      if(w) com_write((unsigned char*)str, w); 			
		}
	  
	 
	 } //and of while (main loop): button state was changed
		
	 goto strt; //button stae was changed: re-enter loop	
		
 }		
		
		
	//sysem initialization
 void SYS_Init(void)
{
    /*---------------------------------------------------------------------------------------------------------*/
    /* Init System Clock                                                                                       */
    /*---------------------------------------------------------------------------------------------------------*/

    /* Enable  XTAL */
    CLK->PWRCTL |= CLK_PWRCTL_HXTEN_Msk;

    /* Select IP clock source */
    /* PCLK divider = 1 (/2) */
    CLK_SetModuleClock(PCLK_MODULE, 0, 1);
    
    /* Update System Core Clock */
    /* Note too high system clock will cause over-spec of SPI Flash read command on running code on SPI Flash. */
    CLK_SetCoreClock(100000000); //set 100MHz CPU clock
    SystemCoreClockUpdate();

	  //Enable clocks of used periferal modules
 	  CLK_EnableModuleClock(UART1_MODULE);
    CLK_EnableModuleClock(TMR0_MODULE);
		CLK_EnableModuleClock(TMR2_MODULE);
		CLK_EnableModuleClock(SPI0_MODULE);
		CLK_EnableModuleClock(SPI1_MODULE);

    /* Select IP clock source */
		CLK_SetModuleClock(UART1_MODULE, CLK_UART1_SRC_EXT, 0);
    CLK_SetModuleClock(TMR0_MODULE, CLK_TMR0_SRC_EXT, 0);
		CLK_SetModuleClock(TMR2_MODULE, CLK_TMR2_SRC_EXT, 0);
    CLK_SetModuleClock(SPI0_MODULE, CLK_SPI0_SRC_EXT, 0);
		CLK_SetModuleClock(SPI1_MODULE, CLK_SPI1_SRC_EXT, 0);

    /* Init I/O multi-function pins */  
    SYS->GPA_MFPH = SYS_GPA_MFPH_PA9MFP_UART1_RXD | SYS_GPA_MFPH_PA8MFP_UART1_TXD;
    SYS->GPA_MFPL = 0x00000000;
    SYS->GPB_MFPH = SYS_GPB_MFPH_PB13MFP_SPI1_MISO | SYS_GPB_MFPH_PB12MFP_SPI1_MOSI | SYS_GPB_MFPH_PB11MFP_SPI1_CLK | SYS_GPB_MFPH_PB10MFP_SPI1_SS;
    SYS->GPB_MFPL = SYS_GPB_MFPL_PB5MFP_SPI0_MISO | SYS_GPB_MFPL_PB4MFP_SPI0_MOSI | SYS_GPB_MFPL_PB3MFP_SPI0_CLK | SYS_GPB_MFPL_PB2MFP_SPI0_SS;
    SYS->GPC_MFPH = 0x00000000;
    SYS->GPC_MFPL = 0x00000000;
    SYS->GPD_MFPL = 0x00000000;
}

//Disable periferal interupts before entering power down mode
void SYS_PD(void)
{
 //disable all used interupts
			 NVIC_DisableIRQ(SPI0_IRQn); //disable SPI0 IRQ
       SPI_DisableInt(SPI0, SPI_FIFO_RXTH_INT_MASK); //enable RX interupt
       SPI_Close(SPI0);
			 
			 NVIC_DisableIRQ(SPI1_IRQn); //disable SPI1 IRQ
       SPI_DisableInt(SPI1, SPI_FIFO_TXTH_INT_MASK); //enable TX interupt
       SPI_Close(SPI1);
			 
	     UART_DisableInt(UART1, (UART_INTEN_RDAIEN_Msk | UART_INTEN_RXTOIEN_Msk));
       NVIC_DisableIRQ(UART1_IRQn);  //disable UART1 IRQ
			 UART_Close(UART1);
	
	
	
	     NVIC_DisableIRQ(I2S_IRQn); ////disable I2S IRQ
       I2S_DisableInt(I2S, (I2S_IEN_RDMATIEN_Msk|I2S_IEN_RDMAEIEN_Msk));
	     I2S_Close(I2S);      
}
 
 