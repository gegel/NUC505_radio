#include "NUC505Series.h"

#include "typedefs.h"
#include "gpc.h"

unsigned char cmt_tout=0;

//=================Timings============================

//1us delay for software SPI by loop
void cmt_spi3_delay(void)
{
    volatile u16 n = SPI_DELAY_COUNT_MS;
    while(n--);
}

//1-65535 us delay by timer
void system_delay_us (unsigned short us)
{	
 unsigned char delay=9; //number of MCU cycles waits for Timer start
	
 TIMER2->CTL = 0;  //clear Timer
 TIMER2->EXTCTL = 0;
 TIMER2->CMP = us; //set delay in uS
 TIMER2->CTL = TIMER_CTL_CNTEN_Msk | TIMER_ONESHOT_MODE | 11; //start timer in oneshot mode with 1Mhz clock (12 Mhz/(11+1))
 for(; delay > 0; delay--) {__NOP();}  //wait a few for Timer start 
 while(TIMER2->CTL & TIMER_CTL_ACTSTS_Msk); //wait wor timer shot
}




void cmt_gpio_pd(void)
{	
	//set all used GPIO as input
	//GPIO_SetMode(PC, BIT14, GPIO_MODE_INPUT);//BTN GND //BTN must work in power-down!
  GPIO_SetMode(PC, BIT13, GPIO_MODE_INPUT);//LED
  
	GPIO_SetMode(PB, BIT8, GPIO_MODE_INPUT); //CS
  GPIO_SetMode(PB, BIT7, GPIO_MODE_INPUT); //CLK
  GPIO_SetMode(PB, BIT6, GPIO_MODE_INPUT);  //DIO
  GPIO_SetMode(PB, BIT9, GPIO_MODE_INPUT); //FB
	
  //GPIO_SetMode(PC, BIT3, GPIO_MODE_INPUT);//TST	 //only for testing: on-board led still active in power-down
  //GPIO_SetMode(PA, BIT10, GPIO_MODE_INPUT); //MIC //Mike control still left enabling in power-down
	//GPIO_SetMode(PA, BIT11, GPIO_MODE_INPUT); //SPC	///Speaker control still left enabling in power-down
	CMT2300A_WriteSPC_0(); //speaker and mike must keep disabled in power-down
	CMT2300A_WriteMIC_0();
}


//=====================GPIO===================================

//initialize GPIO
void cmt_gpio_init(void)
{
  GPIO_SetMode(PC, BIT0, GPIO_MODE_INPUT); //RST (only spike for solve problem: software system reset does not return control to the bootloader in ROM)
	GPIO_SetMode(PC, BIT12, GPIO_MODE_INPUT);//BTN   (fixme: must be replace to PORTA or PORTB !!!)
	GPIO_SetPullMode(PC, BIT12, GPIO_PULL_UP_EN);
	GPIO_SetMode(PC, BIT14, GPIO_MODE_OUTPUT); PC14_DOUT=0; //GND for BTN
  GPIO_SetMode(PC, BIT13, GPIO_MODE_OUTPUT); PC13_DOUT=1;//LED
 
	GPIO_SetMode(PC, BIT5, GPIO_MODE_OUTPUT); PC5_DOUT=0; //GND for TEST
  GPIO_SetMode(PC, BIT6, GPIO_MODE_INPUT);//TEST
	GPIO_SetPullMode(PC, BIT4, GPIO_PULL_UP_EN);
	
	GPIO_SetMode(PB, BIT8, GPIO_MODE_OUTPUT); PB8_DOUT=1; //CS
  GPIO_SetMode(PB, BIT7, GPIO_MODE_OUTPUT); PB7_DOUT=1; //CLK
  GPIO_SetMode(PB, BIT6, GPIO_MODE_OUTPUT); PB6_DOUT=1; //DIO
  GPIO_SetMode(PB, BIT9, GPIO_MODE_OUTPUT); PB9_DOUT=1; //FB
  
	GPIO_SetMode(PC, BIT11, GPIO_MODE_INPUT); //GPIO (fixme: must be replace to PORTA or PORTB !!!)
  GPIO_SetMode(PC, BIT3, GPIO_MODE_OUTPUT); PC3_DOUT=1;//TST
	GPIO_SetMode(PA, BIT10, GPIO_MODE_OUTPUT); PA10_DOUT=0; //MIC
	GPIO_SetMode(PA, BIT11, GPIO_MODE_OUTPUT); PA11_DOUT=0; //SPC
	
	GPIO_SetMode(PB, BIT0, GPIO_MODE_INPUT);//GPIO interupt (only spike for solve problem: GPIOC not handle wakeup interupts in power-down mode)
	GPIO_SetMode(PB, BIT1, GPIO_MODE_INPUT);//BTN interupt (only spike for solve problem: GPIOC not handle wakeup interupts in power-down mode)
}


unsigned char cmt_read_p3(void)
{
	return (unsigned char)CMT2300A_ReadGpio3();	
}








//void procedure for emtu macro
void cmtvoid(void)
{
	
}

/*
////void cmt_spi3_csb_out(void){}//      SET_GPIO_OUT(CMT_CSB_GPIO)
void cmt_spi3_fcsb_out(void){}//     SET_GPIO_OUT(CMT_FCSB_GPIO)
void cmt_spi3_sclk_out(void){}//     SET_GPIO_OUT(CMT_SCLK_GPIO)
void cmt_spi3_sdio_out(void){}//     SET_GPIO_OUT(CMT_SDIO_GPIO)
void cmt_spi3_sdio_in(void){}//      SET_GPIO_IN(CMT_SDIO_GPIO)

void cmt_spi3_csb_1(void){}//        SET_GPIO_H(CMT_CSB_GPIO)
void cmt_spi3_csb_0(void){}//        SET_GPIO_L(CMT_CSB_GPIO)

void cmt_spi3_fcsb_1(void){}//       SET_GPIO_H(CMT_FCSB_GPIO)
void cmt_spi3_fcsb_0(void){}//       SET_GPIO_L(CMT_FCSB_GPIO)
    
void cmt_spi3_sclk_1(void){}//       SET_GPIO_H(CMT_SCLK_GPIO)
void cmt_spi3_sclk_0(void){}//       SET_GPIO_L(CMT_SCLK_GPIO)

void cmt_spi3_sdio_1(void){}//       SET_GPIO_H(CMT_SDIO_GPIO)
void cmt_spi3_sdio_0(void){}//       SET_GPIO_L(CMT_SDIO_GPIO)
unsigned char cmt_spi3_sdio_read(void){return 0;}//    READ_GPIO_PIN(CMT_SDIO_GPIO)


void CMT2300A_SetGpio1In(void) //           SET_GPIO_IN(CMT_GPIO1_GPIO)
{
  //set GPIO to input mode
}

void CMT2300A_SetGpio2In(void) //           SET_GPIO_IN(CMT_GPIO2_GPIO)
{

 //set GPIO to input mode
}


void CMT2300A_SetGpio3In(void)//           SET_GPIO_IN(CMT_GPIO3_GPIO)
{
 //set GPIO to input mode
}

unsigned char CMT2300A_ReadGpio1(void)//            READ_GPIO_PIN(CMT_GPIO1_GPIO)
{
 return 0;
}

unsigned char CMT2300A_ReadGpio2(void)//            READ_GPIO_PIN(CMT_GPIO2_GPIO)
{
 return 0;
}

unsigned char CMT2300A_ReadGpio3(void)//            READ_GPIO_PIN(CMT_GPIO3_GPIO)
{
 return 0;
}

*/

