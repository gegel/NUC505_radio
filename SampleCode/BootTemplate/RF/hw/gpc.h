#ifndef __HW_LAYER_H
#define __HW_LAYER_H

#include "typedefs.h"

#ifdef __cplusplus 
extern "C" { 
#endif



// ************************************************************************
//  The following need to be modified by user
//  ************************************************************************ 




//===========================Timings================================

#define SPI_DELAY_COUNT_MS 150  //loop itteration for ius delay depends by HW speed

void cmt_spi3_delay(void); //delay 1mS by loop
void system_delay_us (unsigned short us); //delay in us by timer


//===============================GPIO===================================

void cmt_gpio_init(void);
void cmt_gpio_pd(void);
unsigned char cmt_read_p3(void);


void cmtvoid(void);

//---------------IO_RADIO------------------------
//VSS P
//VDD P52
//CS PB8 P50 out
//CLK PB7 P48 out
//DIO PB7 P46 in/out
//FB PB9 P51 out
//SPI0 MISO P46

//----------------IO BOARD------------------------
//TST PC3 P11 out  onboard led to vcc
//UART1 RX PA9  P29
//UART1 TX PA8  P28
//MIC PA10 P30
//SPC PA11 P31




//-----------------IO AUDIO-----------------------
//VDD P81
//VSS PC14 P85
//BTN PC12 P83 in
//BTN_GND PC14 P85
//LED PC13 P84 out
//GPIO PE5 P54 in
//SPI1 MOSI P80





//pin direction not changed because startup is used
#define cmt_spi3_csb_out()      cmtvoid()
#define cmt_spi3_csb_out()      cmtvoid()
#define cmt_spi3_fcsb_out()     cmtvoid()
#define cmt_spi3_sclk_out()     cmtvoid()

//pin direction not changed because bidirectional mode is used
//#define cmt_spi3_sdio_out()     PE->MODE = (PE->MODE | (0x1 << (BIT1 <<1)))
//#define cmt_spi3_sdio_in()      PE->MODE = (PE->MODE & ~(0x3 << (BIT1 << 1)))

#define cmt_spi3_sdio_out()     PB->MODE = (PB->MODE | (1<<6))
#define cmt_spi3_sdio_in()      PB->MODE = (PB->MODE & ~(1<<6))

//set spi outputs
#define cmt_spi3_sdio_1()       PB6_DOUT=1
#define cmt_spi3_sdio_0()       PB6_DOUT=0

#define cmt_spi3_sclk_1()       PB7_DOUT=1
#define cmt_spi3_sclk_0()       PB7_DOUT=0

#define cmt_spi3_csb_1()        PB8_DOUT=1
#define cmt_spi3_csb_0()        PB8_DOUT=0

#define cmt_spi3_fcsb_1()       PB9_DOUT=1
#define cmt_spi3_fcsb_0()       PB9_DOUT=0
   

//read SPI data
#define cmt_spi3_sdio_read()    (PB6_PIN)

//read gpion input controls
#define CMT2300A_ReadGpio1()    (PC11_PIN)
#define CMT2300A_ReadGpio2()    (PC11_PIN)
#define CMT2300A_ReadGpio3()    (PC11_PIN)

//pin directing not changed because startup is used
#define CMT2300A_SetGpio1In()   cmtvoid()
#define CMT2300A_SetGpio2In()   cmtvoid()
#define CMT2300A_SetGpio3In()   cmtvoid()

//user control 
#define CMT2300A_ReadBTN()     (PC12_PIN)
#define CMT2300A_ReadGPI()     (PC11_PIN)
#define CMT2300A_ReadTST()     (PC6_PIN)

#define CMT2300A_WriteTST_1()  PC3_DOUT=0
#define CMT2300A_WriteTST_0()  PC3_DOUT=1

#define CMT2300A_WriteLED_1()  PC13_DOUT=0
#define CMT2300A_WriteLED_0()  PC13_DOUT=1

#define CMT2300A_WriteMIC_1()  PA10_DOUT=1
#define CMT2300A_WriteMIC_0()  PA10_DOUT=0

#define CMT2300A_WriteSPC_1()  PA11_DOUT=1
#define CMT2300A_WriteSPC_0()  PA11_DOUT=0






/*
void cmt_spi3_csb_out(void);
void cmt_spi3_fcsb_out(void);
void cmt_spi3_sclk_out(void);
void cmt_spi3_sdio_out(void);
void cmt_spi3_sdio_in(void);

void cmt_spi3_csb_1(void);
void cmt_spi3_csb_0(void);

void cmt_spi3_fcsb_1(void);
void cmt_spi3_fcsb_0(void);
    
void cmt_spi3_sclk_1(void);
void cmt_spi3_sclk_0(void);

void cmt_spi3_sdio_1(void);
void cmt_spi3_sdio_0(void);
unsigned char cmt_spi3_sdio_read(void);
*/

//====================================================

/*
void CMT2300A_SetGpio1In(void); //           SET_GPIO_IN(CMT_GPIO1_GPIO)
void CMT2300A_SetGpio2In(void); //           SET_GPIO_IN(CMT_GPIO2_GPIO)
void CMT2300A_SetGpio3In(void); //           SET_GPIO_IN(CMT_GPIO3_GPIO)
unsigned char CMT2300A_ReadGpio1(void);//            READ_GPIO_PIN(CMT_GPIO1_GPIO)
unsigned char CMT2300A_ReadGpio2(void);//            READ_GPIO_PIN(CMT_GPIO2_GPIO)
unsigned char CMT2300A_ReadGpio3(void);//            READ_GPIO_PIN(CMT_GPIO3_GPIO)
*/


#ifdef __cplusplus
} 
#endif

#endif


