/*
 Entering CPU to power down mode with wake-up by PTT button press or radio activity
 Radio must be preiosly set to periodiv RX mode with muting output (20ms receive/500 ms sleep)
 On carrier with radio appropriated frequency jumps radio outut wiil be pulse to 1
 After wake-up by radio pulse MC still monitore radio output during 1 sec
 If there are less then 3 pulses this FDK casual packet (not a carrier) and MC back to pawer down
 If there is carrier or PTT button  was pressed MC set 0 to pin connected with RESET and restart work
*/


#include "NUC505Series.h"
#include "adc.h"
#include "gpio.h"
#include "i2s.h"
#include "gpc.h" //IO procedures
#include "com.h"
#include "pdc.h"

void DisableIPs(void);
void PowerDownFunction(void);

#define PDC_CNT (1000) //radio output testing loop iterations
#define PDC_FIX 10 //number of pulses in test for carrier present
#define PDC_DEL 75 //elementary delay 

#pragma push // Save existing optimization level
  #pragma O0   // Optimization level now O0 
  
volatile int pdc_cnt=0; //radio out testing  loop
volatile int pdc_fix=0; //radio out pulses counter
volatile short pdc_del=0; //delay for testing loop
int ii;

//entering power down mode (exiting over chip reset)
void pdc_init(void)
{
  //set gpio modes for power down 
	cmt_gpio_pd(); //set most gpio as inputs for low power

	//disable power of system periferal modules
	DisableIPs(); //disable system power, down CPU clock to 750 KHz 

	//setup wake-up source
	GPIO_SetMode(PB, BIT1, GPIO_MODE_INPUT); //set button input
  GPIO_SetPullMode(PB, BIT1, GPIO_PULL_UP_EN); //enable pullup
  GPIO_EnableInt(PB, 1, GPIO_INT_FALLING); //enable interupt on press (level 1->0 )
	
	GPIO_SetMode(PB, BIT0, GPIO_MODE_INPUT); //set radio input (0 in no carrier)  
	//GPIO_SetPullMode(PB, BIT0, GPIO_PULL_UP_EN);
  //GPIO_EnableInt(PB, 0, GPIO_INT_RISING); 
	
	GPIO_ENABLE_WAKE_UP(GPIO_INTCTL_EINT0_WKEN_MASK); //enable wake-up from EINT0
  SYS_ENABLE_WAKEUP(SYS_WAKEUP_GPIOWE_Msk); //enable wakeup from GPIO

	GPIO_SetIntGroup(PB, 0, GPIO_INTSRCGP_EINT0); //set radio interupt as INT0
	GPIO_SetIntGroup(PB, 1, GPIO_INTSRCGP_EINT0); //set button interupt as INT0
  NVIC_EnableIRQ(EINT0_IRQn); //enable IRQ for INT0
	
	
	
pde: //entering power down

  while(PB0_PIN) {;} //wait radio is 0
	GPIO_EnableInt(PB, 0, GPIO_INT_RISING); //enable radio interupts
  
			
	PowerDownFunction(); //entering power down, wait interupt from button or radio
  
 
	GPIO_DisableInt(PB, 0); //disable radio interupts
		
				
pdt: //wake up
		
		//check button press: immediately reset for work	
	if(PB1_PIN==0) goto pdr;	
		
	CMT2300A_WriteTST_1(); //on onboard LED (test only!)
		
  //poll radio out pin in loop		
  pdc_cnt=PDC_CNT; pdc_fix=0; //set time conter and clear test value
	while(pdc_cnt--) 
	{
		pdc_del=PDC_DEL; while(pdc_del--){;}
		if(PB0_PIN) pdc_fix++; //count cases gpio=1 during specified time
	}
 
	CMT2300A_WriteTST_0(); //off onboard led (test only!)
	
	//decice is there carrier or not
  if(!pdc_fix) goto pde; //no more activity on gpio: back to power down
  else if(pdc_fix<PDC_FIX) goto pdt; //middle activity: try check once more
	
pdr:   
	//restart chip work entering work mode
	GPIO_SetMode(PC, BIT0, GPIO_MODE_OUTPUT);
	PC0_DOUT=0; //set 0 on GPIO externally connected with chip's hardware RESET pin 
	while(1){;} //SYS_ResetChip();
	
}

#pragma pop  //Restore original optimization level


//INT0 hundler
void EINT0_IRQHandler(void)
{
    int32_t i;
	  	
	 i = GPIO_GET_INT_FLAG(PB, BIT0); //get radio interupt flag
	 if(i) GPIO_CLR_INT_FLAG(PB, BIT0); //clear radio interupt flag
	 
   i = GPIO_GET_INT_FLAG(PB, BIT1); //get button interupt flag
	 if(i) GPIO_CLR_INT_FLAG(PB, BIT1);	//clearerr button interupt flag
}


//disable USB power
void USB_Device_Phy_Disable(void)
{
    CLK_SetModuleClock(USBD_MODULE, CLK_USBD_SRC_EXT, 0);
    CLK_EnableModuleClock(USBD_MODULE);
    CLK_SysTickDelay(1000);
    USBD_DISABLE_PHY();
    USBD_DISABLE_USB();
    CLK_DisableModuleClock(USBD_MODULE);
}

//disable ADC power
void ADC_Disable(void)
{
    ADC_T* adc=0;
    CLK_SetModuleClock(ADC_MODULE, CLK_ADC_SRC_EXT, 1);
    CLK_EnableModuleClock(ADC_MODULE);
    CLK_SysTickDelay(10);
    ADC_Close(adc);
    CLK_SysTickDelay(10);
    //CLK_DisableModuleClock(ADC_MODULE);
}

//disable I2S and codec power
void I2S_Disable(void)
{
    CLK->APBCLK = CLK->APBCLK | 0x4000;         /* Enable I2S clock */
    CLK_EnableModuleClock(I2S_MODULE);
    CLK_SysTickDelay(10);
    I2S_Open(I2S, I2S_MODE_MASTER, 48000, I2S_DATABIT_32, I2S_STEREO, I2S_FORMAT_I2S, I2S_ENABLE_INTERNAL_CODEC);
    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xFF);            // DAC Power Off
    I2S_SET_INTERNAL_CODEC(I2S, 0x0F, 0xFF);            // ADC Power Off

    //I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xF0);          // DAC Power On
    //I2S_SET_INTERNAL_CODEC(I2S, 0x0F, 0xE0);          // ADC Power On
    CLK_SysTickDelay(10);
    //CLK_DisableModuleClock(I2S_MODULE);
}

//disable LVD detector power
void LVD_Disable(void)
{
    //SYS->LVDCTL &= (~0x4);
    SYS_DISABLE_LVR();
    SYS_DISABLE_POR();
}




//Set Power Conspumtion  
void DisableIPs(void)
{
    volatile uint32_t i;

    /* Power consumption measure */
    CLK_SET_AHBCLK(0xFFFFFFFF);
    CLK_SET_APBCLK(0xFFFFFFFF);

    /* USB Device phy suspend */
    USB_Device_Phy_Disable();

    /* ADC */
    ADC_Disable();

    /* External DAC */
    I2S_Disable();

    /* LVD */
    LVD_Disable();


    /* Core clock from external crystal clock divided by 16  is 750KHz*/
    CLK_SetHCLK(CLK_HCLK_SRC_EXT, 15);
    /* Update System Core Clock */
    SystemCoreClockUpdate();

    /* APLL power down */
    CLK_APLL_ENABLE_POWERDOWN();
    /* PLL power down */
    CLK_PLL_ENABLE_POWERDOWN();

    /* Internal embedded SPI Flash MISO/MOSI pins pull enable */
    SYS_SET_EMBEDDED_SPIFLASH_PULL(0x12);

    /* Pre-scalar counter from 0 ~ 0xFFFF */
    CLK_SET_WAKEUP_PRESCALAR    (0x1000);

    /* Wake up time is about 45ms */
    CLK_ENABLE_WAKEUP_PRESCALAR();
}

//entering CPU to Power Down mode
void PowerDownFunction(void)
{
	
	SPIM_ENABLE_IO_MODE(SPIM, SPIM_CTL0_BITMODE_STAN, 0);
    /* Enter to Power-down mode */
  
    CLK_PowerDown();
    SPIM_ENABLE_DMM_MODE(SPIM, SPIM_CTL0_CMDCODE_FAST_READ, 0);
    // Wake up 

   /// CLK_SetCoreClock(96000000);
    // Update System Core Clock 
   /// SystemCoreClockUpdate();
}
