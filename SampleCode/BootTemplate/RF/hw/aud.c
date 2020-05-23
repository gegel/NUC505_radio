/*\
 Audio recording playing with internal audio codec
 Sampling frequency is 8 KHz, Mono, frame size is 160 short 16 bit PCM
 Internal I2S uses two double buffers (2 frames size each): for recording and for playing
 DMA uses for work with buffers
 Sampling frequency the same for recording and playing so only RX DMA interupt is used
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NUC505Series.h"

#include "aud.h"

//internal procedures
void I2S_Init(void); 
void demo_MIC0(void);


//------------------Audio grabbing with ADC / playing with DAC  using DMA------------------------------

//Audio global variables
int16_t aud_in[2][FRAMELENGTH];   //PCM recording double buffer
int16_t aud_out[2][FRAMELENGTH];  //PCM playing double buffer


volatile unsigned char aud_sp_flg=0; //playing compleet flag
volatile unsigned char aud_sr_flg=0; //recording rady flag
volatile static unsigned char aud_s_flag1=0; //double buffer  part flag



//Initilaize audio/uart module
void aud_init(void)
{
	
	//clear global values
	aud_sp_flg=0;
	
	//clear buffers
	memset(aud_in, 0, sizeof(aud_in));
	memset(aud_out, 0, sizeof(aud_out));
	
	// Init I2S, IP clock and multi-function I/O 
    I2S_Init();
	
	 //setup I2S for work with internal codec
	 I2S_Open(I2S, I2S_MODE_MASTER, 8000, I2S_DATABIT_16, I2S_MONO, I2S_FORMAT_I2S, I2S_ENABLE_INTERNAL_CODEC);


    // Open MCLK
    I2S_EnableMCLK(I2S, 8000*256);

	  //set I2S FOFO's trashold levels
    I2S_SET_TX_TH_LEVEL(I2S, I2S_FIFO_TX_LEVEL_WORD_15);
    I2S_SET_RX_TH_LEVEL(I2S, I2S_FIFO_RX_LEVEL_WORD_16);
		
		//setup recording and playing DMA
		I2S_SET_TXDMA_STADDR(I2S, (uint32_t) &aud_out[0]);                                // Tx Start Address
    I2S_SET_TXDMA_THADDR(I2S, (uint32_t) &aud_out[0][FRAMELENGTH-2]);  // Tx Threshold Address
    I2S_SET_TXDMA_EADDR( I2S, (uint32_t) &aud_out[1][FRAMELENGTH-2]);  // Tx End Address

    I2S_SET_RXDMA_STADDR(I2S, (uint32_t) &aud_in[0]);                                // Rx Start Address
    I2S_SET_RXDMA_THADDR(I2S, (uint32_t) &aud_in[0][FRAMELENGTH-2]);  // Rx Threshold Address
    I2S_SET_RXDMA_EADDR( I2S, (uint32_t) &aud_in[1][FRAMELENGTH-2]);  // Rx End Address


    // Open Rx Dma Enable
    I2S_ENABLE_RXDMA(I2S);

    //setup internal codec for work with Mike
    demo_MIC0();
				
		// Clear Interrupt Status
    I2S_CLR_INT_FLAG(I2S, I2S_STATUS_LZCIF_Msk|I2S_STATUS_RZCIF_Msk|I2S_STATUS_TXOVIF_Msk|I2S_STATUS_TXUDIF_Msk|I2S_STATUS_RXOVIF_Msk|I2S_STATUS_RXUDIF_Msk|I2S_STATUS_TDMATIF_Msk|I2S_STATUS_TDMAEIF_Msk|I2S_STATUS_RDMATIF_Msk|I2S_STATUS_RDMAEIF_Msk);
    // Recording Enable
    I2S_ENABLE_RX(I2S);
    //DMA Interupt enable
    NVIC_EnableIRQ(I2S_IRQn);
    I2S_EnableInt(I2S, (I2S_IEN_RDMATIEN_Msk|I2S_IEN_RDMAEIEN_Msk));
	
}

//setup playing volume (in range 0-31)
void aud_vol(unsigned char vol)
{
	if(vol>0x1F) vol=0x1F; //restrict 0-31
	vol=0x1F-vol; 
	I2S_SET_INTERNAL_CODEC(I2S, 0x08, vol);    // Volume headphone of Left channel
  I2S_SET_INTERNAL_CODEC(I2S, 0x09, vol);    // Volume headphone of Right channel
}


//check we need next frame for playing: return pointer to buffer where frame must be palced or 0
short* aud_play(void)
{
  short* ptr=0;
	unsigned char b;
	
	if(aud_sp_flg) //check frame played, needs new
	{
	 b=aud_sp_flg; //capture flag
   aud_sp_flg=0; //clear flag	
	 aud_sr_flg=b; //set flag for racording ready (synchronized with playing)
   b--;	//getc number of double buffer part is empty	
	 ptr=aud_out[b]; //return pointer to empty playing buffer part 	 ready for new frame for next palying
	}
	return ptr;	//return pointer to buffer for new frame for playing or 0	
}


//chech we have grabbed frame: returns pointer to frame or 0
short* aud_grab(void)
{
	short* ptr=0;
	unsigned char b;
	
	if(aud_sr_flg) //check recording DMA is compleet
	{
	 b=aud_sr_flg; //capure flag
	 aud_sr_flg=0; //clear flag
	 b--;	//getc number of recording buufer part with frame was recently recorded
	 ptr=aud_in[b];	//return pointers with recorded frame
	}
	
	return ptr; //returns pointer to recorded frame or 0 if frame not ready yet
}

//-----------------------------------------------------------------------------------
//Internal procedures
//-----------------------------------------------------------------------------------

//for delay
void SysTick_Handler(void)
{

}

//initialize I2S 
void I2S_Init(void)
{
    // Enable I2S Module clock
    CLK_EnableModuleClock(I2S_MODULE);  
    
//------------------------------------------------------	
	  /* I2S module clock from APLL */
    // APLL = 49152031Hz
    /////////CLK_SET_APLL(CLK_APLL_49152031); //this is not exacly 8Khz!!!!
	 
	
    // I2S = 49152031Hz / (0+1) for 8k, 12k, 16k, 24k, 32k, 48k, and 96k sampling rate

    // APLL = 45158425Hz
    // CLK_SET_APLL(CLK_APLL_45158425);
    // I2S = 45158425Hz / (0+1) for 11025, 22050, and 44100 sampling rate
//----------------------------------------------------------
	
	  //set audio PLL clock 51.2 MHz for exactly 8KHz sampling rate
	  CLK_SET_APLL(0x8FC2);  //APLL_freq=51200000. div must be 0x31   51200000/2*(49+1)= 512000
	
	  //set I2S clock from audio PLL
    CLK_SetModuleClock(I2S_MODULE, CLK_I2S_SRC_APLL, 0);
    /* Reset IP */
    SYS_ResetModule(I2S_RST);
}





/*
void demo_MIC0(void)
{
    uint32_t i;

    // IIC Configure Step without PLL: 
    // Add MCLK(256*Fs) in. 

    I2S_SET_INTERNAL_CODEC(I2S, 0x08, 0x1F);    // Mute headphone of Left channel
    I2S_SET_INTERNAL_CODEC(I2S, 0x09, 0x1F);    // Mute headphone of Right channel
    I2S_SET_INTERNAL_CODEC(I2S, 0x10, 0x0F);    //Mute the ADC Left channel volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x11, 0x0F);    //Mute the ADC Right channel volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x12, 0x0F);    //Mute the ADC Side tone volume

    I2S_SET_INTERNAL_CODEC(I2S, 0x02, 0xC0);    //Set CODEC slave

    I2S_SET_INTERNAL_CODEC(I2S, 0x01, 0x90);    //DAC Digital Part Enable   0x80-stereo, 0x90 left mono
    I2S_SET_INTERNAL_CODEC(I2S, 0x0F, 0xC0);    //Enable Analog Part power
    I2S_SET_INTERNAL_CODEC(I2S, 0x0E, 0x02);    //ADC input select MIC0

    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xF3);    //Headphone power disable, pump enable
	I2S_SET_INTERNAL_CODEC(I2S, 0x0D, 0x31);    //Biasing enable: biggest bias (0x40 for disable)
    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xE3);  //pump disable
    for (i=0; i < 15; i++)  //Delay 1.5s~2.5s
        CLK_SysTickDelay(100000);
				I2S_SET_INTERNAL_CODEC(I2S, 0x0A, 0x09); //set output mixer: DAC signals
    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xD0);   //!!!set headphon 1K impenance
    I2S_SET_INTERNAL_CODEC(I2S, 0x00, 0xC0);    //ADC digital enabled  0xD0-stereo, 0xC0 - left mono !!!!!
    CLK_SysTickDelay(100000);   //Delay 100mS

    I2S_SET_INTERNAL_CODEC(I2S, 0x08, 0x02);    //Un-mute Headphone and set volume (0-0dB 0x1E minus 60dB) 2dB/1
    I2S_SET_INTERNAL_CODEC(I2S, 0x09, 0x02);    //Un-mute Headphone and set volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x10, 0x18);    //Un-Mute the ADC Left channel volume  (0x10 0dB 0x1E +22.4dB), 0x10 is mike boost +20dB
    I2S_SET_INTERNAL_CODEC(I2S, 0x11, 0x08);    //Un-Mute the ADC Right channel volume (0x00 0dB, 0x0E +22.4dB) 1.6dB/1
//      I2S_SET_INTERNAL_CODEC(I2S, 0x12, 0x00);    //Un-Mute the ADC Side tone volume

    // If Fs is changed, please Mute Headphone First and soft reset digital part after MCLK is stable. 
}

*/


//ORIGINAL
void demo_MIC0(void)
{
    uint32_t i;

    // IIC Configure Step without PLL: 
    // Add MCLK(256*Fs) in. 

    I2S_SET_INTERNAL_CODEC(I2S, 0x08, 0x1F);    // Mute headphone of Left channel
    I2S_SET_INTERNAL_CODEC(I2S, 0x09, 0x1F);    // Mute headphone of Right channel
    I2S_SET_INTERNAL_CODEC(I2S, 0x10, 0x0F);    //Mute the ADC Left channel volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x11, 0x0F);    //Mute the ADC Right channel volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x12, 0x0F);    //Mute the ADC Side tone volume

    I2S_SET_INTERNAL_CODEC(I2S, 0x02, 0xC0);    //Set CODEC slave

    I2S_SET_INTERNAL_CODEC(I2S, 0x01, 0x80);    //Digital Part Enable
    I2S_SET_INTERNAL_CODEC(I2S, 0x0F, 0xC0);    //Enable Analog Part
    I2S_SET_INTERNAL_CODEC(I2S, 0x0E, 0x02);    //ADC input select MIC0

    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xF3);    //Analog Part Enable
    I2S_SET_INTERNAL_CODEC(I2S, 0x0D, 0x31);    //Biasing enable
    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xE3);
    for (i=0; i < 15; i++)  //Delay 1.5s~2.5s
        CLK_SysTickDelay(100000);
    I2S_SET_INTERNAL_CODEC(I2S, 0x0A, 0x09);
    I2S_SET_INTERNAL_CODEC(I2S, 0x0B, 0xF0);
    I2S_SET_INTERNAL_CODEC(I2S, 0x00, 0xD0);    //ADC digital enabled
    CLK_SysTickDelay(100000);   //Delay 100mS

    I2S_SET_INTERNAL_CODEC(I2S, 0x08, 0x06);    //Un-mute Headphone and set volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x09, 0x06);    //Un-mute Headphone and set volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x10, 0x18);    //Un-Mute the ADC Left channel volume
    I2S_SET_INTERNAL_CODEC(I2S, 0x11, 0x08);    //Un-Mute the ADC Right channel volume
//      I2S_SET_INTERNAL_CODEC(I2S, 0x12, 0x00);    //Un-Mute the ADC Side tone volume

    // If Fs is changed, please Mute Headphone First and soft reset digital part after MCLK is stable. 
}




//I2S interupts (DMA recorded frame is ready)
void I2S_IRQHandler(void)
{
    uint32_t u32I2SIntFlag;

	  //get interupt flags
    u32I2SIntFlag = I2S_GET_INT_FLAG(I2S, (I2S_STATUS_RDMATIF_Msk | I2S_STATUS_RDMAEIF_Msk));

    //Frame ready in first part of double buffer 
    if (u32I2SIntFlag & I2S_STATUS_RDMATIF_Msk)
    {
        I2S_CLR_INT_FLAG(I2S, I2S_STATUS_RDMATIF_Msk); //vlar interupt falg
			  aud_sp_flg=1; //set flag for frame ready in first part
    }
		 //Frame ready in second part of double buffer
    else if (u32I2SIntFlag & I2S_STATUS_RDMAEIF_Msk)
    {
        I2S_CLR_INT_FLAG(I2S, I2S_STATUS_RDMAEIF_Msk); //clear interupt flag
        aud_sp_flg=2; //set flag of frame ready in second part
			  if ( aud_s_flag1 == 0 ) //once on start enable playing DMA for certain synchronize TX with RX
        {
            aud_s_flag1 = 1; //set once flag
            I2S_ENABLE_TXDMA(I2S); //enable playing DMA
            I2S_ENABLE_TX(I2S); //start playing
        }
    }
}


