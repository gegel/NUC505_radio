/*
 * THE FOLLOWING FIRMWARE IS PROVIDED: (1) "AS IS" WITH NO WARRANTY; AND
 * (2)TO ENABLE ACCESS TO CODING INFORMATION TO GUIDE AND FACILITATE CUSTOMER.
 * CONSEQUENTLY, CMOSTEK SHALL NOT BE HELD LIABLE FOR ANY DIRECT, INDIRECT OR
 * CONSEQUENTIAL DAMAGES WITH RESPECT TO ANY CLAIMS ARISING FROM THE CONTENT
 * OF SUCH FIRMWARE AND/OR THE USE MADE BY CUSTOMERS OF THE CODING INFORMATION
 * CONTAINED HEREIN IN CONNECTION WITH THEIR PRODUCTS.
 *
 * Copyright (C) CMOSTEK SZ.
 */

/*!
 * @file    radio.c
 * @brief   Generic radio handlers
 *
 * @version 1.2
 * @date    Jul 17 2017
 * @author  CMOSTEK R@D
 */
 
#include "radio.h"
#include "cmt2300a_params.h" //register's values for CMT2300A work mode
#include "pd_cmt2300a_params.h" //register's values for CMT2300A power save mode


#include <string.h>

//Initialization of CMT2300A in work mode (on=1) or speep mode (on=0)
void RF_Init(unsigned char on)
{
    u8 tmp;
    
		CMT2300A_InitGpio();
		CMT2300A_Init();
    
    /* Config registers */
	  if(on)
		{  //set registers for work mode (direct mode without CDR and no mute GPIO)
     CMT2300A_ConfigRegBank(CMT2300A_CMT_BANK_ADDR       , g_cmt2300aCmtBank       , CMT2300A_CMT_BANK_SIZE       );
     CMT2300A_ConfigRegBank(CMT2300A_SYSTEM_BANK_ADDR    , g_cmt2300aSystemBank    , CMT2300A_SYSTEM_BANK_SIZE    );
     CMT2300A_ConfigRegBank(CMT2300A_FREQUENCY_BANK_ADDR , g_cmt2300aFrequencyBank , CMT2300A_FREQUENCY_BANK_SIZE );
     CMT2300A_ConfigRegBank(CMT2300A_DATA_RATE_BANK_ADDR , g_cmt2300aDataRateBank  , CMT2300A_DATA_RATE_BANK_SIZE );
     CMT2300A_ConfigRegBank(CMT2300A_BASEBAND_BANK_ADDR  , g_cmt2300aBasebandBank  , CMT2300A_BASEBAND_BANK_SIZE  );
     CMT2300A_ConfigRegBank(CMT2300A_TX_BANK_ADDR        , g_cmt2300aTxBank        , CMT2300A_TX_BANK_SIZE        );
    }
		else
		{ //set registers for power down mode (periodic receiving then sleep mode with muting GPIO)
		 CMT2300A_ConfigRegBank(CMT2300A_CMT_BANK_ADDR       , d_cmt2300aCmtBank       , CMT2300A_CMT_BANK_SIZE       );
     CMT2300A_ConfigRegBank(CMT2300A_SYSTEM_BANK_ADDR    , d_cmt2300aSystemBank    , CMT2300A_SYSTEM_BANK_SIZE    );
     CMT2300A_ConfigRegBank(CMT2300A_FREQUENCY_BANK_ADDR , d_cmt2300aFrequencyBank , CMT2300A_FREQUENCY_BANK_SIZE );
     CMT2300A_ConfigRegBank(CMT2300A_DATA_RATE_BANK_ADDR , d_cmt2300aDataRateBank  , CMT2300A_DATA_RATE_BANK_SIZE );
     CMT2300A_ConfigRegBank(CMT2300A_BASEBAND_BANK_ADDR  , d_cmt2300aBasebandBank  , CMT2300A_BASEBAND_BANK_SIZE  );
     CMT2300A_ConfigRegBank(CMT2300A_TX_BANK_ADDR        , d_cmt2300aTxBank        , CMT2300A_TX_BANK_SIZE        );	
		}	
			
	
    // xosc_aac_code[2:0] = 2
    tmp = (~0x07) & CMT2300A_ReadReg(CMT2300A_CUS_CMT10);
    CMT2300A_WriteReg(CMT2300A_CUS_CMT10, tmp|0x02);
     
}



//configuration of CMT2300A applied after initialization
//cc is number of freq. channel (0-255)
void RF_Config(unsigned char cc)
{

	  CMT2300A_ConfigGpio( CMT2300A_GPIO1_SEL_DOUT |   // INT1 > GPIO1 
	                       CMT2300A_GPIO2_SEL_INT1 |   // INT2 > GPIO3 
                         CMT2300A_GPIO3_SEL_DCLK 	                      
	                      );
   
	  CMT2300A_ConfigInterrupt( CMT2300A_INT_SEL_RSSI_VLD,  // GPIO1 > SYNC_OK 
	                            CMT2300A_INT_SEL_PKT_DONE  // GPIO3 > PKT_DONE
	                          );
      

    // Enable interrupt 
    //  CMT2300A_EnableInterrupt(
		//		CMT2300A_MASK_PKT_DONE_EN | 
    //    CMT2300A_MASK_PREAM_OK_EN |
    //    CMT2300A_MASK_SYNC_OK_EN  |
		//		CMT2300A_MASK_TX_DONE_EN
//        CMT2300A_MASK_NODE_OK_EN  |
//        CMT2300A_MASK_CRC_OK_EN   |
//         
     //   );

    
    /* Disable low frequency OSC calibration */
    CMT2300A_EnableLfosc(FALSE);
    
    /* Use a single 64-byte FIFO for either Tx or Rx */
    //CMT2300A_EnableFifoMerge(TRUE);
    
    //CMT2300A_SetFifoThreshold(16); // FIFO_TH
    
    /* This is optional, only needed when using Rx fast frequency hopping */
    /* See AN142 and AN197 for details */
    //CMT2300A_SetAfcOvfTh(0x27);
    
    /* Go to sleep for configuration to take effect */
		
		CMT2300A_SetFrequencyStep(4); //2.5KHz*n=10 kHz
		CMT2300A_SetFrequencyChannel(cc);
		
    CMT2300A_GoSleep();
}

