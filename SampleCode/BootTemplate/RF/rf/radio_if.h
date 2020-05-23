#include "gpc.h"

void Mcu_Init(void);
unsigned char Radio_On(unsigned char on);

unsigned char Radio_Send_VarLen(unsigned char* pBuf, unsigned char len);
unsigned char Radio_Recv_VarLen(unsigned char* pBuf);

unsigned char Radio_Get_RX_ID(void);
void Radio_Set_TX_ID(unsigned char d);

unsigned char radio_mode(unsigned char mode); //change TX/RX mode, return 1 if OK
void radio_init(unsigned char on); //change work/pd mode
unsigned char radio_ch(short ch); //if ch positive set radio channel< returns carrent radio channel
unsigned char radio_cmd(unsigned char* cmd); //set radio ch by by string command, return length of answer in cmd

