#define TMR_SPI 149 //12MHz timer period for 80000KHz SPI clock rate
#define BLOCKLENGTH  50 //32bit words in SPI block 


//----------------------SPP interface--------------------------------------------


//initialize SPP module on start
void spp_init(void);

//set SPP receiving rate by modem value
void spp_rate(unsigned short nn);

//poll received bits are ready, returns pointer to bits must be process immediately  or 0
unsigned short* spp_rcvd(void);

//poll bits were send, return pointer to required new bits must be set immediately  or 0 
unsigned short* spp_send(void);

void spp_test(signed char r);





