#define COMBDR 115200 //boudrate of uart
#define COMTOUT 100 //number of bits in pause for set end of block
#define COM_RBUF 64  //length of UART RX buffer
#define COM_TBUF 64  //length of UART TX buffer



void com_init(void); //initialize uart
short com_write(unsigned char* data, short len); //write len bytes in data to com port
short com_poll(void); //poll for receiving data
unsigned char* com_read(void); //get pointer to received data

