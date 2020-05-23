#define FRAMELENGTH  160 //number of PCM sampling in AMR frame
#define PKTLENGTH 12     //number of bytes in UART data packet for AMR mode 0
#define FIFOCNT 3       //maximal number of encoded speech frames in FIFO
#define PRC_RC 20       //rate correction coefficient for playing
#define PRC_ONE 16384   //DSP one for resampler

//initinig
void amr_init(void); //initialize amr codec
void prc_init(void); //inittialize processing
void prc_vt(unsigned char on); //activate tx test

//playing rate correction
void prc_set(void); //set internal playing rate correction value at a moment (thread safe)
signed char prc_get(void); //get internal playing rate correction value (-1 to 1)
unsigned char prc_tst(void); //get current fifo difference value (1 to 3) at this moment

//pcm processing
short prc_enc(short* pcm); //encode recorded speech frame to fifo, returns VAD flag
void prc_dec(short* pcm); //decode speech frame for playing from fifo without resampling
void prc_decr(short* pcm); //decode speech frame for playing from fifo with resampling

//bits processing
unsigned short prc_dmd(unsigned short* bits); //demodulate received bits to fifo
void prc_mdm(unsigned short* bits); //modulate bits for sending from fifo



