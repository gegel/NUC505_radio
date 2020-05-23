#define FRAMELENGTH  160 //number of PCM sampling in AMR frame


void aud_init(void);
void aud_vol(unsigned char vol);
void aud_set(short rate);
short* aud_grab(void);
short* aud_play(void);




