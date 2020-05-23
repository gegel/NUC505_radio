
void mdm_init(void);

//modulate 12 bytes in data to 100 words in out
void mdm_modulate(unsigned char* data, unsigned short* out);


//demodulate 100 words to 12 bytes of data,
//frame number (0-15), bit lag (0-7), frame lag(0-99) and sync flag (0 or 1)
//returns rate correction value (-6 to 6)
unsigned short mdm_demodulate(unsigned short* in, unsigned char* data);


