
/*
   -----------------------------------------------------------------------
   Copyright (c) 2001 Dr Brian Gladman <brg@gladman.uk.net>, Worcester, UK

   TERMS

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   This software is provided 'as is' with no guarantees of correctness or
   fitness for purpose.
   -----------------------------------------------------------------------
   */

/* Example of the use of the AES (Rijndael) algorithm for file  */
/* encryption.  Note that this is an example application, it is */
/* not intended for real operational use.  The Command line is: */
/*                                                              */
/* aesxam input_file_name output_file_name [D|E] hexadecimalkey */
/*                                                              */
/* where E gives encryption and D decryption of the input file  */
/* into the output file using the given hexadecimal key string  */
/* The later is a hexadecimal sequence of 32, 48 or 64 digits   */
/* Examples to encrypt or decrypt aes.c into aes.enc are:       */
/*                                                              */
/* aesxam file.c file.enc E 0123456789abcdeffedcba9876543210    */
/*                                                              */
/* aesxam file.enc file2.c D 0123456789abcdeffedcba9876543210   */
/*                                                              */
/* which should return a file 'file2.c' identical to 'file.c'   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <memory.h>
#include <ctype.h>

#include "aes.h"


#include "mem_access.h"

#include "input.h"
#include "config.h"

/*
 * use input data from memory
 */
//extern char _binary_input_asc_start;
//extern char _binary_input_asc_end;
//extern char _binary_input_asc_size;
//char *input=&_binary_input_asc_start;
//int input_size=0;




/* A Pseudo Random Number Generator (PRNG) used for the     */
/* Initialisation Vector. The PRNG is George Marsaglia's    */
/* Multiply-With-Carry (MWC) PRNG that concatenates two     */
/* 16-bit MWC generators:                                   */
/*     x(n)=36969 * x(n-1) + carry mod 2^16                 */ 
/*     y(n)=18000 * y(n-1) + carry mod 2^16                 */
/* to produce a combined PRNG with a period of about 2^60.  */  
/* The Pentium cycle counter is used to initialise it. This */
/* is crude but the IV does not need to be secret.          */

/* void cycles(unsigned long *rtn)     */
/* {                           // read the Pentium Time Stamp Counter */
/*     __asm */
/*     { */
/*     _emit   0x0f            // complete pending operations */
/*     _emit   0xa2 */
/*     _emit   0x0f            // read time stamp counter */
/*     _emit   0x31 */
/*     mov     ebx,rtn */
/*     mov     [ebx],eax */
/*     mov     [ebx+4],edx */
/*     _emit   0x0f            // complete pending operations */
/*     _emit   0xa2 */
/*     } */
/* } */

#define RAND(a,b) (((a = 36969 * (a & 65535) + (a >> 16)) << 16) + (b = 18000 * (b & 65535) + (b >> 16))  )

void fillrand(char *buf, int len)
{   static unsigned long a[2], mt = 1, count = 4;
	static char          r[4];
	int                  i;

	if(mt) { 
		mt = 0; 
		/*cycles(a);*/
		a[0]=0xeaf3;
		a[1]=0x35fe;
	}

	for(i = 0; i < len; ++i)
	{
		if(count == 4)
		{
			*(unsigned long*)r = RAND(a[0], a[1]);
			count = 0;
		}

		buf[i] = r[count++];
	}
}    



//adaption of "encfile(..)"
int encmem(char *dest, aes *ctx, char *src){ 
	char *dPtr = dest;
	char *sPtr = src;
	char            inbuf[16], outbuf[16];
	unsigned long   i=0, l=0;

	fillrand(outbuf, 16);           /* set an IV for CBC mode           */

	/* write the IV to the output       */
	if(!mem_write(outbuf, 16, &dPtr, dest+ INPUT_SIZE +32)){
		printf("couldn't write IV\n");
		return -100;
	}

	fillrand(inbuf, 1);             /* make top 4 bits of a byte random */
	l = 15;                         /* and store the length of the last */
	/* block in the lower 4 bits        */
	inbuf[0] = ((char)INPUT_SIZE & 15) | (inbuf[0] & ~15);

	while(sPtr < input + INPUT_SIZE)               /* loop to encrypt the input file   */
	{                               /* input 1st 16 bytes to buf[1..16] */
		i = mem_read(inbuf+16-l, l, &sPtr, input + INPUT_SIZE);/*  on 1st round byte[0] */
		/* is the length code    */

		if(i < l) break;            /* if end of the input file reached */

		for(i = 0; i < 16; ++i)         /* xor in previous cipher text  */
			inbuf[i] ^= outbuf[i]; 

		encrypt(inbuf, outbuf, ctx);    /* and do the encryption        */

		if(mem_write(outbuf, 16, &dPtr, dest + INPUT_SIZE + 32) != 16){
			printf("encmem | couldn't write '%.16s'\n", outbuf);
			return -100;
		}
		/* in all but first round read 16   */
		l = 16;                     /* bytes into the buffer            */
	}

	/* except for files of length less than two blocks we now have one  */
	/* byte from the previous block and 'i' bytes from the current one  */
	/* to encrypt and 15 - i empty buffer positions. For files of less  */
	/* than two blocks (0 or 1) we have i + 1 bytes and 14 - i empty    */
	/* buffer position to set to zero since the 'count' byte is extra   */

	if(l == 15)                         /* adjust for extra byte in the */
		++i;                            /* in the first block           */

	if(i)                               /* if bytes remain to be output */
	{
		while(i < 16)                   /* clear empty buffer positions */
			inbuf[i++] = 0;

		for(i = 0; i < 16; ++i)         /* xor in previous cipher text  */
			inbuf[i] ^= outbuf[i]; 

		encrypt(inbuf, outbuf, ctx);    /* encrypt and output it        */

		if(mem_write(outbuf, 16, &dPtr, dest + INPUT_SIZE + 32) != 16){
			printf("encmem | couldn't write '%.16s'", outbuf);
			return -100;
		}
	}

	return 0;
}


//adaption of "decfile(..)" with call syntax like "encmem(..)"
//src = fin, dest=fout
//int decmem(char *src, char *dest, aes *ctx, char* ifn, char* ofn)
int decmem(char *src, aes *ctx, char *dest)
{   
	char *dPtr = dest;
	char *sPtr = src;
	char    inbuf1[16], inbuf2[16], outbuf[16], *bp1, *bp2, *tp;
	int     i, l, flen;

	if(mem_read(inbuf1, 16, &sPtr, src + INPUT_SIZE + 32) != 16)  /* read Initialisation Vector   */
	{
		printf("Error reading from input\n");
		return 9;
	}

	i = mem_read(inbuf2, 16, &sPtr, src + INPUT_SIZE + 32);  /* read 1st encrypted file block    */

	if(i && i != 16)
	{
		printf("\nThe input is corrupt");
		return -10;
	}

	decrypt(inbuf2, outbuf, ctx);   /* decrypt it                       */

	for(i = 0; i < 16; ++i)         /* xor with previous input          */
		outbuf[i] ^= inbuf1[i];

	flen = outbuf[0] & 15;  /* recover length of the last block and set */
	l = 15;                 /* the count of valid bytes in block to 15  */                              
	bp1 = inbuf1;           /* set up pointers to two input buffers     */
	bp2 = inbuf2;

	while(1)
	{
		i = mem_read(bp1, 16, &sPtr, src + INPUT_SIZE + 32);     /* read next encrypted block    */

		/* to first input buffer        */
		if(i != 16)         /* no more bytes in input - the decrypted   */
			break;          /* partial final buffer needs to be output  */

		/* if a block has been read the previous block must have been   */
		/* full lnegth so we can now write it out                       */

		if(mem_write(outbuf+16-l, l, &dPtr, dest + INPUT_SIZE) != l){
			printf("decmem | cannout write %.16s\n", outbuf);
			return -100;
		}     

		decrypt(bp1, outbuf, ctx);  /* decrypt the new input block and  */

		for(i = 0; i < 16; ++i)     /* xor it with previous input block */
			outbuf[i] ^= bp2[i];

		/* set byte count to 16 and swap buffer pointers                */

		l = i; tp = bp1, bp1 = bp2, bp2 = tp;
	}

	/* we have now output 16 * n + 15 bytes of the file with any left   */
	/* in outbuf waiting to be output. If x bytes remain to be written, */
	/* we know that (16 * n + x + 15) % 16 = flen, giving x = flen + 1  */
	/* But we must also remember that the first block is offset by one  */
	/* in the buffer - we use the fact that l = 15 rather than 16 here  */  

	l = (l == 15 ? 1 : 0);
	flen += 1 - l;

	if(flen){
		if(mem_write(outbuf+l, flen, &dPtr, dest + INPUT_SIZE) != flen){
			printf("decmem | cannout write %.16s\n", outbuf);
			return -100;
		}     
	}
	return 0;
}



int main(int argc, char *argv[])
{   
	extern void *ECOS_data;
	extern int ECOS_data_size;

	//input_size = (int)&_binary_input_asc_size;

	char *key_hexstr = "1234567890abcdeffedcba09876543211234567890abcdeffedcba0987654321";

	//char *enc_output=(char *)malloc(input_size+32);
	//bzero(enc_output, input_size+32);
	char enc_output[INPUT_SIZE+32] = { 0 };

	//char *dec_output=(char *)malloc(input_size); //to be compared with input
	//bzero(dec_output, input_size);
	char dec_output[INPUT_SIZE] = { 0 };

/*	//safeguards
	if(!enc_output){
		printf("cannot allocate encryption buffer (%dB)", input_size+32);
		goto exit;
	} 
	if(!dec_output){
		printf("cannot allocate decryption buffer (%dB)", input_size);
		goto exit;
	} 
*/
	
	*FI_START = 1;

	char    *cp, ch, key[32];
	int     i=0, by=0, key_len=0, err = 0;
	aes     ctx[1];

	cp = key_hexstr;   /* this is a pointer to the hexadecimal key digits  */
	i = 0;          /* this is a count for the input digits processed   */

	while(i < 64 && *cp)    /* the maximum key length is 32 bytes and   */
	{                       /* hence at most 64 hexadecimal digits      */
		ch = toupper(*cp++);            /* process a hexadecimal */
		if(ch >= '0' && ch <= '9')
			by = (by << 4) + ch - '0';
		else if(ch >= 'A' && ch <= 'F')
			by = (by << 4) + ch - 'A' + 10;
		else                            /* error if not hexadecimal     */
		{
			printf("key must be in hexadecimal notation\n"); 
			err = -2; goto exit;
		}

		/* store a key byte for each pair of hexadecimal digits         */
		if(i++ & 1) 
			key[i / 2 - 1] = by & 0xff; 
	}



	key_len = i / 2;

	//first encrypt into enc_output
	printf("encrypting binary input data... ");
	set_key(key, key_len, enc, ctx);
	err = encmem(enc_output, ctx, input);
	printf("Done.\n");

	//then decrypt into dec_output
	printf("decrypting previously encrypted data... ");
	set_key(key, key_len, dec, ctx);
	err = decmem(enc_output, ctx, dec_output);
	printf("Done.\n");

	//compare input and dec_output
	//printf("diff? -> %d\n", memcmp(input, dec_output, input_size));


	//ECOS_data = dec_output;
	//ECOS_data_size = input_size;
	//ECOS_BENCHMARK_FINISHED();

	*FI_STOP = 1;
	hashResultData(dec_output, INPUT_SIZE);
	*CPU_DONE = 1;

exit:
	;

	return err;
}
