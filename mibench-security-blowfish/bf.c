#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "blowfish.h"

#include "mem_access.h"
#include "config.h"
#include "result.h"

/*
 * use input data from memory
 */

#include "input.h"

	int
main(int argc, char *argv[])
{
	BF_KEY key;
	unsigned char ukey[8];
	unsigned char indata[40],outdata[40],ivec[8];
	int num;
	int by=0,i=0;
	char *cp,ch;

	//	extern void *ECOS_data;
	//	extern int ECOS_data_size;

	char *t_src, *t_dest;
	//input_size = (int)&_binary_input_asc_size;
	char *key_hexstr = "1234567890abcdeffedcba0987654321";

	//char *enc_output = (char *)malloc(input_size);
	//bzero(enc_output, input_size);
	char enc_output[INPUT_SIZE] = { 0 };

//	char *dec_output = (char *)malloc(input_size);
//	bzero(dec_output, input_size);
	char dec_output[INPUT_SIZE] = { 0 };

	//safeguards
/*	if(!enc_output){
		printf("cannot allocate encryption buffer (%dB)", input_size);
		exit(1);
	}
	if(!dec_output){
		printf("cannot allocate decryption buffer (%dB)", input_size);
		exit(1);
	}
*/

	*FI_START = 1;

	/* Read the key */
	cp = key_hexstr;
	while(i < 64 && *cp)    /* the maximum key length is 32 bytes and   */
	{                       /* hence at most 64 hexadecimal digits      */
		ch = toupper(*cp++);            /* process a hexadecimal digit  */
		if(ch >= '0' && ch <= '9')
			by = (by << 4) + ch - '0';
		else if(ch >= 'A' && ch <= 'F')
			by = (by << 4) + ch - 'A' + 10;
		else                            /* error if not hexadecimal     */
		{
			printf("key must be in hexadecimal notation\n");
			exit(-1);
		}

		/* store a key byte for each pair of hexadecimal digits         */
		if(i++ & 1)
			ukey[i / 2 - 1] = by & 0xff;
	}

	BF_set_key(&key,8,ukey);
	if(*cp)
	{
		printf("Bad key value.\n");
		exit(-1);
	}

	//ENCRYPT
	i=0;
	t_src = input;
	t_dest = enc_output;
	bzero(ivec, 8);
	while(t_src < input + INPUT_SIZE)
	{
		i=mem_read(indata, 40,&t_src, input + INPUT_SIZE);
		BF_cfb64_encrypt(indata,outdata,i,&key,ivec,&num,1);
		mem_write(outdata, i, &t_dest, enc_output + INPUT_SIZE);
		i=0;
	}

	//DECRYPT
	i=0;
	t_src = enc_output;
	t_dest = dec_output;
	bzero(ivec, 8);
	while(t_src < enc_output+INPUT_SIZE)
	{
		i=mem_read(indata, 40, &t_src, enc_output + INPUT_SIZE);
		BF_cfb64_encrypt(indata,outdata,i,&key,ivec,&num,0);
		mem_write(outdata, i, &t_dest, dec_output + INPUT_SIZE);
		i=0;
	}

	//compare plaintext with result
	//printf("diff? -> %d\n", memcmp(input, dec_output, input_size));

	*FI_STOP = 1;
	hashResultData(dec_output, INPUT_SIZE);
	*CPU_DONE = 1;

	exit(1);
}



