/* NIST Secure Hash Algorithm */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "sha.h"

#include "mem_access.h"
#include "config.h"
#include "input.h"
#include "result.h"

/*
 * use input data from memory
 */
//extern char _binary_input_asc_start;
//extern char _binary_input_asc_end;
//extern char _binary_input_asc_size;
//char *input=&_binary_input_asc_start;
//int input_size=0;


int main(int argc, char **argv)
{
	//extern void *ECOS_data;
	//extern int ECOS_data_size;

	*FI_START = 1;
	SHA_INFO sha_info;
	unsigned long output[5] = {0,0,0,0,0};

	ecos_sha_stream(&sha_info, input);
	memcpy(output, sha_info.digest, 5*sizeof(unsigned long));

	*FI_STOP = 1;
	hashResultData(output, 5*sizeof(unsigned long));
	*CPU_DONE = 1;


	return(0);
}
