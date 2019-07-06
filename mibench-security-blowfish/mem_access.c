#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* reat count bytes from src into dest and treat src as a stream pointer */
size_t mem_read( void *dest, size_t count, char **src, void *bound ){
	size_t left = (char *)bound - *src;	//bytes left below bound
	size_t num = left>count ? count : left; //bytes to read

	memcpy(dest, *src, num);
	*src = (*src)+num;

	return num;
}

/* write count bytes from src to dest and treat dest as a stream pointer */
size_t mem_write( void *src, size_t count, char **dest, void *bound){
	size_t free = (char*)bound-*dest;		/* free space in byte */
	size_t num=free>count?count:free;		/* count to be written */

	memcpy(*dest, src, num);
	*dest = *dest+num;

	return num;
}

