// borrowed from https://balau82.wordpress.com/2010/12/16/using-newlib-in-arm-bare-metal-programs/
#include <sys/stat.h>
#include "result.h"

int _close(int file) { return -1; }

int _fstat(int file, struct stat *st) {
	st->st_mode = S_IFCHR;
	return 0;
}

int _isatty(int file) { return 1; }

int _lseek(int file, int ptr, int dir) { return 0; }

int _open(const char *name, int flags, int mode) { return -1; }

int _read(int file, char *ptr, int len) {
	return 0;
}

int _write(int file, char *ptr, int len) {
#if 0
	volatile unsigned int * const UART0DR = (unsigned int*)0x10009000;

	//FIXME: should it be char instead of int?
	for (int i = 0; i < len; i++) {
		*UART0DR = (unsigned int)(*ptr);
		ptr++;
	}
#else
	hashIncrResultData(ptr, len);
#endif
	return len;
}

void _exit(int ret) {
	while (1);
}
