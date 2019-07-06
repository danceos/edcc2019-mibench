#include <stdio.h>
#include <sys/times.h>
#include <sys/types.h>

extern char __heap_start;
static char *cur_heap_end = &__heap_start;

caddr_t _sbrk(int nbytes) {
  cur_heap_end += nbytes;
  return cur_heap_end;
}

