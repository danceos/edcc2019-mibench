#include "config.h"
#include "result.h"
#include "selectionsort.h"

void doSort(int* array, int n) {
	selsort(array, n);
}

int main() {
	extern int array[];
	extern int count;

	*FI_START = 1;
	doSort(array, count);
	*FI_STOP = 1;
	hashResultData(array, count * sizeof(int));
	*CPU_DONE = 1;
}
