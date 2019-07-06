#include <inttypes.h>
#include "config.h"
#include "result.h"

volatile unsigned char __attribute__((section (".resultBlock"))) fiResultData[FI_RESULT_DATA_SIZE_KB * 1024] = { 0 };

volatile unsigned char* getResultSectionPointer() {
	return &fiResultData[0];
}

#if 0
#include <stdio.h>
void printHex(void* dataStart, unsigned int dataLen) {
	uint8_t* data = (uint8_t*) dataStart;

	for (int i = 0; i < dataLen; i++) {
		printf("%02x", data[i]);
	}
	fflush(stdout);
}

void hashResultData(void* dataStart, unsigned int dataLen) {
	printHex(dataStart, dataLen);
}

void hashIncrResultData(void* dataStart, unsigned int dataLen) {
	printHex(dataStart, dataLen);
}

#else

/*
 * Generates a XOR-Hash of some data of variable length (esp. larger than the
 * result data section) and stores it in the result section
 */
void hashResultData(void* dataStart, unsigned int dataLen) {
/*	volatile uint32_t* fiResultData32 = (volatile uint32_t*) fiResultData;
	uint32_t* data32 = (uint32_t*) dataStart;
	uint32_t dataLen32 = (dataLen - dataLen%4)/4;
	
	int xorElementsPerDatum = dataLen / (FI_RESULT_DATA_SIZE_KB * 256); // +1 is missing here... //TODO: validate etc..

	uint32_t sum = 0;

	// fiResultData32 length is FI_RESULT_DATA_SIZE_KB*256
	for (int i = 0; i < FI_RESULT_DATA_SIZE_KB*256; i++) {
		
		uint32_t value = 0;
		for (int j = 0; j < xorElementsPerDatum; j++) {

			if (j * FI_RESULT_DATA_SIZE_KB*256 + i < dataLen32) {
				sum += data32[j * FI_RESULT_DATA_SIZE_KB*256 + i];
				value = value ^ data32[j * FI_RESULT_DATA_SIZE_KB*256 + i];
			}
		}
		fiResultData32[i] ^= value;
	}
	
	fiResultData32[(FI_RESULT_DATA_SIZE_KB * 256) - 4] = fiResultData32[(FI_RESULT_DATA_SIZE_KB * 256) - 4] ^ sum; // copy sum somewhere...
	*/

	// TODO: see if there is an error in the above code and fix it. did not work in mibench-auto-susan benchmark, code below is fine.

	uint8_t* data = (uint8_t*) dataStart;

	for (int i = 0; i < dataLen; i++) {
		fiResultData[i%(FI_RESULT_DATA_SIZE_KB * 1024)] ^= data[i];
	}


}


/*
 * This is basically a result printf-replacement. This funtion hashes the data in a 
 * round robin way into the result array, while xor'ing the data with the existing data.
 * dataLen should be divisible by 4.
 */
int hashIncrResultDataIndex = 0;
void hashIncrResultData(void* dataStart, unsigned int dataLen) {
	
	if (dataLen > FI_RESULT_DATA_SIZE_KB * 1024) {
		// data too big, use normal hash function
		hashResultData(dataStart, dataLen);
		return;
	}

	uint8_t* data = (uint8_t*) dataStart;

	for (int i = 0; i < dataLen; i++) {
		fiResultData[hashIncrResultDataIndex] ^= data[i];
		
		hashIncrResultDataIndex++;
		if (hashIncrResultDataIndex >= FI_RESULT_DATA_SIZE_KB * 1024) {
			hashIncrResultDataIndex = 0;
		}
	}
}

#endif
