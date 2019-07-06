#ifndef __RESULT_H__
#define __RESULT_H__

volatile unsigned char* getResultSectionPointer();
void hashResultData(void* dataStart, unsigned int dataLen);
void hashIncrResultData(void* dataStart, unsigned int dataLen);

#endif
