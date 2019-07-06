#ifndef __CONFIG_H__
#define __CONFIG_H__

#define NUM_FI_CPUS 16
#define FI_CPU_RAM_KB 128 // do not forget to set it in .ld file, too!!
#define FI_RESULT_DATA_SIZE_KB 1 // do not forget to set in .ld file, too!!

#define CONTROLLER_CPU_RAM_KB 360


// like the numbers in config.v, increased by one (because count)
#define NUMBER_FI_FLIPFLOPS 4204
#define NUMBER_FI_GATES 14236

#define CPU_DONE ((volatile int*)0x40001000)
#define FAULT_DETECTED ((volatile int*)0x40001004)
#define FI_START ((volatile int*)0x40001008)
#define FI_STOP ((volatile int*)0x4000100c)


#endif
