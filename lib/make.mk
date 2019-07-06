
GCC := arm-none-eabi-gcc
OBJCOPY := arm-none-eabi-objcopy
NM := arm-none-eabi-nm
LD := arm-none-eabi-ld
AS := arm-none-eabi-as
GCC_LIBDIR := /usr/lib/gcc/arm-none-eabi/4.8

TARGET := $(shell basename `pwd`)

ifeq "$(PLATFORM)" "gem5"
	CFLAGS += -g -O3 -mcpu=cortex-m0 -mthumb -nostartfiles -ffunction-sections -fdata-sections -fno-strict-aliasing -Wall -std=gnu11 -mfloat-abi=soft -I../lib/ -DGEM5
	CRT = $(GCC_LIBDIR)/crtbegin.o $(GCC_LIBDIR)/crtend.o $(GCC_LIBDIR)/crti.o
	LDFLAGS += -T../lib/fi-cpu-startup_gem5.ld -Wl,--gc-sections
	LIBS += -L$(GCC_LIBDIR) -lgcc
	ASFLAGS = -g
else
	CFLAGS += -g -O3 -mcpu=cortex-m0 -mthumb -nostartfiles -ffunction-sections -fdata-sections -fno-strict-aliasing -Wall -std=gnu11 -mfloat-abi=soft -I../lib/
	LDFLAGS += -T../lib/fi-cpu-startup.ld -Wl,--gc-sections
endif


GENDEPFLAGS = -MMD -MP -MF .dep/$(@F).d

SRCFILES = $(wildcard *.c) $(wildcard ../lib/*.c)
OBJFILES += $(SRCFILES:%.c=%.o)

ifeq "$(PLATFORM)" "gem5"
	OBJFILES += ../lib/startup_gem5.o
else
	OBJFILES += ../lib/startup.o
endif




all: $(TARGET).bin $(TARGET).hex

$(TARGET).bin: $(TARGET).bin.full
	printf "%d" 0x`$(NM) $(TARGET).elf | sort | grep __heap_start | cut -d' ' -f1` > $(TARGET).len
	head --bytes `cat $(TARGET).len` $< > $@
	rm $(TARGET).len

$(TARGET).bin.full: $(TARGET).elf
	$(OBJCOPY) -O binary $< $@

$(TARGET).hex: $(TARGET).elf
	$(OBJCOPY) -O verilog $< temp.hex
	grep -v "@" temp.hex > $@
	rm temp.hex

$(TARGET).elf: $(OBJFILES)
ifeq "$(PLATFORM)" "gem5"
	$(GCC) $(LDFLAGS) -nostartfiles -nostdinc -o $@ $(CRT) $^ -static $(LIBS)
else
	$(GCC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)
endif

%.o: %.c
	@mkdir -p .dep
	$(GCC) $(CFLAGS) $(GENDEPFLAGS) -c -o $@ $<

%.o: %.s
ifeq "$(PLATFORM)" "gem5"
	$(AS) $(ASFLAGS) -o $@ -c $<
else
	$(GCC) $(CFLAGS) $(GENDEPFLAGS) -c -o $@ $<
endif

clean:
	rm -f *.bin *.elf *.hex *.o *.mem .dep/* */*.o ../lib/*.o *.bin.full
	rm -rf .dep

