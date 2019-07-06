
	.syntax unified

	.section .vectors

	.macro	except label
	.weak	\label
	.set	\label, __unhandled_exception
	.word	\label
	.endm

vectors_start:
	/* Cortex M0-DS exception vectors */
	.word	__stack
	.word	_start                 //  1 - Reset
	except	NMI_Handler            //  2 - NMI
	except	HardFault_Handler      //  3 - HardFault
	.word	__unhandled_exception  //  4 - reserved
	.word	__unhandled_exception  //  5 - reserved
	.word	__unhandled_exception  //  6 - reserved
	.word	__unhandled_exception  //  7 - reserved
	.word	__unhandled_exception  //  8 - reserved
	.word	__unhandled_exception  //  9 - reserved
	.word	__unhandled_exception  // 10 - reserved

	except	SVC_Handler            // 11 - SVCall

	.word	__unhandled_exception  // 12 - reserved
	.word	__unhandled_exception  // 13 - reserved

	except	PendSV_Handler         // 14 - PendSV
	except	SystickHandler        // 15 - SysTick

	except  IRQHandler_UART        // 16 - UART
	except  Injector_Handler       // 17 - common injector IRQ
	except  IRQ2_Handler           // 18 - Interrupt 2
	except  IRQ3_Handler           // 19 - Interrupt 3

	except  IRQ4_Handler           // 20 - Interrupt 4
	except  IRQ5_Handler           // 21 - Interrupt 5
	except  IRQ6_Handler           // 22 - Interrupt 6
	except  IRQ7_Handler           // 23 - Interrupt 7
	
	except  IRQ8_Handler           // 24 - Interrupt 8
	except  IRQ9_Handler           // 25 - Interrupt 9
	except  IRQ10_Handler          // 26 - Interrupt 10
	except  IRQ11_Handler          // 27 - Interrupt 11
	
	except  IRQ12_Handler          // 28 - Interrupt 12
	except  IRQ13_Handler          // 29 - Interrupt 13
	except  IRQ14_Handler          // 30 - Interrupt 14
	except  IRQ15_Handler          // 31 - Interrupt 15
vectors_end:
 
	.section .text

	/* handler for otherwise-unused exceptions */
	/* shows "dead" on the 7-segment display   */
	.weak __unhandled_exception
	.thumb_func
__unhandled_exception:
	ldr	r0, =0
	ldr	r1, =0xdeaddead
	//str	r1, [r0]
failloop:
	b	failloop
	
	
	.global _start
	.thumb_func
_start:
        /* disable interrupts */
        cpsid   i
        
        /* copy vectors to address 0 */
        ldr     r0, =0
        ldr     r1, =0x00000000
        ldr     r2, =(vectors_end - vectors_start) / 4
vecloop:
        ldr     r3, [r1]
        str     r3, [r0]
        adds    r0, r0, #4
        adds    r1, r1, #4
        subs    r2, r2, #1
        bne     vecloop
        
	/* clear BSS */
        ldr     r0, =0
	ldr	r1, =__bss_start__
	ldr	r2, =__bss_end__
bssloop:
        cmp     r1, r2
        beq     bssend
	str	r0, [r1]
	adds	r1, r1, #4
        b       bssloop
bssend:
        
        // FIXME: Set up process stack if used (msp is set by hardware)

        /* call constructors */
        ldr     r4, =__ctors_start
        ldr     r5, =__ctors_end
ctorloop:
        cmp     r4, r5
        beq     ctorend
        ldr     r6, [r4]
        blx     r6
        adds    r4, r4, #4
        b       ctorloop
ctorend:

        /* enable interrupts */
        cpsie   i

	/* start main() */
//        ldr     r0, =main
//        bx      r0
        bl main
loop:
        b loop
