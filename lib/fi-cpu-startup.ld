OUTPUT_FORMAT(elf32-littlearm)
ENTRY(_start)

MEMORY {
  ram (rwx) : ORIGIN = 0x00000000, LENGTH = 128K
}

SECTIONS {
  .text : {
    KEEP(*(.vectors))
    KEEP(*(.init))
    *(.text)
    *(.text.*)
    *(.gnu.linkonce.t.*)
	*(.glue_7t)
	*(.glue_7)

    /* C++ con-/destructors */
    __ctors_start = . ;
    KEEP(*(.init_array))
    __ctors_end = . ;
    __dtors_start = . ;
    KEEP(*(.fini_array))
    __dtors_end = . ;

    __text_end = .;
  } > ram

  /* .ARM.exidx is sorted, so has to go in its own output section.  */
  __exidx_start = .;
  .ARM.exidx :
  {
    *(.ARM.exidx* .gnu.linkonce.armexidx.*)
  }
  __exidx_end = .;

  .rodata : ALIGN(4) {
    __rodata_start = .;
    *(.rodata)
    *(.rodata.*)
    *(.gnu.linkonce.r.*)
    __rodata_end = .;
  } > ram

  .data : ALIGN(4) {
    __data_start = .;
    *(.data)
    *(.data.*)
    *(.gnu.linkonce.d.*)
    __data_end = .;
  } > ram

  .bss (NOLOAD) : ALIGN(4) {
    __bss_start__ = .;
    *(.bss)
    *(.bss.*)
    *(COMMON)
    . = ALIGN(4);
    __bss_end__ = .;
  } > ram

  PROVIDE(__heap_start = ALIGN(. , 4));

  PROVIDE(__stack = ORIGIN(ram) + LENGTH(ram) - 1K);

  .resultBlock LENGTH(ram) - 1K: {
	  KEEP(*(.resultBlock))
  } > ram

  /* Everyone seems to copy the stuff below straight from somewhere else, so I'll do that too */

  /* Stabs debugging sections.  */
  .stab 0 : { *(.stab) }
  .stabstr 0 : { *(.stabstr) }
  .stab.excl 0 : { *(.stab.excl) }
  .stab.exclstr 0 : { *(.stab.exclstr) }
  .stab.index 0 : { *(.stab.index) }
  .stab.indexstr 0 : { *(.stab.indexstr) }
  .comment 0 : { *(.comment) }
  /* DWARF debug sections.
     Symbols in the DWARF debugging sections are relative to the beginning
     of the section so we begin them at 0.  */
  /* DWARF 1 */
  .debug          0 : { *(.debug) }
  .line           0 : { *(.line) }
  /* GNU DWARF 1 extensions */
  .debug_srcinfo  0 : { *(.debug_srcinfo) }
  .debug_sfnames  0 : { *(.debug_sfnames) }
  /* DWARF 1.1 and DWARF 2 */
  .debug_aranges  0 : { *(.debug_aranges) }
  .debug_pubnames 0 : { *(.debug_pubnames) }
  /* DWARF 2 */
  .debug_info     0 : { *(.debug_info) *(.gnu.linkonce.wi.*) }
  .debug_abbrev   0 : { *(.debug_abbrev) }
  .debug_line     0 : { *(.debug_line) }
  .debug_frame    0 : { *(.debug_frame) }
  .debug_str      0 : { *(.debug_str) }
  .debug_loc      0 : { *(.debug_loc) }
  .debug_macinfo  0 : { *(.debug_macinfo) }
}

