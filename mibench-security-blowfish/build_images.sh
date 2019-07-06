#!/bin/bash
set -e

echo removing old ISOs end ELFs \(if any\)...
rm -f bf_small bf_small.iso bf_large bf_large.iso

#__INSTALL_DIR is defined in make_all.sh for mibench security
make INSTALL_DIR=$__INSTALL_DIR inputdata bf_small bf_large 

#create them ISO images
../../build_grub-iso.sh bf_small bf_small.iso
../../build_grub-iso.sh bf_large bf_large.iso

#move ELF files
mv bf_small bf_small.elf
mv bf_large bf_large.elf

#clenaup
make INSTALL_DIR=$__INSTALL_DIR clean > /dev/null 2>&1

echo the die is cast.
