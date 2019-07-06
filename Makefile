TOPTARGETS := all clean

SUBDIRS := mibench-auto-basicmath mibench-auto-bitcount mibench-auto-susan mibench-security-blowfish mibench-security-rijndael mibench-security-sha sort-bubble-unsorted sort-gnomesort-unsorted sort-heapsort-unsorted sort-insertionsort-unsorted sort-mergesort-unsorted sort-quicksort-unsorted sort-selectionsort-unsorted sort-shellsort-unsorted

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)
