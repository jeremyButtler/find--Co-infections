gccArgs=cStrToNumberFun.c \
    cStrFun.c \
    FCIStatsFun.c \
    printErrors.c \
    samEntryStruct.c \
    fqAndFaFun.c \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    fqGetIds.c

# Program to extract fastq entries by a set of read ids
fqGetIds:
	$(CC) -O3 $(gccArgs) -o fqGetIdsFast || gcc -O3 $(gccArgs) -o fqGetIdsFast || egcc -O3 $(gccArgs) -o fqGetids || cc -O3 $(gccArgs) -o fqGetIdsFast

# Program to extract fastq entries by a set of read ids
fqGetIdsMemFast:
	$(CC) -O3 -DMEM $(gccArgs) -o fqGetIdsMemFast || gcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemFast || egcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemFast || cc -O3 -DMEM $(gccArgs) -o fqGetIdsMemFast

linux:
	gcc -O3 $(gccArgs) -o fqGetIdsFast
linuxMem:
	gcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemFast
# Program to extract fastq entries by a set of read ids
fqGetIdsMemFastO0:
	gcc -O0 -g -Wall -DMEM $(gccArgs) -o fqGetIdsMemFast


