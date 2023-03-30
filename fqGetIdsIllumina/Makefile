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
all:
	$(CC) -O3 $(gccArgs) -o fqGetIdsIll || gcc -O3 $(gccArgs) -o fqGetIdsIll || egcc -O3 $(gccArgs) -o fqGetidsIll || cc -O3 $(gccArgs) -o fqGetIdsIll

# Program to extract fastq entries by a set of read ids
mem:
	$(CC) -O3 -DMEM $(gccArgs) -o fqGetIdsMemIll || gcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemIll || egcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemIll || cc -O3 -DMEM $(gccArgs) -o fqGetIdsMemIll

linux:
	gcc -O3 $(gccArgs) -o fqGetIds
linuxMem:
	gcc -O3 -DMEM $(gccArgs) -o fqGetIdsMemIll
# Program to extract fastq entries by a set of read ids
fqGetIdsMemIllO0:
	gcc -O0 -g -Wall -DMEM $(gccArgs) -o fqGetIdsMemIll


