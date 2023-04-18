CC=cc
CFLAGS=\
  -O3
  #-DMEM to uses unsigned longs for big numbers
  # default summs up big numbers in unsigned long and
  # stores big number in an unsigned int
outFLAG=\
  -o fqGetIdsIll

gccArgs=\
  cStrToNumberFun.c \
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
	$(CC) $(CFLAGS) $(gccArgs) $(outFLAG) ||\
    gcc $(CFLAGS) $(gccArgs) $(outFLAG) ||\
    egcc $(CFLAGS) $(gccArgs) $(outFLAG) ||\
    cc $(CFLAGS) $(gccArgs) $(outFLAG)
