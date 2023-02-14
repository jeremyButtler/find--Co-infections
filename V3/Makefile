CC=gcc
gccArgs=-g \
    -Wall \
    cStrToNumberFun.c \
    printErrors.c \
    samEntryStruct.c \
    trimSam.c \
    scoreReadsFun.c \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    findCoInftBinTree.c \
    findCoInft.c \
    -o findCoInft
    
trimGccArgs=-g \
    -Wall \
    cStrToNumberFun.c \
    printErrors.c \
    samEntryStruct.c \
    trimSam.c \
    trimSamFile.c \
    -o trimSamFile

scoreGccArgs= -g \
         -Wall\
         printErrors.c \
         cStrToNumberFun.c \
         samEntryStruct.c \
         scoreReadsFun.c \
         scoreReads.c \
         -o scoreReads

fqGetIdsArgs=-g \
    -Wall \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    fqGetIds.c \
    -o fqGetIds


# Build findCoInfct
all:
	$(CC) -O2 $(gccArgs) || gcc -O2 $(gccArgs) || egcc -O2 $(gccArgs) || cc -O2 $(gccArgs)

# Here if user just wants scoreReads
scoreReads:
	$(CC) -O2 $(scoreGccArgs) || gcc -O2 $(scoreGccArgs) || egcc -O2 $(scoreGccArgs) || cc -O2 $(scoreGccArgs)

# Here if user just wants trimSamFile
trimSam:
	$(CC) -O2 $(trimGccArgs) || gcc -O2 $(trimGccArgs) || egcc -O2 $(trimGccArgs) || cc -O2 $(trimGccArgs)

fqGetIds:
	$(CC) -O3 $(fqGetIdsArgs) || gcc -O3 $(fqGetIdsArgs) || egcc -O3 $(fqGetIdsArgs) || cc -O3 $(fqGetIdsArgs)

# The || is to find the compiler that is on the users system
