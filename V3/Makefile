CC=gcc
gccArgs=-g \
    -Wall \
    cStrToNumberFun.c \
    cStrFun.c \
    minAlnStatsStruct.c \
    FCIStatsFun.c \
    printErrors.c \
    samEntryStruct.c \
    findCoInftChecks.c \
    fqAndFaFun.c \
    trimSam.c \
    scoreReadsFun.c \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    findCoInftBinTree.c \
    readExtract.c \
    binReadsFun.c \
    buildConFun.c \
    findCoInft.c \
    -o findCoInft

gccBuildConArgs=-g \
    -Wall \
    cStrToNumberFun.c \
    cStrFun.c \
	minAlnStatsStruct.c \
    FCIStatsFun.c \
    printErrors.c \
    samEntryStruct.c \
    findCoInftChecks.c \
    fqAndFaFun.c \
    trimSam.c \
    scoreReadsFun.c \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    findCoInftBinTree.c \
    readExtract.c \
	buildConFun.c \
    buildCon.c \
    -o buildCon
    
gccBinReadsArgs=-g \
    -Wall \
    fqAndFaFun.c \
    cStrToNumberFun.c \
    cStrFun.c \
	minAlnStatsStruct.c \
	FCIStatsFun.c \
    printErrors.c \
    samEntryStruct.c \
    findCoInftChecks.c \
    trimSam.c \
    scoreReadsFun.c \
    findCoInftBinTree.c \
	binReadsFun.c \
    binReads.c \
    -o binReads

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
    minAlnStatsStruct.c \
    FCIStatsFun.c \
    samEntryStruct.c \
    fqAndFaFun.c \
    scoreReadsFun.c \
    scoreReads.c \
    -o scoreReads

fqGetIdsArgs=-g \
    -Wall \
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

fqGIdO=-o fqGetIds
fqGIdOMem=-DMEM -o fqGetIdsMem

# To find the compiler I am using || to default to different names for
# gcc when the default gcc does not exist or failes. || is universal
# across all unix OS's and is valid in the bsd make and gnu make, unlike
# the if statements, which have different syntax for bsd make and gnu
# make. However, this unix support comes at the cost of not being able
# to compile this program on windows. That being said, windows needs
# minigw or visual studio, so windows support would be very limited in
# or require an advanced user who could easily edit this file to support
# windows.

# Build findCoInfct
all:
	$(CC) -O2 $(gccArgs) || gcc -O2 $(gccArgs) || egcc -O2 $(gccArgs) || cc -O2 $(gccArgs)

# My consensus building step as  a separate program
buildCon:
	$(CC) -O2 $(gccBuildConArgs) || gcc -O2 $(gccBuildConArgs) || egcc -O2 $(gccBuildConArgs) || cc -O2 $(gccBuildConArgs)

# Program to bin reads using a set of references
binReads:
	$(CC) -O2 $(gccBinReadsArgs) || gcc -O2 $(gccBinReadsArgs) || egcc -O2 $(gccBinReadsArgs) || cc -O2 $(gccBinReadsArgs)

# Program to get stats for each alignment in a sam 
scoreReads:
	$(CC) -O2 $(scoreGccArgs) || gcc -O2 $(scoreGccArgs) || egcc -O2 $(scoreGccArgs) || cc -O2 $(scoreGccArgs)

# Program for trimming softmasking off cigar, seq, and q-score entries
trimSam:
	$(CC) -O2 $(trimGccArgs) || gcc -O2 $(trimGccArgs) || egcc -O2 $(trimGccArgs) || cc -O2 $(trimGccArgs)

# Program to extract fastq entries by a set of read ids
fqGetIds:
	$(CC) -O3 $(fqGetIdsArgs) $(fqGIdO) || gcc -O3 $(fqGetIdsArgs) $(fqGIdO) || egcc -O3 $(fqGetIdsArgs) $(fqGIdO) || cc -O3 $(fqGetIdsArgs) $(fqGIdO)

# Program to extract fastq entries by a set of read ids
fqGetIdsMem:
	$(CC) -O3 $(fqGetIdsArgs) $(fqGIdOMem) || gcc -O3 $(fqGetIdsArgs) $(fqGIdOMem) || egcc -O3 $(fqGetIdsArgs) $(fqGIdOMem) || cc -O3 $(fqGetIdsArgs) $(fqGIdOMem)

# commands for indvidual compilers
openbsd:
	egcc -O2 $(gccArgs)
linux:
	gcc -O2 $(gccArgs)
cc:
	cc -O2 $(gccArgs)
