CC=gcc

gccArgs=-Wall \
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
    trimPrimersStructs.c \
    trimPrimersAvlTree.c \
    trimPrimersHash.c \
    trimPrimersSearch.c \
    findCoInftBinTree.c \
    readExtract.c \
    binReadsFun.c \
    buildConFun.c \
    findCoInft.c \
    -o findCoInft

gccBuildConArgs=-Wall \
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
    
gccBinReadsArgs=-Wall \
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

trimGccArgs=-Wall \
    cStrToNumberFun.c \
    printErrors.c \
    samEntryStruct.c \
    trimSam.c \
    trimSamFile.c \
    -o trimSamFile

scoreGccArgs=-Wall\
    printErrors.c \
    cStrToNumberFun.c \
    minAlnStatsStruct.c \
    FCIStatsFun.c \
    samEntryStruct.c \
    fqAndFaFun.c \
    scoreReadsFun.c \
    scoreReads.c \
    -o scoreReads

fqGetIdsArgs=-Wall \
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

trimPrimGcc=-Wall \
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
    trimPrimersStructs.c \
    trimPrimersAvlTree.c \
    trimPrimersHash.c \
    trimPrimersSearch.c \
    trimPrimers.c

trimPrimO=-o trimPrimers
trimPrimMemO=-DMEM -o trimPrimers

filterReadsGcc=-Wall \
    cStrToNumberFun.c \
    minAlnStatsStruct.c \
    FCIStatsFun.c \
    printErrors.c \
    samEntryStruct.c \
    fqAndFaFun.c \
    filterReads.c \
    -o filterReads

extractReadsGcc=-Wall \
    cStrToNumberFun.c \
    cStrFun.c \
    printErrors.c \
    minAlnStatsStruct.c \
    FCIStatsFun.c \
    samEntryStruct.c \
    trimSam.c \
    findCoInftChecks.c \
    fqAndFaFun.c \
    scoreReadsFun.c \
    fqGetIdsFqFun.c \
    fqGetIdsStructs.c \
    fqGetIdsAVLTree.c \
    fqGetIdsHash.c \
    fqGetIdsSearchFq.c \
    findCoInftBinTree.c \
    readExtract.c \
    extractTopReads.c \
    -o extractTopReads

alignSeqGcc=-Wall\
    sequenceFun.c \
    cStrToNumberFun.c \
    alignmentsFun.c \
    alignSeq.c \
    -o alignSeq

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
	$(CC) -O3 $(gccArgs) || gcc -O3 $(gccArgs) || egcc -O3 $(gccArgs) || cc -O3 $(gccArgs)

# My consensus building step as  a separate program
buildCon:
	$(CC) -O3 $(gccBuildConArgs) || gcc -O3 $(gccBuildConArgs) || egcc -O3 $(gccBuildConArgs) || cc -O3 $(gccBuildConArgs)

# Program to bin reads using a set of references
binReads:
	$(CC) -O3 $(gccBinReadsArgs) || gcc -O3 $(gccBinReadsArgs) || egcc -O3 $(gccBinReadsArgs) || cc -O3 $(gccBinReadsArgs)

# Program to get stats for each alignment in a sam 
scoreReads:
	$(CC) -O3 $(scoreGccArgs) || gcc -O3 $(scoreGccArgs) || egcc -O3 $(scoreGccArgs) || cc -O3 $(scoreGccArgs)

# Program for trimming softmasking off cigar, seq, and q-score entries
trimSam:
	$(CC) -O3 $(trimGccArgs) || gcc -O3 $(trimGccArgs) || egcc -O3 $(trimGccArgs) || cc -O3 $(trimGccArgs)

# Program to extract fastq entries by a set of read ids
fqGetIds:
	$(CC) -O3 $(fqGetIdsArgs) $(fqGIdO) || gcc -O3 $(fqGetIdsArgs) $(fqGIdO) || egcc -O3 $(fqGetIdsArgs) $(fqGIdO) || cc -O3 $(fqGetIdsArgs) $(fqGIdO)

# Program to extract fastq entries by a set of read ids
trimPrimers:
	$(CC) -O3 $(trimPrimGcc) $(trimPrimO) || gcc -O3 $(trimPrimGcc) $(trimPrimO) || egcc -O3 $(trimPrimGcc) $(trimPrimO) || cc -O3 $(trimPrimGcc) $(trimPrimO)

filterReads:
	$(CC) -O3 $(filterReadsGcc) || gcc -O3 $(filterReadsGcc) || egcc -O3 $(filterReadsGcc) || cc -O3 $(filterReadsGcc)

extractReads:
	$(CC) -O3 $(extractReadsGcc) || gcc -O3 $(extractReadsGcc) || egcc -O3 $(extractReadsGcc) || cc -O3 $(extractReadsGcc)

# For my Needleman Wunsch alignment (May add some more stuff later)
alignSeq:
	$(CC) -Ofast $(alignSeqGcc) || gcc -Ofast $(alignSeqGcc) || egcc -Ofast $(alignSeqGcc) || cc -Ofast $(alignSeqGcc)

# commands for indvidual compilers
openbsd:
	egcc -O3 $(gccArgs)
linux:
	gcc -O3 $(gccArgs)
cc:
	cc -O3 $(gccArgs)
