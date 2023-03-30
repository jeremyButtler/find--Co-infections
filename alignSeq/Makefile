CC=gcc

alignSeqGcc=-Wall\
    twoBitArrays.c \
    sequenceFun.c \
    cStrToNumberFun.c \
    alignmentsFun.c \
    alignSeq.c \
    -o alignSeq

# Build findCoInfct
all:
	$(CC) -Ofast $(alignSeqGcc) || gcc -Ofast $(alignSeqGcc) || egcc -Ofast $(alignSeqGcc) || cc -Ofast $(alignSeqGcc)

linux:
	gcc -Ofast $(alignSeqGcc)
openbsd:
	egcc -Ofast $(alignSeqGcc)
cc:
	cc -Ofast $(alignSeqGcc)
