# Using gmp library for big numbers
# For debain: need libgmp-dev
CC=gcc # most likely default on a unix OS
compileCmdStr=-g \
              -Wall \
              -O2 \
              fastqGrepFastqFun.c \
              fastqGrepStructs.c \
              fastqGrepAVLTree.c \
              fastqGrepHash.c \
              fastqGrepSearchFastq.c \
              fastqGrep.c \
              -o fastqGrep

all: #compile # check if system has gcc or egcc
	$(CC) $(compileCmdStr)

egcc: # openbsd install
	egcc $(compileCmdStr)

cc:  # just in case, but likely not needed
	cc $(compileCmdStr)
