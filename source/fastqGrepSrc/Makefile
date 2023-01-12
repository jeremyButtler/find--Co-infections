CC=gcc

gccArgs=-g \
    -Wall \
    fqGrepFqFun.c \
    fqGrepStructs.c \
    fqGrepAVLTree.c \
    fqGrepHash.c \
    fqGrepSearchFq.c \
    fqGrep.c \
    -o fastqGrep

# Build fastqGrep
all:
	$(CC) -O3 $(gccArgs) || gcc -O3 $(gccArgs) || egcc -O3 $(gccArgs) || cc -O3 $(gccArgs)

debug:
	egcc -O0 $(gccArgs)
	bash runGdbFqGrep.sh

openbsd:
	egcc -O3 $(gccArgs)
linux:
	gcc -O3 $(gccArgs)
cc:
	cc -O3 $(gccArgs)

# The || is to find the compiler that is on the users system
