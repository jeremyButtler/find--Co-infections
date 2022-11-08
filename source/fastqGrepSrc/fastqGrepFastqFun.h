/*##############################################################################
# Name: fastqGrepFastqFun
# Use:
#    Functions to read input from fastq, find fastq header, & print out entries
##############################################################################*/

#ifndef FASTQGREPFASTQFUN_H
#define FASTQGREPFASTQFUN_H

#include <stdio.h>
#include "fastqGrepStructs.h"

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: idBigNum to holder head as hex big number (set to 0 for 
#    Modifies: lenBigNumAryUChar when string needs more longs, set to 0 when
#              memory allocatoin failed
#    Returns:
#        8: memory allocation failed
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char parseFastqHeader(
    char *buffCStr,              /*buffer to hold fread input (can have data)*/
    char **startNameCStr,        /*Start of read id*/
    int *lenInInt,               /*Length of input from fread*/
    int lenBuffInt,             /*Size of buffer to work on*/
    struct bigNum *idBigNum,   /*Big number struct to hold read id*/
    unsigned char *lenBigNumAryUChar,  /*Length of U longs in big number*/
    FILE *fqFILE              /*Fastq file to get data from*/
); /*Reads input from file & sets pointer to start of read name*/

/*##############################################################################
# Output:
#    stdout: prints id & the remaing three lines in fastq for read
#    Modifies: bufferCStr to have the next buffer if empty
#    Modifies: incurments pointInBufferCStr to start of next read
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char printFastqEntry(
    char *bufferCStr,            /*buffer to hold fread input (can have data)*/
    char **readStartCStr,         /*Points to start of read name*/
    int buffSizeInt,              /*Size of buffer to work on*/
    int *lenInputInt,             /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
); /*Reads input from file & marks end of read name*/

/*##############################################################################
# Output:
#    Modifies: bufferCStr to have the next buffer if empty
#    Modifies: incurments pointInBufferCStr to start of next read
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char moveToNextFastqEntry(
    char *bufferCStr,            /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr,     /*Points to locatoin working on in buffer*/
    int buffSizeInt,              /*Size of buffer to work on*/
    int *lenInputInt,             /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
); /*Moves to next fastq read, without printing out*/

#endif
