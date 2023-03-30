/*##############################################################################
# Name: fastqGrepFastqFun
# Use:
#    Functions to read input from fastq, find fastq header, & print out entries
##############################################################################*/

#ifndef FQGREPFQFUN_H
#define FQGREPFQFUN_H

#include "fqGetIdsStructs.h" /*<stdlib.h>, <stdio.h>, <stdint.h>*/
#include "fqAndFaFun.h"      /*Move to next fastq entry*/
     /*Has structure functions & hex look up tables*/

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdInt to hold the length of the read id
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
uint8_t parseFastqHeader(
    char *bufferCStr,     /*buffer to hold fread input (can have data)*/
    char **startNameCStr, /*Will hold the start of the name*/
    char **endNameCStr,   /*start of read name, will point to end*/
    uint64_t *lenInputULng,  /*Length of input from fread*/
    uint32_t buffSizeInt,    /*Size of buffer to work on*/
    int32_t *lenIdInt,     /*Length of the read id*/
    struct bigNum *idBigNum,  /*Will hold big number found*/
    FILE *fastqFile          /*Fastq file to get data from*/
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
uint8_t printFastqEntry(
    char *bufferCStr,    /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr, /*locatoin working on in buffer*/
    char **readStartCStr,  /*start of read name*/
    uint32_t buffSizeInt,              /*Size of buffer to work on*/
    uint64_t *lenInputULng,  /*Length of input from fread*/
    FILE *outFILE,                /*File to print reads to*/
    FILE *fastqFile               /*Fastq file to get data from*/
); /*Reads input from file & marks end of read name*/

#endif
