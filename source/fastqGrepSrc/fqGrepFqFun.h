/*##############################################################################
# Name: fastqGrepFastqFun
# Use:
#    Functions to read input from fastq, find fastq header, & print out entries
##############################################################################*/

#ifndef FQGREPFQFUN_H
#define FQGREPFQFUN_H

#include <stdio.h>
#include "fqGrepStructs.h" /*Has structure functions & hex look up tables*/

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char parseFastqHeader(
    unsigned char *bufferCStr,     /*buffer to hold fread input (can have data)*/
    unsigned char **startNameCStr, /*Will hold the start of the name*/
    unsigned char **endNameCStr,   /*start of read name, will point to end*/
    unsigned long *lenInputULng,        /*Length of input from fread*/
    int buffSizeInt,         /*Size of buffer to work on*/
    unsigned long *lenIdULng,/*Lengtho of the read id*/
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
char printFastqEntry(
    unsigned char *bufferCStr,    /*buffer to hold fread input (can have data)*/
    unsigned char **pointInBufferCStr, /*locatoin working on in buffer*/
    unsigned char **readStartCStr,  /*start of read name*/
    int buffSizeInt,              /*Size of buffer to work on*/
    unsigned long *lenInputULng,  /*Length of input from fread*/
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
    unsigned char *bufferCStr,     /*buffer to hold fread input (can have data)*/
    unsigned char **pointInBufferCStr, /*position working on in buffer*/
    int buffSizeInt,              /*Size of buffer to work on*/
    unsigned long *lenInputULng,  /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
); /*Moves to next fastq read, without printing out*/

/*##############################################################################
# Output:
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdULng to hold the length of the read id
#    Modifies: lenInputInt to hold length of input buffer
#    Returns: readInfo struct with bigNum struct having read id converted to hex
#        - 0 if fails or end of file (lenIdULng < buffSizeInt)
##############################################################################*/
struct readInfo * cnvtIdToBigNum(
    unsigned char *bufferCStr, /*buffer to hold fread input (can have data)*/
    int buffSizeInt,         /*Size of buffer to work on*/
    unsigned char **endNameCStr, /*Points to start of id, will point to end*/
    unsigned long *lenInputULng,        /*Length of input from fread*/
    unsigned long *lenIdULng,/*Lengtho of the read id*/
    unsigned char *lenBigNumChar, /*Holds size to make bigNumber*/
    FILE *fastqFile          /*Fastq file to get data from*/
); /*Converts read id in buffer to bigNum read id, will grab new file input*/

#endif
