/*##############################################################################
# Name: fastqGrepSearchFastq
# Use:
#    Extracts target reads from fastq file
# Requires:
#    fastqGrepAVLTree
#    fastqGrepStructs (called by fastqGrepAVLTree)
##############################################################################*/

#ifndef FQGREPSEARCHFQ_H
#define FQGREPSEARCHFQ_H

#include "fqGetIdsFqFun.h" /*includes fqGetIdsStructs.h*/
#include "fqGetIdsHash.h"
    /*Includes:
          - fqGetIdsAVLTree.h:
              - <string.h>
              - "fqGetIdsStructs.h"
                  - <stdlib.h>
                  - <stdio.h>
                  - <stdint.h>
    */

/*##############################################################################
# Output:
#    Stdout: Prints out kept reads
#    Returns: 0 if not a valid fastq file
##############################################################################*/
uint8_t fastqExtract(
    char *filtPathCStr,        /*Path to file with read ids to extract*/
    char *fqPathCStr,          /*Path to fastq file to extract reads*/
    char *outPathCStr,         /*Path to fastq file to to write reads*/
    uint8_t sizeReadStackUChar, /*Number of elements to use in stack*/
    uint32_t buffSizeUInt,      /*Size of buffer to read input with*/
    uint8_t hashSearchChar,     /*1: do hash search, 0: do Tree search*/
    uint8_t printReverseChar    /*1: Keep reads in filter file
                                  0: ingore reads in filter file*/
); /*Searches and extracts reads from a fastq file using read id's*/

/*##############################################################################
# Output:
#    Returns: balanced readInfo tree with all read id's in filterFile
#    Returns: 0 if malloc errored out
##############################################################################*/
struct readInfo * buildAvlTree(
    FILE *filterFile,                 /*File with read ids to keep or ignore*/
    struct readNodeStack *readStack,  /*Stack to use in building AVL tree*/
    char *lineInCStr,        /*Buffer to hold one line from file*/
    uint32_t buffSizeUInt         /*Size of buffer to read each line*/
); /*Builds a readInfo tree with read id's in filterFile*/

/*##############################################################################
# Output:
#    stdout: Prints out reads in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
uint8_t extractReads(
    FILE *fastqFile,            /*fastq file to search through*/
    FILE *outFile,              /*File to write extracted reads to*/
    char *lineInCStr,        /*Buffer to hold input from fastq file*/
    uint32_t buffSizeUInt,      /*Size of lineInCStr*/
    unsigned long majicNumULng, /*Holds majick number for kunths hash*/
    uint8_t digPerKeyUChar,     /*Digits needed to get a key*/
    uint8_t *printNonMatchBool, /*1: print non-match, 0: print match*/
    struct readInfo *readTree,  /*For AVL search (hashTbl == 0)*/
    struct readInfo **hashTbl   /*Hash table to search for ids in*/
); /*Extract target reads from fastq file with hash table or tree*/

#endif
