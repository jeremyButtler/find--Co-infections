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

#include "fqGrepFqFun.h"
    /*Includes: <stdio.h>*/
#include "fqGrepHash.h"
    /*Includes:
          - fastqGrepAVLTree.h:
              - <string.h>
              - "fastqGrepStructs.h"
                  - <stdlib.h>
                  - <stdio.h>
    */

/*##############################################################################
# Output:
#    Stdout: Prints out kept reads
#    Returns: 0 if not a valid fastq file
##############################################################################*/
char fastqExtract(
    FILE *filterFile,   /*FILE object pointing to file with ID's to search for*/
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char sizeReadStackUChar, /*Number of elements to use in stack*/
    unsigned long buffSizeULng,       /*Size of buffer to read input with*/
    char hashSearchChar,        /*1: do a hash search, 0: only do AVL Tree*/
    char printReverseChar       /*1: Keep reads in filter file
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
    unsigned char *lineInCStr,        /*Buffer to hold one line from file*/
    unsigned int buffSizeULng         /*Size of buffer to read each line*/
); /*Builds a readInfo tree with read id's in filterFile*/

/*##############################################################################
# Output:
#    stdout: Prints extracted reads
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsInTree(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo *readTree  /*root of readInfo tree of read id's to search*/
); /*Extract target reads from fastq file*/

/*##############################################################################
# Output:
#    stdout: Prints out reads not in tree
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsNotInTree(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo *readTree  /*root of readInfo tree of read id's to search*/
); /*Extract reads not given as targets from a fastq file*/

/*##############################################################################
# Output:
#    stdout: Prints out reads in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsInHash(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    unsigned long majicNumULng,   /*Holds majick number for kunths hash*/
    unsigned char digPerKeyUChar, /*Digits needed to get a key*/
    struct readInfo **hashTbl  /*root of readInfo tree of read id's to search*/
); /*Extract target reads from fastq file with hash table*/

/*##############################################################################
# Output:
#    stdout: Prints out reads not in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsNotInHash(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    unsigned long majicNumULng,   /*Holds majick number for kunths hash*/
    unsigned char digPerKeyUChar, /*Digits needed to get a key*/
    struct readInfo **hashTbl  /*root of readInfo tree of read id's to search*/
); /*Extract non-target reads from fastq file using a hash table*/

#endif
