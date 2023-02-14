/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepStructs TOC:
#    struct-1: readInfo: Holds the read name
#    struct-2: readNodeStack: Makes a stack (filo) of read nodes
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef FQGREPSTRUCTS_H
#define FQGREPSTRUCTS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h> /*for intx_t & uintx_t variables*/

/*Look up table to use in converting char to hext, (extern use in other files)*/
extern char hexTblCharAry[];

/*##############################################################################
# Struct-3: bigNum
# Use: Stores hex components of string as a big number (array of longs 64 bit)
##############################################################################*/
typedef struct bigNum
{ /*bigNum*/
    unsigned long
       *bigNumAryULng;
    unsigned char
       lenUsedElmChar,  /*Number of elements used in bigNumAryULng used*/
       lenAllElmChar;   /*Number of elemetns in bigNumAryULng*/
}bigNum;

/*##############################################################################
# Struct-1: readInfo
# Use: Stores the start location of a particler read, also connects to the
#      node it holds in the graph
##############################################################################*/
typedef struct readInfo
{ /*readInfo structer*/
    int8_t balanceChar;     /*Tells if the node is balanced*/
    struct bigNum *idBigNum; /*Holds read id as unique big number*/
    struct readInfo *leftChild, 
                    *rightChild;
}readInfo; /*readInfo structure*/

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    char *readIdCStr,         /*c-string with read name to copy*/
    const int32_t *lenCStrUInt /*Length of cString to convert*/
); /*Allocates memomory and makes a readInfo structer (variables set to 0)*/

void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
); /*frees a readInfo structer*/

/*##############################################################################
# Struct-2: readNodeStack
# Use: Used to make a stack of readInfo structers, used in transvering my
#      tree of reads
# Note: Users makes array of nodes they pass around, for a well balanced AVL
#       tree 73 nodes covers 10^18 nodes
##############################################################################*/
typedef struct readNodeStack
{ /*readNodeStack*/
    struct readInfo *readNode;                /*Node in stack*/
}readNodeStack; /*readNodeStack*/

void pushReadNodeStack(
    struct readNodeStack **readStackAry, /*Array of read info nodes*/
    struct readInfo *readNode            /*readInfo to assing to next node*/
); /*pushes a readNodeStack structer onto a readNodeStack stack*/

void popReadNodeStack(
    struct readNodeStack **readStackAry /*readInfo Array (stack) to pop*/
); /*frees readNodeStack & sets readNodeStack to next readInfo node in stack*/

/*##############################################################################
# Output:
#    returns: bigNum structer with the converted big number
#    returns: 0 if memory alloaction failed
##############################################################################*/
struct bigNum * makeBigNumStruct(
    char *cStrToCnvt,/*C-string to convert hex elements to big number*/
    const int32_t *lenCStrUInt /*Length of cString to convert*/
); /*Converts hex characters in c-string to a bitNum struct with a big number*/

/*##############################################################################
# Output:
#    Modifies: idBigNum to hold the converted big number.
#        - Sets idBigNum->lenUsedULng to 0 if memory reallocation failed
##############################################################################*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,  /*Holds the output big number*/
    char *cStrToCnvt,         /*C-string to convert to large number*/
    const int32_t *lenCStrUInt /*Length of cString to convert*/
); /*Flips c-string & converts to big number*/

/*##############################################################################
# Output:
#    frees: idBigNum & sets pointer to 0
##############################################################################*/
void freeBigNumStruct(
    struct bigNum **idBigNum  /*Address to bigNum structer to free*/
); /*Frees a bigNum struct from memory*/

/*##############################################################################
# Output:
#    Returns: 0 if both equal, < 0 if first is smaller, > 0 if first is bigger
##############################################################################*/
unsigned long cmpBigNums(
    struct bigNum *bigNumOne,
    struct bigNum *bigNumTwo
); /*Compares bigNumOne to bigNumTwo to see if equal, >, or <*/

/*##############################################################################
# Output:
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdULng to hold the length of the read id
#    Modifies: lenInputInt to hold length of input buffer
#    Returns: readInfo struct with bigNum struct having read id converted to hex
#        - 0 if fails or end of file (lenIdULng < buffSizeInt)
##############################################################################*/
struct readInfo * cnvtIdToBigNum(
    char *bufferCStr, /*buffer to hold fread input (can have data)*/
    uint32_t buffSizeInt,         /*Size of buffer to work on*/
    char **endNameCStr, /*Points to start of id, will point to end*/
    uint64_t *lenInputULng,        /*Length of input from fread*/
    unsigned char *lenBigNumChar, /*Holds size to make bigNumber*/
    FILE *idFILE          /*Fastq file to get data from*/
); /*Converts read id in buffer to bigNum read id, will grab new file input*/

#endif
