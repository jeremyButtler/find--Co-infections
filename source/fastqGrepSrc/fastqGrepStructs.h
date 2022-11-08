/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepStructs TOC:
#    struct-1: readInfo: Holds the read name
#    struct-2: readNodeStack: Makes a stack (filo) of read nodes
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef FASTQGREPSTRUCTS_H
#define FASTQGREPSTRUCTS_H

#include <stdlib.h>
#include <stdio.h>

/*##############################################################################
# Struct-3: bigNum
# Use: Stores hex components of string as a big number (array of longs 64 bit)
##############################################################################*/
typedef struct bigNum
{ /*bigNum*/
    unsigned long
       *bigNumAryULng;
    char
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
    char balanceChar;           /*Tells if the node is balanced*/
    struct bigNum *idBigNum; /*Holds read id as unique big number*/
    struct readInfo *leftChild, 
                    *rightChild;
}readInfo; /*readInfo structure*/

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    char *readIdCStr,         /*c-string with read name to copy*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
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
    char *cStrToCnvt,        /*C-string to convert hex elements to big number*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
); /*Converts hex characters in c-string to a bitNum struct with a big number*/

/*##############################################################################
# Output:
#    Modifies: idBigNum to hold the converted big number.
#        - Sets idBigNum->lenUsedULng to 0 if memory reallocation failed
##############################################################################*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,  /*Holds the output big number*/
    char *cStrToCnvt,         /*C-string to convert to large number*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
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

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Output:
#    Returns: unsigned char with the number of unsigned longs needed to hold
#             the hex c-string
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
unsigned char cnvtStrLenToNumHexULng(
    const unsigned long *lenHexCStrULng /*Number of characters in hex c-string*/
); /*Converts the number of unsinged longs needed to store the hex characters*/

#endif
