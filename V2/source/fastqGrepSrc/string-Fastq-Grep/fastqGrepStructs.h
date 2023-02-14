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
# Struct-1: readInfo
# Use: Stores the start location of a particler read, also connects to the
#      node it holds in the graph
##############################################################################*/
typedef struct readInfo
{ /*readInfo structer*/
    char *idCStr,               /*Stores the name/id of sequence*/
         balanceChar;           /*Tells if the node is balanced*/
    struct readInfo *leftChild, 
                    *rightChild;
}readInfo; /*readInfo structure*/

struct readInfo * makeReadInfoStruct(
    char *readNameCStr,          /*c-string with read name to copy*/
    unsigned int lenNameUInt,    /*length of readNameCStr*/
    char ignoreChar             /*Character at start to ignore ('' is nothing)*/
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

#endif
