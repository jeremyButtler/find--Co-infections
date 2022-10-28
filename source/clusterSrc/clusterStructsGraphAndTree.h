/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clustGraphAndTreeStructs TOC:
#    struct-1: graphNode
#    struct-2: readInfo
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef CLUSTERSTRUCTSGRAPHANDTREE_H
#define CLUSTERSTRUCTSGRAPHANDTREE_H

#include <stdlib.h>
#include <stdio.h>

struct graphNode;
struct readInfo;

/*##############################################################################
# Struct-1: graphNode
# Use: is a node in a graph of sequences
##############################################################################*/
typedef struct graphNode
{ /*graphNode structure*/
    unsigned long numNodesULng,  /*Number of nodes in cluster*/
                  readCntULng;   /*Number of edges a read shares with cluster*/
                 /*numNodes and readCnt are only incurmented for the head node*/

    struct readInfo *seqIdNode; /*Links to readInfo with seqId (for printing)*/

    struct graphNode *clustNode,  /*Holds head of the cluster read is in*/
                     *edgesList,  /*List of all nodes not connected previously*/
                     *nextNode,   /*Holds neighboring graphNode in parent node*/
                     *prevNode;/*previous neighbor in list (so can move nodes)*/
}graphNode; /*graphNode structure*/

/*##############################################################################
# Struct-2: readInfo
# Use: Stores the start location of a particler read, also connects to the
#      node it holds in the graph
##############################################################################*/
typedef struct readInfo
{ /*readInfo structer*/
    char *idCStr,               /*Stores the name/id of sequence*/
         balanceChar,           /*Tells if the node is balanced*/
         doneChar;              /*Tells if have already done read*/
    struct readInfo *leftChild, 
                    *rightChild;
    struct graphNode *nodeInGraph; /*Location of each alignment for this read*/
}readInfo; /*readInfo structure*/

/*##############################################################################
# Output:
#    Modifes: The edgeList in insertNode to point to the new node
#    Modifes: seqIdNode.idCStr to point to the new node (if idCStr = 0)
#    Incurments: clustNode.numNodesULng in insert node
#    Returns: new graphNode or 0 (if malloc failed to assign memory)
##############################################################################*/
struct graphNode * pushReadIntoGraph(
    struct graphNode **insertNode, /*graphNode to push read onto*/
    struct readInfo *seqIdNode    /*ReadInfo struct with sequence id for read*/
); /*Allocates memomory and makes a graphNode structer (variables set to 0)*/

void freeGraphNodeStruct(
    struct graphNode **graphNodeStruct /*Stucter to free*/
); /*Frees a graphNode structer*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Struct-2 functions: readInfo
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

struct readInfo * makeReadInfoStruct(
    char *readNameCStr, /*c-string with read name to copy*/
    unsigned char lenNameUChar /*length of readNameCStr*/
); /*Allocates memomory and makes a readInfo structer (variables set to 0)*/

void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
); /*frees a readInfo structer*/

#endif
