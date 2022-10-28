/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphStackStructs.h
#    struct-1: graphNodeStack
#    struct-2: readNodeStack
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef CLUSTERSTRUCTSSTACKS_H
#define CLUSTERSTRUCTSSTACKS_H

#include "clusterStructsGraphAndTree.h" /*Includes: <stdlib.h> <stdio.h>*/

/*##############################################################################
# Struct-1: graphNodeStack
# Use: Makes a stack of graph nodes (this is really a que)
##############################################################################*/
typedef struct graphNodeStack
{ /*graphNodeStack*/
    struct graphNode *graphStruct;    /*Graph node pointing to*/
    struct graphNodeStack *prevNode;
}graphNodeStack; /*graphNodeStack*/

struct graphNodeStack * pushGraphNodeStack(
    struct graphNode *graphNodeStruct,       /*graphNode to push into stack*/
    struct graphNodeStack **tailOfGraphStack /*graphNode stack push into*/
); /*Push graphNode structer onto a graphNode stack*/

void popGraphStackStruct(
    struct graphNodeStack **graphStack /*graphNode stack to pop graphNode off*/
); /*Frees a graph node stack structer and returns the next structer*/

/*##############################################################################
# Struct-2: readNodeStack
# Use: Used to make a stack of readInfo structers, used in transvering my
#      tree of reads
##############################################################################*/
typedef struct readNodeStack
{ /*readNodeStack*/
    char whichChildChar;           /*if node is root, left, or right child*/
    struct readInfo *readNode;                /*Node in stack*/
    struct readNodeStack *previousNode; /*Previous node in stack*/
}readNodeStack; /*readNodeStack*/

struct readNodeStack * pushReadNodeStack(
    struct readNodeStack **readStack, /*Stack to push readInfo node onto*/
    struct readInfo *readNode,        /*readInfo node to push into stack*/
    char whichChildChar             /*-1: path is down left child of readNode
                                       1: path is down right child of readNode*/
); /*pushes a readNodeStack structer onto a readNodeStack stack*/

void popReadNodeStack(
    struct readNodeStack **readStack /*readInfo stack to pop top readInfo off*/
); /*frees readNodeStack & sets readNodeStack to next readInfo node in stack*/

#endif
