/*######################################################################
# Name: trimPrimersAvlTree
# Use: Holds functions to make an AVL tree of readPrim structures
# Includes:
#   - "trimPrimersStructs.h"
#   o "fqGetIdsStructs.h"
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#ifndef TRIMPRIMERSAVLTREE_H
#define TRIMPRIMERSAVLTREE_H

#include "trimPrimersStructs.h"
#define defLenStack 256 /*Length of readStack to use in tree searching*/

/*----------------------------------------------------------------------
| Output:
|   - Modifies:
|     o If node is not duplicate, inserts new readPrim into tree.
|     o If readPrim is a duplicate, then the primer coordinates are 
|       transfered to the orignal node in the tree.
|       - The primCord linked list in idToInsert is set to 0 for
|         duplicates
|   - Returns:
|      o 0 if found duplicate
|      o 1 if no duplicates
----------------------------------------------------------------------*/
uint8_t avlInsPrimReadST(
    struct readPrim *idToInsert,   /*node with read to insert in tree*/
    struct readPrim **readTree,    /*readPrim tree to insert read into*/
    struct readPrimStack *readStackAry/*Stack, (as array) forSearching*/
); /*Inserts readPrim structer (idToInsert) into tree an readPrim AVL.*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns: 
|     o Pointer: to node with read
|     o If the read name was not in the tree
\---------------------------------------------------------------------*/
struct readPrim * searchReadPrimTree(
    struct bigNum *queryIdBigNum, /*Read id to search for (big number)*/
    struct readPrim *readTree /*readInfo tree to search for readIdCStr*/
); /*Find a particler read name in the readPrim avl tree*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: balUC in readNode to be depth of its deepest child + 1
\---------------------------------------------------------------------*/
void updateReadPrimDepth(
    struct readPrim *readNode  /*node in tree of read id's to update*/
); /*Updates the depth of a single node by using its child nodes*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: balance of a node, - = left imbalance, + = right imbalance
\---------------------------------------------------------------------*/
int64_t getReadPrimBal(
    struct readPrim *readNode      /*node to get balance from*/
); /*Gets the balance of a single node in a tree*/

/*---------------------------------------------------------------------\
| Output: Modifies: tree to be rebalanced when nodes unbalanced
\---------------------------------------------------------------------*/
void rebalReadPrimTree(
    struct readPrimStack *readStackAry, /*stack of nodes to check*/
    struct readPrim **readTree  /*Root of tree. Changes if root swaped*/
); /*Uses readInfo nodes in readStack to reblance the tree*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readInfo tree by doing a right left reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlRLBal(
    struct readPrim *parNode
); /*Does a right left rebalance on an avl readPrim tree*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readPrim tree by doing a right right reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlRRBal(
    struct readPrim *parNode
); /*Does a right left rebalance on an avl readPrim tree*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readPrim tree by doing a left right reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlLRBal(
    struct readPrim *parNode
); /*Does a right left rebalance on an avl readPrim tree*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (readPrim pointer) in the rebalance
|    Modifies: readPrim tree by doing a left left reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlLLBal(
    struct readPrim *parNode
); /*Does a left left rebalance on an avl readPrim tree*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: the parent node of readOn to point to newPar
|    Pops: readPrim node off of readStack
\---------------------------------------------------------------------*/
void adjReadPrimTreePtr(
    struct readPrim *readOn,          /*Parent node before a rebalance*/
    struct readPrim *newPar,          /*Parent node after rebalance*/
    struct readPrimStack **readStackAry,/*Stack to get parent nodeFrom*/
    struct readPrim **readTree     /*root node to make sure not swaped*/
); /*Addjusts parent poniter after rebalancing a readPrim tree*/

/*---------------------------------------------------------------------\
| Output:
|   - Frees: every node in the readPrim tree & sets root node to 0
|   - Returns:
|     o 0 If malloc failed to allocate stack memory
|     o 1 if suceeded
|     o 4 If thier was nothing to free
\---------------------------------------------------------------------*/
uint8_t freeReadPrimTree(
    struct readPrim **readTree,  /*Root of readPrim tree to free*/
    struct readPrimStack *readStackAry /*Stack to use in freeing*/
); /*Frees a tree of readPrimNodes*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o readPrimList to be the head of a balanced tree of its nodes
|     o Removes any duplicate read ids in the list
| WARNING:
|   - readPrimList list must only use the rightChild pointer.
\---------------------------------------------------------------------*/
void readPrimListToTree(
    char listOnHeapBl,
        /*1 read id list on heap, ok to free duplicates, 0 do not free*/
    struct readPrim **readPrimList/*List of read ids to convert toTree*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: readPrimListToTree
   '   - Converts a list of readPrimt structures with read ids to a
   '     balanced tree. WARNING: This list must only use the rightChild
   '     pointer.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
