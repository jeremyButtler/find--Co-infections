/*##############################################################################
# Name: fastqGrepAVLTree
# Use: Header file for clustGraphReadTree.c
#      Has functions to build a self blancing tree with read info nodes
##############################################################################*/

#ifndef FASTQGREPAVLTREE_C
#define FASTQGREPAVLTREE_C

#include <string.h> /*strcmp*/
#include "fastqGrepStructs.h"
    /*
      Includes: <stdlib.h>
                <stdio.h>
    */

/*##############################################################################
# Output:
#    Returns: pointer to found or new readInfo node (with readnameCStr)
#             0: If the read name was not in the tree
#    Modifes: readTree to have a new node with readNameCStr in it
##############################################################################*/
struct readInfo * findAddNodeToReadTree(
    char * readIdCStr,                 /*c-string with name of read to find*/
    unsigned char *numElmUChar,         /*Number of U longs needed per big num*/
    struct readInfo **readTree,        /*tree to search for readNameCStr*/
    struct readNodeStack *readStackAry /*Stack, (as array) for searching*/
); /*Finds or creates node with input read name in a read info tree*/

/*##############################################################################
# Output:
#    Modifies: readInfo tree by inserting node, if not duplicate
#    Returns: 0 if found duplicate, 1 if no duplicates
##############################################################################*/
char insertNodeIntoReadTree(
    struct readInfo *idToInsert,        /*node with read to insert in tree*/
    struct readInfo **readTree,         /*readInfo tree to insert read into*/
    struct readNodeStack *readStackAry  /*Stack, (as array) for searching*/
); /*Inserts read node in tree, if it is not a duplicate*/

/*##############################################################################
# Output:
#    Returns: 
#      Pointer: to node with read
#      0: If the read name was not in the tree
##############################################################################*/
struct readInfo * searchTree(
    struct bigNum *queryIdBigNum, /*Read id to search for (as big number)*/
    struct readInfo *readTree    /*readInfo tree to search for readNameCStr*/
); /*Find a particler read name in the avl tree*/

/*##############################################################################
# Output:
#    Modifies: balanceChar in readNode to be the depth of its deepest child + 1
##############################################################################*/
void updateDepth(
    struct readInfo *readNode      /*node in tree of read id's to update*/
); /*Updates the depth of a single node by using its child nodes*/

/*##############################################################################
# Output:
#    Returns: long with balance, - means left imbalanced, + right imbalenced
##############################################################################*/
long getBalance(
    struct readInfo *readNode      /*node to get balance from*/
); /*Gets the balance of a single node in a tree*/

/*##############################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
##############################################################################*/
void rebalanceTree(
    struct readNodeStack *readStackAry, /*stack of nodes to check*/
    struct readInfo **readTree          /*Root node of tree.
                                                      Changed when root swaped*/
); /*Uses readInfo nodes in readStack to reblance the tree*/

struct readInfo * readTreeRightLeftRebalance(
    struct readInfo *parNode
); /*Does a right left rebalance on an avl readInfo tree*/

struct readInfo * readTreeRightRightRebalance(
    struct readInfo *parNode
); /*Does a right left rebalance on an avl readInfo tree*/

struct readInfo * readTreeLeftRightRebalance(
    struct readInfo *parNode
); /*Does a right left rebalance on an avl readInfo tree*/

struct readInfo * readTreeLeftLeftRebalance(
    struct readInfo *parNode
); /*Does a left left rebalance on an avl readInfo tree*/

void readTreeAdjustParPtrAfterSwap(
    struct readInfo *readOn,          /*Parent node before a rebalance*/
    struct readInfo *newPar,          /*Parent node after rebalance*/
    struct readNodeStack **readStackAry, /*Stack to get the parent node from*/
    struct readInfo **readTree        /*root node to make sure not swaped*/
); /*Addjusts parent poniter after rebalancing a readInfo tree*/

/*##############################################################################
# Output:
#    Frees: every node in readTree & sets readTreeRoot to 0
#    Returns: 1: if suceeded
#             0: If malloc failed to allocate stack memory
##############################################################################*/
char freeReadTree(
    struct readInfo **readTree,  /*Root of readInfo tree to free*/
    struct readNodeStack *readStackAry /*Stack to use in freeing*/
); /*Frees a tree of readInfoNodes*/

/*##############################################################################
# Output:
#    int: with max depth of the input readInfo tree
##############################################################################*/
int getTreeDepth(
    struct readInfo *readTree
); /*Gets depth of a readInfo tree (just here for testing)*/

#endif
