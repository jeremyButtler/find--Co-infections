/*##############################################################################
# Name: clustGraphReadTree.h
# Use: Header file for clustGraphReadTree.c
#      Has functions to build a self blancing tree with read info nodes
##############################################################################*/

#ifndef CLUSTERREADTREE_H
#define CLUSTERREADTREE_H

#include <string.h> /*strcmp*/
#include "clusterStructsStacks.h"
    /*
      Includes: clustGraphGraphAndTreeStructs.h
                    - readInfo, bridgeNode, graphNode structers and functions
    */

/*##############################################################################
# Output:
#    Returns: pointer to found or new readInfo node (with readnameCStr)
#             0: If malloc failed to allocate memory
#    Modifes: readTree to have a new node with readNameCStr in it
##############################################################################*/
struct readInfo * findAddNodeToReadTree(
    char * readNameCStr,        /*c-string with name of read to find*/
    unsigned char lenNameUChar, /*length of readNameCStr*/
    struct readInfo **readTree  /*readInfo tree to search for readNameCStr*/
); /*Finds or creates node with input read name in a read info tree*/

/*##############################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
##############################################################################*/
void rebalanceTree(struct readNodeStack **readStack, /*stack of nodes to check*/
                   struct readInfo **readTree        /*Root node of tree.
                                                      Changed when root swaped*/
); /*Uses readInfo nodes in readStack to reblance the tree*/

/*##############################################################################
# Output:
#    Frees: every node in readTree & sets readTreeRoot to 0
#    Returns: 1: if suceeded
#             0: If malloc failed to allocate stack memory
##############################################################################*/
char freeReadTree(
    struct readInfo **readTree  /*Root of readInfo tree to free*/
); /*Frees a tree of readInfoNodes*/

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
    struct readNodeStack **readStack, /*Stack to get the parent node from*/
    struct readInfo **readTree        /*root node to make sure not swaped*/
); /*Addjusts parent poniter after rebalancing a readInfo tree*/

#endif
