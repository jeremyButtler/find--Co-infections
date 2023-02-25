/*######################################################################
# Name: findCoInftBinTree
# Use:
#    - Has functions & structers to build an AVL tree with readBin nodes
# Includes:
#    - <string.h>
#    - <stdlib.h>
#    - <stdio.h>
######################################################################*/

#ifndef FINDCOINFTBINTREE_H
#define FINDCOINFTBINTREE_H

#include <string.h> /*strcmp*/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

/*######################################################################
# Struct-1: readBin
# Use: Stores the start location of a particler read, also connects to
#      the node it holds in the graph
######################################################################*/
typedef struct readBin
{ /*readBin structer*/
    char
        refIdCStr[100],        /*Holds reference name*/
        fqPathCStr[100],       /*Holds fastq file name for this ref*/
        statPathCStr[100],     /*Holds stat file name for this ref*/
        bestReadCStr[100],     /*Holds file with the best read*/
        topReadsCStr[100],     /*Holds file name of top reads fastq*/
        consensusCStr[100];    /*Holds file name of consensus fasta*/

    int8_t
        balUChar;              /*Tells if the node is balanced*/

    unsigned long
        numReadsULng;     /*Number of reads in this bin*/ 

    struct readBin
        *leftChild, 
        *rightChild;
}readBin; /*readBin structure*/

/*######################################################################
# Struct-2: readNodeStack
# Use: Used to make a stack of readBin structers, used in transvering my
#      tree of reads
# Note: Users makes array of nodes they pass around, for a well balanced
#       AVL tree 73 nodes covers 10^18 nodes
######################################################################*/
typedef struct readBinStack
{ /*readNodeStack*/
    struct readBin *readNode;                /*Node in stack*/
}readBinStack; /*readNodeStack*/

/*######################################################################
# Output:
#    returns:
#        - pointer to readBin structure if sucessfull
#        - 0: if had a memory allocation error
######################################################################*/
struct readBin * makeReadBin(
    char *refIdCStr,   /*Name of reference, tab or null deliminated*/
    char *fqPathCStr,  /*Name of fastq bin, ends with '\0'*/
    char *statPathCStr /*Name of file holding stats, ends with '\0'*/
); /*Makes a readBin structure*/

/*######################################################################
# Output:
#    Frees: the input readBin structer & all its child nodes
#    Sets: binToFree to null
######################################################################*/
void freeReadBin(
    struct readBin **binToFree
); /*Frees a readBin structure & its children readBin structure*/

/*######################################################################
# Output:
#    Modifes readStack to piont to last element
######################################################################*/
void pushReadBinStack(
    struct readBinStack **readStackAry, /*Array of read info nodes*/
    struct readBin *readBin        /*readBin to assing to next node*/
); /*pushes a readBinStack structer onto a readBinStack stack*/

/*######################################################################
# Output: Modifies: readBinStack to point to next readBin node in stack
######################################################################*/
void popReadBinStack(
    struct readBinStack **readStackAry /*readBin Array (stack) to pop*/
); /*frees readBinStack & sets readBinStack to next node in stack*/

/*######################################################################
# Output:
#    Returns: pointer to found or new readBin node (with readnameCStr)
#             0: If the read name was not in the tree
#    Modifes: binTree to have a new node with refIdCStr in it
######################################################################*/
struct readBin * insBinIntoTree(
    char *refIdCStr,     /*read id to insert ('\0' or '\t' at end)*/
    char *fqPathCStr,     /*Fastq file name of bin ('\0' at end)*/
    char *statPathCStr,   /*Name of stat file for bin ('\0' at end)*/
    struct readBin **binTree, /*tree to search for refIdCStr*/
    struct readBinStack *readStackAry /*Stack, for searching*/
); /*Finds or creates node with input read name in a read info tree*/

/*######################################################################
# Output:
#    Modifies: balUChar in readBin to be the depth of its deepest
#              child + 1
######################################################################*/
void getBinDepth(
    struct readBin *readBin      /*node in tree of read id's to update*/
); /*Updates the depth of a single node by using its child nodes*/

/*######################################################################
# Output:
#    Returns:
#        - negative number: for left imbalanced
#        - positive number: for right imbalenced
#        - zero: for balanced
######################################################################*/
long getBinTreeBal(
    struct readBin *readBin      /*node to get balance from*/
); /*Gets the balance of a single node in a tree*/

/*######################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
######################################################################*/
void rebalBinTree(
    struct readBinStack *readStackAry, /*stack of nodes to check*/
    struct readBin **binTree    /*Root of tree. Changes if root swaped*/
); /*Uses readBin nodes in readStack to reblance the tree*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a right left reblance
######################################################################*/
struct readBin * binRightLeftRebal(
    struct readBin *parNode
); /*Does a right left rebalance on an avl readBin tree*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a right right reblance
######################################################################*/
struct readBin * binRightRightRebal(
    struct readBin *parNode
); /*Does a right left rebalance on an avl readBin tree*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a left right reblance
######################################################################*/
struct readBin * binLeftRightRebal(
    struct readBin *parNode
); /*Does a right left rebalance on an avl readBin tree*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a left left reblance
######################################################################*/
struct readBin * binTreeleftLeftRebal(
    struct readBin *parNode
); /*Does a left left rebalance on an avl readBin tree*/

/*######################################################################
# Output:
#    Modifies: the parent node of readOn to point to newPar
#    Pops: readBin node off of readStack
#######################################################################*/
void binTreeAdjPtrs(
    struct readBin *readOn,          /*Parent node before a rebalance*/
    struct readBin *newPar,          /*Parent node after rebalance*/
    struct readBinStack **readStackAry,/*Stack to get parent node from*/
    struct readBin **binTree       /*root node to make sure not swaped*/
); /*Addjusts parent poniter after rebalancing a readBin tree*/

/*######################################################################
# Output:
#    Frees: every node in binTree & sets binTreeRoot to 0
#    Returns: 4: If nothing to free
#             1: if suceeded
#             0: If malloc failed to allocate stack memory
######################################################################*/
char freeBinTree(
    struct readBin **binTree   /*Root of readBin tree to free*/
); /*Frees a tree of readBinNodes*/

/*######################################################################
# Output:
#    int: with max depth of the input readBin tree
######################################################################*/
int getBinTreeDepth(
    struct readBin *binTree
); /*Gets depth of a readBin tree (just here for testing)*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The last readBin in the list (for recurisive call)
|    Note: Recursive nature will slow down, but speed is not an issue
|          here
\---------------------------------------------------------------------*/
struct readBin * cnvtBinTreeToList(
    struct readBin **binTree /*Tree to convert to a linked list*/
); /*Converts a readBin tree to a linked list*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - numBinsULng to hold the number of bins
|    Note:
|        - This is a debugging function, so speed is not a big deal
\---------------------------------------------------------------------*/
void getNumBins(
    uint64_t *numBinsULng,  /*Will hold the number of bins in the tree*/
    struct readBin *binTree /*Tree of bins to count*/
); /*Get number of bins an a readBin tree*/

/*---------------------------------------------------------------------\
| Output:
|    frees:
|        - binToRm from list and sets next cluster in bin (right child)
|          or next bin (left child) as the new head
| Note:
|    - Each cluster should only have right children (All left children
|      are null (0))
|    - Each bin may have a link to the next bin (left child) & may have
|      a link to a cluster (right child). Empty pointers are null (0).
|
|        bin_1----------->bin_2-----------...>bin_n---------->null
|         |                |                   |
|         v                v                   v
|      cluster_1->null  cluster_1->null     cluster_1->null
|         |                |                   |
|         v                v                   v
|      cluster_2->null  cluster_2->null     cluster_2->null
|         .                .                   .
|         .                .                   .
|         .                .                   .
|      cluster_n->null  cluster_n->null     cluster_n->null
\---------------------------------------------------------------------*/
void rmBinFromList(
    struct readBin **binToRm
); /*Removes a readBin from a leftChild linked list of read bins*/

/*---------------------------------------------------------------------\
| Output:
|    Deletes: All files in the bin
\---------------------------------------------------------------------*/
void binDeleteFiles(
    struct readBin *binToWipe /*Bin to delete files from*/
); /*Deletes all files in a readBin*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - binToKeep->fqPathCStr to hold binToMerge reads
|    Deletes:
|        - File named after binToMerge->fqPathCStr
\---------------------------------------------------------------------*/
void mergeBins(
    struct readBin *binToKeep,  /*Bin to merge into*/
    struct readBin *binToMerge /*Bin to merge into binToKeep*/
); /*Merge two readBins together into on bin*/

#endif
