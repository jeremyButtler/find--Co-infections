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

#include "trimPrimersAvlTree.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimPrimersAvlTree TOC:
'    fun-1 avlInsPrimReadST:
'      o Insert readinto node into readinfo tree
'    fun-2 searchReadPrimTree:
'      o Find if read is in a tree
'    fun-3 updateReadPrimDepth:
'      o Update depth on single node in tree using children
'    fun-4 getReadPrimBal:
'      o Get balance of single node in tree
'    fun-5 rebalanceTree:
'      o Rebalance a tree with a readPrimStack list
'    fun-6 readPrimAvlRRBal:
'      o do a right left rebalance on avl tree
'    fun-7 readPrimAvlRRBal:
'      o do a right rigth rebalance on avl tree
'    fun-8 readPrimAvlLRBal:
'      o do a left right rebalance on avl tree
'    Fun-9 readPrimAvlLLBal:
'      o do a left left rebalance on an avl tree
'    fun-10 adjReadPrimTreePtr:
'      o Set parent nodes ptr swaped node
'    fun-11 freeReadPrimTree:
'      o free a tree of readPrim nodes
'    fun-12 readPrimListToTree:
'      o Converts a list of readPrimt structures with read ids to a
'        balanced tree. WARNING: This list must only use the rightChild
'        pointer.
\---------------------------------------------------------------------*/

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
) /*Inserts readPrim structer (idToInsert) into tree an readPrim AVL.*/
{ /*avlInsPrimReadST*/
   
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ' Fun-1 TOC: alvInsPrimReadST
    '    fun-2 sec-1: variable declerations
    '    fun-2 sec-2: Search tree to see if is a new or duplicate node
    '    fun-2 sec-3: If duplicate, transfer primCord linked list
    '    fun-2 sec-4: Else add node to tree and reblance
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    long matchL = 0;               /*Tells if read names were a match*/
    struct readPrim *refNode = *readTree; /*root node of tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Search tree to see if is a new or duplicate node
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readStackAry->readNode = 0; /*make sure 0 at start*/

    if(*readTree == 0)
    { /*If there is no tree*/
        *readTree = idToInsert;
        return 1;
    } /*If there is no tree*/

    while(refNode != 0)
    { /*Loop till at a leaf node or found node*/
        /*Compare id numbers to see if the same*/
        matchL = cmpBigNums(idToInsert->idBigNum, refNode->idBigNum);

        if(matchL == 0)
            break;         /*Found a match*/

        /*Add node on to the history of visited nodes*/
        pushReadPrimStack(&readStackAry, refNode);

        if(matchL < 0)         /*If new id < older id, move left*/
            refNode = refNode->leftChild;
        else                       /*else, new id > old id, move right*/
            refNode = refNode->rightChild;
    } /*Loop till at a leaf node or found node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: If duplicate, transfer primCord linked list
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     if(matchL == 0)
     { /*If this was a match*/
        if(idToInsert->primCordST == 0)
            return 0;    /*Do not need to transfer a linked list*/
 
        if(refNode->primCordST == 0)
        { /*If the reference has no coordinates, then swap*/
            refNode->primCordST = idToInsert->primCordST;
            /*Prevent user from freeing list for duplicates*/
            idToInsert->primCordST = 0;

            return 0;
        } /*If the reference has no coordinates, then swap*/

        /*Transfer the primer coordinates over*/
        insPrimCordST(idToInsert->primCordST, &refNode->primCordST);

         /*So the user does not free coordiantes when they free the 
           duplicate node*/
        idToInsert->primCordST = 0;
        return 0;
     } /*If this was a match*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-4: Else add node to tree and reblance
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(matchL < 0)
        readStackAry->readNode->leftChild = idToInsert;
    else
        readStackAry->readNode->rightChild = idToInsert;

    rebalReadPrimTree(readStackAry, readTree);
        /*Check if need to rebalance*/
    return 1;
} /*insertNodeIntoReadTree*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns: 
|     o Pointer: to node with read
|     o If the read name was not in the tree
\---------------------------------------------------------------------*/
struct readPrim * searchReadPrimTree(
    struct bigNum *queryIdBigNum, /*Read id to search for (big number)*/
    struct readPrim *readTree /*readPrim tree to search for readIdCStr*/
) /*Find a particler read name in the readPrim avl tree*/
{ /*searchReadPrimTree*/
   
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: Sec-1 Sub-1: searchReadPrimTree
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    long matchL = 0;                   /*Holds if read names were match*/

    while(readTree != 0)
    { /*Loop till at a leaf node or found node*/

        /*Compare id numbers to see if the same*/
        matchL = cmpBigNums(queryIdBigNum, readTree->idBigNum);

        if(matchL == 0)
            return readTree;                 /*Return match*/
        else if(matchL < 0)                  /*query < ref, move left*/
            readTree = readTree->leftChild;
        else                                 /*query > ref, move right*/
            readTree = readTree->rightChild;
    } /*Loop till at a leaf node or found node*/

    return 0;                               /*Read is not in the tree*/
} /*searchReadPrimTree*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: balUC in readNode to be depth of its deepest child + 1
\---------------------------------------------------------------------*/
void updateReadPrimDepth(
    struct readPrim *readNode  /*node in tree of read id's to update*/
) /*Updates the depth of a single node by using its child nodes*/
{ /*updateReadPrimDepth*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-3 TOC: Sec-1 Sub-1: updateReadPrimDepth
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readNode->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readNode->leftChild == 0)
            readNode->balUC = readNode->rightChild->balUC + 1;
        else if(readNode->rightChild->balUC >readNode->leftChild->balUC)
            readNode->balUC = readNode->rightChild->balUC + 1;
        else
            readNode->balUC = readNode->leftChild->balUC + 1;
    } /*if there is a right child, check if there is a left child*/

    else if(readNode->leftChild != 0)                /*only left child*/
        readNode->balUC = readNode->leftChild->balUC + 1;

    else
        readNode->balUC = 0; /*No children, is a leafnode*/

    return;
} /*updateReadPrimDepth*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: balance of a node, - = left imbalance, + = right imbalance
\---------------------------------------------------------------------*/
int64_t getReadPrimBal(
    struct readPrim *readNode      /*node to get balance from*/
) /*Gets the balance of a single node in a tree*/
{ /*getReadPrimBal*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-4 TOC: Sec-1 Sub-1: getReadPrimBal
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(readNode->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readNode->leftChild == 0)
            return readNode->rightChild->balUC + 1;

      return readNode->rightChild->balUC - readNode->leftChild->balUC+1;
    } /*if there is a right child, check if there is a left child*/

    else if(readNode->leftChild != 0)  /*only left child*/
        return -1 * (readNode->leftChild->balUC + 1);

    return 0;
} /*getReadPrimBal*/

/*---------------------------------------------------------------------\
| Output: Modifies: tree to be rebalanced when nodes unbalanced
\---------------------------------------------------------------------*/
void rebalReadPrimTree(
    struct readPrimStack *readStackAry, /*stack of nodes to check*/
    struct readPrim **readTree  /*Root of tree. Changes if root swaped*/
) /*Uses readPrim nodes in readStack to reblance the tree*/
{ /*rebalReadPrimTree*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-5 TOC: rebalReadPrimTree
   '   o fun-5 sec-1: General set up
   '   o fun-5 sec-2: Deal with right imbalences
   '   o fun-5 sec-3: Deal with left imbalences
   '   o fun-5 sec-4: No major imbalance, just update depths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ^ Fun-5 Sec-1: General set up
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   long balanceLng = 0;  /*Balance of tree*/
   struct readPrim *swapNode = 0;

   while(readStackAry->readNode != 0)
   { /*while not at the root of the readPrim tree, adjust balance*/
      /*Get balance of the node*/
      balanceLng = getReadPrimBal(readStackAry->readNode);

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun-5 Sec-2: Deal with right imbalences
      ^   o fun-5 sec-2 sub-1: Do reblances and pointer swaps
      ^   o fun-5 sec-2 sub-2: Reblance the tree
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /****************************************************************\
      * Fun-5 Sec-2 sub-1: Do reblances and pointer swaps
      \****************************************************************/

      if(balanceLng > 1)
      { /*If have more nodes on the right side (right rebalance)*/
          if(getReadPrimBal(readStackAry->readNode->rightChild) < 0)
             swapNode = readPrimAvlRLBal(readStackAry->readNode);

          else swapNode = readPrimAvlRRBal(readStackAry->readNode);

          adjReadPrimTreePtr(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              readTree
          ); /*Adjust pointers after swap*/

          popReadPrimStack(&readStackAry); /*Move to next node*/

          /************************************************************\
          * Fun-5 Sec-2 sub-2: Reblance the tree
          \************************************************************/

          while(readStackAry->readNode != 0)
          { /*While have nodes to update depths for*/
              if(readStackAry->readNode->leftChild == 0)
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->rightChild->balUC +
                      1;
              else if(readStackAry->readNode->rightChild == 0)
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->leftChild->balUC +
                      1;
              else if(
                  readStackAry->readNode->leftChild->balUC <
                  readStackAry->readNode->rightChild->balUC
              )
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->rightChild->balUC +
                      1;
              else
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->leftChild->balUC +
                      1;

              popReadPrimStack(&readStackAry); /*Move to next node*/
          } /*While have nodes to update depths for*/

          return;
      } /*If have more nodes on the right side (right rebalance)*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun-5 Sec-3: Deal with left imbalences
      ^   o fun-5 sec-3 sub-1: Do reblances and pointer swaps
      ^   o fun-5 sec-3 sub-2: Reblance the tree
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /****************************************************************\
      * Fun-5 Sec-2 sub-3: Do reblances and pointer swaps
      \****************************************************************/

      else if(balanceLng < -1)
      { /*else if have more nodes on the left side (left rebalance)*/
          if(getReadPrimBal(readStackAry->readNode->leftChild) > 0)
             swapNode = readPrimAvlLRBal(readStackAry->readNode);
          else
             swapNode = readPrimAvlLLBal(readStackAry->readNode);

          adjReadPrimTreePtr(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              readTree
          ); /*Adjust pointers after swap*/

          popReadPrimStack(&readStackAry); /*Move to next node*/

          /************************************************************\
          * Fun-5 Sec-2 sub-2: Reblance the tree
          \************************************************************/

          while(readStackAry->readNode != 0)
          { /*While have nodes to update depths for*/
              if(readStackAry->readNode->rightChild == 0)
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->leftChild->balUC +
                      1;
              else if(readStackAry->readNode->leftChild == 0)
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->rightChild->balUC +
                      1;
              else if(
                  readStackAry->readNode->leftChild->balUC >
                  readStackAry->readNode->rightChild->balUC
              )
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->leftChild->balUC +
                      1;
              else
                  readStackAry->readNode->balUC = 
                      readStackAry->readNode->rightChild->balUC +
                      1;

              popReadPrimStack(&readStackAry); /*Move to next node*/
          } /*While have nodes to update depths for*/

          return;
      } /*else if have more nodes on the left side (left rebalance)*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun-5 Sec-4: No major imbalance, just update depths
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      else if(balanceLng < 0) /*Left child is deeper*/
         readStackAry->readNode->balUC =
            readStackAry->readNode->leftChild->balUC + 1;

      else if(balanceLng > 0) /*right child is deeper*/
         readStackAry->readNode->balUC =
            readStackAry->readNode->rightChild->balUC + 1;

      else if(readStackAry->readNode->leftChild == 0)
         readStackAry->readNode->balUC = 0; /*Is a leaf node*/

      else /*Both childern are equally balanced*/
         readStackAry->readNode->balUC =
            readStackAry->readNode->rightChild->balUC + 1;

      popReadPrimStack(&readStackAry); /*Move to next node*/
   } /*while not at the root of the readPrim tree, adjust balance*/

   return;
} /*rebalReadPrimTree*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readPrim tree by doing a right left reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlRLBal(
    struct readPrim *parNode
) /*Does a right left rebalance on an avl readPrim tree*/
{ /*readPrimAvlRLBal*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-6 TOC: do a right left rebalance on an avl tree
    '    fun-6 sec-1: variable declerations
    '    fun-6 sec-2: Do swap
    '    fun-6 sec-3: Adjust balance
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-1: Variable decerlations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readPrim *childNode = parNode->rightChild;
     struct readPrim *swapNode = childNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-2: Do swap
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->leftChild = swapNode->rightChild;
     parNode->rightChild = swapNode->leftChild; 
     swapNode->leftChild = parNode;
     swapNode->rightChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-3: Adjust balance
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /* Update depths (swap node was deepest depth, so can use to find
        the updated depths*/
     if(parNode->leftChild != 0)
        parNode->balUC = parNode->leftChild->balUC + 1;
     else
        parNode->balUC = swapNode->balUC;

     childNode->balUC = swapNode->balUC;
     (swapNode->balUC)++;

     return swapNode;
} /*readPrimAvlRLBal*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readPrim tree by doing a right right reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlRRBal(
    struct readPrim *parNode
) /*Does a right left rebalance on an avl readPrim tree*/
{ /*readPrimAvlRRBal*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-7 TOC: do a right rigth rebalance on an avl tree
    '   o fun-7 sec-1: variable declerations
    '   o fun-7 sec-2: Do swap
    '   o fun-7 sec-3: Adjust balance
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-1: Variable decerlations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readPrim *childNode = parNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-2: Do swap
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->rightChild = childNode->leftChild;
     childNode->leftChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-3: Adjust balance
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     if(parNode->rightChild != 0)
     { /*if the paranet had a chlld, child sets balance*/
         parNode->balUC = parNode->rightChild->balUC + 1;
         childNode->balUC = parNode->balUC + 1;
     } /*if the paranet had a chlld, child sets balance*/
     else
         parNode->balUC = childNode->leftChild->balUC + 1;
         /*Childs balance is unchanged*/

     return childNode;
} /*readPrimAvlRRBal*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (pointer) in the rebalance
|    Modifies: readPrim tree by doing a left right reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlLRBal(
    struct readPrim *parNode
) /*Does a right left rebalance on an avl readPrim tree*/
{ /*readPrimAvlRLBal*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-8 TOC: do a left right rebalance on an avl tree
    '   o fun-8 sec-1: variable declerations
    '   o fun-8 sec-2: Do swap
    '   o fun-8 sec-3: Adjust balance
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-8 Sec-1: Variable decerlations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readPrim *childNode = parNode->leftChild;
     struct readPrim *swapNode = childNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-8 Sec-2: Do swap
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->rightChild = swapNode->leftChild;
     parNode->leftChild = swapNode->rightChild; 
     swapNode->rightChild = parNode;
     swapNode->leftChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-8 Sec-3: Adjust balance
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /* Update depths (can do from swap node depth, since swap node was
       deepest node) */

     if(parNode->rightChild != 0)
        parNode->balUC = parNode->rightChild->balUC + 1;
     else
        parNode->balUC = swapNode->balUC;

     childNode->balUC = swapNode->balUC;
     (swapNode->balUC)++;

     return swapNode;
} /*readPrimAvlLRBal*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The new parent node (readPrim pointer) in the rebalance
|    Modifies: readPrim tree by doing a left left reblance
\---------------------------------------------------------------------*/
struct readPrim * readPrimAvlLLBal(
    struct readPrim *parNode
) /*Does a left left rebalance on an avl readPrim tree*/
{ /*readPrimAvlLLBal*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-9 TOC: do a left left rebalance on an avl tree
    '   o fun-9 sec-1: variable declerations
    '   o fun-9 sec-2: Do swap
    '   o fun-9 sec-3: Adjust balance
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-1: Variable decerlations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readPrim *childNode = parNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-2: Do swap
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->leftChild = childNode->rightChild;
     childNode->rightChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-3: Adjust balance
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     if(parNode->rightChild != 0)
     { /*if the paranet had a chlld, child sets balance*/
         parNode->balUC = parNode->rightChild->balUC + 1;
         childNode->balUC = parNode->balUC + 1;
     } /*if the paranet had a chlld, child sets balance*/
     else
         parNode->balUC = childNode->leftChild->balUC + 1;
         /*Childs balance is unchanged*/
     
     return childNode;
} /*readPrimAvlLLBal*/

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
) /*Addjusts parent poniter after rebalancing a readPrim tree*/
{ /*adjReadPrimTreePtr*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-10 TOC: Sec-1 Sub-1: adjReadTreePtr
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(readOn == *readTree) /*Handle swapping the root node*/
        *readTree = newPar;
    else
    { /*Else need to adjust the parent nodes pointer*/
        popReadPrimStack(readStackAry);

        if((*readStackAry)->readNode->rightChild == readOn)
        { /*If was a right shift or right left shift*/
            (*readStackAry)->readNode->rightChild = newPar;

            if((*readStackAry)->readNode->leftChild == 0)
                (*readStackAry)->readNode->balUC = newPar->balUC + 1;
            else if(
                newPar->balUC >
                (*readStackAry)->readNode->leftChild->balUC
            )/*Else if the right hild newPar is greater*/
                (*readStackAry)->readNode->balUC = newPar->balUC + 1;
            else
                (*readStackAry)->readNode->balUC = 
                    (*readStackAry)->readNode->leftChild->balUC + 1;
        } /*If was a right shift or right left shift*/

        else
        { /*Else was a left of left right swap*/
            (*readStackAry)->readNode->leftChild = newPar;

            if((*readStackAry)->readNode->rightChild == 0)
                (*readStackAry)->readNode->balUC = newPar->balUC + 1;
            else if(
                newPar->balUC >
                (*readStackAry)->readNode->rightChild->balUC
            )/*Else if the left hild newPar is greater*/
                (*readStackAry)->readNode->balUC = newPar->balUC + 1;
            else
                (*readStackAry)->readNode->balUC = 
                    (*readStackAry)->readNode->rightChild->balUC + 1;
        } /*Else was a left of left right swap*/

    } /*Else need to adjust the parent nodes pointer*/

    return;
} /*adjReadPrimTreePtr*/

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
) /*Frees a tree of readPrimNodes*/
{ /*freeReadPrimTree*/
   
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-11 TOC: Sec-1 Sub-1: freeReadPrimTree
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readPrim *deleteNode = *readTree;
    readStackAry++; /*Get off the first node*/

    if(*readTree == 0)
        return 4;

    pushReadPrimStack(&readStackAry, 0); /*Make sure their is a 0*/
    pushReadPrimStack(&readStackAry, *readTree);

    while(readStackAry->readNode != 0)
    { /*While there are readPrim nodes in the tree*/
        deleteNode = readStackAry->readNode;
        popReadPrimStack(&readStackAry); /*Free my stack*/ 

        if(deleteNode == 0)
            continue;
        if(deleteNode->leftChild != 0)
            pushReadPrimStack(&readStackAry, deleteNode->leftChild);
        if(deleteNode->rightChild != 0)
            pushReadPrimStack(&readStackAry, deleteNode->rightChild);

        freeReadPrimST(&deleteNode);
    } /*While there are readPrim nodes in the tree*/

    return 1;
} /*freeReadPrimTree*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: readPrimListToTree
   '   - Converts a list of readPrimt structures with read ids to a
   '     balanced tree. WARNING: This list must only use the rightChild
   '     pointer.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   unsigned char errUC = 0;
   struct readPrim *tmpRead = 0;
   struct readPrim *nextRead = 0;
   struct readPrim *rootPrim = 0;
   struct readPrimStack readStack[defLenStack]; /*For freeing*/

   /*Make sure start & end of my stacks are marked*/
   readStack[0].readNode = 0;
   readStack[defLenStack - 1].readNode = 0;

   tmpRead = *readPrimList;
   nextRead = 0;

   while(tmpRead != 0)
   { /*While I have read ids to add to the AVL tree*/
       nextRead = tmpRead->rightChild;
       errUC = avlInsPrimReadST(tmpRead, &rootPrim, readStack);

       if(errUC == 0 && listOnHeapBl & 1) freeReadPrimST(&tmpRead);

       tmpRead = nextRead;
   } /*While I have read ids to add to the AVL tree*/

   *readPrimList = rootPrim; /*Set root of the avl tree*/
   return;
} /*readPrimListToTree*/
