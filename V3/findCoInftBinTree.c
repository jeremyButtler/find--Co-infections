/*######################################################################
# Name: findCoInftBinTree
# Use:
#    - Has functions & structers to build an AVL tree with readBin nodes
# Includes:
#    - <string.h>
#    - <stdlib.h>
#    - <stdio.h>
######################################################################*/

#include "findCoInftBinTree.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' clustGraphBinTree TOC:
'    fun-1 makeReadBin: Makes a readBin struct for a tree
'    fun-2 freeReadBin: Frees a readBin structure
'    fun-3 pushReadBinStack: Push a readBin struct not a stack
'    fun-4 popReadBinStack: Pop the last readBin struct off a stack
'    fun-5 insBinIntoTree: Insert readBin node into a readBin tree
'    fun-6 getBinDepth: Update depth of one node in tree using children
'    fun-7 getBinTreeBal: Get balance of single node in tree
'    fun-8 rebalBinTree: Rebalance a bin tree with a readBinStack list
'    fun-9 binRightRightRebal: do right left rebalance on tree
'    fun-10 binRightRightRebal: do right rigth rebalance on tree
'    fun-11 binLeftRightRebal: do left right rebalance on tree
'    Fun-12 binTreeleftLeftRebal: do left left rebalance on tree
'    fun-13 binTreeAdjPtrs: Set parent nodes pointers
'    fun-14 freeBinTree: free a tree of readBin nodes
'    fun-15 getTreeDepth: gets max depth of readBin tree (for testing)
'    fun-16 cnvtBinTreeToList: converts readBin tree to linked list
'        - rightChild pointer is left open (set to 0)
'    fun-17 getNumBins: Get number of bins in a tree of readBin nodes
'    fun-18 rmBinFromList: Remove bin from a list of bins
'    fun-19 binDeleteFiles: Delete all files in a readBin
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
) /*Makes a readBin structure*/
{ /*makeReadBin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: makeReadBin
    #    fun-1 sec-1: variable declerations
    #    fun-1 Sec-2: Check for memory erros & initalize variables
    #    fun-1 sec-3: Copy reference id
    #    fun-1 sec-4: Copy fastq file path
    #    fun-1 sec-5: Copy stats file path
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-1 Sec-1: variable declerations & check if memory error
    *******************************************************************/

    char
        *tmpCStr = 0;

    struct readBin
        *retBin = malloc(sizeof(struct readBin));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check for memory erros & initalize variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(retBin == 0)
        return 0;

    /*Variables not set by function*/
    retBin->balUChar = '\0';
    retBin->numReadsULng = 1;  /*Their is only one read in this bin*/
    retBin->leftChild = 0;
    retBin->rightChild = 0;

    /*******************************************************************
    # Fun-1 Sec-3: Copy reference id
    *******************************************************************/

    if(refIdCStr != 0)
    { /*If have a reference id*/
        tmpCStr = retBin->refIdCStr;

        while(*refIdCStr > 16)
        { /*While not at end of reference name*/
            *tmpCStr = *refIdCStr;
            ++tmpCStr;
            ++refIdCStr;
        } /*While not at end of reference name*/

        *tmpCStr = '\0';
    } /*If have a reference id*/

    else
        retBin->refIdCStr[0] = '\0';

    /*******************************************************************
    # Fun-1 Sec-4: Copy fastq file path
    *******************************************************************/

    if(fqPathCStr != 0)
    { /*If have a reference id*/
        tmpCStr = retBin->fqPathCStr;

        while(*fqPathCStr > 16)
        { /*While not at end of reference name*/
            *tmpCStr = *fqPathCStr;
            ++tmpCStr;
            ++fqPathCStr;
        } /*While not at end of reference name*/

        *tmpCStr = '\0';
    } /*If have a reference id*/

    else
        retBin->fqPathCStr[0] = '\0';

    /*******************************************************************
    # Fun-1 Sec-5: Copy stat file path
    *******************************************************************/

    if(statPathCStr != 0)
    { /*If have a reference id*/
        tmpCStr = retBin->statPathCStr;

        while(*statPathCStr > 16)
        { /*While not at end of reference name*/
            *tmpCStr = *statPathCStr;
            ++tmpCStr;
            ++statPathCStr;
        } /*While not at end of reference name*/

        *tmpCStr = '\0';
    } /*If have a reference id*/

    else
        retBin->statPathCStr[0] = '\0';

    return retBin;
} /*makeReadBin*/

/*######################################################################
# Output:
#    Frees: the input readBin structer & all its child nodes
#    Sets: binToFree to null
######################################################################*/
void freeReadBin(
    struct readBin **binToFree
) /*Frees a readBin structure & its children readBin structure*/
{ /*freeReadBin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 Sub-1 TOC: freeReadBin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Free the binRead structure and set to null*/
    free(*binToFree);
    binToFree = 0;

    return;
} /*freeReadBin*/

/*######################################################################
# Output:
#    Modifes readStack to piont to last element
######################################################################*/
void pushReadBinStack(
    struct readBinStack **readStackAry, /*Array of read info nodes*/
    struct readBin *readBin        /*readBin to assing to next node*/
) /*pushes a readBinStack structer onto a readBinStack stack*/
{ /*makeReadBinStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: makes a readBinStack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Move to next node*/
    *readStackAry = (*readStackAry) + 1; /*+ sizeof(readBinStack);*/
    (*readStackAry)->readNode = readBin;

    return;
} /*makeReadBinStack*/

/*######################################################################
# Output: Modifies: readBinStack to point to next readBin node in stack
######################################################################*/
void popReadBinStack(
    struct readBinStack **readStackAry /*readBin Array (stack) to pop*/
) /*frees readBinStack & sets readBinStack to next node in stack*/
{ /*popReadBinStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 TOC: frees readBinStack & returns next node in stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *readStackAry = (*readStackAry) - 1;/* - sizeof(readBinStack);*/
    return;
} /*popReadBinStack*/

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
) /*Finds or creates node with input read name in a read info tree*/
{ /*findAddNodeToBinTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 TOC:
    #    fun-5 sec-1: variable declerations
    #    fun-5 sec-2: Set up stack & check if this is the first node
    #    fun-5 sec-3: Search tree, if match return the match (exit)
    #    fun-5 sec-4: Make new node, add to tree, & rebalance tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int
        matchInt = 0;                 /*Holds if read names were match*/

    struct readBin
        *refNode = *binTree,    /*Starting node of tree*/
        *queryNode = 0;          /*new node to add if read not in tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Set up stack & check if this is the first node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readStackAry->readNode = 0; /*Ensure their is a zero at the start*/

    if(*binTree == 0)
    { /*If there is no tree*/
        /*Make readBin node (converts id to number)*/

        queryNode =
            makeReadBin(
               refIdCStr,
               fqPathCStr,
               statPathCStr
        );

        if(queryNode == 0)
            return 0;      /*Memory allocation error*/

        *binTree = queryNode;
        return *binTree;
    } /*If there is no tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-3: Search tree, if match return the match (exit)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(refNode != 0)
    { /*Loop till at a leaf node or found node*/

        /*See if found a match, or if need to move right or left*/
        matchInt = strcmp(refIdCStr, refNode->refIdCStr);

        if(matchInt == 0)
        { /*If found a match*/
            ++refNode->numReadsULng; /*Incurment read count for bin*/
            return refNode;          /*Return match*/
        } /*If found a match*/

        pushReadBinStack(
            &readStackAry, /*Holds path traveled*/
            refNode
        );

        if(matchInt < 0)         /*If new id < older id, move left*/
            refNode = refNode->leftChild;
        else                       /*else, new id > old id, move right*/
            refNode = refNode->rightChild;
    } /*Loop till at a leaf node or found node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-4: Make new node, add to tree, & rebalance tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    queryNode =
            makeReadBin(
               refIdCStr,
               fqPathCStr,
               statPathCStr
    ); /*Make readBin node for tree*/

    if(queryNode == 0) /*If a memory allocation error happened*/
        return 0;

    if(matchInt < 0)
        readStackAry->readNode->leftChild = queryNode;
    else
        readStackAry->readNode->rightChild = queryNode;

    /*Reblance tree account for new node*/
    rebalBinTree(readStackAry, binTree);

    return queryNode;
} /*findAddNodeToBinTree*/

/*######################################################################
# Output:
#    Modifies: balUChar in readBin to be the depth of its deepest
#              child + 1
######################################################################*/
void getBinDepth(
    struct readBin *readBin      /*node in tree of read id's to update*/
) /*Updates the depth of a single node by using its child nodes*/
{ /*getBinDepth*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-6 Sec-1 Sub-1: Update the depth on the new node
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readBin->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readBin->leftChild == 0)
            readBin->balUChar = readBin->rightChild->balUChar + 1;

        else if(
            readBin->rightChild->balUChar >
            readBin->leftChild->balUChar
        ) /*right child is deepest*/
            readBin->balUChar = readBin->rightChild->balUChar + 1;

        else
            readBin->balUChar = readBin->leftChild->balUChar + 1;
    } /*if there is a right child, check if there is a left child*/

    else if(readBin->leftChild != 0)
        readBin->balUChar = readBin->leftChild->balUChar + 1;

    else
        readBin->balUChar = 0; /*No children, is a leafnode*/

    return;
} /*getBinDepth*/

/*######################################################################
# Output:
#    Returns:
#        - negative number: for left imbalanced
#        - positive number: for right imbalenced
#        - zero: for balanced
######################################################################*/
long getBinTreeBal(
    struct readBin *readBin      /*node to get balance from*/
) /*Gets the balance of a single node in a tree*/
{ /*getBinTreeBal*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-7 Sec-1 Sub-1: Get balance of node
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readBin->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readBin->leftChild == 0)
            return readBin->rightChild->balUChar + 1;

        return
            readBin->rightChild->balUChar -
            readBin->leftChild->balUChar +
            1;
    } /*if there is a right child, check if there is a left child*/

    else if(readBin->leftChild != 0)
        return -1 * (readBin->leftChild->balUChar + 1);

    return 0;
} /*getBinTreeBal*/

/*######################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
######################################################################*/
void rebalBinTree(
    struct readBinStack *readStackAry, /*stack of nodes to check*/
    struct readBin **binTree    /*Root of tree. Changes if root swaped*/
) /*Uses readBin nodes in readStack to reblance the tree*/
{ /*reblanceTree*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-8 Sec-1 Sub-1: rebalBinTree
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long balanceLng = 0;  /*Balance of tree*/
   struct readBin *swapNode = 0;

   while(readStackAry->readNode != 0)
   { /*while not at the root of the readBin tree, adjust balance*/
      balanceLng = getBinTreeBal(readStackAry->readNode); /*Get balance*/

      if(balanceLng > 1)
      { /*If have more nodes on the right side (right rebalance)*/
          if(getBinTreeBal(readStackAry->readNode->rightChild) < 0)
             swapNode = binRightLeftRebal(readStackAry->readNode);
          else
             swapNode = binRightRightRebal(readStackAry->readNode);

          binTreeAdjPtrs(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              binTree
          ); /*Adjust pointers after swap*/

          popReadBinStack(&readStackAry); /*Move to next node*/

          while(readStackAry->readNode != 0)
          { /*While have nodes to update depths for*/
              if(readStackAry->readNode->leftChild == 0)
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->rightChild->balUChar + 1;

              else if(readStackAry->readNode->rightChild == 0)
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->leftChild->balUChar + 1;

              else if(
                  readStackAry->readNode->leftChild->balUChar <
                  readStackAry->readNode->rightChild->balUChar
              )
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->rightChild->balUChar + 1;

              else
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->leftChild->balUChar + 1;

              popReadBinStack(&readStackAry); /*Move to next node*/
          } /*While have nodes to update depths for*/

          return;
      } /*If have more nodes on the right side (right rebalance)*/

      else if(balanceLng < -1)
      { /*else if have more nodes on the left side (left rebalance)*/
          if(getBinTreeBal(readStackAry->readNode->leftChild) > 0)
             swapNode = binLeftRightRebal(readStackAry->readNode);

          else
             swapNode = binTreeleftLeftRebal(readStackAry->readNode);

          binTreeAdjPtrs(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              binTree
          ); /*Adjust pointers after swap*/

          popReadBinStack(&readStackAry); /*Move to next node*/

          while(readStackAry->readNode != 0)
          { /*While have nodes to update depths for*/
              if(readStackAry->readNode->rightChild == 0)
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->leftChild->balUChar + 1;

              else if(readStackAry->readNode->leftChild == 0)
                  readStackAry->readNode->balUChar =
                      readStackAry->readNode->rightChild->balUChar + 1;

              else if(
                  readStackAry->readNode->leftChild->balUChar >
                  readStackAry->readNode->rightChild->balUChar
              )
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->leftChild->balUChar + 1;

              else
                  readStackAry->readNode->balUChar = 
                      readStackAry->readNode->rightChild->balUChar + 1;

              popReadBinStack(&readStackAry); /*Move to next node*/
          } /*While have nodes to update depths for*/

          return;
      } /*else if have more nodes on the left side (left rebalance)*/

      else if(balanceLng < 0) /*Left chlid is deeper*/
         readStackAry->readNode->balUChar =
            readStackAry->readNode->leftChild->balUChar + 1;

      else if(balanceLng > 0) /*right child is deeper*/
         readStackAry->readNode->balUChar =
            readStackAry->readNode->rightChild->balUChar + 1;

      else
         readStackAry->readNode->balUChar = 0; /*Is a leaf node*/

      popReadBinStack(&readStackAry); /*Move to next node*/
   } /*while not at the root of the readBin tree, adjust balance*/

   return;
} /*reblanceTree*/ /*Does a single rebalancing run on a readBin tree*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a right left reblance
######################################################################*/
struct readBin * binRightLeftRebal(
    struct readBin *parNode
) /*Does a right left rebalance on an avl readBin tree*/
{ /*binRightLeftRebal*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 TOC: do a right left rebalance on an avl tree
    #    fun-9 sec-1: variable declerations
    #    fun-9 sec-2: Do swap
    #    fun-9 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readBin *childNode = parNode->rightChild,
                     *swapNode = childNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->leftChild = swapNode->rightChild;
     parNode->rightChild = swapNode->leftChild; 
     swapNode->leftChild = parNode;
     swapNode->rightChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*
       Update depths (can do from swap node depth, since swap node was
       deepest node)
     */
     if(parNode->leftChild != 0)
        parNode->balUChar = parNode->leftChild->balUChar + 1;
     else
        parNode->balUChar = swapNode->balUChar;

     childNode->balUChar = swapNode->balUChar;
     (swapNode->balUChar)++;

     return swapNode;
} /*binRightLeftRebal*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a right right reblance
######################################################################*/
struct readBin * binRightRightRebal(
    struct readBin *parNode
) /*Does a right left rebalance on an avl readBin tree*/
{ /*binRightRightRebal*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 TOC: do a right rigth rebalance on an avl tree
    #    fun-10 sec-1: variable declerations
    #    fun-10 sec-2: Do swap
    #    fun-10 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readBin *childNode = parNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->rightChild = childNode->leftChild;
     childNode->leftChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     if(parNode->rightChild != 0)
     { /*if the paranet had a chlld, child sets balance*/
         parNode->balUChar = parNode->rightChild->balUChar + 1;
         childNode->balUChar = parNode->balUChar + 1;
     } /*if the paranet had a chlld, child sets balance*/
     else
         parNode->balUChar = childNode->leftChild->balUChar + 1;
         /*Childs balance is unchanged*/

     return childNode;
} /*binRightRightRebal*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a left right reblance
######################################################################*/
struct readBin * binLeftRightRebal(
    struct readBin *parNode
) /*Does a right left rebalance on an avl readBin tree*/
{ /*binRightLeftRebal*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 TOC: do a left right rebalance on an avl tree
    #    fun-11 sec-1: variable declerations
    #    fun-11 sec-2: Do swap
    #    fun-11 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readBin *childNode = parNode->leftChild,
                     *swapNode = childNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->rightChild = swapNode->leftChild;
     parNode->leftChild = swapNode->rightChild; 
     swapNode->rightChild = parNode;
     swapNode->leftChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*
       Update depths (can do from swap node depth, since swap node was
       deepest node)
     */
     if(parNode->rightChild != 0)
        parNode->balUChar = parNode->rightChild->balUChar + 1;
     else
        parNode->balUChar = swapNode->balUChar;

     childNode->balUChar = swapNode->balUChar;
     (swapNode->balUChar)++;

     return swapNode;
} /*binLeftRightRebal*/

/*######################################################################
# Output:
#    Returns: The new parent node (readBin struct pointer)
#    Modifies: readBin tree by doing a left left reblance
######################################################################*/
struct readBin * binTreeleftLeftRebal(
    struct readBin *parNode
) /*Does a left left rebalance on an avl readBin tree*/
{ /*binTreeleftLeftRebal*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-12 TOC: do a left left rebalance on an avl tree
    #    fun-12 sec-1: variable declerations
    #    fun-12 sec-2: Do swap
    #    fun-12 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-12 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readBin *childNode = parNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-12 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->leftChild = childNode->rightChild;
     childNode->rightChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-12 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


     if(parNode->rightChild != 0)
     { /*if the paranet had a chlld, child sets balance*/
         parNode->balUChar = parNode->rightChild->balUChar + 1;
         childNode->balUChar = parNode->balUChar + 1;
     } /*if the paranet had a chlld, child sets balance*/
     else
         parNode->balUChar = childNode->leftChild->balUChar + 1;
         /*Childs balance is unchanged*/
     
     return childNode;
} /*binTreeleftLeftRebal*/

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
) /*Addjusts parent poniter after rebalancing a readBin tree*/
{ /*binTreeAdjPtrs*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-1 TOC: binTreeAdjPtrs
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readOn == *binTree) /*Handle swapping the root node*/
        *binTree = newPar;
    else
    { /*Else need to adjust the parent nodes pointer*/
        popReadBinStack(readStackAry);

        if((*readStackAry)->readNode->rightChild == readOn)
        { /*If was a right shift or right left shift*/
            (*readStackAry)->readNode->rightChild = newPar;

            if((*readStackAry)->readNode->leftChild == 0)
               (*readStackAry)->readNode->balUChar = newPar->balUChar+1;
            else if(
                newPar->balUChar >
                (*readStackAry)->readNode->leftChild->balUChar
            )/*Else if the right hild newPar is greater*/
               (*readStackAry)->readNode->balUChar = newPar->balUChar+1;
            else
                (*readStackAry)->readNode->balUChar = 
                    (*readStackAry)->readNode->leftChild->balUChar + 1;
        } /*If was a right shift or right left shift*/

        else
        { /*Else was a left of left right swap*/
            (*readStackAry)->readNode->leftChild = newPar;

            if((*readStackAry)->readNode->rightChild == 0)
               (*readStackAry)->readNode->balUChar = newPar->balUChar+1;

            else if(
                newPar->balUChar >
                (*readStackAry)->readNode->rightChild->balUChar
            )/*Else if the left hild newPar is greater*/
               (*readStackAry)->readNode->balUChar = newPar->balUChar+1;

            else
                (*readStackAry)->readNode->balUChar = 
                  (*readStackAry)->readNode->rightChild->balUChar + 1;
        } /*Else was a left of left right swap*/

    } /*Else need to adjust the parent nodes pointer*/

    return;
} /*binTreeAdjPtrs*/

/*######################################################################
# Output:
#    Frees: every node in binTree & sets binTreeRoot to 0
#    Returns: 4: If nothing to free
#             1: if suceeded
#             0: If malloc failed to allocate stack memory
######################################################################*/
char freeBinTree(
    struct readBin **binTree  /*Root of readBin tree to free*/
) /*Frees a tree of readBinNodes*/
{ /*freeBinTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # TOC Fun-14 Sec-1 Sub-1: Free the tree from memory
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*binTree == 0)
        return 4;

    freeBinTree(&(*binTree)->leftChild);
    freeBinTree(&(*binTree)->rightChild);

    free(*binTree);
    binTree = 0;

    return 1;

    /*Doing this in recursive will slow down, but I am not worried
      about this code going fast, since it is used in findCoInft, which
      uses minimap2 & racon, which will always be slower*/
} /*freeBinTree*/

/*######################################################################
# Output:
#    int: with max depth of the input readBin tree
######################################################################*/
int getBinTreeDepth(
    struct readBin *binTree
) /*Gets depth of a readBin tree (just here for testing)*/
{ /*getTreeDepth*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ TOC Fun-15 Sec-1 Sub-1: Test function to get depth of a tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int
        leftDepthInt = 0,
        rightDepthInt = 0;

    if(binTree->leftChild != 0)
        leftDepthInt = getBinTreeDepth(binTree->leftChild);
    if(binTree->rightChild != 0)
        rightDepthInt = getBinTreeDepth(binTree->rightChild);

    if(leftDepthInt > rightDepthInt)
        return leftDepthInt + 1;
    return rightDepthInt + 1;
} /*getTreeDepth*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: The last readBin in the list (for recurisive call)
|    Note: Recursive nature will slow down, but speed is not an issue
|          here
\---------------------------------------------------------------------*/
struct readBin * cnvtBinTreeToList(
    struct readBin **binTree /*Tree to convert to a linked list*/
) /*Converts a readBin tree to a linked list*/
{ /*cnvtBinTreeToList*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    \ Fun-16 TOC: Sec-1 Sub-1: cnvtBinTreeToList
    /
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readBin
        *lastLeftBin = 0,  /*Will hold last left bin in the tree*/ 
        *lastRightBin = 0; /*will hold the last right bin in the tree*/

    if(*binTree == 0)
        return 0;

    lastLeftBin = cnvtBinTreeToList(&(*binTree)->leftChild);
    lastRightBin = cnvtBinTreeToList(&(*binTree)->rightChild);

    if(lastLeftBin == 0)
        lastLeftBin = *binTree;

    if(lastRightBin != 0)
    { /*If had a right child*/
        /*Move currnet bin to end of right children (right > left)*/
        lastRightBin->leftChild = *binTree;

        /*Move next right child to head of list*/
        *binTree = (*binTree)->rightChild;

        /*Set the bin currnetly working on rightChild pointer to 0*/
        lastRightBin->leftChild->rightChild = 0;
    } /*If had a right child*/

    else
        (*binTree)->rightChild = 0; /*Set right bin to 0*/

    return lastLeftBin;
} /*cnvtBinTreeToList*/

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
) /*Get number of bins an a readBin tree*/
{ /*getNumBins*/
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    \ Fun-17 TOC: Sec-1 Sub-1: getNumBins
    /
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(binTree == 0)
        return;

    ++(*numBinsULng);
    getNumBins(numBinsULng, binTree->leftChild);
    getNumBins(numBinsULng, binTree->rightChild);

    return;
} /*getNumBins*/

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
) /*Removes a readBin from a leftChild linked list of read bins*/
{ /*rmBinFromList*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-18 TOC: Sec-1 Sub-1: rmBinFromList
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readBin
        *tmpBin = 0;
    /*some of these files will not always exists, but I want to make
      sure they are removed when they are. So, I am just going to 
      ignore the error. All files have been created by this program
    */

    remove((*binToRm)->fqPathCStr);
    remove((*binToRm)->statPathCStr);
    remove((*binToRm)->bestReadCStr);
    remove((*binToRm)->topReadsCStr);
    remove((*binToRm)->consensusCStr);

    if((*binToRm)->rightChild != 0)
    { /*if have clusters to deal with, set up as next bin head*/
        tmpBin = (*binToRm)->rightChild; /*So can adjust the list*/
        tmpBin->leftChild = (*binToRm)->leftChild;
            /*A cluster should never have a left child (only right)*/
    } /*if have clusters to deal with, set up as next bin head*/

    else
        tmpBin = (*binToRm)->leftChild;

    free(*binToRm);
    *binToRm = tmpBin;

    return;
} /*rmBinFromList*/

/*---------------------------------------------------------------------\
| Output:
|    Deletes: All files in the bin
\---------------------------------------------------------------------*/
void binDeleteFiles(
    struct readBin *binToWipe /*Bin to delete files from*/
) /*Deletes all files in a readBin*/
{ /*binDeleteFiles*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-19 TOC: Sec-1 Sub-1: binDeleteFiles
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    FILE
        *checkFILE = 0; /*Check if file exists*/

    /*Check if the fastq file needs to be deleted*/
    checkFILE = fopen(binToWipe->fqPathCStr, "r");

    if(checkFILE != 0)
    { /* if need to remove the fastq file*/
        fclose(checkFILE);
        remove(binToWipe->fqPathCStr);
    } /* if need to remove the fastq file*/

    /*Check if the best read file needs to be deleted*/
    checkFILE = fopen(binToWipe->bestReadCStr, "r");

    if(checkFILE != 0)
    { /* if need to remove the fastq file*/
        fclose(checkFILE);
        remove(binToWipe->bestReadCStr);
    } /* if need to remove the fastq file*/

    /*Check if the top reads file needs to be delete*/
    checkFILE = fopen(binToWipe->topReadsCStr, "r");

    if(checkFILE != 0)
    { /* if need to remove the fastq file*/
        fclose(checkFILE);
        remove(binToWipe->topReadsCStr);
    } /* if need to remove the fastq file*/

    /*Check if the consensus file needs to be delete*/
    checkFILE = fopen(binToWipe->consensusCStr, "r");

    if(checkFILE != 0)
    { /* if need to remove the fastq file*/
        fclose(checkFILE);
        remove(binToWipe->consensusCStr);
    } /* if need to remove the fastq file*/

    /*Check if the stats file needs to be deleted*/
    checkFILE = fopen(binToWipe->statPathCStr, "r");

    if(checkFILE != 0)
    { /* if need to remove the fastq file*/
        fclose(checkFILE);
        remove(binToWipe->statPathCStr);
    } /* if need to remove the fastq file*/

    binToWipe->numReadsULng = 0;
    binToWipe->balUChar = 0;

    return;
} /*binDeleteFiles*/

