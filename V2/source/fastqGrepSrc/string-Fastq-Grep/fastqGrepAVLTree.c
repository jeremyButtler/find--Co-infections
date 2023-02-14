/*##############################################################################
# Name: fastqGrepAVLTree.c
# Use: Has functions to build a self blancing tree with read info nodes
##############################################################################*/

#include "fastqGrepAVLTree.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clustGraphReadTree TOC:
#    fun-1: findAddNodeToReadTree: Find node in tree or if missing adds node
#    fun-2: insertNodeIntoReadTree: Insert readinto node into readinfo tree
#    fun-3: searchTree: Find if read is in a tree
#    fun-4: updateDepth: Update depth on single node in tree using children
#    fun-5: getBalance: Get balance of single node in tree
#    fun-6: rebalanceTree: Rebalance a tree with a readInfoStack list
#    fun-7: readTreeRightRightRebalance: do a right left rebalance on avl tree
#    fun-8: readTreeRightRightRebalance: do a right rigth rebalance on avl tree
#    fun-9: readTreeLeftRightRebalance: do a left right rebalance on avl tree
#    Fun-10: readTreeLeftLeftRebalance: do a left left rebalance on an avl tree
#    fun-11: readTreeAdjustParPtrAfterSwap: Set parent nodes ptr swaped node
#    fun-12: freeReadTree: free a tree of readInfo nodes
#    fun-13: getTreeDepth: gets max depth of readInfo tree (for testing)
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Output:
#    Returns: pointer to found or new readInfo node (with readnameCStr)
#             0: If the read name was not in the tree
#    Modifes: readTree to have a new node with readNameCStr in it
##############################################################################*/
struct readInfo * findAddNodeToReadTree(
    char * readNameCStr,            /*c-string with name of read to find*/
    unsigned char lenNameUChar, /*length of readNameCStr*/
    struct readInfo **readTree,     /*readInfo tree to search for readNameCStr*/
    struct readNodeStack *readStackAry  /*Stack, (as array) for searching*/
) /*Finds or creates node with input read name in a read info tree*/
{ /*findAddNodeToReadTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: variable declerations
    #    fun-1 sec-2: Search tree, if match return the match (exit)
    #    fun-1 sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *queryIdCStr = 0;

    int
        strMatchInt = 0;                   /*Holds if read names were match*/

    struct readInfo
        *readInfoNode = *readTree,    /*Starting node of tree*/
        *newNode = 0;                 /*new node to add if read not in tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Search tree, if match return the match (exit)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readStackAry++; /*Get off 0 at start*/

    if(*readTree == 0)
    { /*If there is no tree*/
        *readTree =
            makeReadInfoStruct(
                readNameCStr,  /*Name of read in tree*/
                lenNameUChar,  /*Length of read id (so only have to find once*/
                '@'            /*Ingnore all starting @'s (marks fastq header)*/
        ); /*Make new node and add to tree*/

        return *readTree;
    } /*If there is no tree*/

    queryIdCStr = readNameCStr;            /*point to the query's name*/

    if(*readNameCStr == '@')          /*Make sure no @ header*/
        queryIdCStr++;                /*ref handeled when making tree node*/

    while(readInfoNode != 0)
    { /*Loop till at a leaf node or found node*/

        strMatchInt = strcmp(readInfoNode->idCStr, queryIdCStr);

        if(strMatchInt == 0)
            return readInfoNode;                             /*Return match*/

        pushReadNodeStack(
            &readStackAry,   /*Stack wht history of previously visted nodes*/
            readInfoNode     /*Node that I just visited*/
        );

        if(strMatchInt < 0)         /*If new id < older id, move left*/
            readInfoNode = readInfoNode->leftChild;
        else                       /*else, new id > old id, move right*/
            readInfoNode = readInfoNode->rightChild;
    } /*Loop till at a leaf node or found node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newNode =
            makeReadInfoStruct(
                readNameCStr,  /*Name of read in tree*/
                lenNameUChar,  /*Length of read id (so only have to find once*/
                '@'            /*Ingnore all starting @'s (marks fastq header)*/
    ); /*Make new node and add to tree*/

    if(newNode == 0)
    { /*If malloc did not allocate memory for the new readInfo node*/
        fprintf(stderr,"fun-1: findAddNodeToReadTree: clusterReadTree.c:117\n");
        return 0;
    } /*If malloc did not allocate memory for the new readInfo node*/

    if(strMatchInt < 0)
        readStackAry->readNode->leftChild = newNode;
    else
        readStackAry->readNode->rightChild = newNode;

    rebalanceTree(readStackAry, readTree);/*Reblance tree account for new node*/
    return newNode;                      /*Return the new node*/
} /*findAddNodeToReadTree*/

/*##############################################################################
# Output:
#    Modifies: readInfo tree by inserting node, if not duplicate
##############################################################################*/
void insertNodeIntoReadTree(
    struct readInfo *readToInsert,      /*node with read to insert in tree*/
    struct readInfo **readTree,         /*readInfo tree to insert read into*/
    struct readNodeStack *readStackAry  /*Stack, (as array) for searching*/
) /*Inserts read node in tree, if it is not a duplicate*/
{ /*insertNodeIntoReadTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC:
    #    fun-4 sec-1: variable declerations
    #    fun-4 sec-2: Search tree, if match return the match (exit)
    #    fun-4 sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *refIdCStr = 0,                      /*Will hold Id of node on*/
        *queryIdCStr = readToInsert->idCStr; /*Id of read to insert*/

    int strMatchInt = 0;                     /*Holds if read names were match*/

    struct readInfo
        *readInfoNode = *readTree;           /*Starting (root) node of tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Search tree, if match return the match (exit)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readStackAry++; /*Get off 0 at start*/

    if(*readTree == 0)
    { /*If there is no tree*/
        *readTree = readToInsert;
        return;
    } /*If there is no tree*/

    /*Do not need to check for @ header, because it was removed when I made the
      query node*/

    while(readInfoNode != 0)
    { /*Loop till at a leaf node or found node*/
        refIdCStr = readInfoNode->idCStr;      /*Set to reference name*/
        strMatchInt = strcmp(refIdCStr, queryIdCStr);

        if(strMatchInt == 0)
            return;                             /*Return match*/

        pushReadNodeStack(
            &readStackAry,
            readInfoNode   /*Node visited*/
        ); /*Add node on to the history of visited nodes*/

        if(strMatchInt < 0)         /*If new id < older id, move left*/
            readInfoNode = readInfoNode->leftChild;
        else                       /*else, new id > old id, move right*/
            readInfoNode = readInfoNode->rightChild;
    } /*Loop till at a leaf node or found node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(strMatchInt < 0)
        readStackAry->readNode->leftChild = readToInsert;
    else
        readStackAry->readNode->rightChild = readToInsert;

    rebalanceTree(
        readStackAry,  /*Stack holding search history of tree*/
        readTree       /*Tree to rebalance*/
    );/*Reblance tree account for new node (in fastqGrepAVLTree.c)*/

    return;
} /*insertNodeIntoReadTree*/

/*##############################################################################
# Output:
#    Returns: 
#      Pointer: to node with read
#      0: If the read name was not in the tree
##############################################################################*/
struct readInfo * searchTree(
    char * readNameCStr,        /*c-string with name of read to find*/
    struct readInfo *readTree   /*readInfo tree to search for readNameCStr*/
) /*Find a particler read name in the avl tree*/
{ /*findAddNodeToReadTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: variable declerations
    #    fun-3 sec-2: Search tree, if match return the match (exit)
    #    fun-3 sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *refIdCStr = 0,
         *queryIdCStr = 0;

    int strMatchInt = 0;                   /*Holds if read names were match*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Search tree, if match return the match (exit)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(readTree != 0)
    { /*Loop till at a leaf node or found node*/

        refIdCStr = readTree->idCStr;      /*Set to reference name*/
        queryIdCStr = readNameCStr;            /*point to the query's name*/

        if(*readTree->idCStr == '@')
        { /*If the reference has an @ marking start of read name*/
            if(*readNameCStr != '@')
                refIdCStr++;                   /*Move of the header indicator (@)*/
        } /*If the reference has an @ marking start of read name*/

        else if(*readNameCStr == '@')          /*reference is missing @*/
            queryIdCStr++;                     /*Make sure @ not in query*/

        strMatchInt =
          strcmp(
            refIdCStr,                             /*Read at this node*/
            queryIdCStr                            /*Read to find*/
        ); /*See if query is at this node*/

        if(strMatchInt == 0)
            return readTree;                        /*Return match*/

        if(strMatchInt < 0)                         /*query < ref, move left*/
            readTree = readTree->leftChild;
        else                                        /*query > ref, move right*/
            readTree = readTree->rightChild;
    } /*Loop till at a leaf node or found node*/

    return 0;                                       /*Read is not in the tree*/
} /*findAddNodeToReadTree*/

/*##############################################################################
# Output:
#    Modifies: balanceChar in readNode to be the depth of its deepest child + 1
##############################################################################*/
void updateDepth(
    struct readInfo *readNode      /*node in tree of read id's to update*/
) /*Updates the depth of a single node by using its child nodes*/
{ /*updateDepth*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-4 Sec-1 Sub-1: Update the depth on the new node
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readNode == 0) /*FIND AND DELETE ERROR CAUSING THIS*/
        return;
    if(readNode->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readNode->leftChild == 0)
            readNode->balanceChar = readNode->rightChild->balanceChar + 1;
        else if(readNode->rightChild->balanceChar > readNode->leftChild->balanceChar)
            readNode->balanceChar = readNode->rightChild->balanceChar + 1;
        else
            readNode->balanceChar = readNode->leftChild->balanceChar + 1;
    } /*if there is a right child, check if there is a left child*/

    else if(readNode->leftChild != 0)
        readNode->balanceChar = readNode->leftChild->balanceChar + 1; /*only left child*/

    else
        readNode->balanceChar = 0; /*No children, is a leafnode*/

    return;
} /*updateDepth*/

/*##############################################################################
# Output:
#    Returns: long with balance, - means left imbalanced, + right imbalenced
##############################################################################*/
long getBalance(
    struct readInfo *readNode      /*node to get balance from*/
) /*Gets the balance of a single node in a tree*/
{ /*getBalance*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-5 Sec-1 Sub-1: Get balance of node
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readNode->rightChild != 0)
    { /*if there is a right child, check if there is a left child*/
        if(readNode->leftChild == 0)
            return readNode->rightChild->balanceChar;

        return readNode->rightChild->balanceChar -
               readNode->leftChild->balanceChar;
    } /*if there is a right child, check if there is a left child*/

    else if(readNode->leftChild != 0)
        return -1 * readNode->leftChild->balanceChar; /*only left child*/

    return 0;
} /*getBalance*/

/*##############################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
##############################################################################*/
void rebalanceTree(
    struct readNodeStack *readStackAry, /*stack of nodes to check*/
    struct readInfo **readTree        /*Root of tree. Changes if root swaped*/
) /*Uses readInfo nodes in readStack to reblance the tree*/
{ /*reblanceTree*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # TOC Fun-6 Sec-1 Sub-1: rebalanceTree
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   long balanceLng = 0;  /*Balance of tree*/
   struct readInfo *swapNode = 0;

   while(readStackAry->readNode != 0)
   { /*while not at the root of the readInfo tree, adjust balance*/
      balanceLng = getBalance(readStackAry->readNode); /*Get balance of the node*/

      if(balanceLng > 1)
      { /*If have more nodes on the right side (right rebalance)*/
          if(getBalance(readStackAry->readNode->rightChild) < 0)
             swapNode = readTreeRightLeftRebalance(readStackAry->readNode);
          else
              swapNode = readTreeRightRightRebalance(readStackAry->readNode);

          readTreeAdjustParPtrAfterSwap(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              readTree
          ); /*Adjust pointers after swap*/
      } /*If have more nodes on the right side (right rebalance)*/

      else if(balanceLng < -1)
      { /*else if have more nodes on the left side (left rebalance)*/
          if(getBalance(readStackAry->readNode->leftChild) > 0)
             swapNode = readTreeLeftRightRebalance(readStackAry->readNode);
          else
             swapNode = readTreeLeftLeftRebalance(readStackAry->readNode);

          readTreeAdjustParPtrAfterSwap(
              readStackAry->readNode,
              swapNode,
              &readStackAry,
              readTree
          ); /*Adjust pointers after swap*/
      } /*else if have more nodes on the left side (left rebalance)*/

      else /*else not rebalance, just update the depth*/
          updateDepth(readStackAry->readNode);

      popReadNodeStack(&readStackAry); /*Move to next node*/
   } /*while not at the root of the readInfo tree, adjust balance*/

   return;
} /*reblanceTree*/ /*Does a single rebalancing run on a readInfo tree*/


/*##############################################################################
# Output:
#    Returns: The new parent node (readInfo struct pointer) in the rebalance
#    Modifies: readInfo tree by doing a right left reblance
##############################################################################*/
struct readInfo * readTreeRightLeftRebalance(
    struct readInfo *parNode
) /*Does a right left rebalance on an avl readInfo tree*/
{ /*readTreeRightLeftRebalance*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 TOC: do a right left rebalance on an avl tree
    #    fun-7 sec-1: variable declerations
    #    fun-7 sec-2: Do swap
    #    fun-7 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->rightChild,
                     *swapNode = childNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->leftChild = swapNode->rightChild;
     parNode->rightChild = swapNode->leftChild; 
     swapNode->leftChild = parNode;
     swapNode->rightChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     updateDepth(parNode);
     updateDepth(childNode);
     updateDepth(swapNode);/*call last, due to parNode & childNode are children*/

     return swapNode;
} /*readTreeRightLeftRebalance*/

/*##############################################################################
# Output:
#    Returns: The new parent node (readInfo struct pointer) in the rebalance
#    Modifies: readInfo tree by doing a right right reblance
##############################################################################*/
struct readInfo * readTreeRightRightRebalance(
    struct readInfo *parNode
) /*Does a right left rebalance on an avl readInfo tree*/
{ /*readTreeRightRightRebalance*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 TOC: do a right rigth rebalance on an avl tree
    #    fun-8 sec-1: variable declerations
    #    fun-8 sec-2: Do swap
    #    fun-8 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->rightChild = childNode->leftChild;
     childNode->leftChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     updateDepth(parNode);
     updateDepth(childNode); /*update last, because parNode is a child node*/

     return childNode;
} /*readTreeRightRightRebalance*/

/*##############################################################################
# Output:
#    Returns: The new parent node (readInfo struct pointer) in the rebalance
#    Modifies: readInfo tree by doing a left right reblance
##############################################################################*/
struct readInfo * readTreeLeftRightRebalance(
    struct readInfo *parNode
) /*Does a right left rebalance on an avl readInfo tree*/
{ /*readTreeRightLeftRebalance*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 TOC: do a left right rebalance on an avl tree
    #    fun-9 sec-1: variable declerations
    #    fun-9 sec-2: Do swap
    #    fun-9 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->leftChild,
                     *swapNode = childNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->rightChild = swapNode->leftChild;
     parNode->leftChild = swapNode->rightChild; 
     swapNode->rightChild = parNode;
     swapNode->leftChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     updateDepth(parNode);
     updateDepth(childNode);
     updateDepth(swapNode);/*call last, due to parNode & childNode are children*/

     return swapNode;
} /*readTreeLeftRightRebalance*/

/*##############################################################################
# Output:
#    Returns: The new parent node (readInfo struct pointer) in the rebalance
#    Modifies: readInfo tree by doing a left left reblance
##############################################################################*/
struct readInfo * readTreeLeftLeftRebalance(
    struct readInfo *parNode
) /*Does a left left rebalance on an avl readInfo tree*/
{ /*readTreeLeftLeftRebalance*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 TOC: do a left left rebalance on an avl tree
    #    fun-10 sec-1: variable declerations
    #    fun-10 sec-2: Do swap
    #    fun-10 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->leftChild = childNode->rightChild;
     childNode->rightChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


     updateDepth(parNode);
     updateDepth(childNode); /*update last, because parNode is a child node*/
     
     return childNode;
} /*readTreeLeftLeftRebalance*/

/*##############################################################################
# Output:
#    Modifies: the parent node of readOn to point to newPar
#    Pops: readInfo node off of readStack
##############################################################################*/
void readTreeAdjustParPtrAfterSwap(
    struct readInfo *readOn,          /*Parent node before a rebalance*/
    struct readInfo *newPar,          /*Parent node after rebalance*/
    struct readNodeStack **readStackAry, /*Stack to get the parent node from*/
    struct readInfo **readTree        /*root node to make sure not swaped*/
) /*Addjusts parent poniter after rebalancing a readInfo tree*/
{ /*readTreeAdjustParPtrAfterSwap*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-1 TOC: Set parent nodes pointer after a swap to swaped node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readOn == *readTree) /*Handle swapping the root node*/
        *readTree = newPar;
    else
    { /*Else need to adjust the parent nodes pointer*/
        popReadNodeStack(readStackAry);

        if((*readStackAry)->readNode->rightChild == readOn)
            (*readStackAry)->readNode->rightChild = newPar;
        else
            (*readStackAry)->readNode->leftChild = newPar;
    } /*Else need to adjust the parent nodes pointer*/

    return;
} /*readTreeAdjustParPtrAfterSwap*/

/*##############################################################################
# Output:
#    Frees: every node in readTree & sets readTreeRoot to 0
#    Returns: 4: If nothing to free
#             1: if suceeded
#             0: If malloc failed to allocate stack memory
##############################################################################*/
char freeReadTree(
    struct readInfo **readTree,  /*Root of readInfo tree to free*/
    struct readNodeStack *readStackAry /*Stack to use in freeing*/
) /*Frees a tree of readInfoNodes*/
{ /*freeReadTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # TOC Fun-12 Sec-1 Sub-1: Free the tree from memory
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo *deleteNode = *readTree;
    readStackAry++; /*Get off the first node*/

    if(*readTree == 0)
        return 4;

    pushReadNodeStack(&readStackAry, *readTree);

    while(readStackAry->readNode != 0)
    { /*While there are readInfo nodes in the tree*/
        deleteNode = readStackAry->readNode;
        popReadNodeStack(&readStackAry); /*Free my stack*/ 

        if(deleteNode == 0)
            continue;
        if(deleteNode->leftChild != 0)
            pushReadNodeStack(&readStackAry, deleteNode->leftChild);
        if(deleteNode->rightChild != 0)
            pushReadNodeStack(&readStackAry, deleteNode->rightChild);

        freeReadInfoStruct(&deleteNode);
    } /*While there are readInfo nodes in the tree*/

    return 1;
} /*freeReadTree*/

/*##############################################################################
# Output:
#    int: with max depth of the input readInfo tree
##############################################################################*/
int getTreeDepth(
    struct readInfo *readTree
) /*Gets depth of a readInfo tree (just here for testing)*/
{ /*getTreeDepth*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # TOC Fun-13 Sec-1 Sub-1: Test function to get depth of a tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int
        leftDepthInt = 0,
        rightDepthInt = 0;

    if(readTree->leftChild != 0)
        leftDepthInt = getTreeDepth(readTree->leftChild);
    if(readTree->rightChild != 0)
        rightDepthInt = getTreeDepth(readTree->rightChild);

    if(leftDepthInt > rightDepthInt)
        return leftDepthInt + 1;
    return rightDepthInt + 1;
} /*getTreeDepth*/
