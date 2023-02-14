/*##############################################################################
# Name: clustGraphReadTree.c
# Use: Has functions to build a self blancing tree with read info nodes
##############################################################################*/

#include "clusterReadTree.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clustGraphReadTree TOC:
#    fun-1: findAddNodeToReadTree: Find node in tree or if missing adds node
#    fun-2: rebalanceTree: Rebalance a tree with a readInfoStack list
#    fun-3: freeReadTree: free a tree of readInfo nodes
#    fun-4: readTreeRightRightRebalance: do a right left rebalance on avl tree
#    fun-5: readTreeRightRightRebalance: do a right rigth rebalance on avl tree
#    fun-6: readTreeLeftRightRebalance: do a left right rebalance on avl tree
#    Fun-7: readTreeLeftLeftRebalance: do a left left rebalance on an avl tree
#    fun-8: readTreeAdjustParPtrAfterSwap: Set parent nodes ptr swaped node
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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

    int strMatchInt = 0;                   /*Holds if read names were match*/

    struct readInfo *readInfoNode = *readTree,    /*Starting node of tree*/
                    *newNode = 0;        /*new node to add if read not in tree*/

    struct readNodeStack *readStack = 0;      /*current postion in stack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Search tree, if match return the match (exit)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*readTree == 0)
    { /*If there is no tree*/
        *readTree = makeReadInfoStruct(readNameCStr, lenNameUChar);
        return *readTree;
    } /*If there is no tree*/

    while(readInfoNode != 0)
    { /*Loop till at a leaf node or found node*/
        strMatchInt = strcmp(readInfoNode->idCStr, readNameCStr);

        if(strMatchInt == 0)
        { /*If was a match*/
            while(readStack != 0)
                popReadNodeStack(&readStack); /*Free my stack*/ 
            return readInfoNode;                             /*Return match*/
        } /*If was a match*/

        if(strMatchInt < 0)
        { /*If the new string is < the old string*/
            if(pushReadNodeStack(&readStack, readInfoNode, -1) == 0)
            { /*If malloc failed to allocate memory*/
                while(readStack != 0)
                    popReadNodeStack(&readStack);
                fprintf(
                    stderr,
                    "fun-1: findAddNodeToReadTree: clusterReadTree.c:77\n"
                ); /*Print out the location the error is at*/
                return 0;
            } /*If malloc failed to allocate memory*/

            readInfoNode = readInfoNode->leftChild;
        } /*If the new string is < the old string*/

        else
        { /*else the new string is > the old string*/
            if(pushReadNodeStack(&readStack, readInfoNode, 1) == 0)  /*add node*/
            { /*If malloc failed to allocate memory*/
                while(readStack != 0)
                    popReadNodeStack(&readStack);
                fprintf(
                    stderr,
                    "fun-1: findAddNodeToReadTree: clusterReadTree.c:100\n"
                ); /*Print out the location the error is at*/
                return 0;
            } /*If malloc failed to allocate memory*/
               /*1 tells me that the right branch was choosen*/

            readInfoNode = readInfoNode->rightChild;       /*new sting > old*/
        } /*else the new string is > the old string*/
    } /*Loop till at a leaf node or found node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Add node to tree, rebalance tree, return the new node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newNode = makeReadInfoStruct(readNameCStr, lenNameUChar);

    if(newNode == 0)
    { /*If malloc did not allocate memory for the new readInfo node*/
        fprintf(stderr,"fun-1: findAddNodeToReadTree: clusterReadTree.c:117\n");
        return 0;
    } /*If malloc did not allocate memory for the new readInfo node*/

    if(strMatchInt < 0)
        readStack->readNode->leftChild = newNode;
    else
        readStack->readNode->rightChild = newNode;

    rebalanceTree(&readStack, readTree); /*Reblance tree, account for new node*/
    return newNode;                      /*Return the new node*/
} /*findAddNodeToReadTree*/

/*##############################################################################
# Output: Modifies: tree to be rebalanced when nodes unbalanced
##############################################################################*/
void rebalanceTree(
    struct readNodeStack **readStack, /*stack of nodes to check*/
    struct readInfo **readTree        /*Root of tree. Changes if root swaped*/
) /*Uses readInfo nodes in readStack to reblance the tree*/
{ /*reblanceTree*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-2 TOC: rebalanceTree
   #    fun-2 sec-1: declare variables
   #    fun-2 sec-2: get first readInfo node & adjust balance
   #    fun-2 sec-3: Do right reblance (more nodes on right)
   #    fun-2 sec-4: Do left reblance (more nodes on left)
   #    fun-2 sec-5: Move to the parent node and see if needs rebalance
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-2 Sec-1: declare variables
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct readInfo *readOn = 0, 
                   *swapNode = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-2 Sec-2: get first readInfo node & adjust balance
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(*readStack != 0)
   { /*while not at the root of the readInfo tree, adjust balance*/

      readOn = (*readStack)->readNode;
      readOn->balanceChar += (*readStack)->whichChildChar;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # Fun-2 Sec-3: Do right reblance (more nodes on right)
      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

      if(readOn->balanceChar > 1)
      { /*If have more nodes on the right side*/

          if(readOn->rightChild->balanceChar < 0)
             swapNode = readTreeRightLeftRebalance(readOn);
          else
              swapNode = readTreeRightRightRebalance(readOn);

          readTreeAdjustParPtrAfterSwap(readOn, swapNode, readStack, readTree);

          while(*readStack != 0){popReadNodeStack(readStack);}
          return; /*Tree is balanced*/
      } /*If have more nodes on the right side*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # Fun-2 Sec-4: Do left reblance (more nodes on left)
      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      else if(readOn->balanceChar < -1)
      { /*else if have more nodes on the left side*/

          if(readOn->leftChild->balanceChar > 0)
             swapNode = readTreeLeftRightRebalance(readOn);
          else
             swapNode = readTreeLeftLeftRebalance(readOn);

          readTreeAdjustParPtrAfterSwap(readOn, swapNode, readStack, readTree);

          while(*readStack != 0){popReadNodeStack(readStack);}
          return; /*Tree is balanced*/
      } /*else if have more nodes on the left side*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # Fun-2 Sec-5: Move to the parent node and see if needs rebalance
      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      popReadNodeStack(readStack);
   } /*while not at the root of the readInfo tree, adjust balance*/

   return;
} /*reblanceTree*/ /*Does a single rebalancing run on a readInfo tree*/

/*##############################################################################
# Output:
#    Frees: every node in readTree & sets readTreeRoot to 0
#    Returns: 4: If nothing to free
#             1: if suceeded
#             0: If malloc failed to allocate stack memory
##############################################################################*/
char freeReadTree(
    struct readInfo **readTree  /*Root of readInfo tree to free*/
) /*Frees a tree of readInfoNodes*/
{ /*freeReadTree*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: freeReadTree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo *deleteNode = *readTree;
    struct readNodeStack *readStack = 0;

    if(*readTree == 0)
        return 4;

    if(pushReadNodeStack(&readStack, *readTree, 0) == 0) /*init stack*/
    { /*If malloc did not allocate memory for the new readInfo node*/
        fprintf(stderr, "fun-2: freeReadTree: clusterReadTree.c:232\n");
        return 0;
    } /*If malloc did not allocate memory for the new readInfo node*/

    while(readStack != 0)
    { /*While there are readInfo nodes in the tree*/
        deleteNode = readStack->readNode;
        popReadNodeStack(&readStack); /*Free my stack*/ 

        if(deleteNode == 0)
            continue;

        if(deleteNode->leftChild != 0)
        { /*If have a left child to free later*/
            if(pushReadNodeStack(&readStack, deleteNode->leftChild, -1) == 0)
            { /*If malloc did not allocate memory for the new readInfo node*/
                fprintf(stderr, "fun-2: freeReadTree: clusterReadTree.c:250\n");
                while(readStack != 0) popReadNodeStack(&readStack);
                return 0;
            } /*If malloc did not allocate memory for the new readInfo node*/
        } /*If have a left child to free later*/

        if(deleteNode->rightChild != 0)
        { /*If have a right child to free later*/
            if(pushReadNodeStack(&readStack, deleteNode->rightChild, 1) == 0)
            { /*If malloc did not allocate memory for the new readInfo node*/
                fprintf(stderr, "fun-2: freeReadTree: clusterReadTree.c:257\n");
                while(readStack != 0) popReadNodeStack(&readStack);
                return 0;
            } /*If malloc did not allocate memory for the new readInfo node*/
        } /*If have a right child to free later*/

        freeReadInfoStruct(&deleteNode);
    } /*While there are readInfo nodes in the tree*/

    return 1;
} /*freeReadTree*/

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
    # Fun-4 TOC: do a right left rebalance on an avl tree
    #    fun-4 sec-1: variable declerations
    #    fun-4 sec-2: Do swap
    #    fun-4 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->rightChild,
                     *swapNode = childNode->leftChild;


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->leftChild = swapNode->rightChild;
     parNode->rightChild = swapNode->leftChild; 
     swapNode->leftChild = parNode;
     swapNode->rightChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*Adjust balances for reblance*/
     if(swapNode->balanceChar == -1)
     { /*If the swapNode had only a left child (childNode gets no left child)*/
         if(childNode->rightChild != 0)
             childNode->balanceChar = 1; /*If childNode has only right child*/
         else
             childNode->balanceChar = 0; /*If childNode has no children now*/

         if(parNode->leftChild != 0)
             parNode->balanceChar = 0;   /*Parent node has both children now*/
         else
             parNode->balanceChar = 1;   /*Parent node has only right child*/
     } /*If the swapNode had only a left child (childNode gets no left child)*/

     else if(swapNode->balanceChar == 1)
     { /*Else if swapNode had no left child to give parNode as right child*/
         if(parNode->leftChild != 0)
             parNode->balanceChar = -1;  /*If parNode has only left child now*/
         else
             parNode->balanceChar = 0;   /*If parNode has no children now*/

         if(childNode->rightChild != 0)
             childNode->balanceChar = 0; /*If child node has both children*/
         else
             childNode->balanceChar = -1; /*If child node has only left child*/
     } /*Else if swapNode had no left child to give parNode as right child*/

     else
     { /*Else is perfectly balanced (This is impossible if reguarly balanced)*/
         childNode->balanceChar = 0;
         parNode->balanceChar = 0;
     } /*Else is perfectly balanced*/

     swapNode->balanceChar = 0;

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
    # Fun-5 TOC: do a right rigth rebalance on an avl tree
    #    fun-5 sec-1: variable declerations
    #    fun-5 sec-2: Do swap
    #    fun-5 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->rightChild;


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->rightChild = childNode->leftChild;
     childNode->leftChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->balanceChar = 0; /*Child node is perfectly balanced*/

     if(childNode->leftChild == 0)
     { /*If childNode has no left child to give parNode as right child*/
         if(parNode->leftChild == 0)
             parNode->balanceChar = 0; /*parNode has no children*/
         else
             parNode->balanceChar = -1; /*parNode only has right child*/
     } /*If childNode has no left child to give parNode as right child*/

     else
     { /*Else childNode had a left child to give parNode as a right child*/
         if(parNode->leftChild == 0)
             parNode->balanceChar = 1; /*parNode only has a right child*/
         else
             parNode->balanceChar = 0; /*parNode has no children*/
     } /*Else childNode had a left child to give parNode as a right child*/

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
    # Fun-6 TOC: do a left right rebalance on an avl tree
    #    fun-6 sec-1: variable declerations
    #    fun-6 sec-2: Do swap
    #    fun-6 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->leftChild,
                     *swapNode = childNode->rightChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->rightChild = swapNode->leftChild;
     parNode->leftChild = swapNode->rightChild; 
     swapNode->rightChild = parNode;
     swapNode->leftChild = childNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*Adjust balances for reblance*/
     if(swapNode->balanceChar == -1)
     { /*If the swapNode had only a left child (parNode gets no left child)*/
         if(parNode->rightChild != 0)
             parNode->balanceChar = 1;   /*If parNode has only right child now*/
         else
             parNode->balanceChar = 0;   /*If parNode has no children now*/

         if(childNode->leftChild != 0)
             childNode->balanceChar = 0; /*ChildNode has two children*/
         else
             childNode->balanceChar = 1; /*childNode has only right child*/
     } /*If the swapNode had only a left child (parNode gets no left child)*/

     else if(swapNode->balanceChar == 1)
     { /*Else if swapNode had only a right (childNode has no right child)*/
         if(childNode->leftChild != 0)
             childNode->balanceChar = -1; /*If childNode has only left child*/
         else
             childNode->balanceChar = 0;  /*If childNode has no children now*/

         if(parNode->rightChild != 0)
             parNode->balanceChar = 0;    /*parNode is balanced*/
         else
             parNode->balanceChar = -1;    /*parNode has only a left child*/
     } /*Else if swapNode had only a right (childNode has no right child)*/

     else
     { /*Else is perfectly balanced (This is impossible if reguarly balanced)*/
         childNode->balanceChar = 0;
         parNode->balanceChar = 0;
     } /*Else is perfectly balanced*/

     swapNode->balanceChar = 0;

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
    # Fun-7 TOC: do a left left rebalance on an avl tree
    #    fun-7 sec-1: variable declerations
    #    fun-7 sec-2: Do swap
    #    fun-7 sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1: Variable decerlations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     struct readInfo *childNode = parNode->leftChild;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-2: Do swap
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     parNode->leftChild = childNode->rightChild;
     childNode->rightChild = parNode; 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-3: Adjust balance
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     childNode->balanceChar = 0; /*Child node is perfectly balanced*/

     if(childNode->rightChild == 0)
     { /*If childNode has no right child to give parNode as left child*/
         if(parNode->rightChild == 0)
             parNode->balanceChar = 0; /*parNode has no children*/
         else
             parNode->balanceChar = 1; /*parNode only has right child*/
     } /*If childNode has no rigth child to give parNode as left child*/

     else
     { /*Else childNode had a right child to give parNode as a left child*/
         if(parNode->rightChild == 0)
             parNode->balanceChar = -1; /*parNode only has a left child*/
         else
             parNode->balanceChar = 0; /*parNode has no children*/
     } /*Else childNode had a right child to give parNode as a left child*/

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
    struct readNodeStack **readStack, /*Stack to get the parent node from*/
    struct readInfo **readTree        /*root node to make sure not swaped*/
) /*Addjusts parent poniter after rebalancing a readInfo tree*/
{ /*readTreeAdjustParPtrAfterSwap*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1 TOC: Set parent nodes pointer after a swap to swaped node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(readOn == *readTree) /*Handle swapping the root node*/
        *readTree = newPar;
    else
    { /*Else need to adjust the parent nodes pointer*/
        popReadNodeStack(readStack);

        if((*readStack)->readNode->rightChild == readOn)
            (*readStack)->readNode->rightChild = newPar;
        else
            (*readStack)->readNode->leftChild = newPar;
    } /*Else need to adjust the parent nodes pointer*/

    return;
} /*readTreeAdjustParPtrAfterSwap*/
