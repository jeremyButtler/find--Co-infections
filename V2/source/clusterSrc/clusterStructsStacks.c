/*##############################################################################
# Name: clustGraphStackStructs.c
# Use: Holds the functions for structers to make stacks and ques
##############################################################################*/

#include "clusterStructsStacks.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphStackStructs
#    fun-1: pushGraphNodeStack: makes a graphNodeStack structer
#    fun-2: popGraphStackStruct: Pop grahNode off of a graph node stack
#    fun-3: pushReadNodeStack: pushes a readNodeStack structer onto stack
#    fun-4: popReadNodeStack: frees readNodeStack & returns next node in stack
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Output:
#    Returns: new graphNode stack structer (0 if malloc failed)  
#    Modifies: tailOfGraphStack to have next node be new graphNodeStackstruct
##############################################################################*/
struct graphNodeStack * pushGraphNodeStack(
    struct graphNode *graphNodeStruct,       /*graphNode to push into stack*/
    struct graphNodeStack **graphStack /*graphNode stack push into*/
) /*Push graphNode structer onto a graphNode stack*/
{ /*makeGraphNodeStackStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1 TOC: makes a graphNodeStack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNodeStack *retStruct = malloc(sizeof(graphNodeStack));

    if(retStruct == 0)
    { /*If malloc failed to find memory*/
        fprintf(stderr, "fun-1: pushGraphNodeStack: clusterStructsStack.c:");
        fprintf(stderr, "31  malloc failed to allocate memory\n");
        return 0;
    } /*If malloc failed to find memory*/

    retStruct->graphStruct = graphNodeStruct;
    retStruct->prevNode = 0;

    if(*graphStack != 0)                           /*Check if new stack*/
        retStruct->prevNode = *graphStack;          /*Push onto end of stack*/
    else
        *graphStack = retStruct;                   /*Is a new stack*/

    return retStruct;                              /*tells malloc worked*/
} /*makeGraphNodeStackStruct*/

/*##############################################################################
# Output:
#    Free: graphStackStruct and sets pointer to 0
#    Return: next graphNodeStack structer in stack
##############################################################################*/
void popGraphStackStruct(
    struct graphNodeStack **graphStack /*graphNode stack to pop graphNode off*/
) /*Frees a graph node stack structer and returns the next structer*/
{ /*popGraphNodeStackStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: Frees a graph node stack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNodeStack *prevNode = (*graphStack)->prevNode;

    free(*graphStack);
    *graphStack = prevNode;

    return;
} /*popGraphNodeStackStruct*/

/*##############################################################################
# Output:
#    Modifes readStack to piont to last element
#    Returns: newly created node if malloc allocated memory, else 0
##############################################################################*/
struct readNodeStack * pushReadNodeStack(
    struct readNodeStack **readStack, /*Stack to push readInfo node onto*/
    struct readInfo *readNode,        /*readInfo node to push into stack*/
    char whichChildChar             /*-1: path is down left child of readNode
                                       1: path is down right child of readNode*/
) /*pushes a readNodeStack structer onto a readNodeStack stack*/
{ /*makeReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: makes a readNodeStack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readNodeStack *newNode = malloc(sizeof(readNodeStack)),
                         *oldNode = 0;
    if(newNode == 0)
    { /*If malloc failed to find memory*/
        fprintf(stderr, "fun-3: pushReadNodeStack: clusterStructsStack: 90");
        fprintf(stderr, " malloc failed to allocate memory\n");
        return 0;
    } /*If malloc failed to find memory*/

    if(*readStack != 0)       /*Need to make a new stack*/
       oldNode = *readStack;

    newNode->readNode = readNode;
    newNode->whichChildChar = whichChildChar;

    *readStack = newNode; /*Allow user to access the new node*/
    newNode->previousNode = oldNode;
    /*For some odd reason doing:
        newNode->previousNode = *readStack;
        *readStack = newNode;
      Causes the newNode->previous node to point to newNode, not the previous
        node in the stack
    */
 
    return newNode;
} /*makeReadNodeStack*/

/*##############################################################################
# Output: Modifies: readNodeStack to point to next readInfo node in stack
##############################################################################*/
void popReadNodeStack(
    struct readNodeStack **readStack /*readInfo stack to pop top readInfo off*/
) /*frees readNodeStack & sets readNodeStack to next readInfo node in stack*/
{ /*popReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 TOC: frees readNodeStack & returns next node in stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readNodeStack *previousNode = (*readStack)->previousNode;

    (*readStack)->previousNode = 0; /*Make sure free does not free*/
    free(*readStack);
    *readStack = previousNode;

    return;
} /*popReadNodeStack*/
