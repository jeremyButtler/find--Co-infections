/*##############################################################################
# Name:
# Use: Holds structers needed to make a tree of read names
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphGraphAndTreeStructs
#    fun-1: makeReadInfoStruct: makes a readInfo struction with default values
#    fun-2: freeReadInfoStruct: frees a readInfo structer
#    fun-3: pushReadNodeStack: pushes a readNodeStack structer onto stack
#    fun-4: popReadNodeStack: frees readNodeStack & returns next node in stack
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fastqGrepStructs.h"

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    char *readNameCStr,         /*c-string with read name to copy*/
    unsigned int lenNameUInt, /*length of readNameCStr*/
    char ignoreChar             /*Character at start to ignore ('' is nothing)*/
) /*Allocates memomory and makes a readInfo structer (variables set to 0)*/
{ /*initGrapNodeStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: Make a readInfo structer
    #     fun-1 sec-1: Variable declerations
    #     fun-1 sec-2: Copy the input string
    #     fun-1 sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpChar = readNameCStr,
        *tmp2Char = 0;

    struct readInfo
        *readInfoStruct = malloc(sizeof(readInfo));

    if(*tmpChar == ignoreChar)
    { /*If the first char should be ignored (marks header)*/
        tmpChar++;       /*Move past first character in read id*/
        lenNameUInt--;  /*Reduce length by one to account for header char*/
    } /*If the first char should be ignored (marks header)*/

    if(readInfoStruct == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(stderr,"pushReadIntoGraph (fun-2 makeReadInfoStruct:152");
        fprintf(stderr, "  malloc failed to allocate memory\n");
        return 0;
    } /*If mallac failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Copy the input string
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->idCStr = malloc(sizeof(char) * (lenNameUInt + 1));

    if(readInfoStruct->idCStr == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(stderr,"pushReadIntoGraph (fun-2 makeReadInfoStruct:167");
        fprintf(stderr, "  malloc failed to allocate memory\n");
        return 0;
    } /*If mallac failed to allocate memory*/

    tmp2Char = readInfoStruct->idCStr;

    while(*tmpChar != '\0')
    { /*Loop to copy the read name over*/
        *tmp2Char = *tmpChar;
        tmp2Char++;
        tmpChar++;
    } /*Loop to copy the read name over*/

    *tmp2Char = '\0';
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->balanceChar = 0;
    readInfoStruct->leftChild = 0;
    readInfoStruct->rightChild = 0;

    return readInfoStruct;
} /*initGrapNodeStruct*/

/*##############################################################################
# Output: frees readInforStruct and sets pointer to 0
# Note: Does not free nodeInGraph.
##############################################################################*/
void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
) /*frees a readInfo structer*/
{ /*freeReadInfoStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: free readInfoStruct
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*readInfoStruct)->idCStr != 0)
        free((*readInfoStruct)->idCStr);          /*Need to free the read name*/

    free(*readInfoStruct);                          /*User handles nodeInGraph*/
    *readInfoStruct = 0;

    return;
} /*freeReadInfoStruct*/


/*##############################################################################
# Output:
#    Modifes readStack to piont to last element
#    Returns: newly created node if malloc allocated memory, else 0
##############################################################################*/
void pushReadNodeStack(
    struct readNodeStack **readStackAry, /*Array of read info nodes*/
    struct readInfo *readNode            /*readInfo to assing to next node*/
) /*pushes a readNodeStack structer onto a readNodeStack stack*/
{ /*makeReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: makes a readNodeStack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Move to next node*/
    *readStackAry = (*readStackAry) + 1; /*+ sizeof(readNodeStack);*/
    (*readStackAry)->readNode = readNode;

    return;
} /*makeReadNodeStack*/

/*##############################################################################
# Output: Modifies: readNodeStack to point to next readInfo node in stack
##############################################################################*/
void popReadNodeStack(
    struct readNodeStack **readStackAry /*readInfo Array (stack) to pop*/
) /*frees readNodeStack & sets readNodeStack to next readInfo node in stack*/
{ /*popReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 TOC: frees readNodeStack & returns next node in stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    *readStackAry = (*readStackAry) - 1;/* - sizeof(readNodeStack);*/
    return;
} /*popReadNodeStack*/
