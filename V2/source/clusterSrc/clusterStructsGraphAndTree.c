/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphGraphAndTreeStructs
#    fun-1: pushReadIntoGraph: makes a graphNode structer with default values
#    fun-2: freeGraphNodeStruct: frees a graphNode structer
#    fun-3: makeReadInfoStruct: makes a readInfo struction with default values
#    fun-4: freeReadInfoStruct: frees a readInfo structer
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "clusterStructsGraphAndTree.h"

/*##############################################################################
# Output:
#    Modifes: The edgeList in insertNode to point to the new node
#    Modifes: seqIdNode.idCStr to point to the new node (if idCStr = 0)
#    Incurments: clustNode.numNodesULng in insert node
#    Returns: new graphNode or 0 (if malloc failed to assign memory)
##############################################################################*/
struct graphNode * pushReadIntoGraph(
    struct graphNode **insertNode, /*graphNode to push read onto*/
    struct readInfo *seqIdNode    /*ReadInfo struct with sequence id for read*/
) /*Allocates memomory and makes a graphNode structer (variables set to 0)*/
{ /*makeGrapNodeStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: push read into a graph of reads
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: Check if seqIdNode is a valid node
    #    fun-1 sec-3: set common values to 0 & check if will be the root node
    #    fun-1 sec-4: Not root node, insert node into graph
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declerations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
    
    struct graphNode *graphStruct = malloc(sizeof(graphNode));

    if(graphStruct == 0)
    { /*If mallac failed to allocate memory*/
       fprintf(stderr,"fun-1: pushReadIntoGraph clusterStructsGraphAndTree.c:");
       fprintf(stderr, "38  malloc failed to allocate memory\n");
       return 0;
    } /*If mallac failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check if seqIdNode is a valid node
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(seqIdNode != 0 && seqIdNode->nodeInGraph == 0)
    { /*If readInfo struct was given*/
        seqIdNode->nodeInGraph = graphStruct;   /*need to set readInfo ptr*/
        graphStruct->seqIdNode = seqIdNode;
    } /*If readInfo struct was given*/

    else
    { /*Else this is a duplicate node or seqIdNode not provided*/ 
        seqIdNode->nodeInGraph = 0;   /*need to set readInfo ptr*/
        graphStruct->seqIdNode = 0;
    } /*Else this is a duplicate node or seqIdNode not provided*/ 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: set common values to 0 & check if will be the root node
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    graphStruct->edgesList = 0;
    graphStruct->readCntULng = 0;
    graphStruct->edgesList = 0;

    if(*insertNode == 0)
    { /*If this is the root node of the graph*/
        *insertNode = graphStruct;
        graphStruct->clustNode = graphStruct;
        graphStruct->numNodesULng = 1;              /*Is a new cluster*/
        graphStruct->prevNode = 0;                  /*Only node in graph*/
        graphStruct->nextNode = 0;                  /*Only node in graph*/
        return graphStruct;
    } /*If this is the root node of the graph*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Not root node, insert node into graph
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    graphStruct->numNodesULng = 0;              /*Is a new cluster*/
    graphStruct->clustNode = (*insertNode)->clustNode;
    graphStruct->clustNode->numNodesULng++;

    graphStruct->prevNode = *insertNode;        /*Always adding node at head*/
    graphStruct->nextNode = (*insertNode)->edgesList;

    if((*insertNode)->edgesList != 0)
        (*insertNode)->edgesList->prevNode = graphStruct;

    (*insertNode)->edgesList = graphStruct;

    return graphStruct;
} /*makeGrapNodeStruct*/

/*##############################################################################
# Output: frees the graphNode structer and sets pointer to null
# Notes: Does not free the edgeList or neighborList nodes
#        Make sure you free these nodes before runing this struct
##############################################################################*/
void freeGraphNodeStruct(
    struct graphNode **graphNodeStruct /*Stucter to free*/
) /*Frees a graphNode structer*/
{ /*freeGraphNodeStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: free a GraphNode structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo *readNode = (*graphNodeStruct)->seqIdNode;
       /*Otherwise ocassionaly sets readInfo node to 0, like double ptr*/

    if((*graphNodeStruct)->prevNode != 0)
        (*graphNodeStruct)->prevNode->nextNode = (*graphNodeStruct)->nextNode;

    if((*graphNodeStruct)->nextNode != 0)
        (*graphNodeStruct)->nextNode->prevNode = (*graphNodeStruct)->prevNode;

    free(*graphNodeStruct);
    *graphNodeStruct = 0;

    if(readNode != 0)
    { /*If there is a readInfo node linked to this structer*/
        readNode->doneChar = 1; /*So no longer accesed*/
        readNode->nodeInGraph = 0; /*So no longer accesed*/
    } /*If there is a readInfo node linked to this structer*/

    return;
} /*freeGraphNodeStruct*/

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    char *readNameCStr,        /*c-string with read name to copy*/
    unsigned char lenNameUChar /*length of readNameCStr*/
) /*Allocates memomory and makes a readInfo structer (variables set to 0)*/
{ /*initGrapNodeStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC: Make a readInfo structer
    #     fun-3 sec-1: Variable declerations
    #     fun-3 sec-2: Copy the input string
    #     fun-3 sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpChar = readNameCStr,
         *tmp2Char = 0;

    struct readInfo *readInfoStruct = malloc(sizeof(readInfo));

    if(readInfoStruct == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(stderr,"pushReadIntoGraph (fun-2 makeReadInfoStruct:152");
        fprintf(stderr, "  malloc failed to allocate memory\n");
        return 0;
    } /*If mallac failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Copy the input string
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->idCStr = malloc(sizeof(char) * (lenNameUChar + 1));

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
    # Fun-3 Sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->balanceChar = 0;
    readInfoStruct->leftChild = 0;
    readInfoStruct->rightChild = 0;
    readInfoStruct->nodeInGraph = 0;
    readInfoStruct->doneChar = 0;

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
    # Fun-4 Sec-1 TOC: free readInfoStruct
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *graphStruct = (*readInfoStruct)->nodeInGraph;

    if((*readInfoStruct)->idCStr != 0)
        free((*readInfoStruct)->idCStr);          /*Need to free the read name*/


    free(*readInfoStruct);                          /*User handles nodeInGraph*/
    readInfoStruct = 0;

    if(graphStruct != 0)
       graphStruct->seqIdNode = 0; /*prevent leaks*/

    return;
} /*freeReadInfoStruct*/
