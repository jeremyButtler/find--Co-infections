/*######################################################################
# Name: trimPrimersStructs
# Use:
#   o Holds the structures needed for trimPrimers
# Includes:
#   - "fqGetIdsStructs.h"
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: fqGetsIdsStructs
'   fun-1 makeReadAndPrimST:
'     o Makes a readPrim struction with default values & the input 
'       read id converted to a big number
'   fun-2 freeReadPrimST:
'     o Frees a readPrim structer
'   fun-3 freePrimCordList:
'     o Frees a linked list of primer coordinantes (primCord structer)
'   fun-4 freePrimCordST
'     o Frees a single primer coordinante (primCord strucuter)
'   fun-5 pushReadPrimStack:
'     o Pushes a readPrim structer onto a readPrimStack stack
'   fun-6 popReadPrimStack:
'     o Pushes a pops a readPrim structer of a readPrimStack stack
'   fun-7 insPrimCordST:
'     o Insert a primer coordinate structure (primCord) or linked list
'       into a list of primer coordinates
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "trimPrimersStructs.h"

/*---------------------------------------------------------------------\
| Output: Returns: a readPrim structer heap or 0 for memory errors
\---------------------------------------------------------------------*/
struct readPrim * makeReadPrimST()
/*Makes a blank readPrim structer on the heap*/
{ /*makeReadInfoStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Sec-1 Sub-1: makeReadPrimST
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readPrim *readPrimST = malloc(sizeof(struct readPrim));

    if(readPrimST == 0)
        return 0; /*memory allocation error*/

    readPrimST->idBigNum = 0;
    readPrimST->balUC = 0;
    readPrimST->leftChild = 0;
    readPrimST->rightChild = 0;
    readPrimST->primCordST = 0; /*No primer cordiant list*/

    return readPrimST;
} /*makeReadInfoStruct*/

/*---------------------------------------------------------------------\
| Output: frees a readPrim structer and sets its pointer to 0
\---------------------------------------------------------------------*/
void freeReadPrimST(
    struct readPrim **readPrimST /*struct to free*/
) /*frees a readInfo structer*/
{ /*freeReadInfoStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1 TOC: free readInfoStruct
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    freeBigNumStruct(&((*readPrimST)->idBigNum));
    freePrimCordList(&(*readPrimST)->primCordST);

    free(*readPrimST);               /*User handles nodeInGraph*/
    *readPrimST = 0;

    return;
} /*freeReadInfoStruct*/

/*---------------------------------------------------------------------\
| Output: frees a primCord structer linked list & sets root pointer to 0
\---------------------------------------------------------------------*/
void freePrimCordList(
    struct primCord **cordToFree /*Primer cordinate list to free*/
) /*Frees a primCord structer linked list & sets root pointer to 0*/
{ /*freePrimCorList*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: Sec-1 Sub-1: freePrimCordList
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct primCord *tmpCord = *cordToFree;
    struct primCord *nextCord = 0;

    while(tmpCord != 0)
    { /*While I have a linked list of primer cordinates to free*/
        nextCord = tmpCord->nextCord;
        freePrimCordST(&tmpCord);
        tmpCord = nextCord;
    } /*While I have a linked list of primer cordinates to free*/

    *cordToFree = 0;
    return;
} /*freePrimCorList*/

/*---------------------------------------------------------------------\
| Output: Frees a single primCord structure and sets its pointer to 0
\---------------------------------------------------------------------*/
void freePrimCordST(
    struct primCord **cordToFree /*Single primer cordinate to free*/
) /*Frees a single primCord structure and sets its pointer to 0*/
{ /*freePrimCordST*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: Sec-1 Sub-1: freePrimCordST
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(*cordToFree != 0)
        free(*cordToFree);

    *cordToFree = 0;
    return;
} /*freePrimCordST*/

/*---------------------------------------------------------------------\
| Output:
|   o Pushes a readPrim structer onto a readPrim stack
|   o Modifies readStackAry to piont to the most recent node
\---------------------------------------------------------------------*/
void pushReadPrimStack(
    struct readPrimStack **readStackAry, /*Stack as an array*/
    struct readPrim *readPrimST        /*next node to add to the stack*/
) /*pushes a readPrim structer onto a stack of readStack structures*/
{ /*pushReadPrimStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-5 Sec-1 TOC: makes a readNodeStack structer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Move to next node*/
    *readStackAry = (*readStackAry) + 1; /*+ sizeof(readNodeStack);*/
    (*readStackAry)->readNode = readPrimST;

    return;
} /*makeReadNodeStack*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o readStackAry to point to the next readPrim node in the stack
\---------------------------------------------------------------------*/
void popReadPrimStack(
    struct readPrimStack **readStackAry /*readInfo Array (stack) to pop*/
) /*Sets readStackAray to next readPrim stucter in stack*/
{ /*popReadPrimStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 TOC: Sec-1 Sub-1: popReadPrimStck
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *readStackAry = (*readStackAry) - 1;/* - sizeof(readNodeStack);*/
    return;
} /*popReadPrimStack*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o primCordList to start with primCordST if insertion is at root
|     o primCordList to have primCordSt inserted in
|   - Warnings:
|     o This assumes that both of your lists are sorted by starting 
|       coordinate
\---------------------------------------------------------------------*/
void insPrimCordST(
    struct primCord *primCordST,
        /*primer coordinates to insert, this can be as a linked list*/
    struct primCord **primCordList
        /*Linked list to insert into. Use 0 for first node*/
) /*Insert a primer coordinate or list of primer coordinates into a 
   primer coordinate linked list*/
{ /*insPrimCordST*/      

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-7 TOC: insPrimCordST
    '   o fun-7 sec-1: Variable declerations
    '   o fun-7 sec-2: Check if I have an easy task
    '   o fun-7 sec-3: See how much of the list can be insereted at root
    '   o fun-7 sec-4: Find were to insert nodes after root
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct primCord *tmpCord = 0;
    struct primCord *lastCord = 0;
    struct primCord *swapCord = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-2: Check if I have an easy task
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*primCordList == 0)
    { /*If this is the first node in the list*/
        *primCordList = primCordST;
        return;
    } /*If this is the first node in the list*/

    if(primCordST == 0)
        return;        /*Nothing to insert*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-3: See how much of the list can be insereted at root
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    tmpCord = primCordST;

    while(((*primCordList)->startUI > primCordST->startUI))
    { /*While I am inserting a lower coordinate*/
        lastCord = primCordST;
        primCordST = primCordST->nextCord;

        if(primCordST == 0)
            break;
    } /*While I am inserting a lower coordinate*/

    if(lastCord != 0)
    { /*if I need to swap the root node*/
        swapCord = *primCordList;      /*Save pointer to old root*/
        *primCordList = tmpCord;       /*Set node to insert as root*/
        lastCord->nextCord = swapCord; /*Move old root to new position*/
        lastCord = swapCord;          /*so know the last coordinate on*/
        tmpCord = swapCord->nextCord;  /*Move to next start coordinate*/

        if(primCordST == 0)
            return;                    /*If finshed*/
    } /*if I need to swap the root node*/

    else
    { /*Else starting comparison at the second coordinate*/
        lastCord = *primCordList;
        tmpCord = lastCord->nextCord;
    } /*Else starting comparison at the second coordinate*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-4: Find were to insert nodes after root
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(tmpCord != 0)
    { /*While I have nodes to compare*/
        if(primCordST->startUI < tmpCord->startUI)
        { /*If I need to do a swap*/
            lastCord->nextCord = primCordST;

            while(primCordST->startUI < tmpCord->startUI)
            { /*While the insert node has a lower start coordinate*/
                swapCord = primCordST;
                primCordST = primCordST->nextCord;

                if(primCordST == 0)
                    break; /*No more structures to search for*/
            } /*While the insert node has a lower start coordinate*/

            /*Add the pre-insert coordinate to the end*/
            swapCord->nextCord = tmpCord;

            if(primCordST == 0)
                return;              /*Finshed inserting the list*/
        } /*If I need to do a swap*/

        lastCord = tmpCord;
        tmpCord = tmpCord->nextCord;
    } /*While I have nodes to compare*/

    /*Transfer any remaing coordinates to the end*/
    lastCord->nextCord = primCordST;
    return;
} /*insPrimCordST*/      


/*---------------------------------------------------------------------\
| Output: Returns: a blank primCord structer or 0 for memory errors
\---------------------------------------------------------------------*/
struct primCord * makePrimCord()
/*make a blank primCord structure on the heap*/
{ /*makePrimCord*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-8 TOC: Sec-1 Sub-1: makePrimCord
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct primCord *primCordST = malloc(sizeof(struct primCord));

    if(primCordST == 0)
        return 0;     /*memory alloction error*/

    primCordST->startUI = 0;
    primCordST->endUI = 0;
    primCordST->nextCord = 0;

    return primCordST;
} /*makePrimCord*/
