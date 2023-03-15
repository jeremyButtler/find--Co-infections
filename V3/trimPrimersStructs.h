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

#ifndef TRIMPRIMERSSTRUCTS_H
#define TRIMPRIMERSSTRUCTS_H

#include "fqGetIdsStructs.h"

/*---------------------------------------------------------------------\
| Struct-1: primCord
| Use: Stores the start and end cordinates of a primer mapping on a read
\---------------------------------------------------------------------*/
typedef struct primCord
{ /*primCord*/
    uint32_t startUI; /*Starting corrdinate of primer on read*/
    uint32_t endUI;   /*endind cordinate of primer on read*/
    /*Note minimap2 alwasy makes the starting cordinate lesser*/
    struct primCord *nextCord; /*Next cordinate in primer list*/
}primCord;

/*---------------------------------------------------------------------\
| ST-2: readPrim
| Use: Stores the hex id of a read as a big number & a linked list of
|      mappings for each primer start and end
\---------------------------------------------------------------------*/
typedef struct readPrim
{ /*readPrim*/
    int8_t balUC;            /*Tells if the node is balanced*/
    struct bigNum *idBigNum; /*Holds read id as unique big number*/
    struct primCord *primCordST; /*Primer start and end corridinates*/
    struct readPrim *leftChild; 
    struct readPrim *rightChild;
}readPrim;

/*---------------------------------------------------------------------\
| Struct-3: readNodeStack
| Use: Used to make a stack of readPrim structers, which is used to
|      transvering an AVL tree.
| Note:
|   o Users makes array of nodes they pass around
|   o A well balanced AVL tree coverts 10^18 nodes at a depth of 73
\---------------------------------------------------------------------*/
typedef struct readPrimStack
{ /*readNodeStack*/
    struct readPrim *readNode;                /*Node in stack*/
}readPrimStack; /*readNodeStack*/

/*---------------------------------------------------------------------\
| Output: Returns: a readPrim structer heap or 0 for memory errors
\---------------------------------------------------------------------*/
struct readPrim * makeReadPrimST();
/*Makes a blank readPrim structer on the heap*/

/*---------------------------------------------------------------------\
| Output: frees a readPrim structer and sets its pointer to 0
\---------------------------------------------------------------------*/
void freeReadPrimST(
    struct readPrim **readPrimST /*struct to free*/
); /*frees a readInfo structer*/

/*---------------------------------------------------------------------\
| Output: frees a primCord structer linked list & sets root pointer to 0
\---------------------------------------------------------------------*/
void freePrimCordList(
    struct primCord **cordToFree /*Primer cordinate list to free*/
); /*Frees a primCord structer linked list & sets root pointer to 0*/

/*---------------------------------------------------------------------\
| Output: Frees a single primCord structure and sets its pointer to 0
\---------------------------------------------------------------------*/
void freePrimCordST(
    struct primCord **cordToFree /*Single primer cordinate to free*/
); /*Frees a single primCord structure and sets its pointer to 0*/

/*---------------------------------------------------------------------\
| Output:
|   o Pushes a readPrim structer onto a readPrim stack
|   o Modifies readStackAry to piont to the most recent node
\---------------------------------------------------------------------*/
void pushReadPrimStack(
    struct readPrimStack **readStackAry, /*Stack as an array*/
    struct readPrim *readPrimST        /*Next node to add to the stack*/
); /*pushes a readPrim structer onto a stack of readStack structures*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o readStackAry to point to the next readPrim node in the stack
\---------------------------------------------------------------------*/
void popReadPrimStack(
    struct readPrimStack **readStackAry /*readInfo Array (stack) to pop*/
); /*Sets readStackAray to next readPrim stucter in stack*/

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
); /*Insert a primer coordinate or list of primer coordinates into a 
   primer coordinate linked list*/

/*---------------------------------------------------------------------\
| Output: Returns: a blank primCord structer or 0 for memory errors
\---------------------------------------------------------------------*/
struct primCord * makePrimCord();
/*make a blank primCord structure on the heap*/

#endif
