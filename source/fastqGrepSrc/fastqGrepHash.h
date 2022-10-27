/*##############################################################################
# Name: fastqGrepHash.c
#   Use:
#     Uses an hash + AVL tree to find if reads is in user suplied list
##############################################################################*/

#ifndef FASTQGREPHASH_H
#define FASTQGREPHASH_H

#include "fastqGrepAVLTree.h" /*Structs & funs to build an read name AVL tree*/
    /*
      Includes: 
          - <string.h>
          - "fastqGrepStructs.h"
              - <stdlib.h>
              - <stdio.h>
    */

/*##############################################################################
# Output:
#    Returns: pointer to array of readInfo nodes that is the hash
#             0: if failed
#    Modifies: hashSizeUlng to hold hte length of the returned array
#    Modifies: numCharChar to hold the number of characters used in hash
#    Modifies: digPerKeyUChar to hold the hash size (as multiple of two)
##############################################################################*/
struct readInfo ** makeReadHash(
    FILE * filtFILE,                  /*file with read id's to filter by*/
    char * buffCStr,                  /*Buffer to hold input from file*/
    unsigned long lenBuffULng,        /*Length of buffer (buffCStr)*/
    struct readNodeStack *readStackAry,  /*Stack, (as array) for searching*/
    unsigned long *hashSizeULng,      /*Will hold Size of hash table*/
    unsigned char *digPerKeyUChar,    /*Power of two hash size is at*/
    unsigned long *majicNumULng,      /*Holds majick number for kunths hash*/
    char *failedChar                  /*Tells if did not make hash table*/
); /*Makes a read hash array using input read ids*/

/*##############################################################################
# Output:
#   Returns: integer holding the hashed value
##############################################################################*/
unsigned long calcHash(
    struct bigNum *idBigNum,     /*read id converted to number to hash*/
    unsigned long majicNumULng,  /*Majick number to mulitply by*/
    unsigned char digPerKeyUChar  /*Hash table size is 2^digPerKeyUChar power*/
); /*Calculate the hast for a read id*/

/*##############################################################################
# Output:
#    Modifies: Inserts readNode into the hashTbl.
##############################################################################*/
void insertHashEntry(
    struct readInfo **hashTbl,     /*Hash table to insert read into*/
    struct readInfo *readNode,     /*readNode to insert into hash table*/
    unsigned long majicNumULng,            /*Majick number to mulitply by*/
    unsigned char digPerKeyUChar,  /*Hash table size is 2^digPerKeyUChar power*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
); /*Iinserts a read into a hash table*/

/*##############################################################################
# Output:
#   Returns: readInfo node if found read in hash table, otherwise 0
##############################################################################*/
struct readInfo * findReadInHashTbl(
    struct bigNum *idBigNum,       /*read id converted to number to hash*/
    unsigned long majicNumULng,            /*Majick number to mulitply by*/
    unsigned char digPerKeyUChar,  /*Hash table size is 2^digPerKeyUChar power*/
    struct readInfo **hashTbl      /*Hash table to search for read in*/
);  /*find a readInfo node in a hash table using the read id*/

/*##############################################################################
# Output:
#    Frees: hash table & sets hashTblToFree  to 0
##############################################################################*/
void freeHashTbl(
    struct readInfo ***hashTblToFree,
    unsigned long hashSizeULng,       /*Size of hash table*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
); /*Fres a hash table of read trees*/

#endif
