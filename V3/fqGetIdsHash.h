/*##############################################################################
# Name: fastqGrepHash.c
#   Use:
#     Uses an hash + AVL tree to find if reads is in user suplied list
##############################################################################*/

#ifndef FQGREPHASH_H
#define FQGREPHASH_H

#include "fqGetIdsFqFun.h" /*Big number conversion*/
#include "fqGetIdsAVLTree.h" /*Structs & funs to build an read name AVL tree*/
    /*
      Includes: 
          - <string.h>
          - "fqGetIdsStructs.h"
              - <stdlib.h>
              - <stdio.h>
              - <stdint.h>
    */

/*##############################################################################
# Output:
#    Returns: pointer to array of readInfo nodes that is the hash
#             0: if failed
#    Modifies: hashSizeUlng to hold hte length of the returned array
#    Modifies: numCharChar to hold the number of characters used in hash
#    Modifies: digPerKeyUChar to hold the hash size (as multiple of two)
#    Modifies: magickNumULng to hold the magick number for this hash table
##############################################################################*/
struct readInfo ** makeReadHash(
    FILE * filtFILE,               /*file with read id's to filter by*/
    char * buffCStr,            /*Buffer to hold input from file*/
    uint64_t lenBuffULng,          /*Length of buffer (buffCStr)*/
    struct readNodeStack *readStackAry, /*Stack, (as array) for search*/
    uint64_t *hashSizeULng,       /*Will hold Size of hash table*/
    uint8_t *digPerKeyUChar,      /*Power of two hash size is at*/
    unsigned long *majicNumULng,  /*Holds majick number for kunths hash*/
    uint8_t *failedChar           /*Tells if did not make hash table*/
); /*Makes a read hash array using input read ids*/

/*##############################################################################
# Output:
#   Returns: integer holding the hashed value
##############################################################################*/
uint64_t calcHash(
    struct bigNum *idBigNum, /*read id converted to number to hash*/
    const unsigned long *majicNumULng, /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar /*Hash table size 2^digPerKeyUChar*/
); /*Calculate the hast for a read id*/

/*##############################################################################
# Output:
#    Modifies: Inserts readNode into the hashTbl.
##############################################################################*/
void insertHashEntry(
    struct readInfo **hashTbl,   /*Hash table to insert read into*/
    struct readInfo *readNode,   /*readNode to insert into hash table*/
    const unsigned long *majicNumULng,  /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar, /*Hash table size 2^digPerKeyUChar*/
    struct readNodeStack *readStack /*Stack, (as array) for searching*/
); /*Iinserts a read into a hash table*/

/*##############################################################################
# Output:
#   Returns: readInfo node if found read in hash table, otherwise 0
##############################################################################*/
struct readInfo * findReadInHashTbl(
    struct bigNum *idBigNum,       /*read id converted to number to hash*/
    const unsigned long *majicNumULng, /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar,/*Hash table size is 2^digPerKeyUChar*/
    struct readInfo **hashTbl     /*Hash table to search for read in*/
);  /*find a readInfo node in a hash table using the read id*/

/*##############################################################################
# Output:
#    Frees: hash table & sets hashTblToFree  to 0
##############################################################################*/
void freeHashTbl(
    struct readInfo ***hashTblToFree,
    const uint64_t *hashSizeULng,     /*Size of hash table*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
); /*Frees a hash table of read trees*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    Returns:                                                          |
|        - readInfo: hash table (as heap array)                        |
|        - 0: if failed                                                |
|    Modifies:                                                         |
|        - hashSizeUlng to hold length of the returned array           |
|        - numCharChar: to hold the number of characters used in hash  |
|        - digPerKeyUChar: to hold the hash size (as multiple of two)  |
|        - magickNumULng: to hold magick number for this hash table    |
\---------------------------------------------------------------------*/
struct readInfo ** readListToHash(
    struct readInfo * readList,   /*List of id's to build hash from.
                                  only use the rightChild ptr in list*/
    const uint64_t *lenListULng,   /*Number of id's in readList*/
    struct readNodeStack *readStackAry, /*Stack (array) for searching*/
    uint64_t *hashSizeULng,  /*Will hold Size of hash table*/
    uint8_t *digPerKeyUChar,       /*Power of two hash size is at*/
    unsigned long *majicNumULng  /*Holds majick number for kunths hash*/
); /*Makes a read hash array using input read ids*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o Unsigned long with the first 64 bits of the golden number
\---------------------------------------------------------------------*/
unsigned long findMajicNumber();
/*Gets first 64 bits of the goldent ratio*/

#endif
