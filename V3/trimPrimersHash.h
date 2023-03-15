/*######################################################################
# Name: trimPrimersHash
#   - Use:
#     o Sets up the hash table for the trimPrimers progam/function. It
#       also is has functions for freeing or searching the hash table
# Includes:
#   - "fqGetIdsHash.h"
#   - "trimPrimersAVLTree.h"
#   - "defaultSettings.h"
#   - "cStrFun.h"
#   o "fqGetIdsFqFun.h"
#   o "trimPrimersStructs.h"
#   o "fqGetIdsStructs.h"
#   o "cStrToNumberFun.h"
#   o "fqAndFaFun.h"
#   o "FCIStatsFun.h"     (fqAndFqFun.h)
#   o "minAlnStats.h"     (fqAndFaFun.h)
#   o "samEntryStruct.h"  (fqAndFaFun.h->FCIstatsFun.h)
#   o "printError.h"      (fqAndFaFun.h)
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#ifndef TRIMPRIMERSHASH_H
#define TRIMPRIMERSHASH_H

#include "fqGetIdsHash.h"    /*For hashing functions and trees*/
#include "defaultSettings.h" /*For minimap2 command*/
#include "cStrFun.h"         /*For copying strings*/

/*Structs & funs to build an read name AVL tree*/
#include "trimPrimersAvlTree.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimPrimersHash SOH: Start Of Header
'   - st-1 hashTblVar
'     o Holds the universal variables used in my hashing functions
'   o st-2: readPrimHashTbl
'     - Holds all hashing & tree variables used in trimPrimers
'   o fun-1 makeReadPrimHash:
'     - Make a hash table of read id's from a file of read ids
'   o fun-2 readPrimListToHash:
'     - Convert a lined readList to a hash table
'     - Only use readList->rightChild pointer in the linked list
'   o fun-3 insReadPrimSTInHash:
'     - Insert a read id (is big number) into a hash table
'   o fun-4 findReadPrimInHash:
'     -  Search hash table for a particler read id
'   o fun-5 freeReadPrimHash:
'     - Free a hash table
'   o fun-6 initReadPrimHashST:
'     - Initialize a readPrimHash sructer for building a hash table
'   o fun-7 freeReadPrimHashST:
'     - free a readpPrimHash structure
'   o fun-8 initHashTblVarST:
'     - Sets all variables in a hashTblVar to 0
'   o fun-9 freeHashTblVarST:
'     - "free" a hashTblVar structuer (here if needed for future)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| ST-1: hashTblVar
| Use: Holds the universal variabes used in my hashing functions
\---------------------------------------------------------------------*/
typedef struct hashTblVar
{ /*readPrimHashTbl*/
    unsigned long numIdsUL;      /*Number of ids in the hash table*/
    unsigned long majicNumUL;    /*Majic number to use in hashing*/
    unsigned long lenHashUL;     /*Size of the hash table*/
    unsigned char lenLog2HashUC; /*log2(lenHashUL)*/
}hashTblVar; /*readPrimHashTbl*/


/*---------------------------------------------------------------------\
| ST-1: readPrimHashTbl
| Use: Holds all hashing & tree variables used in trimPrimers
\---------------------------------------------------------------------*/
typedef struct readPrimHash
{ /*readPrimHashTree*/
    struct hashTblVar hashVarST;  /*Gereral hash variables*/
    struct readPrim *readTree;    /*AVL tree or list*/
    struct readPrim **hashTbl;    /*Hash table*/
    struct readPrimStack readStack[defLenStack];
      /*For hash and tree search & free functions*/
      /*defLenStack from fqGetIdsHash*/
}readPrimHash; /*readPrimHashTree*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 if succeded
|     o 2 if could not open the fasta file (primer sequences)
|     o 4 if could not open the fastq file (reads)
|     o 64 for memory allocation errors
|   - Modifies:
|     o Initializes hashST, so make sure that thier is no allocated 
|       memory in hashTbl & readTree (both set to 0)
|     o Variables in hashST to hold the readPrim list (readTree) and 
|       length of the readPrim list (hashVarST.numIdsUL)
\---------------------------------------------------------------------*/
unsigned char makeReadPrimList(
    char *primFaFileCStr,     /*Path to fasta file with primers*/
    FILE *pafFILE,
       /*Paf file to get ids & primer coordinates from; skips minimap2*/
    char *fqFileCStr,         /*Path to fastq file with reads*/
    char *threadsCStr,        /*Number of threads to use with minimap2*/
    struct readPrimHash *hashST /*Holds hash variables*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-1 TOC:
   '  - Make a readPrim list of read ids & primer mappings from a fasta
   '    file with primers & a fastq file of reads. Uses minimap2.
   '  o fun-1 sec-1: variable declerations
   '  o fun-1 sec-2: Check if the fasta and fastq file exists
   '  o fun-1 sec-3: Setup minimap2 command & initalize variables
   '  o fun-1 sec-4: Get read ids in file & convert to readPrim list
   '  o fun-1 sec-5: Close file ane return success
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 for succes
|     o 64 for memory allocation error
|   - Modifies:
|     o Variables in hashST to support a hash table
|     o hashST->readTree is set to 0
\---------------------------------------------------------------------*/
unsigned char readPrimListToHash(
    char listOnHeapBl,
        /*1 read id list on heap, ok to free duplicates, 0 do not free*/
    struct readPrimHash *hashST  /*Hash input to build hash table with*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-2 TOC: readPrimListToHash
   '   - Make a read id hash table from a readPrim list of read ids
   '   o fun-2 sec-1: Variable declerations
   '   o fun-2 sec-2: Find the hash table size & make hash table
   '   o fun-2 sec-3: Build hash
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: Inserts readNode into the hashTbl.
\---------------------------------------------------------------------*/
unsigned char insReadPrimSTInHash(
    struct readPrim *readNode,   /*readNode to insert into hash table*/
    struct readPrimHash *hashST  /*has hash table & hashing variables*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-3 TOC: Sec-1 Sub-1: insReadPrimSTInHash
   '   - Inserts readPrim node with big number read id into hash table
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o readPrim node if found read id in the hash table
|     o 0 If read id is not in the hash table
\---------------------------------------------------------------------*/
struct readPrim * findReadPrimInHash(
    struct bigNum *idBigNum,       /*big number read id to find*/
    struct readPrimHash *hashST    /*has hash table & needed variables*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-4 TOC: Sec-1 Sub-1: findReadPrimInHash
   '   - Finds a big number read id in a hash table of read ids
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Frees:
|    o Hash table if hash table has entries
|    o Read tree if hash table is empty (set to 0)
|  - Modifies:
|    o hashST->hashTbl is set to 0  (to mark has been freed)
|    o hashST->readTree is set to 0 (to mark has been freed)
\---------------------------------------------------------------------*/
void freeReadPrimHashTblOrTree(
    struct readPrimHash *hashST
        /*Has hash table to free & other variables for freeing*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-5 TOC: Sec-1 Sub-1: freeHashTbl
   '   - Frees a hash table or readTree (if hash table empty)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o readTree to be 0
|     o hashTbl to be 0
|     o readStack to have first stack set to 0
|     o majicNumUL to hold the majic number for the hash
|     o lenHashUL to be 0
|     o lenLog2HashUL to be 0
\---------------------------------------------------------------------*/
void initReadPrimHashST(
    struct readPrimHash *hashST /*Structer to initialize*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-6 TOC: Sec-1 Sub-1: initReadPrimHashST
   '   - Sets all variables for hashing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: If on heap, frees the structer, else does nothing
| WARNING:
|  o hashST is not set to 0 (you must do this)
|  o This functions assumes that the hash table is on the heap.
|  o If the hash table is not on the heap, then you will need to free
|    things manually
\---------------------------------------------------------------------*/
void freeReadPrimHashST(
    char stOnHeapBl,              /*1: is on heap; 0 on stack*/
    char hashElmOnHeapBl,
        /*1: Nodes in hash table or tree are on the heap; 0 on stack*/
    struct readPrimHash *hashST    /*Structure to free*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-7 TOC: Sec-1 Sub-1: freeReadInfoHashST
   '   - Frees a readInfoHash structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o majicNumUL to hold the majic number for the hash
|     o lenHashUL to be 0
|     o lenLog2HashUL to be 0
\---------------------------------------------------------------------*/
void initHashTblVarST(
    struct hashTblVar *hashTblVarST /*Structure to defaults*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-8 TOC: Sec-1 Sub-1: initHashTblVarST
   '   - Sets all values in hastTblVarST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: If on heap, frees the structer, else does nothing
| WARINGING: the hashTblVarST pointer is not set to 0 (you must do this)
\---------------------------------------------------------------------*/
void freeHashTblVarST(
    char onHeapBl,                  /*1: is on heap; 0 on stack*/
    struct hashTblVar *hashTblVarST /*Structure to free*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-9 TOC: Sec-1 Sub-1: freeHashTblVarST
   '   - Frees a hashTblVarST structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
