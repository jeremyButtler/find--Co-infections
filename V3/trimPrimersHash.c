/*######################################################################
# Name: trimPrimersHash
#   - Use:
#     o Sets up the hash table for the trimPrimers progam/function. It
#       also is has functions for freeing or searching the hash table
# Includes:
#   - "fqGetIdsHash.h"
#   - "trimPrimersAVLTree.h"
#   o "fqGetIdsFqFun.h"
#   o "trimPrimersStructs.h"
#   o "fqGetIdsStructs.h"
#   o "fqAndFaFun.h"      (Only one function is used)
#   o "FCIStatsFun.h"     (fqAndFqFun.h)
#   o "minAlnStats.h"     (fqAndFaFun.h)
#   o "samEntryStruct.h"  (fqAndFaFun.h->FCIstatsFun.h)
#   o "cStrToNumberFun.h" (fqAndFaFun.h)
#   o "printError.h"      (fqAndFaFun.h)
#   o "defaultSettings.h"
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#include "trimPrimersHash.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimPrimersHash SOF: Start Of Functions
'   o fun-1 makeReadPrimHash:
'     - Make a hash table of read id's from a file of read ids
'   o fun-2 readPrimListToHash:
'     - Convert a lined readList to a hash table
'     - Only use readList->rightChild pointer in the linked list
'   o fun-3 insReadPrimSTInHash:
'     - Insert a read id (is big number) into a hash table
'   o fun-4 findReadPrimInHash:
'     -  Search hash table for a particler read id
'   o fun-5 freeReadPrimHashTblOrTree:
'     - Free a hash table or tree in a readPrimHash structer
'   o fun-6 initReadPrimHashST:
'     - Initialize a readPrimHash sructer for building a hash table
'   o fun-7 freeReadPrimHashST:
'     - free a readpPrimHash structure
'   o fun-8 initHashTblVarST:
'     - Sets all variables in a hashTblVar to 0
'   o fun-9 freeHashTblVarST:
'     - "free" a hashTblVar structuer (here if needed for future)
/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-1 TOC:
   '  - Make a readPrim list of read ids & primer mappings from a fasta
   '    file with primers & a fastq file of reads. Uses minimap2.
   '  o fun-1 sec-1: variable declerations
   '  o fun-1 sec-2: Check if the fasta and fastq file exists
   '  o fun-1 sec-3: Setup minimap2 command & initalize variables
   '  o fun-1 sec-4: Get read ids in file & convert to readPrim list
   '  o fun-1 sec-5: Close file ane return success
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-1 Sec-1: variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned short lenBuffUS = 2048;
   char buffCStr[lenBuffUS];
   char minimap2CmdCStr[2048];   /*Command for running minimap2*/
   char *tmpCStr = 0;
   unsigned char maxHexChar = 1; /*Max hex digits to read in*/

   struct bigNum *bigNumST = 0;
   struct primCord *primCordST = 0;/*For holding primer cordinates*/
   struct readPrim *readPrimST = 0;/*readPrim structure to add to list*/

   FILE *stdinFILE = 0;            /*Points to minimap2 output*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-1 Sec-2: Check if the fasta and fastq file exists
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initReadPrimHashST(hashST);

    if(pafFILE == 0)
    { /*If have to call minimap2*/
        stdinFILE = fopen(primFaFileCStr, "r");

        if(stdinFILE == 0) return 2; /*File could not be opened*/

        fclose(stdinFILE); /*No longer need open, minimap2 will handel it*/
        stdinFILE = 0;

        stdinFILE = fopen(fqFileCStr, "r");

        if(stdinFILE == 0) return 4; /*File could not be opened*/

        fclose(stdinFILE);
        stdinFILE = 0;
    } /*If have to call minimap2*/


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-1 Sec-3: Set up command to run minimap2 & initalize variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(pafFILE == 0)
   { /*If running minimap2*/
       tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, defMinimap2PrimCMD);
       tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
       cpParmAndArg(tmpCStr, primFaFileCStr, fqFileCStr);
       stdinFILE = popen(minimap2CmdCStr, "r");
       /*copy the files to the command*/
   } /*If running minimap2*/

   else stdinFILE = pafFILE;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-1 Sec-4: Get read ids in file & convert to readPrim list
   ^   o sec-4 sub-1: Convert read id to a big number
   ^   o sec-4 sub-2: Find starting & ending coordinates of primer map
   ^   o sec-4 sub-3: Check if is a new id or repeated mapping
   ^   o sec-4 sub-4: Make new readPrim structure for new read id
   ^   o sec-4 sub-5: Add the new readPrim structure to the list
   ^   o sec-4 sub-6: Make sure that I grabbed the entire line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Sec-4 Sub-1: Convert read id to a big number
   \*******************************************************************/

   buffCStr[lenBuffUS - 1] = '\0';
   buffCStr[lenBuffUS - 2] = '\0';

   /*2048 bytes will be enough to read in the parts of the paf file I
     care about (read id, query start, and query end*/
   while(fgets(buffCStr, lenBuffUS, stdinFILE))
    { /*While ids to read in*/
       /*Convert read id to big number*/
       bigNumST = buffToBigNum(buffCStr, &tmpCStr, &maxHexChar);
           /*maxHexChar will hold the number of limbs needed to hold
               largest big number read in.
             tmpCStr will point to the end of the read id*/

       if(bigNumST == 0)
       { /*If I falied to convert the string to a big number*/
           freeBigNumStruct(&bigNumST);

           if(pafFILE == 0) pclose(stdinFILE);

           return 64;
       } /*If I falied to convert the string to a big number*/

       /***************************************************************\
       * Sec-4 Sub-2: Find starting & ending coordinates of primer map
       \***************************************************************/

       /*Read in the primer coordinates*/
       primCordST = makePrimCord();/*Make stucture to hold coordinates*/

       ++tmpCStr; /*Get off the tab after the read id*/

       /*Move off the length entry to starting position on the query*/
       while(*tmpCStr > 32) ++tmpCStr;

       ++tmpCStr; /*Get off the tab*/

       /*Capture the starting position of the primer to query map*/
       tmpCStr = cStrToUInt(tmpCStr, &primCordST->startUI);
       ++tmpCStr; /*Get off the tab*/

       /*Capture the ending position of the primer to query map*/
       tmpCStr = cStrToUInt(tmpCStr, &primCordST->endUI);

       /***************************************************************\
       * Sec-4 Sub-3: Check if is a new id or repeated mapping
       \***************************************************************/

       /*Check if this read id has already been used*/
       if(hashST->readTree != 0 &&
          cmpBigNums(bigNumST, hashST->readTree->idBigNum) == 0
       ) { /*If it is another entry for the same id*/
           insPrimCordST(primCordST, &hashST->readTree->primCordST);
           freeBigNumStruct(&bigNumST); /*No longer need*/
       } /*If it is another entry for the same id*/

       else
       { /*Else is a new read id, make a new structer*/

           /***********************************************************\
           * Sec-4 Sub-4: Make new readPrim structure for new read id
           \***********************************************************/

           /*Add the new unique read id to the count*/
           ++(hashST->hashVarST.numIdsUL);

           /*Make a new readPrim structer to add to the list*/
           readPrimST = makeReadPrimST();

           /*Memory allcoation error, need to free stuff*/
           if(readPrimST == 0)
           { /*If have to free structers*/
               /*Free structures and close files*/
               freePrimCordST(&primCordST);
               freeBigNumStruct(&bigNumST);

               if(pafFILE == 0) pclose(stdinFILE);

               return 64;
           } /*If have to free structers*/

           /***********************************************************\
           * Sec-4 Sub-5: Add the new readPrim structure to the list
           \***********************************************************/

           readPrimST->idBigNum = bigNumST;     /*Set read id*/
           readPrimST->primCordST = primCordST;
               /*copy over primer start and end for this id*/

            /*Add the new readPrim to the list*/
            if(hashST->readTree == 0)
                hashST->readTree = readPrimST;    /*First node in list*/
            else
            { /*Else their are other nodes in the list*/
                readPrimST->rightChild = hashST->readTree;
                hashST->readTree = readPrimST;    /*new head of list*/
            } /*Else their are other nodes in the list*/
       } /*Else is a new read id, make a new structer*/

       /***************************************************************\
       * Sec-4 Sub-6: Make sure that I grabbed the entire line
       \***************************************************************/

        while(
            buffCStr[lenBuffUS - 2] != '\n' &&
            buffCStr[lenBuffUS - 2] != '\0'
        ) { /*While not on next entry*/
           /*See if can grab the rest of the line*/
           tmpCStr = fgets(buffCStr, lenBuffUS, stdinFILE);

           if(tmpCStr == 0)
           { /*If at the end of the file*/
               if(pafFILE == 0) pclose(stdinFILE);
              
               return 1; /*end of file*/
           } /*If at the end of the file*/
        } /*While not on next entry*/

       buffCStr[lenBuffUS - 2] = '\0'; /*Reset the line marker*/
   } /*While ids to read in*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-1 Sec-5: Close file ane return success
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(pafFILE == 0) pclose(stdinFILE);
   return 1;  /*Return success*/
} /*makeReadPrimHash function*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-2 TOC: readPrimListToHash
   '   - Make a read id hash table from a readPrim list of read ids
   '   o fun-2 sec-1: Variable declerations
   '   o fun-2 sec-2: Find the hash table size & make hash table
   '   o fun-2 sec-3: Build hash
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-2 Sec-1: variable declerations                                v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char errUC = 0;
   struct readPrim *tmpRead = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-2 Sec-2: Find hash table size & make hash tabel               v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Just using this as a temporary variable*/
   hashST->hashVarST.lenHashUL = hashST->hashVarST.numIdsUL;

   /*Find power of two to use*/
   while(hashST->hashVarST.lenHashUL > 0)
   { /*While need to find the next power of two*/
       hashST->hashVarST.lenHashUL = hashST->hashVarST.lenHashUL >> 1;
       ++hashST->hashVarST.lenLog2HashUC;
   } /*While need to find the next power of two*/

   /*Make hash table sparase by using double the needed size*/
   ++hashST->hashVarST.lenLog2HashUC;

   /*Find the length of the hash table*/
   hashST->hashVarST.lenHashUL = 1 << hashST->hashVarST.lenLog2HashUC;

   hashST->hashTbl =
      calloc(
          (hashST->hashVarST.lenHashUL + 1),
          sizeof(struct readPrim *)
   ); /*Set up the hash table, Calloc intalizes with 0's*/

   if(hashST->hashTbl == 0) return 64; /*Memory alloction error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-2 Sec-3: Build hash                                           v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(hashST->readTree != 0)
   { /*While their are reads to put in hash table*/
       tmpRead = hashST->readTree->rightChild;
       hashST->readTree->rightChild = 0;

       errUC = insReadPrimSTInHash(hashST->readTree, hashST);

       /*Free the read node if it was a duplicate node*/
       if(errUC == 0 && listOnHeapBl & 1)
           freeReadPrimST(&hashST->readTree);

       hashST->readTree = tmpRead;  /*move to the next read*/
   } /*While their are reads to put in hash table*/

   return 1;  /*Return head of hash table*/
} /*makeReadPrimHash function*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies: Inserts readNode into the hashTbl.
|  - Returns:
|    - 0 IF read is a duplicate
|    - 1 If read is unique
\---------------------------------------------------------------------*/
unsigned char insReadPrimSTInHash(
    struct readPrim *readNode,   /*readNode to insert into hash table*/
    struct readPrimHash *hashST  /*has hash table & hashing variables*/
) { /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: Sec-1 Sub-1: insReadPrimSTInHash
    '   - Inserts readPrim node with big number read id into hash table
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    unsigned long hashUL = 0; /*Hold the hash value*/

    /*Find the index in the hash table this id will be/is stored at*/
    hashUL =
       calcHash(
           readNode->idBigNum,
           &hashST->hashVarST.majicNumUL,
           &hashST->hashVarST.lenLog2HashUC
    );

   if(*(hashST->hashTbl + hashUL) == 0) /*If 1st node at hash value*/
   { /*If is the fist node in the tree*/
       *(hashST->hashTbl + hashUL) = readNode;
       return 1;
   } /*If is the fist node in the tree*/

   return
      avlInsPrimReadST(
          readNode,
          hashST->hashTbl + hashUL,
          hashST->readStack
   ); /*Return the result of the insertion (0 duplicate, 1 unique)*/
} /*insReadPrimSTInHash*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o readPrim node if found read id in the hash table
|     o 0 If read id is not in the hash table
\---------------------------------------------------------------------*/
struct readPrim * findReadPrimInHash(
    struct bigNum *idBigNum,       /*big number read id to find*/
    struct readPrimHash *hashST    /*has hash table & needed variables*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-4 TOC: Sec-1 Sub-1: findReadPrimInHash
   '   - Finds a big number read id in a hash table of read ids
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readPrim *tmpTree = 0;
    unsigned long hashUL = 0;

    hashUL =
       calcHash(
           idBigNum,
           &hashST->hashVarST.majicNumUL,
           &hashST->hashVarST.lenLog2HashUC
    ); /*Find hash value and grab readPrim AVL tree at hash value*/

    tmpTree = *(hashST->hashTbl + hashUL);

   /*See if the read id exists in the AVL tree pulled from the table*/
    return searchReadPrimTree(idBigNum, tmpTree);
} /*findReadPrimInHash*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-5 TOC: Sec-1 Sub-1: freeHashTbl
   '   - Frees a hash table or readTree (if hash table empty)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct readPrim *readTree = 0;
    unsigned long numElmInTblUL = 0;

    if(hashST->hashTbl != 0)
    { /*If have a hash table to free*/
        readTree = *(hashST->hashTbl);

        while(numElmInTblUL <= hashST->hashVarST.numIdsUL)
        { /*While there are entries to free in the hash table*/
            if(readTree != 0)
                freeReadPrimTree(&readTree, hashST->readStack);

            /*Find the next tree in the hash table*/
            numElmInTblUL++;
            readTree = *(hashST->hashTbl + numElmInTblUL);
        } /*While there are entries to free in the hash table*/

        free(hashST->hashTbl);
        hashST->hashTbl = 0;  /*So user knows their is nothing here*/
        hashST->readTree = 0; /*All these nodes were in the hash table*/

        return;
    } /*If have a hash table to free*/

    if(readTree != 0)
        freeReadPrimTree(&hashST->readTree, hashST->readStack);

    hashST->readTree = 0; /*So user knows their is nothing here*/
    hashST->hashTbl = 0;  /*Make sure is 0*/
 
    return;
} /*freeHashTblOrTree*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o readTree to be 0
|     o hashTbl to be 0
|     o readStack to have first stack set to 0
|     o majicNumUL to hold the majic number for the hash
|     o lenHashUL to be 0
|     o lenLog2HashUC to be 0
\---------------------------------------------------------------------*/
void initReadPrimHashST(
    struct readPrimHash *hashST /*Structer to initialize*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-6 TOC: Sec-1 Sub-1: initReadPrimHashST
   '   - Sets all variables for hashing
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   initHashTblVarST(&hashST->hashVarST);
   hashST->readTree = 0;
   hashST->hashTbl = 0;

   /*Make sure start & end of my stacks are marked*/
   hashST->readStack[0].readNode = 0;
   hashST->readStack[defLenStack - 1].readNode = 0;

   return;
} /*initReadPrimHashST*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-7 TOC: Sec-1 Sub-1: freeReadInfoHashST
   '   - Frees a readInfoHash structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   unsigned long numElmInTblUL = 0;

   if(hashElmOnHeapBl & 1)
       freeReadPrimHashTblOrTree(hashST); /*Free the hash table/tree*/
   else
   { /*Else only the hash table (not the nodes) are on the heap*/
       if(hashST->hashTbl != 0)
       { /*If have a hash table to free*/

           /*Make sure I only free the hash table*/
           while(numElmInTblUL <= hashST->hashVarST.numIdsUL)
           { /*While their are nodes in the hash table to blank*/
               *(hashST->hashTbl + numElmInTblUL) = 0;
               numElmInTblUL++;
           } /*While their are nodes in the hash table to blank*/

           free(hashST->hashTbl);
       } /*If have a hash table to free*/
   } /*Else only the hash table (not the nodes) are on the heap*/

   hashST->hashTbl = 0;
   hashST->readTree = 0;
   freeHashTblVarST(0, &hashST->hashVarST);

   if(stOnHeapBl & 1)
       free(hashST);

   return;
} /*freeReadInfoHashST*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o majicNumUL to hold the majic number for the hash
|     o lenHashUL to be 0
|     o lenLog2HashUC to be 0
\---------------------------------------------------------------------*/
void initHashTblVarST(
    struct hashTblVar *hashTblVarST /*Structure to defaults*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-8 TOC: Sec-1 Sub-1: initHashTblVarST
   '   - Sets all values in hastTblVarST to 0
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   hashTblVarST->numIdsUL = 0;
   hashTblVarST->majicNumUL = findMajicNumber();  
   hashTblVarST->lenHashUL = 0;  
   hashTblVarST->lenLog2HashUC = 0;  

   return;
} /*initHashTblVarST*/

/*---------------------------------------------------------------------\
| Output: If on heap, frees the structer, else does nothing
| WARINGING: the hashTblVarST pointer is not set to 0 (you must do this)
\---------------------------------------------------------------------*/
void freeHashTblVarST(
    char onHeapBl,                  /*1: is on heap; 0 on stack*/
    struct hashTblVar *hashTblVarST /*Structure to free*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-9 TOC: Sec-1 Sub-1: freeHashTblVarST
   '   - Frees a hashTblVarST structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(onHeapBl & 1)
       free(hashTblVarST);

   return;
} /*freeHashTblVarST*/
