/*######################################################################
# Name: fastqGrepSearchThread
# Use:
#    Extracts target reads from fastq file while using multiple threads
# Requires:
#    fastqGrepAVLTree
#    fastqGrepStructs (called by fastqGrepAVLTree)
#    fqGetIdsSearchFq (single thread version + some needed functions)
######################################################################*/

#ifndef FQGREPSEARCHTHREAD_H
#define FQGREPSEARCHTHREAD_H

#include <pthread.h> /*For multi-threading*/
#include "fqGetIdsSearchFq.h"

/*---------------------------------------------------------------------\
| Struct-1: getFileLenStruct
|   o Here so that I can get the file length when using multiple threads
\---------------------------------------------------------------------*/
typedef struct getFileLenStruct
{ /*getFileLenStruct*/
    unsigned long lenFileUL;
    FILE *inFILE;
}getFileLenStruct;

/*---------------------------------------------------------------------\
| Struct-2: extReadsST
|   o Structer to hold parameters for a multithread read extract
\---------------------------------------------------------------------*/
typedef struct extReadsST
{ /*extReadsST*/
    uint8_t retValUC;        /*Function error out state*/
    unsigned long *endPosUL; /*Point to end read extraction at*/
    unsigned long *startPosUL; /*Point to start read extraction at*/
    uint32_t lenBuffUI;  /*size of buffer to create*/
    unsigned long majicNumUL; /*Magic number for kunths multiply hash*/
    uint8_t digPerKeyUC;   /*Digits needed to get a key*/
    uint8_t printNonMatchBl; /*1: print non-match, 0: print match*/

    struct readInfo *readTree;  /*For AVL search (hashTbl == 0)*/
    struct readInfo **hashTbl;  /*Hash table to search for ids in*/

    FILE *fqFILE;
    FILE *outFILE;
}extReadsST;


/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 2 if invalid filter file
|     o 4 if invalid input fastq file
|     o 8 if could not open the output file
|     o 16 if both filter and fastq file coming from stdin
\---------------------------------------------------------------------*/
uint8_t fastqThreadExtract(
    char *filtPathCStr,        /*Path to file with read ids to extract*/
    char *fqPathCStr,          /*Path to fastq file to extract reads*/
    char *outPathCStr,         /*Path to fastq file to to write reads*/
    unsigned char threadsUC,   /*Number of threads to use*/
    uint8_t sizeReadStackUC,   /*Number of elements to use in stack*/
    uint32_t lenBuffUI,        /*Size of buffer to read input with*/
    uint8_t hashSearchC,    /*1: do hash search, 0: do Tree search*/
    uint8_t printReverseC   /*1: Keep reads in filter file
                                 0: ingore reads in filter file*/
); /*Searches and extracts reads from a fastq file using read id's*/

/*---------------------------------------------------------------------\
| Output:
|    stdout: Prints out reads in hash table
\---------------------------------------------------------------------*/
void * extractReadsThread(
    void *parmST /*extReadsST Structer with parameters*/
); /*Extract target reads from fastq file with hash table or tree*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o lenFileUL in parmStruct to hold the file length
\---------------------------------------------------------------------*/
void * getFileLen(
   void *parmStruct /*Holds file & return varaibles*/ 
); /*Get the length of the file in parmStruct. Here for multi threading*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 if have more file to read into the buffer
|     o 0 if read remaning part of file into buffer
|   - Modifies:
|     o startOfBuffCStr to piont to start of the next read
|     o If neeed to read in more bytes, will modify buffCStr
|     o lenInUL to hold the number of bytes read into buffCStr
|     o extParmST->startPosUL to point to header of the next read
\---------------------------------------------------------------------*/
unsigned char findStartPos(
    char **buffCStr,             /*Pionts to the end of the buffer*/
    char **startOfBuffCStr,      /*Pionts to the start of the buffer*/
    uint64_t *lenInUL,           /*Number bytes read from fread*/
    struct extReadsST *extParmST /*Holds parameters to use*/
); /*Finds the next read after the starting position*/

#endif
