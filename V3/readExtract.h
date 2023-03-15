/*######################################################################
# Use:
#    Holds functions for doing read id extractions for find
#    co-infections
# Includes:
#   o "FCIStatsFun.h"
#   o "minAlnStats.h"
#   o "cStrFun.h"
#   o "fqAndFaFun.h"
#     - "samEntryStruct.h" (See score reads entry)
#   o "trimSam.h"
#     - "samEntryStruct.h" (See score reads entry)
#   o "findCoInftBinTree.h"
#     - <string.h>
#     - <stdlib.h>
#     - <stdio.h>
#     - <stdint.h>
#   o "findCoInftChecks.h
#     - "defaultSettings.h"
#     - "scoreReadsFun.h"
#        o "samEntryStruct.h"
#          - <stdlib.h>
#          - "cStrToNumberFun.h"
#            o <sdtint.h>
#          - "printError.h"
#            o <stdio.h>
#   o fqGetIdsSearchFq
#     - "fqGetIdsFqFun.h"
#        o "fqGetIdsStructs.h"
#     -fqGetIdsHash.h
#        o fqGetIdsAVLTree.h:
#          - <string.h>
#          - "fqGetIdsStructs.h"
#            o <stdlib.h>
#            o <stdio.h>
#            o <stdint.h>
######################################################################*/

#ifndef READEXTRACT_H
#define READEXTRACT_H

#include "trimSam.h"           /*For trimming reads*/
#include "cStrFun.h"           /*copying c-strings*/
#include "findCoInftChecks.h"  /*Read quality checks & mindiff struct*/
#include "cStrFun.h"           /*C-string functions*/
#include "fqAndFaFun.h"        /*Fastq and fasta functions*/
#include "findCoInftBinTree.h" /*for readBin struct*/
#include "fqGetIdsSearchFq.h"  /*For extracting reads by id*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Prints:                                                         |
|        - The read with the highest medain Q-score to the fastq file  |
|          name stored in binIn->bestReadCStr                          |
|    - Modifies:                                                       |
|        - Fastq file binIn-fqPathCStr to not have the best read       |
|        - Stats file binIn-statsPathCStr to not have the best read    |
|    - Returns:                                                        |
|        - 1: If sucessfull                                            |
|        - 2: For blank structer                                       |
|        - 4: For no fastq file                                        |
|        - 8: For no stats file                                        |
|        - 16: For error when extracting stats                         |
|        - 32: If could not open a temporary fastq file                |
|        - 128: If could not copy the fastq reads to their bins        |
\---------------------------------------------------------------------*/
uint8_t extractBestRead(
    struct readBin *binIn        /*Bin to extract best read from*/
); /*Splits bin into read with highest Q-score & other all other reads*/

/*----------------------------------------------------------------------
# Output:
#    Prints: fastq entry to outFILE if finds
#    Returns:
#        - 1: if found and printed the id
#        - 2: If could not find the id (no printing)
#        - 4: If the fqFILE does not exist
#        - 8: If the outFILE does not exist
#        - 16: If an incomplete entry (EOF, but missing lines)
#        - 32: If the idCStr (id looking for) entry is incomplete
----------------------------------------------------------------------*/
uint8_t fqOneIdExtract(
    char *idCStr,  /*'\0' terminated read id to extract*/
    FILE *fqFILE,     /*fastq file to extract read from*/
    FILE *keptFILE,   /*File with the target read*/
    FILE *outFILE     /*fastq file to write read to*/
); /*Extracts one read id from a fastq file*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|     - 1: if succeeded or no reads mapped (numReadsKept=0 for no reads)
|     - 4: if could not read reference
|     - 8: if could not open the fastq file
|     - 16: if minimap2 errored out or returned nothing
| Note:
|     - Score: percMult * (keptSNPs + keptIns + keptDels) / read length
|     - minSimUSht ranges from 1 (0.01%) to precMult (100%)            
\---------------------------------------------------------------------*/
uint8_t findBestXReads(
    const uint64_t *numReadConsULng, /*# reads for bulding a consensus*/
    uint64_t *numReadsKeptULng,  /*Number of reads binned to con*/
    char *threadsCStr,           /*Number threads to use with minimap2*/
    const char *useMapqBl,       /*1: use mapping quality in selection*/
    struct minAlnStats *minStats,/*Min stats to cluster reads together*/
    struct samEntry *samStruct,  /*Struct to use for reading sam file*/
    struct samEntry *refStruct,  /*holds the reference (0 to ignore)*/
    struct readBin *binTree,     /*Bin working on*/
    char noRefBl,                /*Only use fastq file*/
    char makeNameBl
       /*1: Make a name for the file using the fastq file; 0 use 
         binTree->topReadsCStr ('\0' for stdout)*/
); /*Extract the top reads that mapped to the selected best read*/

/*----------------------------------------------------------------------
| Output:
|    Files Created:
|      o A fastq file with only the best read
|    Files Modified:
|      o Removes the best read from the orignal fastq file
|    Returns:
|      o 1: if found and printed the id
|      o 4: If the fqFILE does not exist
|      o 8: If the outFILE does not exist
|      o 64: For memory allocation error
|        - In this case it does not change the old fastq file
----------------------------------------------------------------------*/
uint8_t fqGetBestReadByMedQ(
    struct readBin *clustOn,    /*Has fastq file to extract from*/
    struct samEntry *samStruct, /*Sam struct to use for extraction*/
    struct samEntry *bestRead   /*Holds best read till end*/
); /*Extracts read with best medain Q-score from file.
    It also considers length if integer Q-scores are the same.*/

#endif
