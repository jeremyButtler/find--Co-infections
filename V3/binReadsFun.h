/*######################################################################
# Use:
#   o Holds functions related to read binning.
# Includes:
#   - "cStrFun.h"
#   - "trimSam.h"
#   - "findCoInftBinTree.h"
#   - "findCoInftChecks.h"
#   o "samEntryStruct.h"
#   o "defaultSettings.h"
#   o "scoreReadsFun.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
# C standard libraries:
#     o <string.h>
#     o <stdlib.h>
#     o <stdio.h>
#     o <stdint.h>
######################################################################*/

#ifndef BINREADSFUN_H
#define BINREADSFUN_H

#include "trimSam.h"          /*Trimming soft masks off sam alignemnts*/
#include "findCoInftChecks.h" /*Checking functions for alignments*/
#include "cStrFun.h"          /*C-string manipuplation*/
#include "findCoInftBinTree.h"/*To build the readBin tree*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o errUC to hold error output
|       - 1 for success
|       - 2 for file error
|       - 4 for unable to create a fastq file for a bin
|       - 8 for unable to create a stats file for a bin
|       - 64 for memory allocation error
|   - Creates:
|     o A fastq file with reads for each bin
|     o A stats file with the stats from scoreReads for each bin
|   - Returns:
|     o Tree of bins that the reads mapped into
\---------------------------------------------------------------------*/
struct readBin * binReads(
    char *fqPathCStr,        /*Fastq file to bin*/
    char *refsPathCStr,      /*References to bin with*/
    char *prefixCStr,
    char *threadsCStr,       /*Numbe of threads to use with minimap2*/
    char rmSupAlnBl,         /*Remove supplementary alignments*/
    char trimBl,             /*1: trim reads, 0: do not*/
    struct samEntry *newSam, /*Holds minimap2 output*/
    struct samEntry *oldSam, /*Holds previous line of minimap2 output*/
    struct minAlnStats *minStats,
    unsigned char *errUC     /*Reports any errors*/
); /*Bin reads with a set of references*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 1: if succeded
|    Modifies:
|        - fastq in binClust->fqPathCStr to be the fastq for the cluster
|        - fastq in binTree->fqPathCStr to not have clustered reads
\---------------------------------------------------------------------*/
uint8_t binReadToCon(
    const uint8_t *clustUChar,      /*Cluster on*/
    struct readBin *binTree,        /*Bin working on*/
    struct readBin *binClust,       /*Bin to assign reads & consensus*/
    struct samEntry *samStruct,     /*To hold temporary input*/
    struct minAlnStats *minStats,   /*Min stats needed to keep a read*/
    char *threadsCStr            /*Number threads to use with Minimap2*/
); /*Maps reads to consensus and keeps reads that meet user criteria*/

#endif
