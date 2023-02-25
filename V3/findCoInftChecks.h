/*######################################################################
# Use:
#   o Holds functions and structers for checking read quality in
#     find co-infections
# Includes:
#   o "defaultSettings.h"
#   o "scoreReadsFun.h"
#      - "samEntryStruct.h"
#        o <stdlib.h>
#        o "cStrToNumberFun.h"
#          - <sdtint.h>
#        o "printError.h"
#          - <stdio.h>
######################################################################*/

#ifndef FINDCOINFTCHECKS_H
#define FINDCOINFTCHECKS_H

#include "defaultSettings.h"   /*Has minimap2 command*/
#include "scoreReadsFun.h"   /*Has minimap2 command*/

/*---------------------------------------------------------------------\
| Struct-1: minDiff
| Use:
|    - Holds the min difference (SNPs, indels, or % sim) needed to 
|      consdier a read the same
\---------------------------------------------------------------------*/
typedef struct minDiff
{ /*minDiff*/
    float
        minSNPsFlt,
        minDelsFlt,
        minInssFlt,
        minIndelsFlt,
        minDiffFlt;
}minDiff;

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        -0: If does not meet min stats
|        -1: If meets min stats
\---------------------------------------------------------------------*/
uint8_t checkRead(
    struct minAlnStats *minStats, /*Structer holding min requirments*/
    struct samEntry *samStruct    /*Structer with the alignemtn stats*/
); /*Checks if the sam entry meets the min user requirements*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 0: Read does not meet minimum thresholds in maxDifference
|        - 1: Keep the read (meets minimum thresholds)
\---------------------------------------------------------------------*/
uint8_t checkIfKeepRead(
    struct minAlnStats *maxDifference, /*Has thresholds to keep reads*/
    struct samEntry *samStruct     /*Has stats for read*/
); /*Checks if read does not meet one threshold in maxDifference*/

#endif
