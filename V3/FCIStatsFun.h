/*######################################################################
# Use:
#   o Holds the Median Q-score functions for scoreReads & findCoInft
# Includes:
#   - "samEntryStruct.h"
#   o "cStrToNumberFun.h"
#   o "printErrors.h"
# C standard libraries
#   o <stdlib.h>
#   o <stdio.h>
#   o <string.h>
#   o <stdint.h>
######################################################################*/

#ifndef FCISTATSFUN_H
#define FCISTATSFUN_H

#include "samEntryStruct.h"

/*######################################################################
# output:
#    returns: The read length (unsigned long)
# Note: 
#    - samStruct->seqQHistUInt must have all values initialized to 0
#    - Requires: qHistToMedian from scoreReadsSamFileFunctions.c
######################################################################*/
void findQScores(
    struct samEntry *samStruct /*Sam entry to find Q-scores for*/
); /*Finds Q-scores for input sam entry*/

/*######################################################################
# Output:
#    returns: double with the median Q-score [or 0 if nothing]
######################################################################*/
float qHistToMed(
    uint32_t qHistUInt[], /*Histogram of Q-scores*/
    uint32_t readLenUInt     /*Number of bases in the read*/
); /*converts histogram of q-scores into samStruct into median Q-score*/

#endif
