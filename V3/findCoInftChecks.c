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

#include "findCoInftChecks.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: findCoInftChecks
'   o fun-1 checkRead:
'     - Checks if the sam entry meets the min user requirements
'   o fun-2 checkIfKeepRead:
'     - Checks if read does not meet one threshold in maxDifference
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        -0: If does not meet min stats
|        -1: If meets min stats
\---------------------------------------------------------------------*/
uint8_t checkRead(
    struct minAlnStats *minStats, /*Structer holding min requirments*/
    struct samEntry *samStruct    /*Structer with the alignemtn stats*/
) /*Checks if the sam entry meets the min user requirements*/
{ /*checkRead*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-4 TOC: Fun-4 Sec-1 Sub-1: checRead
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(samStruct->mapqUChar < minStats->minMapqUInt)
       return 0; /*move onto the next entry (has to low mapq)*/
           /*Minimap2 only gives MAPQ for best match*/

   if(
      samStruct->medianQFlt < minStats->minMedianQFlt &&
      *samStruct->qCStr !='*'
   )
       return 0; /*move onto the next entry (median Q-score to low)*/

   if(
      samStruct->meanQFlt < minStats->minMeanQFlt &&
      *samStruct->qCStr !='*'
   )
       return 0; /*move onto the next entry (mean Q-score to low)*/

   if(
       samStruct->readLenUInt > minStats->maxReadLenULng &&
       minStats->maxReadLenULng > 0
   )
            return 0; /*move onto the next entry (read is to long)*/

    if(
        samStruct->medianAligQFlt < minStats->minAlignedMedianQFlt &&
        *samStruct->qCStr !='*'
    )
        return 0; /*aligned median Q-score low*/

    if(
        samStruct->meanAligQFlt < minStats->minAlignedMeanQFlt &&
        *samStruct->qCStr != '*'
    )
        return 0; /*aligned mean Q-score to low*/

    if(samStruct->readAligLenUInt < minStats->minReadLenULng)
        return 0; /*Read to to short*/

    return 1;
} /*checkRead*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 0: Read does not meet minimum thresholds in maxDifference
|        - 1: Keep the read (meets minimum thresholds)
\---------------------------------------------------------------------*/
uint8_t checkIfKeepRead(
    struct minAlnStats *maxDifference, /*Has thresholds to keep reads*/
    struct samEntry *samStruct     /*Has stats for read*/
) /*Checks if read does not meet one threshold in maxDifference*/
{ /*checkIfKeepRead*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-16 TOC: Sec-1 Sub-1: checkIfKeepRead
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    float
        percSnpFlt = 0,
        percInsFlt = 0,
        percDelFlt = 0;

    /*This might be better if I found these separately and had if 
          statements between each calcuation, however this is clearer
      Need * 1.0 to ensure it does not truncate*/
    percSnpFlt =
        samStruct->numKeptSNPUInt / (samStruct->readLenUInt * 1.0);

    percInsFlt =
        samStruct->numKeptInsUInt / (samStruct->readLenUInt * 1.0);

    percDelFlt =
        samStruct->numKeptDelUInt / (samStruct->readLenUInt * 1.0);

    if(maxDifference->minSNPsFlt < percSnpFlt) /*percent SNPs to high*/
        return 0;
    /*The sequences in the sam alignment are to different*/
    if(maxDifference->minDiffFlt < percSnpFlt + percDelFlt + percInsFlt)
        return 0;
    /*To many indels is to hight*/
    if(maxDifference->minIndelsFlt < percInsFlt + percDelFlt)
        return 0;
    if(maxDifference->minInssFlt < percInsFlt) /*To many insertions*/
        return 0;
    if(maxDifference->minDelsFlt < percDelFlt) /*To many deletions*/
        return 0;

    return 1;
} /*checkIfKeepRead*/
