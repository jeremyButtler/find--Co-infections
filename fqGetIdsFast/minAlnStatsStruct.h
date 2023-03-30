/*######################################################################
# Use:
#   o Holds the minAlnStats structer and its supporting functions.
# Includes:
#   - defaultSettings.h
# C standard libraries
#   - <stdint.h>
######################################################################*/

#ifndef MINALNSTATSSTRUCT_H
#define MINALNSTATSSTRUCT_H

#include "defaultSettings.h" /*Has devault settings for minAlnStats*/
#include <stdint.h>

/*######################################################################
# Name: minAlnStats
# Use: Stores the min stats to keep an read, mismatch, or indel
######################################################################*/
typedef struct minAlnStats
{ /*minAlnStats structer*/
    uint8_t minQChar; /*Min Q-score to keep a mismatch or indel*/
    uint32_t minMapqUInt;     /*default min mapping quality*/
    
    float minMedianQFlt;         /*default median Q-score*/
    float minMeanQFlt;           /*default mean Q-score*/
    float minAlignedMedianQFlt;  /*default median aligend Q*/
    float minAlignedMeanQFlt;    /*default mean aligend Q-score*/

    uint32_t
        maxReadLenULng,    /*default max read length (entire read)*/
        minReadLenULng,    /*min read length (looks at aligned)*/
        maxHomoInsAry[16], /*Array of limits for deletions*/
        maxHomoDelAry[16]; /*Array of limits for deletions*/
        /*A->0, C->1, G->3, T/U->10;
          find with: (baseChar & ~(1+32+64+138)) >> 1*/

     /*These are not used in scoreReads, but are used in findCoInft*/
     float minSNPsFlt;
     float minDelsFlt;
     float minInssFlt;
     float minIndelsFlt;
     float minDiffFlt;
}minAlnStats; /*minAlnStats structer*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
######################################################################*/
void blankMinStats(
    struct minAlnStats *minStats
); /*Sets minStats minimum requirements for sam alingemtns to defaults*/

/*THESE FUNCTIONS ARE HERE FOR findCoInft AND ARE NOT NEEDED FOR 
  scoreReads*/
/*######################################################################
# Output:
#    modifes minStats to have default entries
# Uses default read to read mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsReadRead(
    struct minAlnStats *minStats
); /*Sets minStats minimum requirements for sam alingemtns to defaults*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# Uses default read to consensus mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsReadCon(
    struct minAlnStats *minStats
); /*Sets minStats minimum requirements for sam alingemtns to defaults*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# default consensus to consensus mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsConCon(
    struct minAlnStats *minStats
); /*Sets minStats minimum requirements for sam alingemtns to defaults*/

#endif
