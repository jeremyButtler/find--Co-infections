/*##############################################################################
# TOC:
#    void blankSamFileEntry: sets variables in sameFileEntry struct to 0
#        fun-1 sec-1: No sections
#    void blankSamFileEntryVariables: Sets all non-pointers to 0
#        fun-2 sec-1: No sections
#    void blankMinStats: sets variables in minStats structure to default
#        fun-3 sec-1: No sections
#    void deepCpSamFileEntry: copy same file entry (struct and c-string)
#        fun-4 sec-1: No sectoins (Adjust Q-scores)
#    void cpSamFileEntry: Copy variblles from one samEntry structer to another
#        fun-5 sec-1: No sections
#    void makeSamEntryPosStructAry
#        fun-6 sec-1: No sections
#    void initSamEntryPosStruct: Intailize a samEntryPos structer
#        fun-7 sec-1: No sections
#    void freeSamEntryPosStruct: Frees alocated memory in samEntryPos
#        fun-8 sec-1: No sections
#    void delSamEntryPosStructAry: deletes a samEntryPos array
#        fun-9 sec-1: no sections
#    fun-10: makeSamEntryReadSetAry: Make array of samEntryReadSet structures
#    fun-11: initSamEntryReadSetStruct: Set variables in samEntryReadSet
#                                       structure to 0
#    fun-12: freeSamEntryReadSetStruct: free variables in samEntryReadSet struct
#    fun-13: delSamEntryReadSetAry: Free samEntryReadSet structer array
##############################################################################*/

#include "scoreReadsStructers.h"

/*##############################################################################
# Name: blankSamFileEntry
# Use: Sets all variables in samEntryStruct to 0
# Input: samEntryStruct: samEntry struct to blank
# Output: modifes samEntryStruct to have every entry = 0
##############################################################################*/
void blankSamFileEntry(samEntry *samEntryStruct)
{ /*blankSamFileEntry*/
    /*Fun-1 Sec-1*/
    (*samEntryStruct).queryCStr = 0;
    (*samEntryStruct).refCStr = 0;
    (*samEntryStruct).cigarCStr = 0;
    (*samEntryStruct).seqCStr = 0;
    (*samEntryStruct).qCStr = 0; 

    (*samEntryStruct).flagUInt = 0;
    (*samEntryStruct).mapqUInt = 0;

    (*samEntryStruct).medianQDbl = 0;
    (*samEntryStruct).medianAligQDbl = 0;
    (*samEntryStruct).meanQDbl = 0;
    (*samEntryStruct).meanAligQDbl = 0;

    (*samEntryStruct).readLenULng = 0;
    (*samEntryStruct).readAligLenULng = 0;

    (*samEntryStruct).numMatchULng = 0;
    (*samEntryStruct).numKeptMatchULng = 0;
    (*samEntryStruct).numMisULng = 0;
    (*samEntryStruct).numIgnoreMisULng = 0;
    (*samEntryStruct).numInsULng = 0;
    (*samEntryStruct).numIgnoreInsULng = 0;
    (*samEntryStruct).numDelULng = 0;
    (*samEntryStruct).numIgnoreDelULng = 0;

    return;
} /*blankSamFileEntry*/

/*##############################################################################
# Name: blankSamFileEntryVariables
# Use: Sets all non pointer variables in samEntryStruct to 0
# Input: samEntryStruct: samEntry struct to blank
# Output: modifes samEntryStruct to have every non-pointer = 0
##############################################################################*/
void blankSamFileEntryVariables(samEntry *samEntryStruct)
{ /*blankSamFileEntry*/
    /*Fun-2 Sec-1*/
    (*samEntryStruct).flagUInt = 0;
    (*samEntryStruct).mapqUInt = 0;

    (*samEntryStruct).medianQDbl = 0;
    (*samEntryStruct).medianAligQDbl = 0;
    (*samEntryStruct).meanQDbl = 0;
    (*samEntryStruct).meanAligQDbl = 0;

    (*samEntryStruct).readLenULng = 0;
    (*samEntryStruct).readAligLenULng = 0;

    (*samEntryStruct).numMatchULng = 0;
    (*samEntryStruct).numKeptMatchULng = 0;
    (*samEntryStruct).numMisULng = 0;
    (*samEntryStruct).numIgnoreMisULng = 0;
    (*samEntryStruct).numInsULng = 0;
    (*samEntryStruct).numIgnoreInsULng = 0;
    (*samEntryStruct).numDelULng = 0;
    (*samEntryStruct).numIgnoreDelULng = 0;

    return;
} /*blankSamFileEntry*/


/*##############################################################################
# Name: blankMinStats
# Use: Sets all variables in minStats to defaults
# Input: minReadStats: minStats struct to blank
# Output: modifes minReadStats to have default entries
##############################################################################*/
void blankMinStats(minStats *minReadStats)
{ /*blankMinStats*/
    /*Fun-3 Sec-1*/
    (*minReadStats).minMapqUInt = 20;          /*default min mapping quality*/
    (*minReadStats).minQChar = 13;             /*default min Q-score*/
    (*minReadStats).minMedianQDbl = 13;        /*default median Q-score*/
    (*minReadStats).minMeanQDbl = 13;          /*default mean Q-score*/
    (*minReadStats).minAlignedMedianQDbl = 13; /*default median aligend Q*/
    (*minReadStats).minAlignedMeanQDbl = 13;   /*default mean aligend Q-score*/
    (*minReadStats).maxReadLenULng = 1000;     /*default max read length*/
    (*minReadStats).minReadLenULng = 600;      /*default max read length*/
    
    (*minReadStats).maxInsAHomoULng = 2;    /*max A insertion homopolymer size*/ 
    (*minReadStats).maxInsTHomoULng = 2;    /*max T insertion homopolymer size*/ 
    (*minReadStats).maxInsCHomoULng = 1;    /*max C insertion homopolymer size*/ 
    (*minReadStats).maxInsGHomoULng = 1;    /*max G insertion homopolymer size*/ 

    (*minReadStats).maxDelAHomoULng = 0;    /*max A deletion homopolymer size*/ 
    (*minReadStats).maxDelTHomoULng = 0;    /*max T deletion homopolymer size*/ 
    (*minReadStats).maxDelCHomoULng = 0;    /*max C deletion homopolymer size*/ 
    (*minReadStats).maxDelGHomoULng = 0;    /*max G deletion homopolymer size*/ 
        /*Setting default to 0, because have not quality scores to looks at*/

    return;
} /*blankMinStats*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SamEntryPos structer functions
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Name: deepCpSamFileEntry
# Use: Copies one sam file entry struct into another and the c-string with
#      the sam file entry into another provided c-string
# Input:
#    samEntryStruct: pionter to samEntry to copy
#    cpSamEntryStruct: pointer to samEntry to copy stats to
#    entryCStr: C-string with the same file line (entry) to copy
#    cpEntryCStr: C-string to copy the sam file line (entry) into
#        - entryCStr needs to be the same size as oldEntryCStr
# Output: 
#    Modifies: meanQDbl, medianQDbl, mapqUInt, and readLenULng variables 
#              of samEntryStruct to cpSamEntryStruct values
#    Modifies: cpEntryCStr to be a copy of entryCStr
#    Modifies: seqCStr and qCStr in oldSamEntry to piont to the sequence and
#              Q-score enteries in cpEntryCStr
##############################################################################*/
void deepCpSamFileEntry(samEntry *samEntryStruct,
                        samEntry *cpSamEntryStruct,
                        char *entryCStr,
                        char *cpEntryCStr)
{ /*cpSamFileEntry*/
    /*Fun-4 Sec-1 (Copy stats)*/

    char *lineCStr = entryCStr, *cpLineCStr = cpEntryCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Copy the sam file line (entry) and stats for entry
    #    fun-4 sec-2 sub-1: Copy c-string with sam file line to new buffer
    #    fun-4 sec-2 sub-2: Copy shared query stats in old entry to new
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-4 Sec-2 Sub-1: Copy c-string with sam file line to new buffer
    ***************************************************************************/

    while(*lineCStr != '\0') /*So do not overwrite the sequence entry*/
    { /*Loop and copy the samfile line*/
        *cpLineCStr = *lineCStr; /*On a new sequence line*/
        cpLineCStr++;
        lineCStr++;
    } /*Loop and copy the samfile line*/

    *(cpLineCStr + 1) = '\0'; /*Make into a c-string*/

    /***************************************************************************
    # Fun-4 Sec-2 Sub-2: Copy shared query stats in old entry to new
    ***************************************************************************/

    /*Set pointers to sequence and Q-score line in copied sam entry*/
    (*cpSamEntryStruct).seqCStr = cpEntryCStr +
                                  ((*samEntryStruct).seqCStr -
                                   (*samEntryStruct).queryCStr);
    (*cpSamEntryStruct).qCStr = cpEntryCStr +
                                ((*samEntryStruct).qCStr -
                                 (*samEntryStruct).queryCStr);
    (*cpSamEntryStruct).mapqUInt = (*samEntryStruct).mapqUInt;

    /*Record shared stats between same querys*/
    (*cpSamEntryStruct).readLenULng = (*samEntryStruct).readLenULng;
    (*cpSamEntryStruct).meanQDbl = (*samEntryStruct).meanQDbl;
    (*cpSamEntryStruct).medianQDbl = (*samEntryStruct).medianQDbl;

    return;
} /*cpSamFileEntry*/

/*##############################################################################
# Name: cpSamFileEntry
# Use: Copies shared query variables of one sam file entry struct into another
# Input:
#    samEntryStruct: pionter to samEntry to copy
#    cpSamEntryStruct: pointer to samEntry to copy stats to
# Output: 
#    Modifies: meanQDbl, medianQDbl, readLenULng, seqCStr, and qCStr variablesi
#              of samEntryStruct to cpSamEntryStruct values
##############################################################################*/
void cpSamFileEntry(samEntry *samEntryStruct,
                    samEntry *cpSamEntryStruct)
{ /*cpSamFileEntry*/
    /*Fun-5 Sec-1 (Copy stats)*/

    /*Set pointers to sequence and Q-score line in copied sam entry*/
    (*cpSamEntryStruct).seqCStr = (*samEntryStruct).seqCStr;
    (*cpSamEntryStruct).qCStr = (*samEntryStruct).qCStr;

    /*Record shared stats between same querys*/
    (*cpSamEntryStruct).readLenULng = (*samEntryStruct).readLenULng;
    (*cpSamEntryStruct).meanQDbl = (*samEntryStruct).meanQDbl;
    (*cpSamEntryStruct).medianQDbl = (*samEntryStruct).medianQDbl;

    return;
} /*cpSamFileEntry*/

/*##############################################################################
# Name: initScoreReadStruct
# Use: Intalizes a score read structer
# Input:
#    structIn: scoreReadStruct to intialize (scoreReadStruct pointer)
#    samEntryStruct: samEntry sturcter with sequence data and q-score to grab
##############################################################################*/

void initScoreReadStruct(scoreReadStruct *structIn, samEntry *samEntryStruct)
{ /*initScoreReadsStruct*/
    (*structIn).sequenceCStr = (*samEntryStruct).seqCStr; 
    (*structIn).qScoreCStr = (*samEntryStruct).qCStr; 

    (*structIn).cigarEntryULng = 0;
    (*structIn).totalQScoreULng = 0;

    /*Intalize the Q-score array*/
    for(int intQ = 0; intQ < MAX_Q_SCORE_CHAR; intQ++)
        (*structIn).seqQHistULng[intQ] = 0;
} /*initScoreReadsStruct*/
