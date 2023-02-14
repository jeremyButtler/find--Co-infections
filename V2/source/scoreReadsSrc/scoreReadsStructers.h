/*##############################################################################
# Name: samEntry
# Use: Stores stats about a sam file entry
##############################################################################*/

#ifndef SCOREREADSSTRUCTERS_H
#define SCOREREADSSTRUCTERS_H

#include "scoreReadsGlobalVar.h" /*Has my global variables values*/

typedef struct samEntry
{ /*samEntry*/
    char *queryCStr,        /*points to the query name in the sam entry tab 0*/
         *refCStr,          /*points ot the reference in the sam entry  tab 1*/
         *cigarCStr,        /*points to the cigar in the sam entry      tab 6*/
         *seqCStr,          /*points to the seqeuence in the sam entry  tab 10*/
         *qCStr;            /*points to the q-score in the same entry   tab 11*/
    unsigned int flagUInt,  /*Holds the flag in the sam entry*/
                 mapqUInt;  /*Hols the mapping quality in the same entry*/
    double medianQDbl,      /*holds the median read Q-score*/
           medianAligQDbl,  /*holds the aligned read median Q-score*/
           meanQDbl,        /*holds the mean read Q-score*/
           meanAligQDbl;    /*holds the mean aligned read Q-score (no low Q)*/
    unsigned long
        readLenULng,        /*Holds the read length of sam entry*/
        readAligLenULng,    /*Aligned read length of sam entry*/
        numMatchULng,       /*Holds number of matches*/
        numKeptMatchULng,   /*Holds number of matches with Q-score > min Q*/
        numMisULng,         /*number of mismatches in sam entry*/
        numIgnoreMisULng,   /*number ingored mismatches in sam entry*/
        numDelULng,         /*number of indels in sam entry*/
        numIgnoreDelULng,   /*number of ignored indels in sam entry*/
        numInsULng,         /*number of indels in sam entry*/
        numIgnoreInsULng;   /*number of ignored indels in sam entry*/
}samEntry; /*samFileEntry*/

/*##############################################################################
# Name: minStats
# Use: Stores the min stats to keep an read, mismatch, or indel
##############################################################################*/
typedef struct minStats
{ /*minStats structer*/
    char minQChar; /*Min Q-score to keep a mismatch or indel*/

    unsigned int minMapqUInt;     /*default min mapping quality*/
    
    double minMedianQDbl,         /*default median Q-score*/
           minMeanQDbl,           /*default mean Q-score*/
           minAlignedMedianQDbl,  /*default median aligend Q*/
           minAlignedMeanQDbl;    /*default mean aligend Q-score*/

    unsigned long maxReadLenULng,    /*default max read length (entire read)*/
                  minReadLenULng,    /*min read length (looks at aligned)*/
                  maxInsAHomoULng,   /*max A insertion homopolymer size*/ 
                  maxInsTHomoULng,   /*max T insertion homopolymer size*/ 
                  maxInsCHomoULng,   /*max C insertion homopolymer size*/ 
                  maxInsGHomoULng,   /*max C insertion homopolymer size*/ 
                  maxDelAHomoULng,   /*max A deletion homopolymer size*/ 
                  maxDelTHomoULng,   /*max T deletion homopolymer size*/ 
                  maxDelCHomoULng,   /*max C deletion homopolymer size*/ 
                  maxDelGHomoULng;   /*max C deletion homopolymer size*/ 
}minStats; /*minStats structer*/


/*##############################################################################
# Name: scoreReadStruct
# Use: Used in scoreReads functions to store the small bits of information, 
#      such as sequence pionters, Q-score pointers, total bases, ect...
##############################################################################*/
typedef struct scoreReadStruct
{ /*scoreReadStruct structure*/
    char *sequenceCStr,
         *qScoreCStr;

    unsigned long cigarEntryULng,
                  totalQScoreULng,
                  seqQHistULng[94];
}scoreReadStruct;

/*Reset all data in the sam file entry to 0*/
void blankSamFileEntry(samEntry *samEntryStruct);

/*Only reset non-pointer variables in the sam file entry to 0*/
void blankSamFileEntryVariables(samEntry *samEntryStruct);

/*reset the variables in the minStats sturcter to defaults*/
void blankMinStats(minStats *minReadStats);

/*Makes copy of samEntry and c-string with sam entry*/
void deepCpSamFileEntry(samEntry *samEntryStruct,
                        samEntry *cpSamEntryStruct,
                        char *entryCStr,
                        char *cpEntryCStr);

/*Makes copy of samEntry*/
void cpSamFileEntry(samEntry *samEntryStruct,
                    samEntry *cpSamEntryStruct);

/*Intalize a scoreReadStruct*/

void initScoreReadStruct(scoreReadStruct *structIn, samEntry *samEntryStruct);
#endif
