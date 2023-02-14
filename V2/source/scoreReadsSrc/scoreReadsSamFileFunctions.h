/*##############################################################################
# Name: ScoreReadsSamFileFunctions
# Use: Has functions to score reads, including median Q-score functions
##############################################################################*/

#ifndef SCOREREADSSAMFILEFUNCTIONS_H
#define SCOREREADSSAMFILEFUNCTIONS_H

#include "scoreReadsStructers.h" /*Holds my structers for samfile data storage*/
#include "scoreReadsConversion.h" /*conversion functions*/
#include "scoreReadsCheckCigarEntry.h" /*functions to check cigar entries*/
#include "scoreReadsRefCheckCigarEntry.h" /*functions to check cigar entries*/
#include "scoreReadsGlobalVar.h" /*Has my global variables values*/

/*Scores reads*/
void scoreRead(minStats *minReadStats, samEntry *samEntryStruct);

/*Index the entries in a single sam file entry*/
char indexSamEntry(char *entryCStr, samEntry *samEntryStruct);

/*##############################################################################
# Output: Modifies: samEntryStruct to hold aligned length, Q-scores,
#                   mismatches, & indels
##############################################################################*/
void refScoreRead(
    struct minStats *minReadStats,   /*Min requirments to keep aligment*/
    struct samEntry *samEntryStruct, /*sequence & Q-score line & holds stats*/
    struct samEntry *refEntry,       /*sequence & Q-score line for reference*/
    char refForDelChar              /*Marks if only using ref for deltions*/
); /*Find similarity between two reads, aligned read length, & aligned Q-scores*/

/*get median values from a Q-score histogram*/
double qHistToMedian(unsigned long qHistULng[],
                     unsigned long readLenULng);

/*Finds the mean and median Q-scores (calls qHistToMedian)*/
unsigned long findQScores(char *qCStr, double *meanQDbl, double *medianQDbl);

/*##############################################################################
# Output:
#    returns: unsigned long with number of bases
#    modifies: cigarCStr to point to the entry type (mismatch, indel, ect...)
##############################################################################*/
unsigned long readCigEntry(
    char **cigarCStr,       /*C-string with cigar to read and incurment*/
    unsigned int flagUInt   /*Flag telling if is reverse sequence (16)*/
); /*Reads a single entry from a eqx cigar line*/

#endif
