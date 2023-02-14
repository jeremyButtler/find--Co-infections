/*##############################################################################
# Name: scoreReadsCheckCigarEntry
# Use: Functions to checks cigar entries and incurments pointers for samfile
##############################################################################*/

#ifndef SCOREREADSCHECKCIGARENTRY_H
#define SCOREREADSCHECKCIGARENTRY_H

#include <stdlib.h> /*for strtoul*/
#include<stddef.h> /*for NULL (use with strtoul)*/

#include "scoreReadsStructers.h"
#include "scoreReadsGlobalVar.h" /*Has my global variables values*/

/*check mismatch entries in cigar*/
void checkMismatches(scoreReadStruct *seqReadStruct,
                     minStats *minReadStats,
                     samEntry *samEntryStruct);

/*check matches entries in cigar*/
void checkMatches(scoreReadStruct *seqReadStruct,
                  minStats *minReadStats,
                  samEntry *samEntryStruct);

/*check insertions entries in cigar*/
void checkInsertion(scoreReadStruct *seqReadStruct,
                    minStats *minReadStats,
                    samEntry *samEntryStruct);

/*check deletion entries in cigar*/
void checkDeletions(scoreReadStruct *seqReadStruct,
                    minStats *minReadStats,
                    samEntry *samEntryStruct);

void checkSoftMasks(scoreReadStruct *seqReadStruct);

/*check if bases are the same (for indels checks) in cigar*/
char checkIfBasesMatch(char baseOneChar, char baseTwoChar);

#endif
