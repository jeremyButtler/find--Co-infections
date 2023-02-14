/*##############################################################################
# Name: scoreReadsRefCheckCigarEntry
# Use: Checks if bases should be kept for mismatches, matches, & deletions
#      using a reference
##############################################################################*/

#ifndef SCOREREADSREFCHECKCIGARENTRY_H
#define SCOREREADSREFCHECKCIGARENTRY_H

#include <stdlib.h> /*for strtoul*/
#include<stddef.h> /*for NULL (use with strtoul)*/

#include "scoreReadsStructers.h"
#include "scoreReadsGlobalVar.h" /*Has my global variables values*/
#include "scoreReadsCheckCigarEntry.h" /*has checkIfBasesMatch*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept mismatches, ignored mismatches
##############################################################################*/
void refCheckMismatches(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep mismatch*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
); /*Checks if mismatch is valid or not & updates histogram if valid*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept matches, total matches
##############################################################################*/
void refCheckMatches(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep mismatch*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
); /*Counts total matches, checks if matches are kept, & incurments pointers*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept mismatches, ignored insertions
# Note:
#    This function acts a little oddly, in that it only keeps a deletion if 
#      the reference has to little support. This is tricky, since I have no
#      Q-score for the read here
##############################################################################*/
void refCheckDeletions(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep insertion*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
); /*Uses refernce to checks if deletion is supported*/

#endif
