#ifndef SCOREREADSSCORINGWRAPPERS_H
#define SCOREREADSSCORINGWRAPPERS_H

#include <stdlib.h>
#include <stdio.h>

#include "scoreReadsStructers.h" /*Holds my structers for samfile data storage*/
#include "scoreReadsFileIO.h" /*Input and output functions for this code*/
#include "scoreReadsSamFileFunctions.h" /*samfile functions to score reads*/
#include "scoreReadsGlobalVar.h" /*Has my global variables values*/

void scoreReadsFile(FILE *samFile,
                    unsigned long maxLineLenULng,
                    minStats *minReadStats
); /*Get read scores using input from a file*/

/*Get read scores using input from stdin*/
void scoreReadsStdin(unsigned long maxLineLenULng, minStats *minReadStats);

/*##############################################################################
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Frees: seqCStr & qCStr in refEntry struct
#    Returns: 4: If reference beneath min quality
#             2: if sucess
#             0: If malloc failed
##############################################################################*/
char refScoreReadsStdin(
    unsigned long maxLineLenULng,   /*Max line length to read in sam file*/
    struct minStats *minReadStats,  /*Min thresholds to keep an aligment*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    struct samEntry *refEntry       /*Refence to compare read to*/
); /*Score each alignment using the provided reference & cigar*/

/*##############################################################################
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Frees: seqCStr & qCStr in refEntry struct
#    Closes: samFile
#    Returns: 4: If reference beneath min quality
#             2: if sucess
#             0: If malloc failed
##############################################################################*/
char refScoreReadsFile(
    FILE *samFile,                  /*Sam file with aligments to score*/
    unsigned long maxLineLenULng,   /*Max line length to read in sam file*/
    struct minStats *minReadStats,  /*Min thresholds to keep an aligment*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    struct samEntry *refEntry       /*Refence to compare read to*/
); /*Score each alignment using the provided reference & cigar*/

void initForScoring(samEntry *samStruct, 
                    char **samLineCStr,
                    unsigned long maxLineLenULng
); /*Prepare and intalize scoreRead* variables*/

void scoreAligment(char *samLineCStr,     /*line with sam entry (alignment)*/
                  samEntry *samStruct,    /*will have index of samLineCStr*/
                  char *oldSamLineCStr,   /*line with last sam entry*/
                  samEntry *oldSamStruct, /*holds index of last sam entry*/
                  minStats *minReadStats, /*min stats to keep an alignment*/
                  char *printHeadChar     /*1: print header, 0 do not*/
);

/*##############################################################################
# Output: 
#    stdout: if alignment in sam file is kept, prints aligment stats to stdout
#    Modifies: When move to new set of alignments, modifies oldSamEntry and 
#        oldSamLineCStr to have the new entry
##############################################################################*/
char refScoreAligment(
    struct samEntry *refEntry,     /*Holds refence sequence & Q-score lines*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    char *samLineCStr,             /*line with sam entry (alignment)*/
    struct samEntry *samStruct,    /*will have index of samLineCStr*/
    char *oldSamLineCStr,          /*line with last sam entry*/
    struct samEntry *oldSamStruct, /*holds index of last sam entry*/
    struct minStats *minReadStats, /*min stats to keep an alignment*/
    char *printHeadChar     /*1: print header, 0 do not*/
); /*Scores alignments & after scoring determins if should keep alignment*/

/*##############################################################################
# Output:
#    Modifes: refEntry: to have read length and Q-scores of reference
#    Returns: 4: If reference is beneath min stats
#             2: If reference is ok
#             0: If reference does not have Q-score line to check
##############################################################################*/
char checkRefEntry(
    struct samEntry *refEntry,
    struct minStats *minReadStats
); /*Checks if reference meets the min requirments provided by the user*/


#endif
