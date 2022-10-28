#ifndef SCOREREADSFILEIO_H
#define SCOREREADSFILEIO_H

#include <stdio.h> /*file IO library*/
#include <stdlib.h> /*malloc*/
#include "scoreReadsStructers.h" /*Holds my structers for samfile data storage*/

/*Open sam file and check length of the longest line*/
FILE * openFile(char *fileInCStr,
                unsigned long *maxLineLenULng,
                int buffSizeInt
);

/*Output samfile entry stats to stdout*/
void stdoutReadScore(samEntry *samEntryStruct, char *printHeaderChar);

/*
 Output:
    Allocates: seqCStr to have the first sequence from the fastq file
        - Set to 0 & freed if no sequence line
    Allocates: qCStr to have the first sequences Q-score line from the fastq
        - Set to 0 & freed if no q-score line
    Returns: 1: if no malloc errors
             0: If malloc failed to find memory
*/
char readFirstSeqInFastq(
    char *fileInCStr, /*C-string with path & name of file to open*/
    char **seqCStr,   /*points to c-string with sequence, c-string is on heap*/
    char **qCStr,     /*points to c-string with q-score line, c-string on heap*/
    int buffSizeInt   /*Size of buffer to read file in with*/
); /*Gets the frist reads sequence & q-score line from a fastq file*/

#endif
