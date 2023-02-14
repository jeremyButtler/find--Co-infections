/*######################################################################
# Name: trimSam
# Use:
#    - trims soft mask regions off all alignments with sequences in a
#      sam file. Aligments without sequences are ignored & not printed
#      out.
# Includes:
#    - "samEntryStruct.h"
#        - <stdlib.h>
#        - "cStrToNumberFun.h"
#            - <stdint.h>
#        - "printError.h"
#            - <stdio.h>
######################################################################*/

#ifndef TRIMSAM_H
#define TRIMSAM_H

#include "samEntryStruct.h"

/*######################################################################
# Output:
#    Prints: Trimmed sam entries with sequences to outFILE, but ignores
#            sam entries without sequences
######################################################################*/
void trimSamReads(
    FILE *samFILE,               /*Sam file to convert*/
    FILE *outFILE                /*File to store output*/
); /*Prints trimmed sam entries with sequences to file*/

/*######################################################################
# Output:
#    Returns:
#        - 0 if suceeded
#        - 2 if header (invalid and ignored)
#        - 4 if an unmapped read (no reference)
#        - 8 if no sequence line
#    Modifies:
#        - Trims cigar, sequence, & q-score entries in samStruct.
######################################################################*/
uint8_t trimSamEntry(
    struct samEntry *samStruct   /*has sam line to trim softmasks*/
); /*Trims soft masked regions at start & end of sam entry*/

#endif
