/*######################################################################
# Name: scoreReadsFun
# Use:
#  - Holds unique functions needed for scoreReads.
# Includes:
#   - "minAlnStatsStruct.h"
#   - "fqAndFqFun.h"
#   - "FCIStatsFun.h"
#   o "defaultSettings.h"
#   o "samEntryStruct.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
# C standard libraries
#   o <stdlib.h>
#   o <sdtint.h>
#   o <stdio.h>
######################################################################*/

#ifndef SCOREREADSFUN_H
#define SCOREREADSFUN_H

#include "minAlnStatsStruct.h" /*For checking thresholds*/
#include "FCIStatsFun.h"       /*For median Q-scores*/
#include "fqAndFaFun.h"

/*######################################################################
# Output:
#    File: Prints the scores of the reads to the input file
#    Returns:
#        -1 If their were no problems
#        -64 If memory allocation failed
######################################################################*/
uint8_t scoreReads(
    struct minAlnStats *minStats, /*Min stats to keep an alignment*/
    const uint8_t *useRefForDelBool,/*use reference only for deletions*/
    FILE *samFILE,                /*Sam file with alignments to score*/
    FILE *refFILE,                /*reference to use in scoring*/
    FILE *outFILE                 /*File to ouput kept alignments to*/
); /*Scores all alignments in a sam file*/

/*######################################################################
# Output: 
#    Modifies: 
#        - stat variabls in samStruct to have alignment scores
#        - Swaps samSruct & oldSamStruct address when starting new read
#    Returns:
#        -1: If alignment meets min stats
#        -2: If is a header
#        -4: If alignment does not meer the min stats
#        -8: If is not mapped to any reference
######################################################################*/
uint8_t setUpScoreAln(
    struct samEntry **samStruct,   /*will have index of samLineCStr*/
    struct samEntry **oldSamStruct,/*holds index of last sam entry*/
    struct samEntry *refStruct,    /*Holds reference sequence &q-score*/
    struct minAlnStats *minStats,  /*min stats to keep an alignment*/
    const uint8_t *useRefForDelBool,/*use reference only for deletions*/
    const uint8_t *refQBool        /*0 no q-score entry; 1 present*/
); /*scores & finds if should keep an sam aligment*/

/*######################################################################
# Name: scoreAln
# Use: Finds the similarity between two reads, aligned read length.
# Output:
#    modifies: samEntry struct to have kept matches, SNPs, & indels
######################################################################*/
void scoreAln(
    struct minAlnStats *minStats,
    struct samEntry *samStruct,
    struct samEntry *refStruct,
    const uint8_t *refQBool,
    const uint8_t *useRefForDelBool /*use reference only for deletions*/
); /*Scores a signle alignment in a sam file*/

/*######################################################################
# Output:
#    modifies: samStruct Q-score histogram, running Q-score total
#    incurments: samStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
######################################################################*/
void checkSNPs(
    const uint32_t *cigEntryUInt,  /*Cigar entry with number of SNPs*/
    struct minAlnStats *minStats, /*Has min q-score to keep SNP*/
    struct samEntry *samStruct,    /*Holds stats & sam entry*/
    struct samEntry *refStruct,    /*Holds reference sequence &q-score*/
    const uint8_t *qLineBool,      /*0 no q-score entry; 1 present*/
    const uint8_t *refQBool,       /*0 no ref q-score entry; 1 present*/
    const uint8_t *useRefBool,     /*0 do not use reference, 1 use ref*/
    const int32_t *incInt          /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
); /*Checks if should keep SNPs or discard*/

/*######################################################################
# Output:
#    modifies: samStruct Q-score histogram, running Q-score total
#    incurments: samStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
######################################################################*/
void checkMatches(
    const uint32_t *cigEntryUInt, /*Cigar entry with number of matches*/
    struct minAlnStats *minStats, /*Has min q-score to keep matches*/
    struct samEntry *samStruct,    /*Holds stats & sam entry*/
    struct samEntry *refStruct,    /*Holds reference sequence &q-score*/
    const uint8_t *qLineBool,      /*0 no q-score entry; 1 present*/
    const uint8_t *refQBool,       /*0 no ref q-score entry; 1 present*/
    const uint8_t *useRefBool,     /*0 do not use reference, 1 use ref*/
    const int32_t *incInt          /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
); /*Checks if should keep match or not*/

/*######################################################################
# Output:
#    modifies: samStruct Q-score histogram, running Q-score total
#    incurments: samStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
# Note:
#    - No reference provided, because a insertion means a deletion in
#      the reference (so reference provides no additional information)
######################################################################*/
void checkInss(
    const uint32_t *cigEntryUInt, /*Cigar entry with number of ins's*/
    struct minAlnStats *minStats, /*Has min q-score to keep insertion*/
    struct samEntry *samStruct,    /*Holds stats & sam entry*/
    const uint8_t *qLineBool,      /*0 no q-score entry; 1 present*/
    const uint8_t *useRefForDelBool, /*Handeling deletion like ins*/
    const int32_t *incInt          /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
); /*Check if should keep insertions*/

/*######################################################################
# Output:
#    modifies: samStruct with stats
######################################################################*/
void checkDels(
    const uint32_t *cigEntryUInt,  /*Cigar with number deletions*/
    struct minAlnStats *minStats,  /*Has min q-score to keep deletions*/
    struct samEntry *samStruct     /*Holds stats & sam entry*/
); /*Checks if should keep the deletion*/

/*######################################################################
# Output:
#    modifies: seqCStr to point at bast after soft mask
#    modifies: qCStr to point at q-score after soft mask
######################################################################*/
void checkSoftMasks(
    const uint32_t *cigUInt,    /*Number of bases soft masked*/
    struct samEntry *samStruct, /*Holds sequence and Q-score entry*/
    const uint8_t *qLineBool,   /*0 no q-score entry; 1 present*/
    const int32_t *incInt          /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
); /*Moves sequence & Q-score pointers past soft mask*/

/*######################################################################
# Output:
#    Modifies: refStruct to hold the reference sequence
#    Returns:
#        -0: if no reference file provided
#        -1 if the reference is good
#        -2 if refFILE was not a fastq
#        -4 if the reference was beneath the min requirements for an 
#           alignment
# Note:
#    - This is here instead of in the main fucntion for flexability.
######################################################################*/
uint8_t readAndCheckRef(
    FILE *refFILE,                /*Reference file to read in*/
    struct samEntry *refStruct,   /*Structer to hold reference*/
    struct minAlnStats *minStats  /*Minimum stats to keep an alignment*/
); /*Reads reference from fastq & check that it meets min requirments*/

/*######################################################################
# Output:
#    modifies: retULng to hold number of bases in cigar entry
#    modifies: cigarUCStr to point to entry type (SNP, indel, ect...)
# Warning:
#    - This does not check for long overflows, however, there should be
#      no sequence with more bases than an unsigned long can count
######################################################################*/
void readCigEntry(
    char **cigarUCStr,   /*c-string cigar to read & incurment*/
    uint32_t *retUInt             /*Holds returned long*/
); /*Reads a single entry from a eqx cigar line*/

/*######################################################################
# Output:
#    modifies: retULng to hold number of bases in cigar entry
#    modifies: cigarUCStr to point to entry type (SNP, indel, ect...)
# Warning:
#    - This does not check for long overflows, however, there should be
#      no sequence with more bases than an unsigned long can count
#    - You should have a pointer to get the entry type, since it comes
#      first
######################################################################*/
void readReverseCigEntry(
    char **cigarUCStr,     /*c-string cigar to read & incurment*/
    uint32_t *retUInt         /*Holds returned long*/
); /*Reads a single entry from a eqx cigar line*/

#endif
