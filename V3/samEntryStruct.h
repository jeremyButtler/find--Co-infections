/*######################################################################
# Name: samEntryStruct
# Use:
#    - Holds structer to hold a sam file entry. This also includes
#      the functions needed to support this structer.
# Includes:
#    - <string.h>
#    - <stdlib.h>
#    - "cStrToNumberFun.h"
#        - <sdtint.h>
#    - "printError.h"
#        - <stdio.h>
# Note: End of file has sam file format & flags
######################################################################*/

#ifndef SAMENTRYSTRUCT_H
#define SAMENTRYSTRUCT_H

#include <stdlib.h>          /*memory allocation*/
#include <string.h>
#include "cStrToNumberFun.h"
#include "printErrors.h"

#define Q_ADJUST 33 /*offest to get q-score of 0*/
#define MAX_Q_SCORE 94 /*highest possible Q-score*/

/*---------------------------------------------------------------------\
| Struct-1: samEntry
\---------------------------------------------------------------------*/
typedef struct samEntry
{ /*samEntry*/
    char *samEntryCStr;/*Holds the c-string for the sam entry*/

    char
        *queryCStr,    /*points to the query name in samEntryCStr*/
        *refCStr,      /*points ot the reference in samEntryCStr*/
        *cigarCStr,    /*points to the cigar in samEntryCStr*/
        *seqCStr,      /*points to the seqeuence in samEntryCStr*/
        *qCStr;        /*points to the q-score in samEntryCStr*/

    uint8_t
        mapqUChar;     /*Holds the mapping quality in samEntryCStr*/

    uint16_t flagUSht; /*Holds the flag in samEntryCStr*/

    float
       medianQFlt,     /*holds the median read Q-score*/
       medianAligQFlt, /*holds the aligned read median Q-score*/
       meanQFlt,       /*holds the mean read Q-score*/
       meanAligQFlt;   /*holds mean aligned read Q-score (no low Q)*/

    uint32_t
        posOnRefUInt,      /*Points to starting position on reference*/
        unTrimReadLenUInt, /*Holds length of untrimmed read*/
        readLenUInt,       /*Holds the read length of sam entry*/
        readAligLenUInt,   /*Aligned read length of sam entry*/
        numMatchUInt,      /*Holds number of matches*/
        numKeptMatchUInt,  /*Holds number matches with Q-score > min Q*/
        numKeptSNPUInt,    /*number of kept mismatches in sam entry*/
        numSNPUInt,        /*total number of mismatches in sam entry*/
        numKeptDelUInt,    /*number of kept deletions in sam entry*/
        numDelUInt,        /*total number of deletions in sam entry*/
        numKeptInsUInt,    /*number of kept insertions in sam entry*/
        numInsUInt,        /*total number of insertions in sam entry*/
        seqQHistUInt[MAX_Q_SCORE],  /*Histogram of base Q-scores*/
        seqQAlnHistUInt[MAX_Q_SCORE]; /*Histogram of kept base Q-score*/


    uint64_t
        totalQScoreULng, /*Q-score of all bases added together*/
        totalAlnQScoreULng; /*Q-score of kept bases added together*/

    unsigned long lenBuffULng; /*# bytes allocated to samEntryCStr*/
}samEntry; /*samFileEntry*/

/*---------------------------------------------------------------------\
| Struct-2: readStat
| Use:
|    - Holds query id, reference id, & stats for a single read
\---------------------------------------------------------------------*/
typedef struct readStat
{ /*readStat*/
    char
        queryIdCStr[100],
        refIdCStr[100];

    unsigned char
        mapqUChar;

    float
       medianQFlt,     /*holds the median read Q-score*/
       medianAligQFlt, /*holds the aligned read median Q-score*/
       meanQFlt,       /*holds the mean read Q-score*/
       meanAligQFlt;   /*holds mean aligned read Q-score (no low Q)*/

    uint32_t
        readLenUInt,       /*Holds the read length of sam entry*/
        readAligLenUInt,   /*Aligned read length of sam entry*/
        numMatchUInt,      /*Holds number of matches*/
        numKeptMatchUInt,  /*Holds number matches with Q-score > min Q*/
        numKeptSNPUInt,    /*number of kept mismatches in sam entry*/
        numSNPUInt,        /*total number of mismatches in sam entry*/
        numKeptDelUInt,    /*number of kept deletions in sam entry*/
        numDelUInt,        /*total number of deletions in sam entry*/
        numKeptInsUInt,    /*number of kept insertions in sam entry*/
        numInsUInt;        /*total number of insertions in sam entry*/
}readStat;


/*######################################################################
# Output: Modifies: Sets every variable but samEntryCStr to 0
######################################################################*/
void blankSamEntry(
    struct samEntry *samEntryStruct /*samEntryStruct to blank*/
); /*Sets all non-alloacted variables in samEntryStruct to 0*/

/*######################################################################
# Output:
#    - Modifies: Sets every variable but samEntryCStr to 0
#    - Modifies: the first character in samEntryCStr to be '\0'
######################################################################*/
void blankReadStats(
    struct samEntry *samEntryStruct /*Structure to blank*/
); /*modifes samEntryStruct to have every non-pointer variable set to 0*/

/*######################################################################
# Output: 
#    Modifies: newSamEntry to hold variables in oldSamEntry
#    Returns: 0 if failed (memory allowcation error) or pointer to new
#             sam entry if suceeded. Also prints warning to stdout.
######################################################################*/
uint8_t deepCpSamReadStats(
    struct samEntry *oldSamEntry,   /*sam entry to copy*/
    struct samEntry *newSamEntry    /*sam entry to copy to*/
); /*Copies read stats and sam file line from on samEntry to another.
    This ignores alignment specific items, such as cigars and number of
    matches.*/

/*######################################################################
# Output: 
#    Modifies: newSamEntry to use oldSamEntry q-score & sequence
#              pointers
######################################################################*/
void cpSamEntry(
    samEntry *oldSamEntry, /*Copy q-score & sequence pointers from*/
    samEntry *newSamEntry  /*Copy q-score & sequence pointers to*/
); /*Copies address of q-score & sequence pointers from old*/

/*######################################################################
# Output: Modifies: samEntry to have all variables set to 0
# Warning: This sets samEntryCStr to 0 without freeing.
######################################################################*/
void initSamEntry(
    struct samEntry *samEntry
); /*Initalize a samEntry struct to 0's*/

void freeStackSamEntry(
    struct samEntry *samEntry
); /*Frees heap allocations in a stack allocated sam entry struct*/

/*######################################################################
# Output: Frees: samEntry & sets to 0
######################################################################*/
void freeHeapSamEntry(
    struct samEntry **samEntry
); /*Frees the samEntry structer*/

/*######################################################################
# Output:
#     Modifies: read id, reference id, q-score, cigar, & sequence
#               pointers to piont to their entires in the sam entry.
#               This also grabs the flag & mapq.
######################################################################*/
void processSamEntry(
    struct samEntry *samEntry /*Has sam file entry to find data in*/
); /*Sets Q-score, cigar, & sequence pionters in samEntry. Also finds
    the mapping quality, sam flag, & sequence length*/

/*######################################################################
# Output:
#    Returns:
#        - 1 if succeded
#        - 2 if end of file
#        - 64 if memory allocation error
#    Modifies: All pointers in samStruct & the variables for read Lenth
######################################################################*/
uint8_t readSamLine(
    struct samEntry *samStruct,
    FILE *inFILE
); /*Reads in sam entry & sets pointers in samStruct*/

/*######################################################################
# Output:
#     Returns:
#        - 1 if succeded
#        - 2 if end of file
#        - 4 if end of file, but not a complete sam entry
#        - 64 if memory allocation error
#     Modifies: samStruct->samEntryCStr, samStruct->lenBuffULng, &
#               samStruct->unTrimedLenULng
# Note:
#    - This requires that samStruct->samEntryCStr has read in at least
#      part of the sam entry
######################################################################*/
uint8_t findSamCig(
    FILE *inFILE,               /*File to search for cigar in*/
    struct samEntry *samStruct
); /*Finds the cigar in same entry (uses fgets if needs more buffer) &
    also finds the number of bases in query from cigar*/

/*######################################################################
# Output:
#    Returns: pointer to unsigned char with fgets return value
#        - 0 for eof or memory allocation error (no data read in)
#        - pointer to start of buffer with characters read in
#    Sets: the second to last entry (samStruct->lenBuffULng - 2) in
#          samStruct->samEntryCStr to '\0'. This allows the user to 
#          check if fgets read in the entire line.
#          This position will be '\0' or '\n' if entire line read in.
######################################################################*/
char * readNextPartOfLine(
    FILE *inFILE,                   /*File to read from*/
    struct samEntry *samStruct,
    uint64_t buffChangeULng  /*How much to increase buffer by*/
); /*Reads in next part of line when fgets did not get a full line*/

/*######################################################################
# Output:
#     Prints: line with stats from samEntryStruct to file
#     Modifies: printHeaderChar to be 0 if set to 1
######################################################################*/
void printSamStats(
    struct samEntry *samEntryStruct, /*Sam entry to print stats for*/
    uint8_t *printHeaderChar,/*1: print header, 0: do not print header*/
    FILE *outFILE             /*File to output to*/
);/*Prints stats in samEntry structer for a sam file entry to tsv file*/

/*######################################################################
# Output:
#    Prints: sam file entry in samStrcut to a file (outFILE).
######################################################################*/
void printSamEntry(
    struct samEntry *samStruct, /*Has sam entry line to print out*/
    FILE *outFILE               /*File to print sam entry to*/
); /*Prints the sam entry in samStruct to a file (outFILE)*/

/*---------------------------------------------------------------------\
| Output:
|   o Prints stat file header to the input file
\---------------------------------------------------------------------*/
void printStatHeader(
    FILE *statFILE /*File to print the hader to*/
); /*Prints the stat file header made using a sam entry struct*/

/*######################################################################
# Output:
#    File: fastq with the read id, sequence, & q-score entry of the sam
#                entry
######################################################################*/
void samToFq(
    struct samEntry *samStruct, /*Has sam entry to print as fastq*/
    FILE *outFILE               /*file to print fastq entry to*/
); /*Converts sam entry into a fastq entry & prints to a fastq file*/

/*######################################################################
# Output:
#    Prints: Sam file entry in samStruct to outFILE.
######################################################################*/
void printSamLine(
    struct samEntry *samStruct, /*Has sam entry to print out*/
    FILE *outFILE               /*File to print sam entry to*/
); /*Prints out the sam entry in the samStruct, does not print stats*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - readToBlank to have all stats set to 0 & c-strins to start 
|          with '\0'
\---------------------------------------------------------------------*/
void blankReadStat(
    struct readStat *readToBlank
);

/*---------------------------------------------------------------------\
| Output:
|    Modifies: newReadList to have the same values as oldReadList
\---------------------------------------------------------------------*/
void cpReadStat(
    struct readStat *newReadStat, /*Read to copy stats to*/
    struct readStat *oldReadStat /*Read to copy stats from*/
); /*Copies stats from one read list to another*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Modifies:                                                       |
|        - readStruct to have stats from the next line in statsFILE    |
\---------------------------------------------------------------------*/
uint8_t readStatsFileLine(
    FILE *statsFILE,            /*File with line to grab*/
    uint8_t *onHeaderBool,      /*1: skip one line, 0: grab first line*/
    struct readStat *readStruct /*Holds the stats from the stats file*/
); /*Reads single line from printSamStats function in samEntryStruct.c*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|   Prints: Prints out variables in readToPrint structer               |
\---------------------------------------------------------------------*/
void printReadStat(
    struct readStat *readToPrint, /*Read to print stats out for*/
    FILE *outFILE                  /*File to print read to*/
); /*Prints out the stats in a readStat structure to a file*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies:
|        - newBin to have stats in samStruct
\---------------------------------------------------------------------*/
void samEntryToReadStat(
    struct readStat *newBin,   /*Read bin to hold stats from samStruct*/
    struct samEntry *samStruct /*copy stats from this struct*/
); /*Copies stats from a samEntry struct to a readStat struct*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies
|    - refStruct: To hold the sequence (no header)
|  - Returns
|    - 1 if succeeded
|    - 2 if file does not exist
|    - 4 invalid file
|    - 64 memory allocation error
| Note:
|  - Fasta file should only have one sequence
\---------------------------------------------------------------------*/
unsigned char readInConFa(
    char *conFaToReadCStr, /*Name of fasta file with the consensus*/
    struct samEntry *refStruct /*Sam struct to hold consensus*/
); /*Reads in reference sequence in fasta file*/

#endif

/*
Sam file table for first 11 columns (all sam files have)
| Col | Field |  Type  |        Brief description              |
|:---:|:-----:|:------:|:-------------------------------------:|
|  1  | QNAME | String |       Query template NAME             |
|  2  | FLAG  |  Int   |          bitwise FLAG                 |
|  3  | RNAME | String |     Reference sequence NAME           |
|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
|  5  | MAPQ  |  Int   |          MAPping Quality              |
|  6  | CIGAR | String |            CIGAR string               |
|  7  | RNEXT | String | Reference name of the mate/next read  |
|  8  | PNEXT |  Int   |   Position of the mate/next read      |
|  9  | TLEN  |  Int   |      observed Template LENgth         |
| 10  |  SEQ  | String |          segment SEQuence             |
| 11  | QUAL  | String | ASCII of Phred-scaled base Quality+33 |
*/

/*eqx cigar entry
    matches are '=' or blank (no symbol) at end or if only matches
    mimsatches are 'X'
    deletions are 'D'
    insertions are 'I'
    soft masks are 'S'
    Numbers come before entry (=, X, D, I, or S) & show number of times
        the entry is repeated
    Everything else is a hard mask (was removed)

    EX: 10S5=1X701 (10 soft masked bases, 5 matches, 1 SNP, 701 matches)
*/

/* Sam file flag values tables
| Bit  | FLAG  |                        Description                                 |
|:----:|:-----:|:------------------------------------------------------------------:|
| 1    | 0x1   | template having multiple segments in sequencing                    |
| 2    | 0x2   | each segment properly aligned according to the aligner             |
| 4    | 0x4   | segment unmapped                                                   |
| 8    | 0x8   | next segment in the template unmapped                              |
| 16   | 0x10  | SEQ being reverse complemented                                     |
| 32   | 0x20  | SEQ of the next segment in the template being reverse complemented |
| 64   | 0x40  | the first segment in the template                                  |
| 128  | 0x80  | the last segment in the template                                   |
| 256  | 0x100 | secondary alignment                                                |
| 512  | 0x200 | not passing filters, such as platform/vendor quality controls      |
| 1024 | 0x400 | PCR or optical duplicate                                           |
| 2048 | 0x800 | supplementary alignment                                            |
*/
