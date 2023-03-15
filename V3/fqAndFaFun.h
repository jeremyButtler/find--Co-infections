/*######################################################################
# Use:
#   o Holds functions for manipulating fastq and fasta files.
# Non c-standard includes:
#   - "minAlnStatsStruct.h"
#   - "FCIStatsFun.h"
#   o "defaultSettings.h"
#   o "samEntryStruct.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
# C standard Includes (all though non c-standard includes):
#   o <stdlib.h>
#   o <sdtint.h>
#   o <stdio.h>
######################################################################*/

#ifndef FQANDFAFUN_H
#define FQANDFAFUN_H

#include "minAlnStatsStruct.h"
#include "FCIStatsFun.h"

/*---------------------------------------------------------------------\
| Output:
|   o Creates:
|      - Fasta file with reads from the fastq file
|   o Returns:
|      - 1: Success
|      - 2: No fastq file (or invalid fastq file)
|      - 4: No fasta file
|      - 64: memory allocation error
\---------------------------------------------------------------------*/
unsigned char fqToFa(
    char *fqToCnvtCStr,         /*File name of fastq to convert*/
    char *outFaCStr,            /*File name of new fasta file*/
    struct samEntry *samStruct  /*Holds fastq file entries*/
); /*Converts a fastq file to a fasta file*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: refStruct to hold the read in fastq entry & sets its
|              pointers
|    Returns:
|        - 0: if no reference file provided
|        - 1: if succeded
|        - 2: If file was not a fastq file
|        - 64: If malloc failed to find memory
\---------------------------------------------------------------------*/
uint8_t readRefFqSeq(
    FILE *refFILE,       /*Pointer to fastq file to grab reference from*/
    struct samEntry *refStruct,/*Sam entry struct to hold reference*/
    char quickRunBl
      /*1: do not remove new lines, 0: remove any extra new lines*/
); /*Gets the frist reads sequence & q-score line from a fastq file*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o buffCStr to hold the next line.
|       - buffCStr is resizied if it is to small to hold the next line.
|       - buffCStr + lenBuffUL - 2 will be '\0' or '\n'
|       - buffCStr will be 0 if had a memory allocation error
|     o curBuffUL: To hold the number of characters read into the buffer
|     o lenBuffUL: To hold resized buffer size if buffCStr is resized
|     o inFILE: To point to the next line (fgets does this automaticly)
|   - Returns:
|     o 0 if was end of file (EOF)
|     o 1 if read in the next line
|     o 64 if had a memory allocation error
| Note:
|   - Yes, I am aware of getline, but getline will overwrite my buffer 
|     for each new line and is unix dependent
\---------------------------------------------------------------------*/
unsigned char addLineToBuff(
    char **buffCStr,          /*Buffer to add data to*/
    unsigned long *lenBuffUL, /*Size of the buffer*/
    unsigned long *curBuffUL, /*Length buffer with valid data*/
    unsigned long resBuffUL,  /*Amount to resize buffer by if full*/
    FILE *inFILE              /*File to grab data from*/
); /*Add characters from file to buffer, if needed resize.
    This will only read in till the end of the line.*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 0 if the fastq file could not be opened
|     o The number of reads in the fastq file
| Note:
|   - fqFILE will be set back to its starting position at the end
\---------------------------------------------------------------------*/
unsigned long getNumReadsInFq(
    char *fqFileCStr /*Path to fastq file to get number of reads in*/
); /*Find the number of reads in a fastq file*/

/*---------------------------------------------------------------------\
| Output:
|   o Writes a copy of the original file (orgFqCStr):
|   o Returns:
|     - 1 if no problems happened
|     - 2 if could not open orgFqCStr
|     - 4 if could not open newFqCStr
\---------------------------------------------------------------------*/
unsigned char copyFile(
    char *orgFqCStr,  /*Path to fastq file to copy (orignal)*/
    char *newFqCStr   /*Path to the new duplicate Fastq file*/
); /*Copies a fastq file to a new file*/

/*---------------------------------------------------------------------\
| Output:
|    Creates File: from outFqPathCStr with filtered reads
|    Returns:
|        - 1: if succeded
|        - 2: If file was not a fastq file
|        - 130: If file had an invalide entry
|            - This error uses to flags, first it uses 2 to specify that
|              it is not a fastq file (invalid file).
|            - 2nd it uses 128 to specifty that it is not an blank file
|        - 64: If malloc failed to find memory
\---------------------------------------------------------------------*/
unsigned char filterReads(
    char *fqCStr,/*Fastq file with reads to filter (null for stdin)*/
    char *outCStr,/*Name of fastq to write reads to (null for stdout)*/
    struct samEntry *samST, /*For reading in lines for the fastq file*/
    struct minAlnStats *minStats
        /*Has min/max lengths & mean/median Q-score*/
); /*Filters reads in a fastq file by length and mean/median Q-score*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: bufferCStr to have the next buffer if empty
|    Modifies: incurments pointInBufferCStr to start of next read
|    Returns:
|        4: If the end of the file
|        2: if nothing went wrong
|        0: If ran out of file
\---------------------------------------------------------------------*/
uint8_t moveToNextFastqEntry(
    char *bufferCStr,     /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr,     /*position working on in buffer*/
    uint32_t buffSizeInt,     /*Size of buffer to work on*/
    uint64_t *lenInputULng,  /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
); /*Moves to next fastq read, without printing out*/

#endif
