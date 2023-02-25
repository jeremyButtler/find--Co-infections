/*######################################################################
# Use:
#   o Holds functions for manipulating fastq and fasta files.
# Includes:
#   o "samEntryStruct.h"
#     - <stdlib.h>
#     - "cStrToNumberFun.h"
#       o <sdtint.h>
#     - "printError.h"
#       o <stdio.h>
######################################################################*/

#ifndef FQANDFAFUN_H
#define FQANDFAFUN_H

#include "samEntryStruct.h"

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


#endif
