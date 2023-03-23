/*######################################################################
# Name: sequenceFun
# Use:
#  o Holds functions for reading in or manipulating sequences
#  o All functions in this file will only use c standard libraries
# C standard libraries
#  - <stdlib.h>
#  - <stdio.h>
#  - <stdint.h>
######################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 reverseComplementSeq:
'     o Reverse complement a sequence
'  - fun-02 complementBase:
'     o Returns the complement of a base
'  - fun-03 readFqSeq:
'    o Reads a fastq sequence from a fastq file
'  - fun-04 readFaSeq:
'     o Grabs the next read in the fasta file
'  - fun-05 addLineToBuff:
'     o Add characters from file to buffer, if needed resize. This
'       will only read in till the end of the line
'  - fun-06 reverseCStr;
'     o Reverse a c-string to be backwards (here for Q-score entries)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef SEQUENCEFUN_H
#define SEQUENCEFUN_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*---------------------------------------------------------------------\
| Output: Modfies seqCStr to be the reverse complement sequence
\---------------------------------------------------------------------*/
void reverseComplementSeq(
    char *seqCStr,    // Sequence to refeverse complement
    uint32_t lenSeqUI // Length of input sequence (index 1)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: reverseComplementSeq
   '  - Reverse complement a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Returns the complement of the input base (0 if invalid base)
\---------------------------------------------------------------------*/
char complementBase(
    const char *baseC
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: complementBase
   '  - Return the complement of a base
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets its pointers
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: If file was not a fastq file
|     o 130: If file had an invalide entry
|         - This error uses to flags, first it uses 2 to specify that
|           it is not a fastq file (invalid file).
|         - 2nd it uses 128 to specifty that it is not an blank file
|     o 64: If malloc failed to find memory
\---------------------------------------------------------------------*/
uint8_t readFqSeq(
    FILE *fqFILE,      // Pointer to fastq file to grab sequence from
    char **headerCStr,  // Will hold the sequence header
    uint32_t *lenHeadUI, // Holds the size of the header array
    char **seqCStr,     // Will Hold the new sequence; resized if needed
    uint32_t *lenSeqUI, // Length of the sequence c-string
    char **qCStr,       // Will hold the q-score entry;resized if needed
    uint32_t *lenQUI,   // Length of the q-score c-string
    uint32_t *basesInSeqUI  // Will hOld the number of bases in the seq
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: readRefFqSeq
   '  -  Grabs the next read in the fastq file
   '    fun-2 sec-1: Variable declarations
   '    fun-2 sec-2: Check if need to allocate memory for buffer
   '    fun-2 sec-3: Read in the first data
   '    fun-2 sec-4: If not at file end, see if have the full entry
   '    fun-2 sec-5: Read till end of file, check if valid fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o refStruct to hold the read in a fasta entry
|  - Returns:
|     o 0: if EOF
|     o 1: if succeded
|     o 2: for an invalid fasta entry
|     o 64: If malloc failed to find memory
| Note:
|   - This will remove new lines from the sequence.
|   - This will only remove spaces or tabs at the end of each sequence
|     entry, so "atg atc \n" will got to "atc atc".
\---------------------------------------------------------------------*/
uint8_t readFaSeq(
    FILE *faFILE,      // Pointer to fastq file to grab sequence from
    char **headerCStr,  // Will hold the sequence header
    uint32_t *lenHeadUI, // Holds the size of the header array
    char **seqCStr,     // Will Hold the new sequence; resized if needed
    uint32_t *lenSeqUI, // Length of the sequence c-string
    uint32_t *basesInSeqUI  // Will hOld the number of bases in the seq
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: readFaSeq
   '  -  Grabs the next read in the fasta file
   '    fun-04 sec-1: Variable declarations
   '    fun-04 sec-2: Check if need to allocate memory for buffer
   '    fun-04 sec-3: Read in the sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
\---------------------------------------------------------------------*/
unsigned char addLineToBuff(
    char **buffCStr,          /*Buffer to add data to*/
    uint32_t *lenBuffUL, /*Size of the buffer*/
    unsigned long *curBuffUL, /*Length buffer with valid data*/
    unsigned long resBuffUL,  /*Amount to resize buffer by if full*/
    FILE *inFILE              /*File to grab data from*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: addLineToBuff
   '  - Add characters from file to buffer, if needed resize. This
   '    will only read in till the end of the line
   '   o fun-3 sec-1: variable declerations
   '   o fun-3 sec-2: Check if need to resize the buffer
   '   o fun-3 sec-3: Read in the next line in the buffer
   '   o fun-3 sec-4: If at end of file, update read in lengths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*---------------------------------------------------------------------\
| Output: Modifies inCStr to be backwards (end at start, start at end)
\---------------------------------------------------------------------*/
void reverseCStr(
    char *inCStr,       // C-string to refeverse
    uint32_t lenCStrUI // Length of input string (index 1)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: reverseCStr
   '  - Reverse a c-string to be backwards (here for Q-score entries)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
