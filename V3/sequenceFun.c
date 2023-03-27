/*######################################################################
# Name: sequenceFun
# Use:
#  o Holds functions for reading in or manipulating sequences
#  o All functions in this file will only really on c standard libraries
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
'     o Reads a fastq sequence from a fastq file
'  - fun-04 readFaSeq:
'     o Grabs the next read in the fasta file
'  - fun-05 addLineToBuff:
'     o Add characters from file to buffer, if needed resize. This
'       will only read in till the end of the line
'  - fun-06 reverseCStr;
'     o Reverse a c-string to be backwards (here for Q-score entries)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "sequenceFun.h"

/*---------------------------------------------------------------------\
| Output: Modfies seqCStr to be the reverse complement sequence
\---------------------------------------------------------------------*/
void reverseComplementSeq(
    char *seqCStr,    // Sequence to refeverse complement
    uint32_t lenSeqUI // Length of input sequence (index 1)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: reverseComplementSeq
   '  - Reverse complement a sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *endCStr = seqCStr + lenSeqUI - 1;
    char swapC = 0;

    while(endCStr > seqCStr)
    { // While I have bases to reverse complement
        swapC = *seqCStr;

        *seqCStr = complementBase(endCStr);
        *endCStr = complementBase(&swapC);
        ++seqCStr;
        --endCStr;
    } // While I have bases to reverse complement

    // Check if ended on the same base (if so complement the base)
    if(endCStr == seqCStr)
        *seqCStr = complementBase(seqCStr);

    return;
} // reverseComplementSeq

/*---------------------------------------------------------------------\
| Output: Returns the complement of the input base (0 if invalid base)
\---------------------------------------------------------------------*/
char complementBase(
    const char *baseC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: Sec-1 Sub-1: complementBase
   '  - Returns the complement of a base
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    switch(*baseC & (~32)) // remove 32 to convert to upper case
    { // switch, reverse complement
        case 'A': return 'C';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'U': return 'A';
        case 'W': return 'W'; // W is its reverse complement (AT=TA)
        case 'S': return 'S'; // S is its reverse complement (CG=GC)
        case 'M': return 'K'; // A/C to T/G
        case 'K': return 'M'; // T/G to A/C
        case 'R': return 'Y'; // A/G to T/C
        case 'Y': return 'R'; // T/C to A/G
        case 'B': return 'V'; // C/G/T to G/C/A
        case 'D': return 'H'; // G/T/A to C/A/T
        case 'H': return 'D'; // C/A/T to G/T/A
        case 'V': return 'B'; // A/C/G to T/G/C
        case 'N': return 'N'; // A/C/G/T to A/C/G/T
    } // switch, reverse complement

    return 0; // Invalid enry
} // complementBase

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: readRefFqSeq
   '  -  Grabs the next read in the fastq file
   '    fun-03 sec-1: Variable declarations
   '    fun-03 sec-2: Check if need to allocate memory for buffer
   '    fun-03 sec-3: Read in the first data
   '    fun-03 sec-4: If not at file end, see if have the full entry
   '    fun-03 sec-5: Read till end of file, check if valid fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-1: Variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char tmpC = 'C';               // To initilize the sequence loop
    unsigned char numLinesUC = 0; /*How many lines in sequence entry*/
    uint16_t extraBuffUS = 1024;
    unsigned long bytesInBuffUL = 0;
    unsigned long tmpBuffUL = 0;

    uint8_t errUC = 0;
    char *oldIterCStr = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-2: Read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(fqFILE == 0)
        return 2;    /*No file provided*/

    errUC =
        addLineToBuff(
            headerCStr,            // C-string to hold the header
            lenHeadUI,             // Length of the header buffer
            &bytesInBuffUL,        // Number of bytes in the buffer
            extraBuffUS,           // Amount to resize buffer by if full
            fqFILE                // Fastq file to get header from
    ); // Get the header (will resize as needed)

    // Check if have EOF (0) or memory allocation error (64)
    if(!(errUC & 1)) return errUC;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-3: Read in the sequence & spacer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    bytesInBuffUL = 0; // No bytes in the sequence c-string
    /*need to set this up so the loop does not error out*/
    oldIterCStr = &tmpC;

    while(*oldIterCStr != '+')
    { /*While I have not reached the spacer entry*/
        tmpBuffUL = bytesInBuffUL; // So can get to this position later

        errUC =
            addLineToBuff(
                seqCStr,          // Buffer to hold the sequence entry
                lenSeqUI,         // Size of the sequence buffer
                &bytesInBuffUL,   // Number of bytes in the buffer
                extraBuffUS,      // Amount to resize buffer by if full
                fqFILE           // Fastq file to get sequence from
        ); /*Get the header*/

        // Make sure new lines are removed on the next read in
        oldIterCStr = *seqCStr + bytesInBuffUL - 1;

        while(*oldIterCStr < 33) // Will catch some white space
        { // While removing new lines
            --bytesInBuffUL;
            --oldIterCStr;
        } // While removing new lines

        if(errUC & 64) return errUC;   // Memory allocation error (64)
        if(errUC == 0) return 2 + 128; // Invalid fastq entry

        /*Get on first character in the new buffer*/
        oldIterCStr = *seqCStr + tmpBuffUL;
        ++numLinesUC; /*Count number of new lines in sequence entry*/
    } /*While I have not reached the spacer entry*/

    --numLinesUC; /*Account for the overcounting*/

    oldIterCStr = *seqCStr + bytesInBuffUL;

    while(*oldIterCStr != '+')
    { // While not at start of spacer entry
        --oldIterCStr;
        --bytesInBuffUL; // Acount for white space
    } // While not at start of spacer entry

    *oldIterCStr = '\0';  // Set spacer to null
    --oldIterCStr;
    --bytesInBuffUL; // Acount for white space

    while(*oldIterCStr < 33)
    { // While have white space at end to remove
        *oldIterCStr = '\0';
        --oldIterCStr;
        --bytesInBuffUL; // Acount for white space
    } // While have white space at end to remove
    
    *basesInSeqUI = bytesInBuffUL + 1; // Index 0

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-03 Sec-4: Read in the Q-score entry
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    bytesInBuffUL = 0;

    while(numLinesUC > 0)
    { // While I need to read in the Q-score entry
        errUC =
            addLineToBuff(
                qCStr,             // Buffer to hold the Q-score entry
                lenQUI,            // Size of the q-score entry buffer
                &bytesInBuffUL,    // Number of bytes in q-score buffer
                extraBuffUS,       // Amount to resize buffer by if full
                fqFILE            // Fastq file to move past spacer
        ); /*Get the header*/

        // Make sure new lines are removed on the next read in
        oldIterCStr = *qCStr + bytesInBuffUL - 1;

        while(*oldIterCStr < 33) // Will catch some white space
        { // While removing new lines
            --bytesInBuffUL;
            --oldIterCStr;
        } // While removing new lines

        if(errUC & 64) return errUC;   /*memory allocation error (64)*/
        if(errUC == 0) return 2 + 128; /*Invalid fastq entry*/

        --numLinesUC; /*Count number of new lines in sequence entry*/
    } // While I need to read in the Q-score entry

    // Remove any white space at the end
    oldIterCStr = *qCStr + bytesInBuffUL;

    while(*oldIterCStr < 33)
    { // While have white space at end to remove
        *oldIterCStr = '\0';
        --oldIterCStr;
        --bytesInBuffUL; // Acount for white space
    } // While have white space at end to remove

    return errUC; /*Is 1 for another entry or 0 for EOF*/
} /*readFqSeq*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: readFaSeq
   '  -  Grabs the next read in the fasta file
   '    fun-04 sec-1: Variable declarations
   '    fun-04 sec-2: Check if need to allocate memory for buffer
   '    fun-04 sec-3: Read in the sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-04 Sec-1: Variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char tmpC = 'C';               // To initilize the sequence loop
    uint16_t extraBuffUS = 1024;
    unsigned long bytesInBuffUL = 0;
    unsigned long tmpBuffUL = 0;

    uint8_t errUC = 0;
    char *oldIterCStr = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-04 Sec-2: Read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(faFILE == 0)
        return 2;    /*No file provided*/

    errUC =
        addLineToBuff(
            headerCStr,            // C-string to hold the header
            lenHeadUI,             // Length of the header buffer
            &bytesInBuffUL,        // Number of bytes in the buffer
            extraBuffUS,           // Amount to resize buffer by if full
            faFILE                // Fastq file to get header from
    ); // Get the header (will resize as needed)

    // Check if have EOF (0) or memory allocation error (64)
    if(!(errUC & 1)) return errUC;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-04 Sec-3: Read in the sequence
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    bytesInBuffUL = 0; // No bytes in the sequence c-string
    /*need to set this up so the loop does not error out*/
    oldIterCStr = &tmpC;
    if(seqCStr != 0) *seqCStr = '\0'; // So know if missed a sequence

    while(*oldIterCStr != '>')
    { /*While I have not reached the spacer entry*/
        tmpBuffUL = bytesInBuffUL; // So can get to this position later

        errUC =
            addLineToBuff(
                seqCStr,          // Buffer to hold the sequence entry
                lenSeqUI,         // Size of the sequence buffer
                &bytesInBuffUL,   // Number of bytes in the buffer
                extraBuffUS,      // Amount to resize buffer by if full
                faFILE           // Fastq file to get sequence from
        ); /*Get the header*/

        // Make sure new lines are removed on the next read in
        oldIterCStr = *seqCStr + bytesInBuffUL - 1;

        while(*oldIterCStr < 33) // Will catch some white space
        { // While removing new lines
            --bytesInBuffUL;
            --oldIterCStr;
        } // While removing new lines

        if(errUC & 64) return errUC;   // Memory allocation error (64)
        if(errUC == 0) break;          // End of file (finsh processing)

        /*Get on first character in the new buffer*/
        oldIterCStr = *seqCStr + tmpBuffUL;
    } /*While I have not reached the spacer entry*/

    if(*seqCStr == 0 || **seqCStr == '\0') return 2; //Not a valid fasta

    // Make sure the new line at the end is removed
    oldIterCStr = *seqCStr + bytesInBuffUL;

    while(*oldIterCStr < 33) // Will catch some white space
    { // While I have white space at the end
        *oldIterCStr = '\0';
        --oldIterCStr;
        --bytesInBuffUL;
    } // While I have white space at the end

    *basesInSeqUI = bytesInBuffUL + 1;
    return errUC; /*Is 1 for another entry or 0 for EOF*/
} /*readFaSeq*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: addLineToBuff
   '  - Add characters from file to buffer, if needed resize. This
   '    will only read in till the end of the line
   '   o fun-05 sec-1: variable declerations
   '   o fun-05 sec-2: Check if need to resize the buffer
   '   o fun-05 sec-3: Read in the next line in the buffer
   '   o fun-05 sec-4: If at end of file, update read in lengths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-05 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    unsigned long spareBuffUL = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-05 Sec-2: Check if need to resize the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*curBuffUL == 0 && *lenBuffUL > 0)
        spareBuffUL = *lenBuffUL;

    else if(*lenBuffUL == 0 || *curBuffUL - 1 >= *lenBuffUL)
    { /*If need to resize the buffer (-1 for 0 index)*/
        *lenBuffUL += resBuffUL - 1;
            /*-1 to account for adding two index one items*/
        *buffCStr = realloc(*buffCStr, sizeof(char) * *lenBuffUL);

        if(*buffCStr == 0)
            return 64; /*Memory allocation error*/

        spareBuffUL = resBuffUL; /*Amount of extra space in the buffer*/
    } /*If need to resize the buffer*/

    else
        spareBuffUL = *lenBuffUL - *curBuffUL;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-05 Sec-3: Read in the next line in the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /*Set up marker to mark when the entire line read was in*/
    *(*buffCStr + *lenBuffUL - 2) = '\0';
    tmpCStr = *buffCStr + *curBuffUL;

    while(fgets(tmpCStr, spareBuffUL, inFILE))
    { /*While I have lines to read*/

        if(*(*buffCStr + *lenBuffUL - 2) == '\0' ||
           *(*buffCStr + *lenBuffUL - 2) == '\n'
        ) { /*If read in the line*/

            if(*(*buffCStr + *lenBuffUL - 2) == '\n')
                *curBuffUL = *lenBuffUL; /*used entire buffer*/
            else
            { /*Else only read in part of the buffer*/
                while(*tmpCStr != '\0')
                { /*While have characters in the buffer*/
                    ++tmpCStr;
                    ++(*curBuffUL);
                } /*While have characters in the buffer*/
            } /*Else only read in part of the buffer*/

            return 1; /*Read in entire line, return end of buff*/
        } /*If read in the line*/

        /*Else resize the buffer*/
        *curBuffUL = *lenBuffUL - 1; /*Filled up the buffer*/
            /*-1 to account for 0 index*/
        *lenBuffUL += resBuffUL - 1;
           /*-1 to account for adding two index 1 items*/
        *buffCStr = realloc(*buffCStr, sizeof(char) * *lenBuffUL);

        if(*buffCStr == 0)
            return 64; /*Memory allocation error*/

        /*Reset my maker for entire line read in*/
        *(*buffCStr + *lenBuffUL - 2) = '\0';
        spareBuffUL = resBuffUL; /*Amount of extra space in the buffer*/
        tmpCStr = *buffCStr + *curBuffUL;
    } /*While I have lines to read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-05 Sec-4: If at end of file, update read in lengths
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*(*buffCStr + *lenBuffUL - 2) == '\n')
        *curBuffUL = *lenBuffUL; /*used entire buffer*/
    else
    { /*Else only read in part of the buffer*/
        while(*tmpCStr != '\0')
        { /*While have characters in the buffer*/
            ++tmpCStr;
            ++(*curBuffUL);
        } /*While have characters in the buffer*/
    } /*Else only read in part of the buffer*/

    return 0; /*End of file*/
} /*addLineToBuff*/

/*---------------------------------------------------------------------\
| Output: Modifies inCStr to be backwards (end at start, start at end)
\---------------------------------------------------------------------*/
void reverseCStr(
    char *inCStr,       // C-string to refeverse
    uint32_t lenCStrUI // Length of input string (index 1)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: reverseCStr
   '  - Reverse a c-string to be backwards (here for Q-score entries)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *endCStr = inCStr + lenCStrUI - 1;
    char swapC = 0;

    while(endCStr > inCStr)
    { // While I have bases to reverse complement
        swapC = *inCStr;

        *inCStr = *endCStr;
        *endCStr = swapC;
        ++inCStr;
        --endCStr;
    } // While I have bases to reverse complement

    return;
} // reverseCStr

