/*######################################################################
# Name: trimSam
# Use:
#    - trims soft mask regions off all alignments with sequences in a
#      sam file. Aligments without sequences are ignored & not printed
#      out.
# Includes:
#   - "samEntryStruct.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
# C standard library includes:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>
######################################################################*/

#include "trimSam.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
| SOF: start of trimSam functions
|    fun-1 trimSamReads: Trims soft mask regions for all reads with
|                        a sequence in a sam file
|    fun-2 trimSamEntry: Trim soft mask regions off end of sam entry
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|    Prints: Trimmed sam entries with sequences to outFILE, but ignores
|            sam entries without sequences
\---------------------------------------------------------------------*/
void trimSamReads(
    FILE *samFILE,               /*Sam file to convert*/
    FILE *outFILE,               /*File to store output*/
    char keepUnmappedReadsBl     /*1: keep unmapped reads, 0: do not*/
) /*Prints trimmed sam entries with sequences to file*/
{ /*trimSamReads*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 Sec-1 Sub-1 TOC: trimSamReads
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    unsigned char
        errorFlagUChar = 0;/*Tells me if memory allocation error*/
    struct samEntry samStruct;

    initSamEntry(&samStruct);        /*Blank the structure for reading*/
    errorFlagUChar = readSamLine(&samStruct, samFILE);

    while(errorFlagUChar & 1)
    { /*While there are lines in the same file to convert*/

        if(*(samStruct.samEntryCStr) == '@')
        { /*If was a header*/
            printSamEntry(&samStruct, outFILE);
            blankSamEntry(&samStruct);
            errorFlagUChar = readSamLine(&samStruct, samFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        /*Convert & print out sam file entry*/
        errorFlagUChar = trimSamEntry(&samStruct);

        /*Print out the converted entry*/
        if(!(errorFlagUChar >> 2))
            printSamEntry(&samStruct, outFILE);
            /*If was header or had sequence, then print it out*/
        else if(errorFlagUChar & 4 && keepUnmappedReadsBl & 1)
            printSamEntry(&samStruct, outFILE);
            /*Else if printing umapped reads as well as mapped*/

        blankSamEntry(&samStruct);
        errorFlagUChar = readSamLine(&samStruct, samFILE);
    } /*While there are lines in the same file to convert*/

    freeStackSamEntry(&samStruct);
    return;
} /*trimSamReads*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 0 if suceeded
|        - 2 if header (invalid and ignored)
|        - 4 if an unmapped read (no reference)
|        - 8 if no sequence line
|    Modifies:
|        - Trims cigar, sequence, & q-score entries in samStruct.
\---------------------------------------------------------------------*/
uint8_t trimSamEntry(
    struct samEntry *samStruct   /*has sam line to trim softmasks*/
) /*Trims soft masked regions at start & end of sam entry*/
{ /*trimSamEntry*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: trimSamEntry
    '    fun-1 sec-1: Variable declerations
    '    fun-2 sec-2: Find how much to trim & trim cigar entry
    '    fun-2 sec-3: Trim the sequence entry
    '    fun-2 sec-4: Trim the q-score entry
    '    fun-2 sec-5: Shift characters for other parts of the cigar
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *incSamUCStr = samStruct->cigarCStr,
        *delUCStr = samStruct->cigarCStr;    /*Bases to remove*/

    unsigned char
        uCharTabCnt = 0;                     /*Counts number of tabs*/

    uint32_t
        lenStartTrimUInt = 0, /*Number of bases soft masked at start*/
        lenEndTrimUInt = 0;   /*Number of bases soft masked at end*/
        
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Find how much to trim & trim cigar entry
    ^    fun-2 sec-2 sub-1: Check start of cigar & trim if needed
    ^    fun-2 sec-2 sub-2: Check end of cigar & trim if needed
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-2 Sec-2 Sub-1: Check start of cigar & trim if needed
    \******************************************************************/

    if(*(samStruct->samEntryCStr) == '@')
        return 2;             /*Is a header line, invalid entry*/

    if(samStruct->flagUSht & 4)
        return 4;             /*Unmapped read*/
    if(*(samStruct->seqCStr) == '*')
        return 8;             /*No seqence provided for this alignment*/

    incSamUCStr = cStrToUInt(incSamUCStr, &lenStartTrimUInt);

    if(*incSamUCStr != 'S')
    { /*If no softmasking at start, set varaibles for no trim*/
        lenStartTrimUInt = 0; /*If no soft masking at start*/

        while(*incSamUCStr > 32)
            ++incSamUCStr;        /*Move to end of cigar entry*/
        delUCStr = incSamUCStr;
    } /*If no softmasking at start, set varaibles for no trim*/

    else
    { /*Else is a soft mask and need to remove*/
        ++incSamUCStr;            /*Move off soft mask marker*/

        while(*incSamUCStr > 32)
        { /*While not at the end of the cigar*/
            *delUCStr = *incSamUCStr;
            ++incSamUCStr;        /*Move to end of cigar entry*/
            ++delUCStr;           /*Move to next base to replace*/
        } /*While not at the end of the cigar*/
    } /*Else is a soft mask and need to remove*/

    /******************************************************************\
    * Fun-2 Sec-2 Sub-2: Check end of cigar & trim if needed
    \******************************************************************/

    --delUCStr;        /*Move back to last cigar entry*/

    if(*delUCStr != 'S')
    { /*If their is not softmasking at the end*/
        lenEndTrimUInt = 0; /*If no soft masking at start*/

        if(lenStartTrimUInt == 0)
            return 0; /*Nothing to trim*/
        ++delUCStr;             /*Account for the minus one*/
    } /*If their is not softmasking at the end*/

    else
    { /*Else I am trimming bases off the end*/
        /*Get number of bases to trim*/
        --delUCStr;       /*Move of 'S' marker for soft masking*/
        delUCStr = backwarsCStrToUInt(delUCStr, &lenEndTrimUInt);
        ++delUCStr;       /*Move off the entry before the soft mask*/
    } /*Else I am trimming bases off the end*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-3: Trim the sequence entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(uCharTabCnt < 4)
    { /*If need to find the sequence entry*/
        /*Find the sequence entry*/
        for(uint8_t uCharCnt = uCharTabCnt; uCharCnt < 4; ++uCharCnt)
        { /*Loop past RNEXT, PNEXT, & TLEN*/
            while(*incSamUCStr > 32)        /*Move to next entry*/
            { /*While not at the end of the entry*/
                *delUCStr = *incSamUCStr;
                ++delUCStr;
                ++incSamUCStr;
            } /*While not at the end of the entry*/

            /*Move leader of the tab*/
            *delUCStr = *incSamUCStr;
            ++delUCStr;
            ++incSamUCStr;
        } /*Loop past RNEXT, PNEXT, & TLEN*/
    } /*If need to find the sequence entry*/

    else
        incSamUCStr = samStruct->seqCStr;

    samStruct->seqCStr = delUCStr;   /*Record new sequence start*/
    incSamUCStr += lenStartTrimUInt; /*skip the trim region at start*/

    /*Find how long the sequence will be after triming*/
    samStruct->readLenUInt =
        samStruct->unTrimReadLenUInt - lenStartTrimUInt -lenEndTrimUInt;

    for(uint32_t uIntSeq=0; uIntSeq < samStruct->readLenUInt; ++uIntSeq)
    { /*Loop to shift bases back in sequence*/
        *delUCStr = *incSamUCStr;
        ++delUCStr;
        ++incSamUCStr;
    } /*Loop to shift bases back in sequence*/

    incSamUCStr += lenEndTrimUInt + 1; /*Move off trimed bases & tab*/
    *delUCStr = *(incSamUCStr - 1);    /*Save the tab*/
    ++delUCStr;                        /*Move past tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-4: Trim the q-score entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*(samStruct->qCStr) != '*' || *(samStruct->qCStr + 1) != '\t')
    { /*If their is a Q-score entry to trim*/
        samStruct->qCStr = delUCStr;    /*new start of q-score entry*/
        incSamUCStr += lenStartTrimUInt; /*skip starting trim region*/

        for(
            uint32_t uIntSeq=0;
            uIntSeq < samStruct->readLenUInt;
            ++uIntSeq
        ) { /*Loop to shift bases back in q-score entry*/
            *delUCStr = *incSamUCStr;
            ++delUCStr;
            ++incSamUCStr;
        } /*Loop to shift bases back in q-score entry*/

        incSamUCStr += lenEndTrimUInt; /*Move off end trim to tab*/
    } /*If their is a Q-score entry to trim*/

    else samStruct->qCStr = delUCStr;    /*new start of q-score entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-5: Shift characters for other parts of the cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*incSamUCStr != '\0' && *incSamUCStr != '\n')
    { /*While their are other elements in the q-score entry to shift*/
        *delUCStr = *incSamUCStr;
        ++delUCStr;
        ++incSamUCStr;
    } /*While their are other elements in the q-score entry to shift*/

    *delUCStr = *incSamUCStr;

    if(*delUCStr == '\n')
        *(delUCStr + 1) = '\0'; /*Make sure null terminated*/

    return 0;
} /*trimSamEntry*/
