/*######################################################################
# Name: samEntryStruct
# Use:
#    - Holds structer to hold a sam file entry. This also includes
#      the functions needed to support this structer.
# Includes:
#    - <stdlib.h>
#    - "cStrToNumberFun.h"
#        - <stdint.h>
#    - "printError.h"
#        - <stdio.h>
######################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    fun-1 blankSamEntry: Sets samEntry structer variables to 0
#        - Does not set samEntryCStr or lenBuffULng to 0
#    fun-2 blankReadStats: Sets all non-pointers to 0 in a samEntry
#    fun-3 deepCpSamReadStats: cpSamEntry, but also copies samEntryCStr
#    fun-4 cpSamEntry: Copies q-score & sequence pointer addresss
#    fun-5 initSamEntry: Initalize a new sam entry struct
#    fun-6 freeStackSamEntry: Free heap allocated variables in samEntry
#    fun-7 freeHeapSamEntry: Fee samEntry struct
#    fun-8 processSamEntry: Extract data from a sam file entry
#    fun-9: readSamLine: Read in a line form a sam file,
#        - calls processSamEntry
#    fun-10 findSamCig: find cigar in buffer with sam file entry. Will
#                       resize buffer if an incomplete entry.
#    fun-11 readNextPartOfLine: read part of a line with fgets
#    fun-12 printSamStats: Print stats in samEntry struct to tsv filE
#    fun-13 samToFq: Prints out sam entry as fastq file
#        - This does not print out the stats printed with printSamStats
#    fun-14 printSamEntry: Prints sam entry in a samEntry structer
#        - This does not print out the stats printed with printSamStats
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "samEntryStruct.h"

/*######################################################################
# Output:
#    - Modifies: Sets every variable but samEntryCStr to 0
#    - Modifies: the first character in samEntryCStr to be '\0'
######################################################################*/
void blankSamEntry(
    struct samEntry *samEntryStruct /*samEntryStruct to blank*/
) /*Sets all non-alloacted variables in samEntryStruct to 0*/
{ /*blankSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1 Sub-1 TOC: blankSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samEntryStruct->queryCStr = 0;
    samEntryStruct->refCStr = 0;
    samEntryStruct->cigarCStr = 0;
    samEntryStruct->seqCStr = 0;
    samEntryStruct->qCStr = 0; 

    if(samEntryStruct->samEntryCStr != 0)
        *samEntryStruct->samEntryCStr = '\0'; /*So user knows blanked*/

    samEntryStruct->flagUSht = 0;
    samEntryStruct->mapqUChar = 0;

    samEntryStruct->medianQFlt = 0;
    samEntryStruct->medianAligQFlt = 0;
    samEntryStruct->meanQFlt = 0;
    samEntryStruct->meanAligQFlt = 0;

    samEntryStruct->posOnRefUInt = 0;

    samEntryStruct->unTrimReadLenUInt = 0;
    samEntryStruct->readLenUInt = 0;
    samEntryStruct->readAligLenUInt = 0;

    samEntryStruct->numMatchUInt = 0;
    samEntryStruct->numKeptMatchUInt = 0;
    samEntryStruct->numKeptSNPUInt = 0;
    samEntryStruct->numSNPUInt = 0;
    samEntryStruct->numKeptInsUInt = 0;
    samEntryStruct->numInsUInt = 0;
    samEntryStruct->numKeptDelUInt = 0;
    samEntryStruct->numDelUInt = 0;

    samEntryStruct->totalQScoreULng = 0;
    samEntryStruct->totalAlnQScoreULng = 0;

    /*Reset the q-score histogram*/
    for(uint16_t uShtCnt = 0; uShtCnt < MAX_Q_SCORE; ++uShtCnt)
    { /*Set all values in histograms to 0*/
        samEntryStruct->seqQHistUInt[uShtCnt] = 0;
        samEntryStruct->seqQAlnHistUInt[uShtCnt] = 0;
    } /*Set all values in histograms to 0*/

    return;
} /*blankSamEntry*/

/*######################################################################
# Output: Modifies: samEntryStruct
######################################################################*/
void blankReadStats(
    struct samEntry *samEntryStruct /*Structure to blank*/
) /*modifes samEntryStruct to have every non-pointer variable set to 0*/
{ /*blankReadStats*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 Sub-1 TOC: blankSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*samEntryStruct->flagUSht = 0;
    samEntryStruct->mapqUChar = 0; These are taken from the sam file*/

    samEntryStruct->medianQFlt = 0;
    samEntryStruct->medianAligQFlt = 0;
    samEntryStruct->meanQFlt = 0;
    samEntryStruct->meanAligQFlt = 0;

    /*samEntryStruct->readLenUInt = 0;
    samEntryStruct->unTrimReadLenUInt = 0;
        found when adding read entry to a samEnty struct*/
    samEntryStruct->readAligLenUInt = 0;

    samEntryStruct->numMatchUInt = 0;
    samEntryStruct->numKeptMatchUInt = 0;
    samEntryStruct->numKeptSNPUInt = 0;
    samEntryStruct->numSNPUInt = 0;
    samEntryStruct->numKeptInsUInt = 0;
    samEntryStruct->numInsUInt = 0;
    samEntryStruct->numKeptDelUInt = 0;
    samEntryStruct->numDelUInt = 0;

    samEntryStruct->totalQScoreULng = 0;
    samEntryStruct->totalAlnQScoreULng = 0;

    /*Reset the q-score histogram*/
    for(uint16_t uShtCnt = 0; uShtCnt < MAX_Q_SCORE; ++uShtCnt)
    { /*Set all values in histograms to 0*/
        samEntryStruct->seqQHistUInt[uShtCnt] = 0;
        samEntryStruct->seqQAlnHistUInt[uShtCnt] = 0;
    } /*Set all values in histograms to 0*/

    return;
} /*blankReadStats*/

/*######################################################################
# Output: 
#    Modifies: newSamEntry to hold variables in oldSamEntry
#    Returns: 64 if failed (memory allowcation error)
#             1 if succeded
######################################################################*/
uint8_t deepCpSamReadStats(
    struct samEntry *oldSamEntry,   /*sam entry to copy*/
    struct samEntry *newSamEntry    /*sam entry to copy to*/
) /*Copies read stats and sam file line from on samEntry to another.
    This ignores alignment specific items, such as cigars and number of
    matches.*/
{ /*deepCpSamReadStats*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC: deepCpSamReadStats
    #    fun-3 sec-1: Variable declerations
    #    fun-3 sec-2: Check if need to resize arrays in newSamEntry
    #    fun-3 sec-3: Copy the sam file line (entry) and stats for entry
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char
        *lineCStr = oldSamEntry->samEntryCStr,
        *cpLineCStr = newSamEntry->samEntryCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Check if need to resize arrays in newSamEntry
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(oldSamEntry->lenBuffULng > newSamEntry->lenBuffULng)
    { /*If need to resize the new array*/
        if(newSamEntry->samEntryCStr != 0)
            free(newSamEntry->samEntryCStr);

        newSamEntry->samEntryCStr =
            malloc(sizeof(char) * oldSamEntry->lenBuffULng);

        if(newSamEntry->samEntryCStr == 0)
        { /*If memory allocation failed*/
            printMemAlocErr(
                "samFileStructs.c",
                "deepCpSamFileEntry",
                3,
                172
            ); /*Let user know of memory allocation failure*/

            return 64;
        } /*If memory allocation failed*/

        cpLineCStr = newSamEntry->samEntryCStr;
    } /*If need to resize the new array*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Copy the sam file line (entry) and stats for entry
    #    fun-3 sec-2 sub-1: Copy c-string with sam entry
    #    fun-3 sec-2 sub-2: Copy stats in from old entry to new
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-3 Sec-2 Sub-1: Copy c-string with sam file line to new buffer
    *******************************************************************/

    while(*lineCStr != '\0') /*So do not overwrite the sequence entry*/
    { /*Loop and copy the samfile line*/
        *cpLineCStr = *lineCStr; /*On a new sequence line*/
        cpLineCStr++;
        lineCStr++;
    } /*Loop and copy the samfile line*/

    *(cpLineCStr + 1) = '\0'; /*Make into a c-string*/

    /*******************************************************************
    # Fun-3 Sec-2 Sub-2: Copy shared query stats in old entry to new
    *******************************************************************/

    /*Set pointers to sequence and Q-score line in copied sam entry*/
    newSamEntry->seqCStr =
        newSamEntry->samEntryCStr +
        (oldSamEntry->seqCStr - oldSamEntry->queryCStr);

    (*newSamEntry).qCStr =
        newSamEntry->samEntryCStr +
        (oldSamEntry->qCStr - oldSamEntry->queryCStr);

    /*Record shared stats between oldSame querys*/
    newSamEntry->unTrimReadLenUInt = 0;
    newSamEntry->readLenUInt = oldSamEntry->readLenUInt;
    newSamEntry->meanQFlt = oldSamEntry->meanQFlt;
    newSamEntry->medianQFlt = oldSamEntry->medianQFlt;
    newSamEntry->mapqUChar = oldSamEntry->mapqUChar;
    newSamEntry->posOnRefUInt = 0;

    return 1;
} /*deepCpSamReadStats*/

/*######################################################################
# Output: 
#    Modifies: newSamEntry to use oldSamEntry q-score & sequence
#              pointers
######################################################################*/
void cpSamEntry(
    samEntry *oldSamEntry, /*Copy q-score & sequence pointers from*/
    samEntry *newSamEntry  /*Copy q-score & sequence pointers to*/
) /*Copies address of q-score & sequence pointers from old*/
{ /*cpSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 Sub-1: cpSamEntry
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Set pointers to sequence and Q-score line in copied sam entry*/
    newSamEntry->seqCStr = oldSamEntry->seqCStr;
    newSamEntry->qCStr = oldSamEntry->qCStr;

    /*Record shared stats between oldSame querys*/
    newSamEntry->unTrimReadLenUInt = 0;
    newSamEntry->readLenUInt = oldSamEntry->readLenUInt;
    newSamEntry->meanQFlt = oldSamEntry->meanQFlt;
    newSamEntry->medianQFlt = oldSamEntry->medianQFlt;

    return;
} /*cpSamEntry*/

/*######################################################################
# Output: Modifies: samEntry to have all variables set to 0
# Warning: This sets samEntryCStr to 0 without freeing.
######################################################################*/
void initSamEntry(
    struct samEntry *samEntry
) /*Initalize a samEntry struct to 0's*/
{ /*initSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 Sub-1 TOC: initSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samEntry->samEntryCStr = 0;
    samEntry->lenBuffULng = 0;
    blankSamEntry(samEntry);
} /*initSamEntry*/

/*######################################################################
# Output: Frees: allocated memory in samEntry
######################################################################*/
void freeStackSamEntry(
    struct samEntry *samEntry
) /*Frees heap allocations in a stack allocated sam entry struct*/
{ /*freeStackSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1 Sub-1 TOC: freeStackSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samEntry->samEntryCStr != 0)
        free(samEntry->samEntryCStr);

    return;
} /*freeStackSamEntry*/

/*######################################################################
# Output: Frees: samEntry & sets to 0
######################################################################*/
void freeHeapSamEntry(
    struct samEntry **samEntry
) /*Frees the samEntry structer*/
{ /*freeHeapSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1 Sub-1 TOC: freeHeapSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*samEntry)->samEntryCStr != 0)
        free((*samEntry)->samEntryCStr);

    free(*samEntry);
    *samEntry = 0;

    return;
} /*freeStackSamEntry*/

/*######################################################################
# Output:
#     Modifies: read id, reference id, q-score, cigar, & sequence
#               pointers to piont to their entires in the sam entry.
#               This also grabs the flag & mapq.
######################################################################*/
void processSamEntry(
    struct samEntry *samEntry /*Has sam file entry to find data in*/
) /*Sets Q-score, cigar, & sequence pionters in samEntry. Also finds
    the mapping quality, sam flag, & sequence length*/
{ /*processSamEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 TOC: processSamEntry
    # Fun-8 Sec-1: Variable declerations
    # Fun-8 Sec-2: Check if their is a sam entry
    # Fun-8 Sec-3: Extract data from the sam entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *samIterUChar = samEntry->samEntryCStr; /*iterator*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-2: Check if their is a sam entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samIterUChar == 0 || *samIterUChar == '@')
    { /*If their is nothing to do or is header entry*/
        samEntry->queryCStr = 0;
        samEntry->refCStr = 0;
        samEntry->cigarCStr = 0;
        samEntry->seqCStr = 0;
        samEntry->mapqUChar = 0;
        samEntry->qCStr = 0;
        samEntry->flagUSht = 0;

        return;
    } /*If their is nothing to do or is header entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-3: Extract data from the sam entry
    #    fun-8 sec-3 sub-1: Find non-sequence & q-score entries
    #    fun-8 sec-3 sub-2: sequence length, sequence entry & q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-8 Sec-3 Sub-1: Find non-sequence & q-score entries
    *******************************************************************/

    samEntry->queryCStr = samIterUChar; /*Story query id address*/

    /*Get the flag*/
    while(*samIterUChar > 32)            /*Move to flag*/
        samIterUChar++;
    samIterUChar++;                     /*Move of tab after query id*/
    samIterUChar = cStrToUSht(samIterUChar, &samEntry->flagUSht);

    /*Find the reference id*/
    samIterUChar++;                    /*Move off tab after flag*/
    samEntry->refCStr = samIterUChar;  /*Set pointer to reference id*/

    while(*samIterUChar > 32)         /*Move to tab after reference id*/
        samIterUChar++;
    samIterUChar++;                   /*Move of tab after reference id*/

    /*Find the starting position on the reference*/
    samIterUChar = cStrToUInt(samIterUChar, &samEntry->posOnRefUInt);
    samIterUChar++;                     /*Move of tab after position*/

    /*Find the mapq entry*/
    samIterUChar = cStrToUChar(samIterUChar, &samEntry->mapqUChar);

    /*Find the cigar entry*/
    samIterUChar++;                     /*Move of tab after mapq entry*/
    samEntry->cigarCStr = samIterUChar;

    /*******************************************************************
    # Fun-8 Sec-3 Sub-2: Find sequence length, sequence entry & q-score
    *******************************************************************/

    /*Find the sequence entry*/
    for(uint8_t uCharCnt = 0; uCharCnt < 4; ++uCharCnt)
    { /*Loop past cigar, RNEXT, PNEXT, & TLEN*/
        while(*samIterUChar > 32)        /*Move to next entry*/
            samIterUChar++;
        samIterUChar++;                 /*Move of tab after last entry*/
    } /*Loop past cigar, RNEXT, PNEXT, & TLEN*/

    samEntry->seqCStr = samIterUChar;   /*Recored sequence address*/
    samEntry->unTrimReadLenUInt = 0;
    samEntry->readLenUInt = 0;

    /*Find the sequence length & find q-score entry*/
    if(*samIterUChar != '*')
    { /*If this entry has a sequence*/
        while(*samIterUChar > 32)        /*Move to next entry*/
        { /*While I am not at the end of the sequence*/
            samIterUChar++;
            ++(samEntry->unTrimReadLenUInt);
        } /*While I am not at the end of the sequence*/

        samEntry->readLenUInt = samEntry->unTrimReadLenUInt;
    } /*If this entry has a sequence*/
    else
        ++samIterUChar;                 /*Move onto tab after entry*/

    samIterUChar++;                 /*Move of tab after last entry*/
    samEntry->qCStr = samIterUChar; /*Set q-score pionter*/

    return;
} /*processSamEntry*/

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
) /*Reads in sam entry & sets pointers in samStruct*/
{ /*readSamLine*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 TOC: readSamLine
    #     fun-9 sec-1: Variable declerations
    #     fun-9 sec-2: Check if have buffer to read into & prepare buffer
    #     fun-9 sec-3: Check if can read line in one go
    #     fun-9 sec-4: Check if is header & need to read in more
    #     fun-9 sec-5: If not read in cigar & find sequence length
    #     fun-9 sec-6: Check if have full sequence from cigar read in
    #     fun-9 sec-7: If not find sequence length
    #     fun-9 sec-8: Use sequence length to get sequence & q-score
    #     fun-9 sec-9: If not, read rest of cigar by adding 1kb chunks
    #     fun-9 sec-10: Set sequence & q-score pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *oldBufAddUChar = 0,
        *samIterUChar = 0;  /*Tells if fget error or end of file*/
    uint64_t
        newBuffSizeULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Check if have buffer to read into & prepare buffer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samStruct->samEntryCStr == 0)
    { /*If have no buffer*/
        samStruct->lenBuffULng = 1400;
        samStruct->samEntryCStr = malloc(sizeof(char) * 1401);
    } /*If have no buffer*/

    /*Set second to last character in buffer to null, so can tell if
      read entire line*/
    *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2) = '\0';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: Check if can read line in one go
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar =
        fgets(
            samStruct->samEntryCStr,
            samStruct->lenBuffULng,
            inFILE
    ); /*See if can read in a single line*/

    if(samIterUChar == 0)
        return 2;           /*At end of the file*/

    if(
        *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2)=='\0' ||
        *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2) =='\n'
    ) { /*If grabed the entire line*/
        processSamEntry(samStruct);
        return 1;
    } /*If grabed the entire line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-4: Check if is header & need to read in more
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*samStruct->samEntryCStr == '@')
    { /*If is a header line*/
        newBuffSizeULng = 1000; /*Should be close, read in small chunks*/
    
        while(
          !(*(samStruct->samEntryCStr+samStruct->lenBuffULng-2)=='\0' ||
           *(samStruct->samEntryCStr + samStruct->lenBuffULng -2)=='\n')
        ) { /*While have not read in the full line*/
            samIterUChar =
                readNextPartOfLine(inFILE, samStruct, newBuffSizeULng);
    
            if(samStruct->samEntryCStr == 0)
            { /*If memory allocation errory*/
                printMemAlocErr(
                    "samFileStructs.c",
                    "readSamLine",
                    10,
                    533
                ); /*Let user know of memory allocation failure*/
    
                return 64;
            } /*If memory allocation errory*/
        } /*While have not read in the full line*/

        return 1;
    } /*If is a header line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-5: If not read in cigar & get sequence length from cigar
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(findSamCig(inFILE, samStruct) == 4)
    { /*If memory allocation error or end of file*/
        if(samStruct->samEntryCStr == 0)
        { /*If memory allocation errory*/
            printMemAlocErr(
                "samFileStructs.c",
                "readSamLine",
                10,
                454
            ); /*Let user know of memory allocation failure*/

            return 64;
        } /*If memory allocation errory*/
    } /*If memory allocation error or end of file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-6: Check if have full sequence from cigar read in
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if have full sequence and finish setting pionters*/
    if(
        *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2)=='\0' ||
        *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2) =='\n'
    ) { /*If now have the entire line*/
        samIterUChar = samStruct->cigarCStr;

        for(uint32_t intTab = 0; intTab < 4; ++intTab)
        { /*Loop through cigar, RNEXT, PNEXT, & TLEN entires*/
            while(*samIterUChar != '\t')
                ++samIterUChar;
            ++samIterUChar; /*get off the tab*/
        } /*Loop through cigar, RNEXT, PNEXT, & TLEN entires*/

        samStruct->seqCStr = samIterUChar;
        samStruct->qCStr = samIterUChar + samStruct->readLenUInt + 1;
            /*reaLenUint to move past sequence; +1 to get off the tab*/
        return 1;
    } /*If now have the entire line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-7: If not find sequence length
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar = samStruct->cigarCStr;

    while(*samIterUChar != '\t')
    { /*While their are entries in the cigar*/
        ++samIterUChar;
        ++newBuffSizeULng; /*Find length of the cigar (avoid guessing)*/
    } /*While their are entries in the cigar*/

    /*MD entries will be bit longer than cigars, but can double cigar
      to get a rough size. This way I only have to read in once or
      twice*/
    newBuffSizeULng = newBuffSizeULng << 1;

    /*Figure out buffer size with sequence & q-score*/
    newBuffSizeULng +=                    /*newBuff has Cigar length*/
        (samStruct->readLenUInt << 1) +   /*Times by 2: seq + q-score*/
        samStruct->cigarCStr - samStruct->samEntryCStr + /*pre cigar*/
        2402;                             /*Extra buffer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-7: Use sequence length to get sequence & q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(newBuffSizeULng > samStruct->lenBuffULng)
    { /*If the new buffer size with the sequence is bigger*/
        newBuffSizeULng -= samStruct->lenBuffULng;
        oldBufAddUChar = samStruct->samEntryCStr; /*For resting point*/

        samIterUChar =
            readNextPartOfLine(inFILE, samStruct, newBuffSizeULng);

        if(samStruct->samEntryCStr == 0)
        { /*If memory allocation errory*/
            printMemAlocErr(
                "samFileStructs.c",
                "readSamLine",
                10,
                510
            ); /*Let user know of memory allocation failure*/

            return 64;
        } /*If memory allocation errory*/

        /*Readjust previous pointers*/
        samStruct->queryCStr = samStruct->samEntryCStr;
    
        samStruct->refCStr =
            samStruct->samEntryCStr +
            (samStruct->refCStr - oldBufAddUChar);

        samStruct->cigarCStr =
            samStruct->samEntryCStr +
            (samStruct->cigarCStr - oldBufAddUChar);
    } /*If the new buffer size with the sequence is bigger*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-9: If not, read rest of cigar by adding 1kb chunks
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newBuffSizeULng = 1000; /*Should be close, so read in small chunks*/

    while(
      !(*(samStruct->samEntryCStr + samStruct->lenBuffULng - 2)=='\0' ||
        *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2) =='\n')
    ) { /*While have not read in the full line*/
        oldBufAddUChar = samStruct->samEntryCStr;
        samIterUChar =
            readNextPartOfLine(inFILE, samStruct, newBuffSizeULng);

        if(samStruct->samEntryCStr == 0)
        { /*If memory allocation errory*/
            printMemAlocErr(
                "samFileStructs.c",
                "readSamLine",
                10,
                533
            ); /*Let user know of memory allocation failure*/

            return 64;
        } /*If memory allocation errory*/

        /*Readjust previous pointers*/
        samStruct->queryCStr = samStruct->samEntryCStr;
    
        samStruct->refCStr =
            samStruct->samEntryCStr +
            (samStruct->refCStr - oldBufAddUChar);

        samStruct->cigarCStr =
            samStruct->samEntryCStr +
            (samStruct->cigarCStr - oldBufAddUChar);
    } /*While have not read in the full line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-10: Set sequence & q-score pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Have entire entry, just need to find additional pointers*/
    samIterUChar = samStruct->cigarCStr;

    /*SEQ IS ON Q-SCORE, FIX*/
    for(uint32_t intTab = 0; intTab < 4; ++intTab)
    { /*Loop through cigar, RNEXT, PNEXT, & TLEN entires*/
        while(*samIterUChar != '\t')
            ++samIterUChar;
        ++samIterUChar; /*get off the tab*/
    } /*Loop through cigar, RNEXT, PNEXT, & TLEN entires*/

    samStruct->seqCStr = samIterUChar;

    if(*samStruct->seqCStr != '*')  /*if their is a sequence*/
    { /*If have a sequence entry*/
        if(*samStruct->cigarCStr != '*')
            samStruct->qCStr = samIterUChar + samStruct->readLenUInt +1;
        else
        { /*Else have no cigar entry, need to find length manually*/
            while(*samIterUChar > 32)
            { /*While have bases in the sequence*/
                ++samStruct->unTrimReadLenUInt;
                ++samIterUChar;
            } /*While have bases in the sequence*/

            samStruct->readLenUInt = samStruct->unTrimReadLenUInt;
            ++samIterUChar;  /*Get off tab after sequence*/
            samStruct->qCStr = samIterUChar;
        } /*Else have no cigar entry, need to find length manually*/
        /*reaLenUint to move past sequence; +1 to get off the tab*/
    } /*If have a sequence entry*/

    else
        samStruct->qCStr = samIterUChar + 2; /*no sequence present*/
   
    return 1;                    /*Read in full line, nothing to do*/
} /*readSamLine*/

/*######################################################################
# Output:
#     Returns:
#        - 1 if succeded
#        - 4 if end of file or memory allocation error
#     Modifies: samStruct->samEntryCStr, samStruct->lenBuffULng, &
#               samStruct->unTrimedLenULng
# Note:
#    - This requires that samStruct->samEntryCStr has read in at least
#      part of the sam entry
######################################################################*/
uint8_t findSamCig(
    FILE *inFILE,               /*File to search for cigar in*/
    struct samEntry *samStruct
) /*Finds the cigar in same entry (uses fgets if needs more buffer) &
    also finds the number of bases in query from cigar*/
{ /*findSamCig*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 TOC: findSamCig
    #    fun-10 sec-1: Variable declerations
    #    fun-10 sec-2: Set query pointer & read to start of flag
    #    fun-10 sec-3: Read in flag & move to start of reference id
    #    fun-10 sec-4: Set reference id pointer & move to position
    #    fun-10 sec-5: Read position & move to start of mapq
    #    fun-10 sec-6: Read in mapq & move to start of cigar
    #    fun-10 sec-7: Read in cigar
    #    fun-10 sec-8: Find length of sequence in the cigar
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *oldBufAddUChar = 0,    /*Used for refinding pointers*/
        *samIterUChar = samStruct->samEntryCStr;
    uint32_t
        cigEntryUInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-2: Set query pointer & read to start of flag
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->queryCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*While on the query id*/
        if(*samIterUChar == '\0')
        { /*If at end of buffer*/
            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
        } /*If at end of buffer*/

        ++samIterUChar;
    } /*While on the query id*/

    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-3: Read in flag & move to start of reference id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get how many characters the mapq is in the sam entry*/
    cigEntryUInt = samIterUChar - samStruct->samEntryCStr;

    /*Need to make sure entire flag is read in*/ 
    for(uint32_t intIter = 0; *(samIterUChar+intIter) != '\t';++intIter)
    { /*loop till at end of number*/
        if(*(samIterUChar + intIter) == '\0')
        { /*If at end of buffer*/
            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Make sure pionter is set start of flag*/
            samIterUChar = samStruct->samEntryCStr + cigEntryUInt;
            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
        } /*If at end of buffer*/
    } /*loop till at end of number*/

    samIterUChar = cStrToUSht(samIterUChar, &samStruct->flagUSht);
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-4: Set reference id pointer & move to position
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->refCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*While on the reference id*/
        if(*samIterUChar == '\0')
        { /*If at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If at end of buffer*/

        ++samIterUChar;
    } /*While on the reference id*/

    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-5: Get off position & move to start of mapq
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(uint32_t intIter = 0; *(samIterUChar+intIter) != '\t';++intIter)
    { /*Loop while not sure if read in entire positoin entry*/
        if(*samIterUChar == '\0')
        { /*If at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If at end of buffer*/
    } /*Loop while not sure if read in entire positoin entry*/

    /*Record the position entry*/
    samIterUChar = cStrToUInt(samIterUChar, &samStruct->posOnRefUInt);
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-6: Read in mapq & move to start of cigar
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get how many characters the mapq is in the sam entry*/
    cigEntryUInt = samIterUChar - samStruct->samEntryCStr;

    /*Need to make sure entire mapq is read in*/ 
    for(uint32_t intIter = 0; *(samIterUChar+intIter) != '\t';++intIter)
    { /*loop till at end of number*/
        if(*(samIterUChar + intIter) == '\0')
        { /*If at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Make sure pionter points to start of mapq*/
            samIterUChar = samStruct->samEntryCStr + cigEntryUInt;

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);
        } /*If at end of buffer*/
    } /*loop till at end of number*/

    samIterUChar = cStrToUChar(samIterUChar, &(samStruct->mapqUChar));
    ++samIterUChar;   /*Get off tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-7: Read in cigar
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->cigarCStr = samIterUChar;

    while(*samIterUChar != '\t')
    { /*While checking if cigar is compelete*/
        if(*samIterUChar == '\0')
        { /*If at end of buffer*/
            oldBufAddUChar = samStruct->samEntryCStr;

            samIterUChar =
                readNextPartOfLine(
                    inFILE,
                    samStruct,
                    samStruct->lenBuffULng /*Doubles buffer size*/
            ); /*Read in next section of line & double buffer size*/

            if(samIterUChar == 0)
                return 4;         /*incomplete sam aligment*/

            /*Readjust previous pointers*/
            samStruct->queryCStr = samStruct->samEntryCStr;
    
            samStruct->refCStr =
                samStruct->samEntryCStr +
                (samStruct->refCStr - oldBufAddUChar);

            samStruct->cigarCStr =
                samStruct->samEntryCStr +
                (samStruct->cigarCStr - oldBufAddUChar);
        } /*If at end of buffer*/

        ++samIterUChar;
    } /*While checking if cigar is compelete*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-8: Find length of sequence in the cigar
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samIterUChar = samStruct->cigarCStr;

    while(*samIterUChar != '\t')
    { /*While have cigar entries to read in*/
        samIterUChar = cStrToUInt(samIterUChar, &cigEntryUInt);

        if(
            *samIterUChar == 'X' || /*Mismatch --eqx format*/
            *samIterUChar == '=' || /*Match --eqx format*/
            *samIterUChar == 'M' || /*Mismatch or match (no --eqx)*/
            *samIterUChar == 'I' || /*Indel*/
            *samIterUChar == 'S' || /*Soft mask*/
            *samIterUChar == '\t'   /*Last number is a match (--eqx)*/
        )
            samStruct->unTrimReadLenUInt += cigEntryUInt;
        ++samIterUChar; /*Get off entry character type*/
    } /*While have cigar entries to read in*/

    samStruct->readLenUInt = samStruct->unTrimReadLenUInt;
    return 1;
} /*findSamCig*/

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
    uint64_t buffChangeULng         /*How much to increase buffer by*/
) /*Reads in next part of line when fgets did not get a full line*/
{ /*readNextPartOfLine*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-1 Sub-1 TOC: readNextPartOfLine
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpUCStr = 0;

    samStruct->samEntryCStr =
        realloc(
            samStruct->samEntryCStr, 
            sizeof(uint8_t) * (samStruct->lenBuffULng + buffChangeULng)
        ); /*Rellocate memory*/

    tmpUCStr = samStruct->samEntryCStr + samStruct->lenBuffULng - 1;
    samStruct->lenBuffULng += buffChangeULng;

    while(*tmpUCStr == '\0') /*Make sure on the first null*/
    { /*While on a null, loop till not on null so can find first null*/
        ++buffChangeULng;
        --tmpUCStr;
    } /*While on a null, loop till not on null so can find first null*/

    ++tmpUCStr;
    /*Not adjusting buff change because fgets reads in one off*/

    if(samStruct->samEntryCStr == 0)
        return 0;

    /*Set my marker for non-entire line read in*/
    *(samStruct->samEntryCStr + samStruct->lenBuffULng - 2) = '\0';

    return fgets(tmpUCStr, buffChangeULng, inFILE);
} /*readNextPartOfLine*/

/*######################################################################
# Output:
#     Prints: line with stats from samEntryStruct to file
#     Modifies: printHeaderChar to be 0 if set to 1
######################################################################*/
void printSamStats(
    struct samEntry *samEntryStruct, /*Sam entry to print stats for*/
    uint8_t *printHeaderChar,/*1: print header, 0: do not print header*/
    FILE *outFILE             /*File to output to*/
) /*Prints stats in samEntry structer for a sam file entry to tsv file*/
{ /*printSamStats*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-12 TOC: printSamStats
    #    fun-12 sec-1: Variable declerations
    #    fun-12 sec-2: Print out the header
    #    fun-12 sec-3: Convert query & reference ids to c-strings
    #    fun-12 sec-4: Print out stats from the samEntry structer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-12 Sec-1: Variable declerations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char
        *tmpQueryCStr = (*samEntryStruct).queryCStr,
        *tmpRefCStr = (*samEntryStruct).refCStr;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-12 Sec-2: Print out the header
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(*printHeaderChar == 1)
    { /*If wanted to print the header*/
        fprintf(
            outFILE,
            "Read\tRef\tMAPQ\treadLength\talignedLength\tmatches"
        );
        fprintf(
            outFILE,
            "\tkeptMatches\tmismatches\tinsertions\tdeletions"
        ); 
        fprintf(
            outFILE,
            "\tmedianQ\tmeanQ\talignedMedianQ\talignedMeanQ"
        ); 
        fprintf(
            outFILE,
            "\tTotalMismatches\tTotalInsertions\tTotalDeletions\n"
        ); 
        *printHeaderChar = 0;
    } /*If wanted to print the header*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-12 Sec-3: Convert query & reference ids to c-strings
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    while(*tmpQueryCStr != '\t')
        tmpQueryCStr++; /*move to end of query name*/
    *tmpQueryCStr = '\0'; /*turn into c-string*/

    while(*tmpRefCStr != '\t')
        tmpRefCStr++; /*move to end of reference name*/
    *tmpRefCStr = '\0'; /*turn into c-string*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-12 Sec-4: Print out stats from the samEntry structer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Print out the entry stats*/
    fprintf(
        outFILE,
        "%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%f\t%f\t%f",
        samEntryStruct->queryCStr, 
        samEntryStruct->refCStr,
        samEntryStruct->mapqUChar,
        samEntryStruct->readLenUInt,
        samEntryStruct->readAligLenUInt,
        samEntryStruct->numMatchUInt,
        samEntryStruct->numKeptMatchUInt,
        samEntryStruct->numKeptSNPUInt,
        samEntryStruct->numKeptInsUInt,
        samEntryStruct->numKeptDelUInt,
        samEntryStruct->medianQFlt,
        samEntryStruct->meanQFlt,
        samEntryStruct->medianAligQFlt,
        samEntryStruct->meanAligQFlt
    ); /*1st printf: print out stats*/

    fprintf(
        outFILE,
        "\t%u\t%u\t%u\n",
        samEntryStruct->numSNPUInt,
        samEntryStruct->numInsUInt,
        samEntryStruct->numDelUInt
    ); /*2nd fprintf: print out stats*/

    *tmpQueryCStr = '\t'; /*turn back into single string*/
    *tmpRefCStr = '\t'; /*turn back into singel string*/

    return;
} /*printSamStats*/

/*######################################################################
# Output:
#    File: fastq with the read id, sequence, & q-score entry of the sam
#                entry
######################################################################*/
void samToFq(
    struct samEntry *samStruct, /*Has sam entry to print as fastq*/
    FILE *outFILE               /*file to print fastq entry to*/
) /*Converts sam entry into a fastq entry & prints to a fastq file*/
{ /*samToFq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 TOC: samToFq
    #    fun-13 sec-1: Variable declerations
    #    fun-13 sec-2: Find length of query id
    #    fun-13 sec-3: Print sam entry as fastq entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpQueryCStr = samStruct->queryCStr;
    uint16_t lenQueryUSht = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-2: Find length of query id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpQueryCStr > 32)
    { /*wile not at end of query id*/
        tmpQueryCStr++; /*move to end of query name*/
        lenQueryUSht++; /*Find length of query id*/
    } /*wile not at end of query id*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-3: Print sam entry as fastq entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fwrite("@", sizeof(uint8_t), 1, outFILE); /*Print out @ header*/
    /*print out the query id*/
    fwrite(samStruct->queryCStr, sizeof(uint8_t), lenQueryUSht,outFILE);

    /*Write new line after query id*/
    fwrite("\n", sizeof(uint8_t), 1, outFILE);

    fwrite(
        samStruct->seqCStr,
        sizeof(uint8_t),
        samStruct->readLenUInt,
        outFILE
    ); /*Write sequence entry*/

    fwrite("\n+\n", sizeof(uint8_t), 3, outFILE); /*Write spacer entry*/

    fwrite(
        samStruct->qCStr,
        sizeof(uint8_t),
        samStruct->readLenUInt,
        outFILE
    ); /*Write Q-score entry*/

    /*Write new line for next entry*/
    fwrite("\n", sizeof(uint8_t), 1, outFILE);

    return;
} /*samToFq*/

/*######################################################################
# Output:
#    Prints: Sam file entry in samStruct to outFILE.
######################################################################*/
void printSamEntry(
    struct samEntry *samStruct, /*Has sam entry to print out*/
    FILE *outFILE               /*File to print sam entry to*/
) /*Prints out the sam entry in the samStruct, does not print stats*/
{ /*printSamEntry*/
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-14 Sec-1 Sub-1 TOC: printSamEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*printing out query name since is start of sam entry. This allows
      me to store the sam entry in a different samStruct, but still
      have a handle on the entire entry with the query. I do this when
      their might be multiple entries for a single read. In this case
      only the first entry has the sequence.*/
    if(samStruct->queryCStr != 0)
        fprintf(outFILE, "%s\n", samStruct->queryCStr);
    else
        fprintf(outFILE, "%s\n", samStruct->samEntryCStr); /*header*/
} /*printSamEntry*/
