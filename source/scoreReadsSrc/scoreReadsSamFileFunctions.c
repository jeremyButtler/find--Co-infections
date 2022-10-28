/*##############################################################################
# Name: ScoreReadsSamFileFunctions
# Use: Has functions to score reads, including median Q-score functions
##############################################################################*/

#include "scoreReadsSamFileFunctions.h"

/*##############################################################################
# TOC:
#   fun-1: indexSamEntry: finds mapq, flag, and pointer locations in sam entry
#   fun-2: scoreReads: Scores a single read based on sam file entry
#   fun-3: refScoreReads: Scores a single read based using sam entry & refernce
#   fun-4: qHistToMedian: Convert histogram of Q-scores to median
#   fun-5: findQScores: Find mean and median Q-scores
#   fun-6: readCigEntry: Reads a single entry form a cigar
##############################################################################*/

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
| 11  | QUAL  | String | ASCII of Phred-scaled base QUALity+33 |
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

/*##############################################################################
# Name: indexSamEntry
# Use: sets pointers to the query, reference, cigar, sequence, and q-score lines.
#      Also extracts the flag and mapq.
# Input:
#    entryCStr: start of the same file entry (c-string)
#    samEntryStruct: Structer to store flags, mapq, and pionters in (samEntry)
# Output:
#    modfies: samEntry to have the starting locations, mapq, and flag
#    Returns: 1 if indexed, 0 or -1 otherwise
# Note: will not set Q-score and sequence pointers when '*' (nothing)
##############################################################################*/
char indexSamEntry(char *entryCStr, samEntry *samEntryStruct)
{ /*indexSamEntry*/
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declarations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *incCStr = entryCStr; /*loop through C-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Find the entry lines
    #     fun-1 sec-2: Find the flag line (tab 2)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-2: Find the flag line (tab 2)
    ***************************************************************************/

    (*samEntryStruct).queryCStr = entryCStr; /*query first entry in samfile entry*/

    /*find the flag header*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;

    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the flag*/

    if(cStrToUInt(incCStr, &((*samEntryStruct).flagUInt), '\t') < 1)
        return -1; /*No flag, so likely not a sam file*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-2: Find the reference line (tab 3)
    ***************************************************************************/

    /*move to the reference line*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;

    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the reference*/
    (*samEntryStruct).refCStr = incCStr; /*move off the tab to the reference*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-3: Find the mapq line (tab 5)
    ***************************************************************************/

    /*Find the position line*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the position entry*/

    /*Find the mapq entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the mapq entry*/

    if(cStrToUInt(incCStr, &((*samEntryStruct).mapqUInt),'\t')<1 && *incCStr!='*')
        return -1; /*No flag, so likely not a sam file*/
    else if(*incCStr == '*')
        (*samEntryStruct).mapqUInt = 0; /*Not present, likely an unmaped read*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-4: Find the cigar line (tab 6)
    ***************************************************************************/

    /*Find the cigar entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the cigar entry*/
    (*samEntryStruct).cigarCStr = incCStr; /*Recored the cigar location*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-5: Find the sequence line (tab 10)
    ***************************************************************************/

    /*Find the RNEXT (ref name of mate read) entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the RNEXt entry*/

    /*Find the PNEXT (Position of of mate read) entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the PNEXt entry*/

    /*Find the TLEN (Observed template length) entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/
    incCStr++; /*move off the tab to the TLEN entry*/

    /*Find the sequence entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/

    incCStr++;  /*move off the tab to the sequence entry*/
    (*samEntryStruct).seqCStr = incCStr; /*recored the sequence start*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-5: Find the Q-score line (tab 11)
    ***************************************************************************/

    /*Find the Q-score entry*/
    while(*incCStr != '\t' && *incCStr != '\0')
        incCStr++;
    if(*incCStr == '\0')
        return -1; /*at end of file*/

    incCStr++; /*move off the tab to the q-score entry*/
    (*samEntryStruct).qCStr = incCStr; /*recored the q-score start*/

    return 1; /*Noting else to do*/
} /*indexSamEntry*/

/*##############################################################################
# Name: scoreRead
# Use: Finds the similarity between two reads, aligned read length.
#      If Q-score entry returns the aligned mean and median Q-score and
# Input:
#    minReadStats.minQChar:
#        - min Q-Score to use (should be offset by 30) (character)
#    samEntry:
#        - Holds the pointers for the cigar, sequence, and q-score entries
#        - Will have mismatches and indels input 
# Output:
#    modifies: samEntry struct to have the mismatches and indels
##############################################################################*/
void scoreRead(minStats *minReadStats, samEntry *samEntryStruct)
{ /*scoreReads*/
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: declare variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *cigarCStr = (*samEntryStruct).cigarCStr;    /*start of cigar entry*/
    scoreReadStruct seqReadStruct;          /*holds sequence information*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Score reads and find aligned read length
    #    fun-2 sec-2 sub-1: Blank the Q-score array
    #    fun-2 sec-2 sub-2: Get the number bases from the cigar entry
    #    fun-2 sec-2 sub-3: Check cigar entry type (match, mismatch, indel, ...)
    #    fun-2 sec-2 sub-4: Find aligned mean and median Q-scores and return
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-2 Sec-2 Sub-1: Blank the Q-score array
    ***************************************************************************/

    initScoreReadStruct(&seqReadStruct, samEntryStruct);

    /***************************************************************************
    # Fun-2 Sec-2 Sub-2: Get the number of bases from the cigar entry
    ***************************************************************************/

    while(*cigarCStr != '\t')
    { /*loop through all entries in the Cigar*/

        seqReadStruct.cigarEntryULng = readCigEntry(&cigarCStr, 0);

        /***********************************************************************
        # Fun-2 Sec-2 Sub-3: Check cigar entry type (match,mismatch, indel, ...)
        ***********************************************************************/

        if(*cigarCStr == 'X')
            checkMismatches(&seqReadStruct, minReadStats, samEntryStruct);

        else if(*cigarCStr == 'I')
            checkInsertion(&seqReadStruct, minReadStats, samEntryStruct);

        else if(*cigarCStr == 'D')
            checkDeletions(&seqReadStruct, minReadStats, samEntryStruct);

        else if(*cigarCStr == '=')
            checkMatches(&seqReadStruct, minReadStats, samEntryStruct);

        else if(*cigarCStr == 'S')  /*was soft mask*/
            checkSoftMasks(&seqReadStruct);

        else if(*cigarCStr == '\t') /*End of cigar is a match*/
            checkMatches(&seqReadStruct, minReadStats, samEntryStruct);

        /*If did not hit if statment is a Hard mask or unknown entry, ignore*/
        cigarCStr++; /*move to the next cigar eqx entry*/
    } /*loop through all entries in the Cigar*/

    /***************************************************************************
    # Fun-2 Sec-2 Sub-4: Find aligned mean and median Q-scores and return
    ***************************************************************************/

    (*samEntryStruct).meanAligQDbl = seqReadStruct.totalQScoreULng /
                                     (*samEntryStruct).readAligLenULng;
    (*samEntryStruct).medianAligQDbl = qHistToMedian(seqReadStruct.seqQHistULng,
                                             (*samEntryStruct).readAligLenULng);
    return;
} /*scoreReads*/

/*##############################################################################
# Output: Modifies: samEntryStruct to hold aligned length, Q-scores,
#                   mismatches, & indels
##############################################################################*/
void refScoreRead(
    struct minStats *minReadStats,   /*Min requirments to keep aligment*/
    struct samEntry *samEntryStruct, /*sequence & Q-score line & holds stats*/
    struct samEntry *refEntry,       /*sequence & Q-score line for reference*/
    char refForDelChar              /*Marks if only using ref for deltions*/
) /*Find similarity between two reads, aligned read length, & aligned Q-scores*/
{ /*refScoreReads*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: declare variables
    #    fun-3 sec-2: Blank the Q-score array & handle reverse complement
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: declare variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *cigarCStr = samEntryStruct->cigarCStr,  /*start of cigar entry*/
        *refSeqCStr = refEntry->seqCStr,         /*Start of reference sequence*/
        *refQCStr = refEntry->qCStr;             /*Start of reference Q-score*/
    struct scoreReadStruct seqReadStruct;        /*holds sequence information*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Blank the Q-score array & handle reverse complement
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initScoreReadStruct(&seqReadStruct, samEntryStruct);

    if(samEntryStruct->flagUInt & 16)
    { /*If is a reverse complement alignment (need to work backwards)*/
        while(*refSeqCStr != '\0')              /*go to end of sequence line*/
            refSeqCStr++;

        while(*refQCStr != '\0')                /*go to end of Q-score line*/
            refQCStr++;

        refSeqCStr--;                           /*Set on last base (not null)*/
        refQCStr--;                             /*Set on last base (not null)*/
    } /*If is a reverse complement alignment (need to work backwards)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Score reads and find aligned read length
    #    fun-3 sec-2 sub-1: Read in cigar entry & check mismatches
    #    fun-3 sec-3 sub-2: Check insertions
    #    fun-3 sec-3 sub-3: Check deletions
    #    fun-3 sec-3 sub-4: Check matches
    #    fun-3 sec-3 sub-5: Check softmasks
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-3 Sec-3 Sub-1: Read in cigar entry & check mismatches
    ***************************************************************************/

    while(*cigarCStr != '\t')
    { /*loop through all entries in the Cigar*/

        seqReadStruct.cigarEntryULng = readCigEntry(&cigarCStr, 0);

        if(*cigarCStr == 'X')
        { /*If the cigar entry was a mismatch*/

            if(refForDelChar == 1)
            { /*If not using refrence to check for mismatches*/
                checkMismatches(
                    &seqReadStruct,
                    minReadStats,
                    samEntryStruct
                ); /*Check if should keep mismatches*/
            } /*If not using refrence to check for mismatches*/
 
            else
            { /*Else using reference to check for mismatches*/
                refCheckMismatches(
                    &seqReadStruct,
                    &refSeqCStr,
                    &refQCStr,
                    minReadStats,
                    samEntryStruct
                ); /*Check if should keep mismatches*/
            } /*Else using reference to check for mismatches*/
        } /*If the cigar entry was a mismatch*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-2: Check insertions
        ***********************************************************************/

        else if(*cigarCStr == 'I')
            checkInsertion(
                &seqReadStruct,
                minReadStats,
                samEntryStruct
            ); /*Check if should keep insertions*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-3: Check deletions
        ***********************************************************************/

        else if(*cigarCStr == 'D')
            refCheckDeletions(
                &seqReadStruct,
                &refSeqCStr,
                &refQCStr,
                minReadStats,
                samEntryStruct
            ); /*Check if should keep deletion*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-4: Check mathces
        ***********************************************************************/

        else if(*cigarCStr == '=' || *cigarCStr == '\t')
        { /*Else if was a match (\t is end of sequence*/

            if(refForDelChar == 1)
            { /*If not using reference to check for matches*/
                checkMatches(
                    &seqReadStruct,
                    minReadStats,
                    samEntryStruct
                ); /*Check if should keep the match*/
            } /*If not using reference to check for matches*/

            else
            { /*Else using reference to check for matches*/
                refCheckMatches(
                    &seqReadStruct,
                    &refSeqCStr,
                    &refQCStr,
                    minReadStats,
                    samEntryStruct
                ); /*Check if should keep the match*/
            } /*Else using reference to check for matches*/
        } /*Else if was a match (\t is end of sequence*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-5: Check softmatches
        ***********************************************************************/

        else if(*cigarCStr == 'S')  /*was soft mask (no mapped bases in ref)*/
            checkSoftMasks(&seqReadStruct);

        /*If did not hit if statment is a Hard mask or unknown entry, ignore*/
        cigarCStr++; /*move to the next cigar eqx entry*/
    } /*loop through all entries in the Cigar*/

    /***************************************************************************
    # Fun-3 Sec-2 Sub-4: Find aligned mean and median Q-scores and return
    ***************************************************************************/

    samEntryStruct->meanAligQDbl = seqReadStruct.totalQScoreULng /
                                   samEntryStruct->readAligLenULng;
    samEntryStruct->medianAligQDbl = qHistToMedian(seqReadStruct.seqQHistULng,
                                                 samEntryStruct->readAligLenULng
    ); /*Get median Q-score*/

    return;
} /*scoreReads*/

/*##############################################################################
# Name: qHistToMedian
# Use: Takes in a histogram of Q-scores and finds the median
# Input:
#    qHistULng: Histogram holding q-scores (unsigned long array)
#    readLenULng: Length of the read [find midpoint with] (unsinged long)
# Output:
#    returns: double with the median Q-score [or 0 if nothing]
##############################################################################*/
double qHistToMedian(unsigned long qHistULng[],
                     unsigned long readLenULng)
{ /*qHistToMedian*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long *incULng = 0, curPosULng = 0, midPointULng = readLenULng / 2;
    double medianQDbl;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Find the median Q-score from the histogram
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(int intQ = 0; intQ < MAX_Q_SCORE_CHAR; intQ++)
    { /*Loop through the Q-score histogram to find mean and median Q-score*/
        curPosULng = curPosULng + qHistULng[intQ];

        if(curPosULng >= midPointULng)
        { /*if found the midpoint, then find the median*/

            if(curPosULng > midPointULng)
                medianQDbl = intQ;  /*2 postions have same Q past midpoint*/
            else
            { /*else if no clear midpiont*/
                incULng = qHistULng + intQ + 1; /*Next Q-score in histo*/
                midPointULng = 1; /*Difference between two Q-scores*/

                while(*incULng == 0)
                { /*loop till found the next Q-score*/
                    incULng++;
                    midPointULng++; /*How much greater is the next Q-score*/
                } /*loop till found the next Q-score*/

                if(readLenULng % 2 == 0) /*not even break, find mean*/
                    medianQDbl = intQ + (midPointULng / 2); /*difference / 2*/
                else  /*Even break (odd #), answer floored*/
                    medianQDbl = intQ + midPointULng; 
            } /*else if no clear midpiont*/

            return medianQDbl; /*Found median Q-score*/
        } /*if found the midpoint, then find the median*/
    } /*Loop through the Q-score histogram to find mean and median Q-score*/

    return 0; /*Just in case was all 0's in the array*/
} /*qHistToMedian*/

/*##############################################################################
# Name: findQScores
# Use: Finds the mean and median q-score for a q-score entry
# Input:
#    qCStr: q-score entry to find the mean and median q-score of (c-string)
#    meanQDbl: Holds the mean Q-score (double pointer)
#    medianQDbl: Holds the median Q-score (double pointer)
# output:
#    returns: The read length (unsigned long)
##############################################################################*/
unsigned long findQScores(char *qCStr, double *meanQDbl, double *medianQDbl)
{ /*findQScores*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *seqQCStr = qCStr;

    unsigned long qTotalULng = 0,
                  qScoreHistULng[MAX_Q_SCORE_CHAR],
                  readLenULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Find median and mean Q-scores
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*seqQCStr == 0)
        return 0;       /*No Q-score entry*/
    if(*seqQCStr == '*' && *(seqQCStr + 1) == '\t')
        return 0;       /*No Q-score entry*/

    /*Fill my Q-score array with 0's first*/
    for(int intQ = 0; intQ < MAX_Q_SCORE_CHAR; intQ++)
        qScoreHistULng[intQ] = 0;
    
    while(*seqQCStr != '\t' && *seqQCStr != '\n' && *seqQCStr != '\0')
    { /*loop through the Q-score entry*/
        qScoreHistULng[*seqQCStr - Q_ADJUST_CHAR]++; /*Build histogram*/
        qTotalULng = qTotalULng + (*seqQCStr - Q_ADJUST_CHAR); /*total score*/
        readLenULng++;
        seqQCStr++;
    } /*loop through the Q-score entry*/
    
    /*Find the mean and median*/
    *meanQDbl = qTotalULng / readLenULng;
    *medianQDbl = qHistToMedian(qScoreHistULng, readLenULng);

    return readLenULng;
} /*findQScores*/

/*##############################################################################
# Output:
#    returns: unsigned long with number of bases
#    modifies: cigarCStr to point to the entry type (mismatch, indel, ect...)
##############################################################################*/
unsigned long readCigEntry(
    char **cigarCStr,       /*C-string with cigar to read and incurment*/
    unsigned int flagUInt   /*Flag telling if is reverse sequence (16)*/
) /*Reads a single entry from a eqx cigar line*/
{ /*readCigEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6: readCigEntry
    #    fun-6 sec-1: variable declerations
    #    fun-6 sec-2: Get a single cigar entry for reverse complement entry
    #    fun-6 sec-3: Get a single cigar entry for regular
    #    fun-6 sec-4: Retrun the number of digits in the cigar entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char cigarEntryCStr[19],        /*holds one cigar entry*/
         holdChar;                  /*Holds cigar digit when reversing number*/
    int numDigitsInt = 0,           /*Number of digits in cigar entry*/
        stopInt = 0;                /*When to stop a loop*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Get a single cigar entry for reverse complement entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(flagUInt & 16)
    { /* if is a reverse complement alignment, read ref cigar backwards*/

        while(**cigarCStr > 47 && **cigarCStr < 59)
        { /*Loop till of number entry of cigar*/
            cigarEntryCStr[numDigitsInt] = **cigarCStr;
            (*cigarCStr)--;
            numDigitsInt++;
        } /*Loop till of number entry of cigar*/

        cigarEntryCStr[numDigitsInt] = '\0'; /*Make into c-string*/

        numDigitsInt--; /*get of null*/
        stopInt = numDigitsInt / 2; /*When to stopp swapping*/

        for(int intPow = 0; intPow < stopInt; intPow++)
        { /*Loop number is backward, flip to forwards*/
            holdChar = cigarEntryCStr[numDigitsInt];
            cigarEntryCStr[numDigitsInt] = cigarEntryCStr[intPow];
            cigarEntryCStr[intPow] = holdChar;
            numDigitsInt--; 
        } /*Loop number is backward, flip to forwards*/

    } /* if is a reverse complement alignment, read ref cigar backwards*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Get a single cigar entry for regular
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    else
    { /*else is normal alignment, read forward*/

        while(**cigarCStr > 47 && **cigarCStr < 59)
        { /*Loop till of number entry of cigar*/
            cigarEntryCStr[numDigitsInt] = **cigarCStr;
            (*cigarCStr)++;
            numDigitsInt++;
        } /*Loop till of number entry of cigar*/

        cigarEntryCStr[numDigitsInt] = '\0'; /*Make into c-string*/
    } /*else is normal alignment, read forward*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-4: Retrun the number of digits in the cigar entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return strtoul(cigarEntryCStr, NULL, 10); /*cStr->Ulong*/
} /*readCigEntry*/
