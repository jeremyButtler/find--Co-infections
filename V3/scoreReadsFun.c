/*######################################################################
# Name: scoreReadsFun
# Use:
#    - Holds unique functions needed for scoreReads.
#    - Use scoreAln to score a single alignment
# Includes:
#    - "samEntryStruct.h"
#        - <stdlib.h>
#        - "cStrToNumberFun.h"
#            - <sdtint.h>
#        - "printError.h"
#            - <stdio.h>
######################################################################*/

#include "scoreReadsFun.h" /*Holds structers to hold samfile stats*/

/*######################################################################
# TOC:
#    fun-1 scoreReads: Scores all alignments in a sam file
#        - Reads alignment and calls setUpScoreAln
#    fun-2 setUpScoreAln: Checks to see if should keep sam alignment and
#        - calls scoreAln
#    fun-3 scoreAln: Scores a single alignment in a sam file
#    fun-4 checkSNPs: Check if want to keep or ignore a mismatch
#    fun-5 checkMatches: Check if matches should be kept
#    fun-6 checkInss:
#        - Checks if insertion should be kept or if specified checks
#          if deletion is valid
#    fun-7 checkDels: Checks if deletions should be kept
#    fun-8 checkSoftMasks: Increments pointer when soft mask is present
#    fun-9 readAndCheckRef: Read in reference (calls readRefFqSeq)and
#                            check if reference meets requirements
#    fun-10 readRefFqSeq: Read in refence from fastq file
#    fun-11 findQScores: find Q-score of read using q-score entry
#    fun-12 qHistToMed: Find median Q-score from a q-score histogram
#    fun-13 readCigEntry: Get number of bases in a single cigar entry
#    fun-14 readReverseCigEntry: readCigEntry, but works backwards
#    fun-15 blankMinStats: set a minAlnStats structer to defaults
######################################################################*/

static uint32_t MAX_UINT = 0xFFFFFFFF; /*max value of 32 bit number*/

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
) /*Scores all alignments in a sam file*/
{ /*scoreReads*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1: scoreReads
    #    fun-1 sec-1: Declare variables & initalize structers
    #    fun-1 sec-2: Initalize the samEntry structers & read first line
    #    fun-1 sec-3: Read in reference from fastq file if provided
    #    fun-1 sec-4: Call the scoring function & clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Declare variables & initalize structers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t
        errUChar = 0,              /*Holds error output from functions*/
        printHeadChar = 1,           /*Tells to print out header*/
        refQBool = 0,
        keepAlnUChar = 0;

    /*Structers are from scoreReadsStructers.h*/
    struct samEntry
        *samStruct = malloc(sizeof(struct samEntry)),
        *refStruct =  malloc(sizeof(struct samEntry)),
        *oldSamStruct = malloc(sizeof(struct samEntry));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Initalize the samEntry structers & read first line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initSamEntry(samStruct);    /*Blank variables*/
    initSamEntry(oldSamStruct); /*Blank variables*/
    initSamEntry(refStruct); /*Blank variables*/

    keepAlnUChar = readSamLine(samStruct, samFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Read in reference from fastq file if provided
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errUChar = readAndCheckRef(refFILE, refStruct, minStats);

    switch(errUChar)
    { /*check if reference was provided*/
        case 0:                              /*No reference provided*/
             refQBool = 0;
             freeHeapSamEntry(&refStruct);   /*Pointer Is set to 0*/
             break;
        case 1:                             /*Valid reference*/
            refQBool = 1; /*This case will be a q-score entry*/
            break;
         default:                           /*Issue with reference*/
            freeHeapSamEntry(&refStruct);
            return errUChar; /*2 = under min quality, 4 = invalid file*/
    } /*check if reference was provided*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Call the scoring function & clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(keepAlnUChar & 1)
    { /*While I have lines to read in the sam file*/
         if(
             setUpScoreAln(
                 &samStruct,
                 &oldSamStruct,
                 refStruct,    /*pointer to reference*/
                 minStats,
                 useRefForDelBool, /*Use reference for only deletions*/
                 &refQBool        /*telling if ref has Q-score*/
            ) & 1)
                printSamStats(samStruct, &printHeadChar, outFILE);
                /*Logic is >> 2 removes first two flags*/
         keepAlnUChar = readSamLine(samStruct, samFILE);
    } /*While I have lines to read in the sam file*/

    freeHeapSamEntry(&samStruct);
    freeHeapSamEntry(&oldSamStruct);

    if(refStruct != 0)
        freeHeapSamEntry(&refStruct);

    if(keepAlnUChar & 64)
        return keepAlnUChar;
    return 1;
} /*scoreReads*/

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
)/*scores & finds if should keep an sam aligment*/
{ /*setUpScoreAln*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC: scoreAln
    #    fun-2 sec-1: Check if is a valid aligment & if new read
    #    fun-2 sec-2: check if alignment meets min thresholds
    #    fun-2 sec-3: Score alignment & see if meets min aligned stats
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Check if is a valid aligment & if new read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct samEntry *swapStruct = 0;

    if(*(*samStruct)->samEntryCStr == '@')
        return 2; /*This is a header line*/

    if(*(*samStruct)->refCStr == '*')
        return 8;

    blankReadStats(*samStruct); /*Remove old stats in sam file*/

    if(*(*samStruct)->seqCStr != '*')
    { /*If starting a new sequence*/
        findQScores(*samStruct);

        /*Swap pointers (this avoids doing a deep copy*/
        swapStruct = *samStruct;
        *samStruct = *oldSamStruct;
        *oldSamStruct = swapStruct;

        /*Copy over some variables not copied in my copy function*/
        (*samStruct)->flagUSht = (*oldSamStruct)->flagUSht;
        (*samStruct)->cigarCStr = (*oldSamStruct)->cigarCStr;
        (*samStruct)->mapqUChar = (*oldSamStruct)->mapqUChar;
        (*samStruct)->refCStr = (*oldSamStruct)->refCStr;
        (*samStruct)->queryCStr = (*oldSamStruct)->queryCStr;
        (*samStruct)->posOnRefUInt = (*oldSamStruct)->posOnRefUInt;

        /*Not setting samEntryCStr, because that would result in memory
          loss or require a deep copy. Instead working with query, ref,
          cigar, sequence (set in cpSamEntry), & q-score (set in 
          cpSamEntry) pointers.*/
    } /*If starting a new sequence*/

    cpSamEntry(*oldSamStruct, *samStruct);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: check if alignment meets min thresholds
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((*samStruct)->flagUSht & 4)
         return 4; /*move onto the next entry (not mapped)*/

    if((*oldSamStruct)->mapqUChar < minStats->minMapqUInt)
        return 4; /*move onto the next entry (has to low mapq)*/
            /*Minimap2 only gives MAPQ for best match*/

    if((*samStruct)->medianQFlt < minStats->minMedianQFlt &&
       *(*samStruct)->qCStr !='*')
        return 4; /*move onto the next entry (median Q-score to low)*/

    if((*samStruct)->meanQFlt < minStats->minMeanQFlt &&
       *(*samStruct)->qCStr !='*')
        return 4; /*move onto the next entry (mean Q-score to low)*/

    if((*samStruct)->readLenUInt > minStats->maxReadLenULng &&
       minStats->maxReadLenULng > 0)
        return 4; /*move onto the next entry (read is to long)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Score alignment & see if meets min aligned thresholds
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    scoreAln(
        minStats,
        *samStruct,          /*Holds alignment entry*/
        refStruct,          /*Holds reference information*/
        refQBool,           /*0 no q-score entry; 1 present*/
        useRefForDelBool    /*use reference only for deletions*/
    );

    if((*samStruct)->medianAligQFlt < minStats->minAlignedMedianQFlt &&
       *(*samStruct)->qCStr !='*')
        return 4; /*aligned median Q-score low*/

    if((*samStruct)->meanAligQFlt < minStats->minAlignedMeanQFlt &&
       *(*samStruct)->qCStr != '*')
        return 4; /*aligned mean Q-score to low*/

    if((*samStruct)->readAligLenUInt < minStats->minReadLenULng)
        return 4; /*Read to to short*/

    return 1;
} /*setUpScoreAln*/

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
) /*Scores a signle alignment in a sam file*/
{ /*scoreAln*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: declare variables
    #    fun-3 sec-2: Check if reading backwards & if have q-score entry
    #    fun-3 sec-3: Process each cigar entry
    #    fun-3 sec-4: Find aligned median & means
    #    fun-3 sec-5: Clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: declare variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t
        qLineBool = 0,
        nullUChar = 0,
        oneUChar = 1,
        useRefBool = !*useRefForDelBool; /*If using ref for non del*/
        /*reverseBool = 0, *//*Uncomment for reverse complement*/

    char
        *refStartCStr = 0,
        *refQCStr = 0,
        *tmpCigCStr = 0,                  /*Holds errpr type*/
        *cigarCStr = samStruct->cigarCStr;

    int32_t
        intOne = 1,  /*for checkInss with deletions*/
        intCnt = 0;  /*Incurment sequence & qscore by*/
    uint32_t cigEntryUInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Check if reading backwards & if have q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refStruct != 0)
    { /*If have a reference sequence*/
        /*recored starting sequence positions so can reset at end*/
        refStartCStr = refStruct->seqCStr;
        refQCStr = refStruct->qCStr;

        /*Adjust reference position for read position on reference*/
        refStruct->seqCStr += samStruct->posOnRefUInt - 1;
        refStruct->qCStr += samStruct->posOnRefUInt - 1;
        /*-1 to account for base 0 being base 1*/
    } /*If have a reference sequence*/

    /*Check if read has a q-score entry*/
    qLineBool =
       !!((*samStruct->qCStr ^ '*') | (*(samStruct->qCStr + 1) ^ '\t'));
        /* qCStr ^ '*' sets qCStr to 0 only when qCStr == '*'
           qCStr ^ '\t' sets qCStr to 0 only when qCStr == '\t'
           | ensures that qLineBool is 0 only if qCStr == "*\t"
           !! ensures that numbers greater than 1 are converted to 1
        */
 
    /*Note I need to set up the sam system to handle reverse complement
      for reads without a sequence (alignment with seq is normal, while
      other alignments are reverense complemnt (16), or vise versa, 
      this code is a start, but not quite their. It needs to check for
      the initial seq alignment. Also, minimap2 automatically corrects
    */
    /*if(samStruct->flagUSht & 16)
    {*/ /*if working on a reverse complement read*/
    /*    while(*cigarCStr != '\t')
            ++cigarCStr;

        --cigarCStr; *//*Get off the tab*/
        /*Move to last base*/
        /*samStruct->seqCStr += samStruct->readLenUInt - 1;
        samStruct->qCStr += samStruct->readLenUInt - 1;
    }*/ /*if working on a reverse complement read*/

    /*-1 (for incurmenting backwards) or 1 for incurmenting forwards*/
    /*if(samStruct->flagUSht & 16)
    {*/ /*If the sequqnce is reverse complement*/
        /*intCnt = -1;
        reverseBool = 1;
    }*/ /*If the sequqnce is reverse complement*/

    /*else
    {*/ /*Else the sequence is not reverse complement*/
        intCnt = 1;
        /*reverseBool = 0;*/ /*Uncomment if need to reverse complent*/
    /*}*/ /*Else the sequence is not reverse complement*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Process each cigar entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*cigarCStr != '\t')
    { /*loop through all entries in the Cigar*/

        /*switch(reverseBool)
        {*/ /*Switch: check if reading cigar backwards or forwards*/
          /*  case 0:*/
                readCigEntry(&cigarCStr, &cigEntryUInt);
                tmpCigCStr = cigarCStr;
                ++cigarCStr; /*move to the next cigar eqx entry*/
                /*break;*/
            /*case 1:*/                   /*Is reverse complement*/
                /*tmpCigCStr = cigarCStr;
                readReverseCigEntry(&cigarCStr, &cigEntryUInt);*/
                /*readReveres puts on next entry*/
                /*break;
        }*/ /*Switch: check if reading cigar backwards or forwards*/

        switch(*tmpCigCStr)
        { /*switch, check what kind of cigar entry*/
            case '=':
                checkMatches(
                    &cigEntryUInt,
                    minStats,
                    samStruct,
                    refStruct,
                    &qLineBool,
                    refQBool,
                    &useRefBool, /*Tells if using reference for check*/
                    &intCnt
                ); /*Check if keeping the matches*/
                break;
            case 'X':
                checkSNPs(
                    &cigEntryUInt,
                    minStats,
                    samStruct,
                    refStruct,
                    &qLineBool,
                    refQBool,
                    &useRefBool, /*Tells if using reference for check*/
                    &intCnt
                ); /*Check if keeping the SNPs*/
                break;
            case 'I':
                checkInss(
                    &cigEntryUInt,
                    minStats,
                    samStruct,
                    &qLineBool,
                    &nullUChar,   /*Note that this is Not a deletion*/
                    &intCnt
                ); /*Check if should keep the insertions*/
                break;
            case 'D':
            { /*Case: deletion*/
                if(refStruct == 0)
                    checkDels( &cigEntryUInt, minStats, samStruct );
                else
                    checkInss(         /*Using ref to check dels*/
                        &cigEntryUInt,
                        minStats,
                        refStruct,/*Holds stats for reference*/
                        refQBool, /*0 no q-score entry; 1 present*/
                        &oneUChar, /*deletion instead of insertion*/
                        &intOne    /*Incurment as forwards*/
                    ); /*Checks if should keep the deletion*/

                break;
            } /*Case: deletion*/

            case '\t':
                checkMatches(
                    &cigEntryUInt,
                    minStats,
                    samStruct,
                    refStruct,
                    &qLineBool,
                    refQBool,
                    &useRefBool, /*Tells if using reference for check*/
                    &intCnt
                ); /*Check if keeping the matches*/
                break;
            case 'S':  /*was soft mask*/
                checkSoftMasks(
                    &cigEntryUInt,
                    samStruct,
                    &qLineBool,
                    &intCnt
                );
                break;
        } /*switch, check what kind of cigar entry*/

        /*else is a hard mask, which has been removed. Ignore*/
    } /*loop through all entries in the Cigar*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-4: Find aligned median & means
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(qLineBool != 0)
    { /*If need to find Q-scores*/
        samStruct->meanAligQFlt =
            samStruct->totalAlnQScoreULng /
            ((float) samStruct->readAligLenUInt);

        samStruct->medianAligQFlt =
            qHistToMed(
                samStruct->seqQAlnHistUInt,
                samStruct->readAligLenUInt
        ); /*Get median Q-score*/
    } /*If need to find Q-scores*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-5: clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Minimap2 reverse complements reads, so always working in forward
      cigar. However, if this is not the case, uncomment this as well*/
    /*Need to move back to my origanl positions*/
    /*if(!(samStruct->flagUSht & 16))
    {*/ /*If was a forward sequence*/
        samStruct->seqCStr -= samStruct->readLenUInt;

       if(qLineBool & 1) /*If had a Q-score entry*/
           samStruct->qCStr -= samStruct->readLenUInt;
    /*}*/ /*If was a forward sequence*/

    /*else
    {*/ /*Else will be on the tab before the sequence*/
        /*++samStruct->seqCStr;

       if(qLineBool & 1)*/ /*If had a Q-score entry*/
            /*++samStruct->qCStr;*/

    /*}*/ /*Else will be on the tab before the sequence*/

    if(refStruct != 0)
    { /*If have a reference sequence, need to reset to start position*/
        /*recored starting sequence positions so can reset at end*/
        refStruct->seqCStr = refStartCStr;
        refStruct->qCStr = refQCStr;

        /*Deletions need to be copied over (were counted in ref struct*/
        samStruct->numDelUInt = refStruct->numDelUInt;
        samStruct->numKeptDelUInt = refStruct->numKeptDelUInt;

        /*Reset deletions for next alignment*/
        refStruct->numDelUInt = 0;
        refStruct->numKeptDelUInt = 0;
    } /*If have a reference sequence, need to reset to start position*/

    return;
} /*scoreAln*/

/*######################################################################
# Output:
#    modifies: samStruct Q-score histogram, running Q-score total
#    incurments: samStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
######################################################################*/
void checkSNPs(
    const uint32_t *cigEntryUInt,  /*Cigar entry with number of SNPs*/
    struct minAlnStats *minStats,  /*Has min q-score to keep SNP*/
    struct samEntry *samStruct,    /*Holds stats & sam entry*/
    struct samEntry *refStruct,    /*Holds reference sequence &q-score*/
    const uint8_t *qLineBool,      /*0 no q-score entry; 1 present*/
    const uint8_t *refQBool,       /*0 no ref q-score entry; 1 present*/
    const uint8_t *useRefBool,     /*0 do not use reference, 1 use ref*/
    const int32_t *incInt          /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
) /*Checks if should keep SNPs or discard*/
{ /*checkSNPs*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 Sub-1 TOC: checkMismatches
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t keepBool = 0; /*Keep or ignore mismatch*/
    int32_t
        refQInt = 0,
        tmpQInt = 0; /*Want negatives*/


    for(uint32_t uIntBase = 0; uIntBase < *cigEntryUInt; uIntBase++)
    { /*Loop through all bases that are mismatches*/
        tmpQInt = (*samStruct->qCStr-Q_ADJUST) & (!*qLineBool+MAX_UINT);

        /*Here to avlid an if statement*/
        keepBool =
          ((uint32_t)(minStats->minQChar -tmpQInt) >> 31) | !*qLineBool;
            /*minQChar - tmpQInt: is negative if tmpQInt is > minQChar
              (uint32_t) converts minQChar - tmpQInt from to unsigned
              >> 31 grabs the old negative flag so is 1 when
                   minQChar - tmpQInt is negative
              | !qlineBoolBool sets to one if have no Q-score entry*/

        if(refStruct != 0)
        { /*switch: Check if reference was provided*/
            refQInt = (*refStruct->qCStr - Q_ADJUST + 1);
                /*-1 = minQChar + 1, note the next line (keepBool &=)
                    handles the no ref q-score case*/

            keepBool &=
                ((
                     ((uint32_t) (refQInt - minStats->minQChar) >> 31) |
                     !*refQBool
                 ) | /*Check if a valid Q-score or have Q-score entry*/
                 !*useRefBool /*Check if even using the refrence*/
                );/*Sets keepBool to 1 if not using ref*/
             /*(minQChar - refQInt): is negative if minQChar > refQInt
               uint32_t: Converts to unisgned, so flag not negative
               >> 31: only 1 if refQInt > minQChar (grabs negative flag)
               | !refQBool: Sets to one if their was no q-score line
               !useRefBool: Sets to one if not using the reference
               keepBool &=: Sets keep bool to zero if previous calcs
                            are 0 (reference is high quality and used)*/

            refStruct->qCStr += *refQBool;
            ++refStruct->seqCStr;
        } /*switch: Check if reference was provided*/

        ++samStruct->numSNPUInt; /*! flips bit*/
        samStruct->numKeptSNPUInt += keepBool;

        /*Get the stats*/
        ++(samStruct->seqQAlnHistUInt[tmpQInt]);    /*aligned median*/
        samStruct->totalAlnQScoreULng += tmpQInt; /*For aligned mean*/
        ++(samStruct->readAligLenUInt);     /*incurment aligned length*/

        /*Move to the next base in the sequence*/
        switch(*qLineBool){case 1: samStruct->qCStr += *incInt;}
        samStruct->seqCStr += *incInt;     /*move to next/previous base*/
    } /*Loop through all bases that are mismatches*/

    return;
} /*checkSNPs*/

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
) /*Checks if should keep match or not*/
{ /*checkMatches*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 Sub-1 TOC: checkMatches
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t keepBool = 0; /*Keep or ignore mismatch*/

    int32_t
        refQInt = 0,
        tmpQInt = 0; /*Want negatives*/

    for(uint32_t uIntBase = 0; uIntBase < *cigEntryUInt; ++uIntBase)
    { /*Loop through all bases that are matches*/
        tmpQInt =
            (*samStruct->qCStr - Q_ADJUST) & (!*qLineBool+MAX_UINT);

        keepBool =    /*Find if keeping the match*/
          ((uint32_t)(minStats->minQChar -tmpQInt) >> 31) | !*qLineBool;
            /*minQChar - tmpQInt: is negative if tmpQInt is > minQChar
              (uint32_t) converts minQChar - tmpQInt from to unsigned
              >> 31 grabs the old negative flag so is 1 when
                   minQChar - tmpQInt is negative
              | !qlineBoolBool sets to one if have no Q-score entry*/

        if(refStruct != 0)
        { /*switch: Check if reference was provided*/
            refQInt = (*refStruct->qCStr - Q_ADJUST - 1);
                /*-1 = minQChar + 1, note the next line (keepBool &=)
                    handles the no ref q-score case*/
            keepBool &=
                ((
                     ((uint32_t) (minStats->minQChar - refQInt) >> 31) |
                     !*refQBool
                 ) | /*Check if a valid Q-score or have Q-score entry*/
                 !*useRefBool /*Check if even using the refrence*/
                );/*Sets keepBool to 1 if not using ref*/
             /*(minQChar - refQInt): is negative if refQInt > minQChar
               uint32_t: Converts to unisgned, so flag not negative
               >> 31: only 1 if refQInt > minQChar (grabs negative flag)
               | !refQBool: Sets to one if their was no q-score line
               !useRefBool: Sets to one if not using the reference
               keepBool &=: Sets keep bool to zero if previous calcs
                            are 0 (reference is low quality and used)*/

            refStruct->qCStr += *refQBool;
            ++refStruct->seqCStr;
        } /*switch: Check if reference was provided*/

        ++samStruct->numMatchUInt; /*Total matches*/
        samStruct->numKeptMatchUInt += keepBool;

        /*Get the stats*/
        ++(samStruct->seqQAlnHistUInt[tmpQInt]);    /*aligned median*/
        samStruct->totalAlnQScoreULng += tmpQInt; /*For aligned mean*/
        ++(samStruct->readAligLenUInt);     /*incurment aligned length*/

        switch(*qLineBool){case 1: samStruct->qCStr += *incInt;}
        samStruct->seqCStr += *incInt;     /*move to next/previous base*/
    } /*Loop through all bases that are mismatches*/

    return;
} /*checkMatches*/

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
)/*Check if should keep insertions*/
{ /*checkInsertions*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6: checkInsertions
    #     fun-6 sec-1: variable declerations
    #     fun-6 sec-2: Find if is a valid indel or not
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t 
        keepBool = 0; /*Keep or ignore mismatch*/

    char
        *tmpBaseCStr = 0;

    int32_t tmpQInt = 0; /*Want negatives*/

    uint32_t lenHomoUInt = 0;    /*Size of homopolymer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Find if is a valid indel or not
    #    fun-6 sec-2 sub-1: Check if Q-scores eliminate insertion
    #    fun-6 sec-2 sub-2: Find the homopolymer length
    #    fun-6 sec-2 sub-3: Check if insertion > max homopolymer length
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-6 Sec-2 Sub-1: Check if Q-scores eliminate insertion
    *******************************************************************/

    for(uint32_t intIndel = 0; intIndel < *cigEntryUInt; intIndel++)
    { /*Loop through all insertions*/
        tmpQInt = (*samStruct->qCStr-Q_ADJUST) & (!*qLineBool+MAX_UINT);

        keepBool =
          ((uint32_t)(minStats->minQChar -tmpQInt) >> 31) | !*qLineBool;
            /*minQChar - tmpQInt: is negative if tmpQInt is > minQChar
              (uint32_t) converts minQChar - tmpQInt from to unsigned
              >> 31 grabs the old negative flag so is 1 when
                   minQChar - tmpQInt is negative
              | !qlineBoolBool sets to one if have no Q-score entry*/

        /*adds 0 when useRefForDelBool = 1*/

        switch(*useRefForDelBool)
        { /*Switch: check if checing for deletions*/
            case 0: /*Looking at insertions*/
                samStruct->readAligLenUInt += 1 & !(*useRefForDelBool);

                samStruct->seqQHistUInt[tmpQInt] +=
                    1 & !*useRefForDelBool;
                samStruct->totalAlnQScoreULng +=
                    tmpQInt & !*useRefForDelBool;

                switch(keepBool)
                { /*switch: decided if base is good quality*/
                    case 0:
                    { /*case: low base quality*/
                        ++samStruct->numInsUInt;
                        switch(*qLineBool)
                            {case 1: samStruct->qCStr += *incInt;}
                        samStruct->seqCStr += *incInt;

                        continue; /*Restart at top of loop*/
                    } /*case: low base quality*/
                } /*switch: decided if base is good quality*/

                break;
            /*Case 0:*/

            case 1: /*Looking at deletions*/
                ++samStruct->numDelUInt;
                switch(*qLineBool){case 1: samStruct->qCStr += *incInt;}
                samStruct->seqCStr += *incInt;
                continue;
        } /*Switch: check if checing for deletions*/

        /***************************************************************
        # Fun-6 Sec-2 Sub-2: Find the homopolymer length
        ***************************************************************/

        lenHomoUInt = 1; /*Start homopolymer at base*/
        tmpBaseCStr = samStruct->seqCStr - 1; /*base before indel*/

        /*& ~(33) converts to uppercase and clears first bit
          T is 84 and U is 85, so clearing first bits puts these
          together. This will convert C to B (not A), so their are some
          problems, but should be good enough for this code*/
        while((*samStruct->seqCStr & ~(33))-(*tmpBaseCStr & ~(33)) == 0)
        { /*loop backwards though the homopolymer*/
            ++lenHomoUInt;                /*Count bases in homopolymer*/
            --tmpBaseCStr;                /*Move to next base*/
        } /*loop backwards though the homopolymer*/

        tmpBaseCStr = samStruct->seqCStr + 1; /*base after insert*/

        while((*samStruct->seqCStr & ~(33))-(*tmpBaseCStr & ~(33)) == 0)
        { /*loop forwards though the homopolymer*/
            ++lenHomoUInt;             /*Number of forward bases*/
            ++tmpBaseCStr;             /*Move to next base*/
        } /*loop forwards though the homopolymer*/

        /***************************************************************
        # Fun-6 Sec-2 Sub-3: Check if insertion > max homopolymer length
        ***************************************************************/

        switch(*useRefForDelBool)
        { /*Check if checking for deletions or insertions*/
            case 0: /*Checking for insertions*/
                keepBool =
                    ((uint32_t)
                      (
                         lenHomoUInt - 1 -
                         minStats->maxHomoInsAry[
                             (*samStruct->seqCStr & ~(1+32+64+128)) >> 1
                         ] /*See if am keeping*/
                      )
                    ) >> 31;

                ++samStruct->numInsUInt;
                samStruct->numKeptInsUInt += keepBool;
                break;

            case 1:     /*Is a deletion*/
                keepBool =
                    ((uint32_t)
                      (
                          lenHomoUInt - 1 -
                          minStats->maxHomoDelAry[
                             (*samStruct->seqCStr & ~(1+32+64+128)) >> 1
                          ] /*See if am keeping*/
                      )
                    ) >> 31;

                    ++samStruct->numDelUInt;
                    samStruct->numKeptDelUInt += keepBool;
                    break;
        } /*Check if checking for deletions or insertions*/
        /* Logic behind keepBool =
              & ~(1 + 32 + 64 + 128) clears 1st bit and last 3 bits
                  - 1: Merges T and U
                  - cearing 32: convert to uppler
                  - cearing 64 & 128: Ensures array size is under 32
                  - >> 1 removes first bit (reduces to 16 bytes)
              minStats->maxInsHomoUIntAry has the min insertions
                  homopolymer size
                  - minStats->maxDelHomoUIntAry for deletions
              lenHomoUInt - : makes negative if keeping hompolmer
                  - The -1 is to account for being equal to the min
                    threshold
              (uint32_t): Converts int32_t to unsigned, which ensures
                  that the bit value is positive
              >> 31 keeps negative flag (1 if under, 0 if over)

              Result: 1 if keeping, 0 if discarding
              */

        switch(*qLineBool){case 1: samStruct->qCStr += *incInt;}
        samStruct->seqCStr += *incInt;
    } /*Loop through all insertions*/

    return;
} /*checkInsertions*/

/*######################################################################
# Output:
#    modifies: samStruct with stats
######################################################################*/
void checkDels(
    const uint32_t *cigEntryUInt,  /*Cigar with number deletions*/
    struct minAlnStats *minStats, /*Has min q-score to keep deletions*/
    struct samEntry *samStruct     /*Holds stats & sam entry*/
) /*Checks if should keep the deletion*/
{ /*checkDeletions*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7: checkDeletions
    #     fun-7 sec-1: Variable declerations
    #     fun-7 sec-2: Find if is a valid deletion
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpBaseCStr = 0,           /*Loop through the homopolymers*/
        *backBaseCStr = samStruct->seqCStr - 1;

    uint8_t
        keepBool = 0,
        keepBackBool = 0;           /*Back homopolymer is longer*/
           /*backBaseCStr decleration here is to avoid compile warnings*/

    int32_t
         lenBackHomoInt = 0,
         lenForHomoInt = 0,
         lenHomoInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-2: Find if is a valid deletion
    #    fun-7 sec-2 sub-1: Find hompolymer lengths around deletion
    #    fun-7 sec-2 sub-2: Find total homopolymer length
    #    fun-7 sec-2 sub-3: Check if should keep deltion
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-7 Sec-2 Sub-1: Find hompolymer lengths around deletion
    *******************************************************************/

    tmpBaseCStr = backBaseCStr - 1; /*move to base before del*/
    lenBackHomoInt = 1;  /*Homopolymer Is at least one base long*/

    /*& ~(33) converts to uppercase and clears first bit
      T is 84 and U is 85, so clearing first bits puts these
      together. This will convert C to B (not A), so their are some
      problems, but should be good enough for this code*/
    while((*backBaseCStr & ~(33)) - (*tmpBaseCStr & ~(33)) == 0)
    { /*loop backwards though the homopolymer*/
        ++lenBackHomoInt;  /*At least 1 more base in back homopolymer*/
        --tmpBaseCStr;     /*Move to the previous base*/
    } /*loop backwards though the homopolymer*/

    /*samStruct->seqCStr points to the first base after the deletion*/
    tmpBaseCStr = samStruct->seqCStr + 1; /*base after current base*/
    lenForHomoInt = 1;   /*Has at least one base*/

    while((*samStruct->seqCStr & ~(33)) - (*tmpBaseCStr & ~(33)) == 0)
    { /*loop backwards though the homopolymer*/
        ++lenForHomoInt; /*At least 1 more base in forward homopolymer*/
        ++tmpBaseCStr;   /*Move to the previous base*/
    } /*loop backwards though the homopolymer*/

    /*******************************************************************
    # Fun-7 Sec-2 Sub-2: Find total homopolymer length
    *******************************************************************/

    keepBackBool =
        ((int32_t) -
             (!!
               ((*samStruct->seqCStr & ~(33)) - (*backBaseCStr & ~(33)))
             ) &
        (1 | (lenForHomoInt - lenBackHomoInt) >> 31)) + 1;
        /* The first line determines if the bases are the same & saves
           me several if checks (convert case, check if same, check if
           T or U present (alternative combination for same base))
               *seqCStr & ~(33) sets a=A, c=C, g=G, and t=T, u, or U
               seqCStr - backBaseCStr: gives 0 if bases are same
               (!!: keeps 0 & converts anything else to 1
               (int32_t) -: converts to signed int32_t & makes negative
                   - Need to keep the negative when using & with the 
                     first and second line
             Result: 0 if bases are same, 1 if are differnt
          The second line needs three operations to determine which
          homopolymer is greater (if(x == y) else if (x > y) else)
              lenForHomoInt - lenBackHomoInt: is negative if back > for
              >> 31: Removes all bits except for negative bit
                  - I am working with int32_t, so >> 31 gives a -1 when
                    a number is negative (lenFor & lenBack are int32_t)
              - So 1 if negative, 0 if positive
           Result: -1 if back homoplymer is > forward hompolymer
                    1 if forward homopolymer is > back homoplymer
          First line & second lines sets 0 if the bases are the same &
            -1 or 1 if the bases are differnt.
          + 1 removes the negative
       Result: 1 if bases are the same, 0 if backwards homopolymer is
               larger, & 2 if forwards hompolymer is larger
       These values are used in the switch statment, which will likely
         be converted to a look up table to determine the size.
       */

    switch(keepBackBool)
    { /*Switch: Get my hompolymer size*/
        case 0:
            lenHomoInt = lenBackHomoInt;
            break;
        case 1:
            lenHomoInt = lenBackHomoInt + lenForHomoInt;
            break;
        case 2:
            lenHomoInt = lenForHomoInt;
            break;
    } /*Switch: Get my hompolymer size*/
    
    /*******************************************************************
    # Fun-7 Sec-2 Sub-3: Check if should keep deltion
    *******************************************************************/

    keepBool =
        ((uint32_t)
         (minStats->maxHomoDelAry[
             (*samStruct->seqCStr & ~(1+32+64+128)) >> 1
         ] - lenHomoInt
        )) >> 31;
        /*& ~(1 + 32 + 64 + 128) clears 1st bit and last 3 bits
              - 1: Merges T and U
              - cearing 32: convert to uppler
              - cearing 64 & 128: Ensures array size is under 32
              - >> 1 removes first bit (reduces to 16 bytes)
          minStats->maxInsHomoUIntAry has value for homopolymers
          - lenHomoUInt: makes negative if discarding homopolymer
          (uint32_t): Make sure is an unsigned value.
              ensures negative bit does not make negative
          >> 31 keeps negative flag (0 if under, 1 if over)
         Result: 0 if keeping, 1 if discarding
        */

    /*Record how many deletions were present*/
    samStruct->numDelUInt += *cigEntryUInt;
    samStruct->numKeptDelUInt += *cigEntryUInt & (keepBool + MAX_UINT);
        /*MAX_UINT + 1 (discarding) = 0, but all bits in uint32_t set to
          1 if keeping (MAX_UINT + 0 (keeping) = 0xFFFFFFFF)*/

    return; /*Do not need to incurment pointers, since no base present*/
} /*checkDeletions*/

/*######################################################################
# Output:
#    modifies: seqCStr to point at bast after soft mask
#    modifies: qCStr to point at q-score after soft mask
######################################################################*/
void checkSoftMasks(
    const uint32_t *cigUInt,    /*Number of bases soft masked*/
    struct samEntry *samStruct, /*Holds sequence and Q-score entry*/
    const uint8_t *qLineBool,   /*0 no q-score entry; 1 present*/
    const int32_t *incInt       /*-1 to incurment sequence backwards
                                     1: to incurement sequence forwards
                                   */
) /*Moves sequence & Q-score pointers past soft mask*/
{ /*checkSoftMasks*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1 Sub-1 TOC: checkSoftMasks
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    samStruct->seqCStr += (*incInt * (int32_t) *cigUInt);
    switch(*qLineBool)
        {case 1: samStruct->qCStr += (*incInt * (int32_t) *cigUInt);}

    return;
} /*checkSoftMask*/

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
) /*Reads reference from fastq & check that it meets min requirments*/
{ /*readAndCheckRef*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1 Sub-1 TOC: readAndCheckRef
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t errUChar = readRefFqSeq(refFILE, refStruct);
    
    switch(errUChar)
    { /*switch, check if valid file provided*/
        case 0:
            return 0;
        case 2:
            return 2;
    } /*switch, check if valid file provided*/

    findQScores(refStruct);

    if(
        refStruct->meanQFlt < minStats->minMedianQFlt ||   /*Medain q*/
        refStruct->meanQFlt < minStats->minMeanQFlt ||     /*Mean Q*/
        refStruct->readLenUInt < minStats->minReadLenULng  /*length*/
    ) /*If the reference is beneath the min stats for a reference*/
        return 4;

    return 1;
} /*readAndCheckRef*/

/*######################################################################
# Output:
#    Modifies: refStruct to hold the read in fastq entry & sets its
#              pointers
#    Returns:
#        - 0: if no reference file provided
#        - 1: if succeded
#        - 2: If file was not a fastq file
#        - 64: If malloc failed to find memory
######################################################################*/
uint8_t readRefFqSeq(
    FILE *refFILE,       /*Pointer to fastq file to grab reference from*/
    struct samEntry *refStruct /*Sam entry struct to hold reference*/
) /*Gets the frist reads sequence & q-score line from a fastq file*/
{ /*readRefFqSeq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 TOC: readRefFqSeq
    #    fun-10 sec-1: Variable declarations
    #    fun-10 sec-2: Check if need to allocate memory for buffer
    #    fun-10 sec-3: Read in the first data
    #    fun-10 sec-4: If not at file end, see if have the full entry
    #    fun-10 sec-5: Read till end of file, check if valid fastq entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-10 Sec-1: Variable declarations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    uint8_t
        numNewLineUChar = 0;      /*Tells me which fastq entry I am on*/
    char
        *fqIterUCStr = 0,         /*Marks spot working on in fastq*/
        *oldBuffUCStr = 0;        /*For reallocs*/
    uint32_t numCharInLineUInt = 0;
    uint64_t bytesReadInULng = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-10 Sec-2: Check if need to allocate memory for buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(refFILE == 0)
        return 0;    /*No file provided*/

    if(refStruct->samEntryCStr == 0)
    { /*If have no buffer*/
        refStruct->lenBuffULng = 1400;
        refStruct->samEntryCStr = malloc(sizeof(char) * 1401);

        if(refStruct->samEntryCStr == 0)
        { /*If memory allocation errory*/
            printMemAlocErr(
                "scoreReadsScoringWrappers.c",
                "readRefFqSeq",
                5,
                559
            ); /*Let user know of memory allocation failure*/
    
            return 64;
        } /*If memory allocation errory*/
    } /*If have no buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-10 Sec-3: Read in the first data
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    bytesReadInULng =
        fread(
            refStruct->samEntryCStr,
            sizeof(uint8_t),
            refStruct->lenBuffULng,
            refFILE
    ); /*Do my first read*/

    if(refStruct->samEntryCStr == 0)
        return 2;                          /*Nothing in the file*/

    /*Set up the query id pointer (+ 1 to get off the @)*/
    refStruct->queryCStr = refStruct->samEntryCStr + 1;
    fqIterUCStr = refStruct->samEntryCStr;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-10 Sec-4: If not at file end, see if have the full entry
    #    fun-10 sec-4 sub-1: Check if have read in the entire entry
    #    fun-10 sec-4 sub-2: Check if need to read in more buffer
    #    fun-10 sec-4 sub-3: Read in more buffer
    #    fun-10 sec-4 sub-4: Reset pointers for new buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    while(fqIterUCStr != 0)
    { /*While not at end of buffer*/

        /***************************************************************
        # Fun-10 Sec-4 Sub-1: Check if have read in the entire entry
        ***************************************************************/

        if(*fqIterUCStr == '\n')
        { /*If found the end of an entry*/
            ++numNewLineUChar;

            switch(numNewLineUChar)
            { /*Switch detect witch new line I am on*/
                case 1:
                    numCharInLineUInt = 0;
                    refStruct->seqCStr = fqIterUCStr + 1;
                    break;
                case 3:
                    --numCharInLineUInt; /*Will be one base off*/
                    refStruct->readLenUInt = numCharInLineUInt;
                    refStruct->qCStr = fqIterUCStr + 1;
                    break;
                case 4:
                    return 1;
            } /*Switch detect witch new line I am on*/
        } /*If found the end of an entry*/

        /***************************************************************
        # Fun-10 Sec-4 Sub-2: Check if need to read in more buffer
        ***************************************************************/

        if(*fqIterUCStr == '\0')
        { /*If at the end of my buffer*/
            oldBuffUCStr = refStruct->samEntryCStr;
            bytesReadInULng = refStruct->lenBuffULng;
            refStruct->lenBuffULng = refStruct->lenBuffULng << 1;

            refStruct->samEntryCStr =
                realloc(
                    refStruct->samEntryCStr,
                    sizeof(uint8_t) * refStruct->lenBuffULng
            ); /*double the buffer size*/

            if(refStruct->samEntryCStr == 0)
            { /*If have no buffer*/
                refStruct->lenBuffULng = 1400;
                refStruct->samEntryCStr = malloc(sizeof(char) * 1401);

                if(refStruct->samEntryCStr == 0)
                { /*If memory allocation errory*/
                    printMemAlocErr(
                        "scoreReadsScoringWrappers.c",
                        "readRefFqSeq",
                        5,
                        559
                    ); /*Let user know of memory allocation failure*/
            
                    return 64;
                } /*If memory allocation errory*/
            } /*If have no buffer*/

            /***********************************************************
            # Fun-10 Sec-4 Sub-3: Read in more buffer
            ***********************************************************/

            /*Reset the iterator for the new buffer*/
            fqIterUCStr = refStruct->samEntryCStr + bytesReadInULng;

            bytesReadInULng =  /*Read in the next buffer*/
                fread(
                    fqIterUCStr,
                    sizeof(uint8_t),
                    bytesReadInULng,
                    refFILE
            ); /*Do my first read*/

            /***********************************************************
            # Fun-10 Sec-4 Sub-4: Reset pointers for new buffer
            ***********************************************************/

            refStruct->refCStr = refStruct->samEntryCStr;

            if(refStruct->seqCStr != 0)
            { /*If need to update the sequence position*/
                refStruct->seqCStr =
                    refStruct->samEntryCStr +
                    (refStruct->seqCStr - oldBuffUCStr);
            } /*If need to update the sequence position*/

            if(refStruct->qCStr != 0)
            { /*If need to update the q-score position*/
                refStruct->qCStr =
                    refStruct->samEntryCStr +
                    (refStruct->qCStr - oldBuffUCStr);
            } /*If need to update the q-score position*/
        } /*If at the end of my buffer*/

        else
            ++numCharInLineUInt; /*Count number of characters in entry*/

        ++fqIterUCStr;
    } /*While not at end of buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-10 Sec-5: Read till end of file, check if is valid fastq entry
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(numNewLineUChar == 3 && *refStruct->qCStr > 32)
        return 1;                  /*No new line after q-score entry*/

    return 2; /*End of file, but did not have four lines*/
} /*readRefFqSeq*/

/*######################################################################
# output:
#    returns: The read length (unsigned long)
# Note: 
#    - samStruct->seqQHistUInt must have all values initialized to 0
#    - Requires: qHistToMedian from scoreReadsSamFileFunctions.c
######################################################################*/
void findQScores(
    struct samEntry *samStruct /*Sam entry to find Q-scores for*/
) /*Finds Q-scores for input sam entry*/
{ /*findQScores*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-1 Sub-1 TOC: findQScores
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *seqQCStr = samStruct->qCStr;
    uint64_t
        lenReadULng = 0;

    samStruct->totalQScoreULng = 0;
    samStruct->meanQFlt = 0;
    samStruct->medianQFlt = 0;

    if(seqQCStr == 0)
        return;       /*No Q-score entry*/

           /* *qCStr ^ '*' is 0 when qCStr == '*'
              *(qCStr + 1) ^ '\t' is 0 when qCStr == '\t'
              | ensures value is only 0 when qCStr = "*\t"*/
    if(((*samStruct->qCStr ^ '*') | (*(samStruct->qCStr+1) ^ '\t')) ==0)
        return;       /*No Q-score entry*/

    while(*seqQCStr > 32)
    { /*loop through the Q-score entry*/
        ++(samStruct->seqQHistUInt[*seqQCStr - Q_ADJUST]);
        samStruct->totalQScoreULng += (*seqQCStr - Q_ADJUST);
        ++seqQCStr;
        ++lenReadULng;
    } /*loop through the Q-score entry*/
    
    /*Find the mean and median*/
    samStruct->meanQFlt =
        samStruct->totalQScoreULng / ((float) lenReadULng);

    samStruct->medianQFlt =
        qHistToMed(samStruct->seqQHistUInt, lenReadULng);

    return;
} /*findQScores*/

/*######################################################################
# Output:
#    returns: double with the median Q-score [or 0 if nothing]
######################################################################*/
float qHistToMed(
    uint32_t qHistUInt[], /*Histogram of Q-scores*/
    uint32_t readLenUInt     /*Number of bases in the read*/
) /*converts histogram of q-scores into samStruct into median Q-score*/
{ /*qHistToMedian*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-12 Sec-1 Sub-1 TOC: qHistToMedian
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint32_t
        numBasesUInt = 0,
        midPointULng = readLenUInt / 2;

    for(uint32_t intQ = 0; intQ < MAX_Q_SCORE; intQ++)
    { /*Loop through Q-score histogram to find mean and median Q-score*/
        numBasesUInt += qHistUInt[intQ];

        if(numBasesUInt >= midPointULng)
        { /*if found the midpoint, then find the median*/

            if(numBasesUInt > midPointULng)
                return intQ;  /*at least 2 bases past mid have same Q*/

            else if(numBasesUInt % 2 == 1)
                return intQ;  /*Is odd, so Their is a middle number*/

            else
            { /*Else is even, so no middle number (two numbers at mid)*/
                numBasesUInt = intQ;
                ++intQ;

                while(qHistUInt[intQ] == 0)
                    ++intQ;

                return (numBasesUInt + qHistUInt[intQ]) / ((float) 2);
            } /*Else is even, so no middle number (two numbers at mid)*/
        } /*if found the midpoint, then find the median*/
    } /*Loop through Q-score histogram to find mean and median Q-score*/

    return 0; /*Just in case was all 0's in the array*/
} /*qHistToMedian*/

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
) /*Reads a single entry from a eqx cigar line*/
{ /*readCigEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-1 Sub-1 TOC: readCigEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *retUInt = 0;

    if((**cigarUCStr) > 57 || (**cigarUCStr) < 48)
        ++(*cigarUCStr);               /*Need to move off old entry*/

    while((**cigarUCStr) < 58 && (**cigarUCStr) > 47)
    { /*Loop till of number entry of cigar*/
        *retUInt = (*retUInt * 10) + (**cigarUCStr - 48);
        ++(*cigarUCStr);
    } /*Loop till of number entry of cigar*/

    return;
} /*readCigEntry*/

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
) /*Reads a single entry from a eqx cigar line*/
{ /*readReverseCigEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-14 Sec-1 Sub-1 TOC: readReverseCigEntry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t tensPosULng = 1;

    *retUInt = 0;

    if(**cigarUCStr > 57 || **cigarUCStr < 48)
        --(*cigarUCStr);                    /*move off entry type*/
     
    while(**cigarUCStr < 58 && **cigarUCStr > 47)
    { /*Loop till of number entry of cigar*/
        *retUInt = (*retUInt) + ((**cigarUCStr - 48) * tensPosULng);
            /*Building number backwards, so a bit odd*/
        tensPosULng *= 10;
        --(*cigarUCStr);
    } /*Loop till of number entry of cigar*/

    return;
} /*readReverseCigEntry*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
######################################################################*/
void blankMinStats(
    struct minAlnStats *minStats
) /*Sets minStats minimum requirements for sam alingemtns to defaults*/
{ /*blankMinStats*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-15 Sec-1 Sub-1 TOC: blankMinStats
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    minStats->minMapqUInt = 20;          /*min mapping quality*/
    minStats->minQChar = 10;             /*default min Q-score*/
    minStats->minMedianQFlt = 13;        /*default median Q-score*/
    minStats->minMeanQFlt = 13;          /*default mean Q-score*/
    minStats->minAlignedMedianQFlt = 13; /*median aligend Q*/
    minStats->minAlignedMeanQFlt = 13;   /*mean alinged Q-score*/
    minStats->maxReadLenULng = 1000;     /*default max read length*/
    minStats->minReadLenULng = 600;      /*default max read length*/
    
    for(uint8_t uCharCnt = 0; uCharCnt < 16; ++uCharCnt)
    { /*loop till have initialized the deltion and insertion arrays*/
        minStats->maxHomoInsAry[uCharCnt] = 0;
        minStats->maxHomoDelAry[uCharCnt] = 0;
    } /*loop till have initialized the deltion and insertion arrays*/

    minStats->maxHomoInsAry[0] = 1; /*max A insertion homopolymer*/ 
    minStats->maxHomoInsAry[10] = 1; /*max T insertion homopolymer*/ 
    minStats->maxHomoInsAry[1] = 1; /*max C insertion homopolymer*/ 
    minStats->maxHomoInsAry[3] = 1; /*max G insertion homopolymer*/ 

    minStats->maxHomoDelAry[0] = 0; /*max A deletion homopolymer*/ 
    minStats->maxHomoDelAry[10] = 0; /*max T deletion homopolymer*/ 
    minStats->maxHomoDelAry[1] = 0; /*max C deletion homopolymer*/ 
    minStats->maxHomoDelAry[3] = 0; /*max G deletion homopolymer*/ 

    return;
} /*blankMinStats*/
