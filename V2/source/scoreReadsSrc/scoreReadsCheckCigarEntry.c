/*##############################################################################
# Name: scoreReadsCheckCigarEntry
# Use: Functions to checks cigar entries and incurments pointers for samfile
##############################################################################*/

#include "scoreReadsCheckCigarEntry.h" /*Holds structers to hold samfile stats*/

/*##############################################################################
# TOC:
#    fun-1: checkMIsmatches: Check if want to keep or ignore a mismatch
#    fun-2: checkMatches: Incurment seq and q-score entry pointers 
#    fun-3: checkInsertion: Checks if insertion should be kept
#    fun-4: checkDeletions: Checks if deletions should be kept
#    fun-5 checkIfBasesMatch: Check if two bases are matches
#    fun-6 checkSoftMasks: Increments pointer when soft mask is present
#    fun-7: readCigEntry: Reads in cigar entry and returns unsigned long
##############################################################################*/

/*##############################################################################
# Name: checkMismatches
# Use: Checks is a mismatche is valid or not. Also, adds value to q-score
#      histogram and tally.
# Input: 
#    scoreStruct: scoreReadStruct with sequence, q-score, number of bases,
#                   running Q-score, and Q-score histogram (scoreReadStruct ptr)
#    minReadStats: Min q-score to keep a mismatch (character)
#    samStruct: Stores the stats for the read (samEntry struct)
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: scoreReadStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
##############################################################################*/
void checkMismatches(scoreReadStruct *scoreStruct,
                     minStats *minReadStats,
                     samEntry *samStruct)
{ /*checkMismatches*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1: checkMismatches
    #    fun-1 sec-1: variable declerations
    #    fun-1 sec-2: Check if there is a Q-score line
    #    fun-1 sec-3: If Q-score line, determine if should keep mismatch
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int tmpQInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check if there is a Q-score line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*scoreStruct->qScoreCStr == '*' && *(scoreStruct->qScoreCStr++) == '\t')
    { /*if no Q-score line*/
        (*samStruct).numMisULng = (*samStruct).numMisULng + 
                                       (*scoreStruct).cigarEntryULng;

        /*Move to the bases after the mismatch*/
        (*scoreStruct).sequenceCStr += (*scoreStruct).cigarEntryULng;

        (*samStruct).readAligLenULng++;    /*Add base to aligned length*/
        return;
    } /*if no Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: If Q-score line, determine if should keep mismatch
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(unsigned long uLngBase = 0;
        uLngBase < (*scoreStruct).cigarEntryULng;
        uLngBase++)
    { /*Loop through all bases that are mismatches*/
        tmpQInt = *(*scoreStruct).qScoreCStr - Q_ADJUST_CHAR; /*Get Q-score*/

        if(tmpQInt < (*minReadStats).minQChar)
            (*samStruct).numIgnoreMisULng++;    /*Ignore, < min Q-score*/
        else
            (*samStruct).numMisULng++;          /*Keep, > min Q-score*/

        /*Get the stats*/
        (*scoreStruct).seqQHistULng[tmpQInt]++;     /*For aligned median*/
        (*scoreStruct).totalQScoreULng += tmpQInt; /*For aligned mean*/
        (*samStruct).readAligLenULng++;       /*incurment aligned length*/

        /*Move to the next base in the sequence*/
        (*scoreStruct).qScoreCStr++;
        (*scoreStruct).sequenceCStr++;
    } /*Loop through all bases that are mismatches*/

    return;
} /*checkMismatches*/

/*##############################################################################
# Name: checkMatches
# Use: Pretends to check mathces, but really incurments the sequence and q-score
#      pointers. Also, gets the aligned Q-score histogram and q-score total
# Input: 
#    scoreStruct: scoreReadStruct with sequence, q-score, number of bases,
#                   running Q-score, and Q-score histogram (scoreReadStruct ptr)
#    minReadStats: Min q-score to keep a mismatch (character)
#    samStruct: Stores the aligned length (samStruct)
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: scoreReadStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
##############################################################################*/
void checkMatches(scoreReadStruct *scoreStruct,
                  minStats *minReadStats,
                  samEntry *samStruct)
{ /*checkMatches*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2: checkMatches
    #     fun-2 sec-1: Find if is a valid match or not
    #     fun-2 sec-2: Check if has a Q-score entry
    #     fun-2 sec-3: Check if bases have a valid Q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Find if is a valid match or not
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int tmpQInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Check if has a Q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*scoreStruct->qScoreCStr == '*' && *(scoreStruct->qScoreCStr++) == '\t')
    { /*If there is no Q-score line, increment pointers*/

       /*Move to the bases after the mismatch*/
       scoreStruct->sequenceCStr += scoreStruct->cigarEntryULng;

       samStruct->numMatchULng++;
       samStruct->numKeptMatchULng++;
       samStruct->readAligLenULng++; /*Add base to aligned length*/

       return;
    } /*If there is no Q-score line, increment pointers*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Check if bases have a valid Q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(unsigned long uLngBase = 0;
        uLngBase < scoreStruct->cigarEntryULng;
        uLngBase++
    ) { /*Loop through all bases that are matches*/
        tmpQInt = *scoreStruct->qScoreCStr - Q_ADJUST_CHAR;  /*get Q-score*/
        scoreStruct->sequenceCStr++;              /*Move to the next base*/
        scoreStruct->qScoreCStr++;                /*Move to the next base*/
        scoreStruct->seqQHistULng[tmpQInt]++;     /*For aligned median*/
        scoreStruct->totalQScoreULng += tmpQInt;  /*For aligned mean*/

        samStruct->numMatchULng++;

        if(tmpQInt > minReadStats->minQChar)
            samStruct->numKeptMatchULng++;

        samStruct->readAligLenULng++;          /*aligned length*/
    } /*Loop through all bases that are mismatches*/

    return;
} /*checkMatches*/

/*##############################################################################
# Name: checkInsertion
# Use: Checks is a insertion is valid or not. Also, adds value to q-score
#      histogram and tally if the insertion is valid.
# Input: 
#    scoreStruct: scoreReadStruct with sequence, q-score, number of bases,
#                   running Q-score, and Q-score histogram (scoreReadStruct ptr)
#    minReadStats: Holds the min requirements to keep a base (minStats struct)
#    samStruct: Stores the number of kept and ignored indels (samEntry)
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: scoreReadStruct sequence and Q-score c-strings
#    modifies: samStruct with stats
##############################################################################*/
void checkInsertion(scoreReadStruct *scoreStruct,
                    minStats *minReadStats,
                    samEntry *samStruct)
{ /*checkInsertions*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3: checkInsertions
    #     fun-3 sec-1: variable declerations
    #     fun-3 sec-2: Find if is a valid indel or not
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpBaseCStr = 0;
    int tmpQInt = 0;
    unsigned long lenHomoULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Find if is a valid indel or not
    #    fun-3 sec-2 sub-1: Check Q-scores eliminate insertion
    #    fun-3 sec-2 sub-2: Find the homopolymer length
    #    fun-3 sec-2 sub-3: Check if insertion is > max homopolymer length
    #    fun-3 sec-2 sub-4: Increment pointers and record Q-scores
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-3 Sec-2 Sub-1: Check Q-scores eliminate insertion
    ***************************************************************************/

    for(int intIndel = 0;
        intIndel < (*scoreStruct).cigarEntryULng;
        intIndel++
       )
    { /*Loop through all insertions*/
        lenHomoULng = 1;             /*Have at least one base (the insertion)*/
        tmpQInt = *(*scoreStruct).qScoreCStr - Q_ADJUST_CHAR;  /*get Q-score*/

        if(*(*scoreStruct).qScoreCStr != '*' &&
           *(scoreStruct->qScoreCStr++) == '\t'&&
           tmpQInt < (*minReadStats).minQChar
          )
        { /*if insertion under min quality, ignore*/
            ((*samStruct).numIgnoreInsULng)++;       /*Insertion > min Q*/
            ((*samStruct).readAligLenULng)++;        /*aligned length*/

            (*scoreStruct).seqQHistULng[tmpQInt]++;     /*For aligned median*/
            (*scoreStruct).totalQScoreULng += tmpQInt;  /*For aligned mean*/
            (*scoreStruct).qScoreCStr++;
            (*scoreStruct).sequenceCStr++;
            continue;
        } /*if insertion under min quality, ignore*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-2: Find the homopolymer length
        ***********************************************************************/

        tmpBaseCStr = (*scoreStruct).sequenceCStr - 1; /*base before insert*/

        while(checkIfBasesMatch(*(*scoreStruct).sequenceCStr,
                                *tmpBaseCStr) > 0)
        { /*loop backwards though the homopolymer*/
            lenHomoULng++;                       /*Count bases in homopolymer*/
            tmpBaseCStr--;                       /*Move to next base*/
        } /*loop backwards though the homopolymer*/

        tmpBaseCStr = (*scoreStruct).sequenceCStr + 1;  /*base after insert*/

        while(checkIfBasesMatch(*(*scoreStruct).sequenceCStr,
                                *tmpBaseCStr) > 0)
        { /*loop forwards though the homopolymer*/
            lenHomoULng++;                             /*bases in homopolymer*/
            tmpBaseCStr++;                             /*Move to next base*/
        } /*loop forwards though the homopolymer*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-3: Check if insertion is > max homopolymer length
        ***********************************************************************/

        if(checkIfBasesMatch(*(*scoreStruct).sequenceCStr, 'A') == 1)
        { /*If was an Hompolymer of A's*/
           if(lenHomoULng > (*minReadStats).maxInsAHomoULng)
                ((*samStruct).numIgnoreInsULng)++;
            else
                ((*samStruct).numInsULng)++;
        } /*If was an Hompolymer of A's*/
    
        else if(checkIfBasesMatch(*(*scoreStruct).sequenceCStr, 'T') == 1)
        { /*else If was an Hompolymer of T's*/
            if(lenHomoULng > (*minReadStats).maxInsTHomoULng)
                ((*samStruct).numIgnoreInsULng)++;
            else
                ((*samStruct).numInsULng)++;
        } /*else If was an Hompolymer of T's*/
    
        else if(checkIfBasesMatch(*(*scoreStruct).sequenceCStr, 'C') == 1)
        { /*else If was an Hompolymer of C's*/
            if(lenHomoULng > (*minReadStats).maxInsCHomoULng)
                ((*samStruct).numIgnoreInsULng)++;
            else
                ((*samStruct).numInsULng)++;
        } /*else If was an Hompolymer of C's*/
    
        else if(checkIfBasesMatch(*(*scoreStruct).sequenceCStr, 'G') == 1)
        { /*else If was an Hompolymer of G's*/
            if(lenHomoULng > (*minReadStats).maxInsGHomoULng)
                ((*samStruct).numIgnoreInsULng)++;
            else
                ((*samStruct).numInsULng)++;
        } /*else If was an Hompolymer of G's*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-4: increment pointers and record Q-score
        ***********************************************************************/

        ((*samStruct).readAligLenULng)++; /*Add base to aligned length*/

        if(*(*scoreStruct).qScoreCStr != '*' &&
           *(scoreStruct->qScoreCStr++) == '\t'
        ) { /*If we have Q-scores*/
            (*scoreStruct).seqQHistULng[tmpQInt]++;     /*For aligned median*/
            (*scoreStruct).totalQScoreULng += tmpQInt;   /*For aligned mean*/
        } /*If we have Q-scores*/

        (*scoreStruct).qScoreCStr++; /*Move to the next Q-score*/
        (*scoreStruct).sequenceCStr++; /*Move to the next base*/
    } /*Loop through all insertions*/

    return;
} /*checkInsertions*/

/*##############################################################################
# Name: checkDeletions
# Use: Checks is a deletion is valid or not and increments the Q-score entry
# Input: 
#    scoreStruct: scoreReadStruct with sequence, q-score, number of bases,
#                   running Q-score, and Q-score histogram (scoreReadStruct ptr)
#    minReadStats: Holds the min requirements to keep a base (minStats struct)
#    samStruct: Stores the number of kept and ignored indels (samEntry)
# Output:
#    modifies: samStruct with stats
# Note: Without reference, this is a bit tricky. So, I am just going to check
#       the largest hompolymer around the deletion. Also assuming seqCStr starts
#       on the base before the deletion
##############################################################################*/
void checkDeletions(scoreReadStruct *scoreStruct,
                    minStats *minReadStats,
                    samEntry *samStruct)
{ /*checkDeletions*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4: checkDeletions
    #     fun-4 sec-1: Variable declerations
    #     fun-4 sec-2: Find if is a valid indel or not
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpBaseCStr = 0,                    /*Loop through the homopolymers*/
         baseCStr = 0,                       /*Holds longest homopolymer base*/
         *backBaseCStr = (*scoreStruct).sequenceCStr - 1;
           /*backBaseCStr decleration here is to avoid compile warnings*/
    unsigned long lenBackHomoULng = 0,
                  lenForHomoULng = 0,
                  lenHomoULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Find if is a valid indel or not
    #    fun-4 sec-2 sub-1: Find the homopolymer length behind the deletion
    #    fun-4 sec-2 sub-2: Find the homopolymer length after the deletion
    #    fun-4 sec-2 sub-3: Check if deletion in homoplymer > max length
    #    fun-4 sec-2 sub-4: Increment pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # fun-4 sec-2 sub-1: Find homopolymer length before the deletion
    ***************************************************************************/

    tmpBaseCStr = backBaseCStr - 1; /*move to base behind deletions*/
    lenBackHomoULng++;  /*Homopolymer Is at least one base long*/

    while(checkIfBasesMatch(*backBaseCStr, *tmpBaseCStr) > 0)
    { /*loop backwards though the homopolymer ('\t' will output as non-match)*/
        lenBackHomoULng++;  /*Count number of bases in homopolymer*/
        tmpBaseCStr--; /*Move to the next base*/
    } /*loop backwards though the homopolymer*/

    /***************************************************************************
    # fun-4 sec-2 sub-2: Find the homopolymer length after the deletion
    ***************************************************************************/

    tmpBaseCStr = (*scoreStruct).sequenceCStr + 1; /*base after del*/

    if(checkIfBasesMatch(*backBaseCStr, *(*scoreStruct).sequenceCStr) < 1)
    { /*If the bases around the homopolymer are different*/
        lenForHomoULng++;  /*Homopolymer is at least one base long*/

        while(checkIfBasesMatch(*(*scoreStruct).sequenceCStr,
                                *tmpBaseCStr) > 0)
        { /*loop forwards though the homopolymer*/
            lenForHomoULng++;  /*Count number of bases in homopolymer*/
            tmpBaseCStr++; /*Move to the next base*/
        } /*loop forwards though the homopolymer*/

        if(lenForHomoULng > lenBackHomoULng) /*Find the longest homopolymer*/
        { /* if the bast after the deletion had a longer homopolymer*/
            lenHomoULng = lenForHomoULng;
            baseCStr = *(*scoreStruct).sequenceCStr;
        } /* if the bast after the deletion had a longer homopolymer*/

        else
        { /* if the bast before the deletion had a longer homopolymer*/
            lenHomoULng = lenBackHomoULng;
            baseCStr = *backBaseCStr;
        } /* if the bast before the deletion had a longer homopolymer*/
    } /*If the bases around the homopolymer are different*/

    else
    { /*else, bases are same and need to find length of forward bases*/
        lenHomoULng = lenBackHomoULng + 1;
        baseCStr = *(*scoreStruct).sequenceCStr;         /*Homopolymer base*/

        while(checkIfBasesMatch(*(*scoreStruct).sequenceCStr,
                                *tmpBaseCStr) > 0)
        { /*loop forwards though the homopolymer*/
            lenHomoULng++;             /*Count number of bases in homopolymer*/
            tmpBaseCStr++;             /*Move to the next base*/
        } /*loop forwards though the homopolymer*/
    } /*else, bases are same and need to find length of forward bases*/
    
    /***************************************************************************
    # Fun-4 Sec-2 Sub-3: Check if deletion in hompolymer > max length
    ***************************************************************************/

    if(checkIfBasesMatch(baseCStr, 'A') == 1)
    { /*If was an Hompolymer of A's*/
       if(lenHomoULng > (*minReadStats).maxDelAHomoULng)
          ((*samStruct).numIgnoreDelULng)+=(*scoreStruct).cigarEntryULng;
        else
            ((*samStruct).numDelULng) += (*scoreStruct).cigarEntryULng;
    } /*If was an Hompolymer of A's*/

    else if(checkIfBasesMatch(baseCStr, 'T') == 1)
    { /*else If was an Hompolymer of T's*/
        if(lenHomoULng > (*minReadStats).maxDelTHomoULng)
          ((*samStruct).numIgnoreDelULng)+=(*scoreStruct).cigarEntryULng;
        else
            ((*samStruct).numDelULng) += (*scoreStruct).cigarEntryULng;
    } /*else If was an Hompolymer of T's*/

    else if(checkIfBasesMatch(baseCStr, 'C') == 1)
    { /*else If was an Hompolymer of C's*/
        if(lenHomoULng > (*minReadStats).maxDelCHomoULng)
          ((*samStruct).numIgnoreDelULng)+=(*scoreStruct).cigarEntryULng;
        else
            ((*samStruct).numDelULng) += (*scoreStruct).cigarEntryULng;
    } /*else If was an Hompolymer of C's*/

    else if(checkIfBasesMatch(baseCStr, 'G') == 1)
    { /*else If was an Hompolymer of G's*/
        if(lenHomoULng > (*minReadStats).maxDelGHomoULng)
          ((*samStruct).numIgnoreDelULng)+=(*scoreStruct).cigarEntryULng;
        else
            ((*samStruct).numDelULng) += (*scoreStruct).cigarEntryULng;
    } /*else If was an Hompolymer of G's*/

    /***************************************************************************
    # Fun-4 Sec-2 Sub-4: increment pointers
    ***************************************************************************/
    /*Do not need to increment pointers, since deletion not in the sequence*/

    return;
} /*checkDeletions*/


/*##############################################################################
# Name: CheckIfBasesMatch
# Use: checks if two bases are matches (ignores case)
# Input:
#    baseOneChar: first base to check (character)
#    baseTwoChar: base to compare the first base to (character)
# Output:
#    returns: 0 if no match
#             1 if is a match
# Note: Right know pretty simple, This is were you would add in anonymous base
#       checks
##############################################################################*/
char checkIfBasesMatch(char baseOneChar, char baseTwoChar)
{ /*checkIfBasesMatch*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5: checkIfBasesMatch
    #    fun-5 sec-1: Convert lower to upper case and replace U with T
    #    fun-5 sec-2: Check if bases are a match
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: Convert lower to upper case and replace U with T
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Convert to lower case*/
   if(baseOneChar > 90)
       baseOneChar = baseOneChar - 32;
   if(baseTwoChar > 90)
       baseTwoChar = baseTwoChar - 32;

   /*Convert U to T*/
   if(baseOneChar == 'U')
       baseOneChar = 'T';
   if(baseTwoChar == 'U')
       baseTwoChar = 'T';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Check if bases are a match
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(baseOneChar != baseTwoChar)
       return 0; /*No match*/

   return 1; /*match*/
} /*checkIfBasesMatch*/

/*##############################################################################
# Name: checkSoftMasks
# Use: Pretends to check soft masks, but really just increments the points 
# Input: 
#    scoreStruct: scoreReadStruct with sequence line, q-score line,
#                   and number of bases, (scoreReadStruct pointer)
# Output:
#    modifies: seqCStr to point at bast after soft mask
#    modifies: qCStr to point at q-score after soft mask
##############################################################################*/
void checkSoftMasks(scoreReadStruct *scoreStruct)
{ /*checkSoftMasks*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6: checkSoftMasks
    #    fun-6 sec-1: move off the soft mask
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #    fun-6 sec-1: move off the soft mask
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    (*scoreStruct).sequenceCStr += (*scoreStruct).cigarEntryULng;
    (*scoreStruct).qScoreCStr += (*scoreStruct).cigarEntryULng;

    return;
} /*checkSoftMask*/
