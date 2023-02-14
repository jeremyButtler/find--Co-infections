/*##############################################################################
# Name: scoreReadsRefCheckCigarEntry
# Use: Checks if bases should be kept for mismatches, matches, & deletions
#      using a reference
##############################################################################*/

#include "scoreReadsRefCheckCigarEntry.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    fun-1: refCheckMismatches: Check if mismatches should be kept
#    fun-2: refCheckMatches: Check if matches should be kept
#    fun-3: refCheckDeletions: Check if deletions should be kept
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept mismatches, ignored mismatches
##############################################################################*/
void refCheckMismatches(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep mismatch*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
) /*Checks if mismatch is valid or not & updates histogram if valid*/
{ /*refCheckMismatches*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: variable declerations
    #    fun-1 sec-2: Check if reference & sequence is missing the Q-score line
    #    fun-1 sec-3: Check if only the read has a Q-score line
    #    fun-1 sec-4: Check if only the reference has a Q-score line
    #    fun-1 sec-5: Check if both sequence & reference support mismatch
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Check if reference & sequence is missing the Q-score line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int
        seqQInt = 0,
        refQInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check if there is a Q-score line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(
       *scoreStruct->qScoreCStr == '*' &&
       *(scoreStruct->qScoreCStr++) == '\t' &&
       **refQCStr == 0
    ) { /*if no Q-score line*/
        samStruct->numMisULng=samStruct->numMisULng+scoreStruct->cigarEntryULng;
        scoreStruct->sequenceCStr += scoreStruct->cigarEntryULng; /*mv of miss*/
        samStruct->readAligLenULng++;    /*Add base to aligned length*/

        /*Incurement the reference sequence (check if reverse complement (16))*/
        if(samStruct->flagUInt & 16)
            refSeqCStr -= scoreStruct->cigarEntryULng;
        else
            refSeqCStr += scoreStruct->cigarEntryULng;

        return;
    } /*if no Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Check if only the read has a Q-score line
    #    fun-1 sec-3 sub-1: Check mismatches in the sequence (incurment pointer)
    #    fun-1 sec-3 sub-2: Incurment reference sequence pointer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-3 Sub-1: Check mismatches in the sequence (incurment pointer)
    ***************************************************************************/

    if(*refQCStr == 0)
    { /*if the reference does not have a Q-score line*/
        for(unsigned long uLngBase = 0;
            uLngBase < scoreStruct->cigarEntryULng;
            uLngBase++
        ) { /*for bases that have mismatches*/
            seqQInt = *scoreStruct->qScoreCStr - Q_ADJUST_CHAR; /*Get Q-score*/
    
            if(seqQInt < minReadStats->minQChar)
                samStruct->numIgnoreMisULng++;    /*Ignore, < min Q-score*/
            else
                samStruct->numMisULng++;          /*Keep, > min Q-score*/
    
            /*Get the stats*/
            scoreStruct->seqQHistULng[seqQInt]++;     /*For aligned median*/
            scoreStruct->totalQScoreULng += seqQInt; /*For aligned mean*/
            samStruct->readAligLenULng++;       /*incurment aligned length*/
    
            /*Move to the next base in the sequence*/
            scoreStruct->qScoreCStr++;

            scoreStruct->sequenceCStr++;

            /*******************************************************************
            # Fun-1 Sec-3 Sub-2: Incurment reference sequence pointer
            *******************************************************************/
    
            if(samStruct->flagUInt & 16) /*Reference is reverse complement*/
                refSeqCStr -= scoreStruct->cigarEntryULng;
            else
                refSeqCStr += scoreStruct->cigarEntryULng;
        } /*for bases that have mismatches*/

        return;
    } /*if the reference does not have a Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Check if only the reference Q-score line
    #    fun-1 sec-4 sub-1: Check if keep mismtach
    #    fun-1 sec-4 sub-2: Move onto the next base
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-4 Sub-1: Check if keep mismatche
    ***************************************************************************/

    if(*scoreStruct->qScoreCStr == '*' && *(scoreStruct->qScoreCStr++) == '\t')
    { /*if read does not have a Q-score line*/
        for(unsigned long uLngBase = 0;
            uLngBase < scoreStruct->cigarEntryULng;
            uLngBase++
        ) { /*for bases that have mismatches*/
            refQInt = **refQCStr - Q_ADJUST_CHAR; /*Get Q-score*/
    
            if(refQInt < minReadStats->minQChar)
                samStruct->numIgnoreMisULng++;      /*Ignore, < min Q-score*/
            else
                samStruct->numMisULng++;            /*Keep, > min Q-score*/
    
            /*Get the stats*/
            scoreStruct->seqQHistULng[refQInt]++;   /*For aligned median*/
            scoreStruct->totalQScoreULng += refQInt;/*For aligned mean*/
            samStruct->readAligLenULng++;           /*incurment aligned length*/

            /*******************************************************************
            # Fun-1 Sec-4 Sub-2: Move onto the next base
            *******************************************************************/

            if(samStruct->flagUInt & 16)
            { /*If reference is reverse complement of read*/
                (*refSeqCStr)--;
                (*refQCStr)--;
            } /*If reference is reverse complement of read*/
    
            else
            { /*Else refernce and read are in same direction*/
                (*refSeqCStr)++;
                (*refQCStr)++;
            } /*Else refernce and read are in same direction*/
    
            scoreStruct->sequenceCStr++;       /*Move to next base in sequence*/
        } /*for bases that have mismatches*/

        return;
    } /*if read does not have a Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-5: Check if both sequence & reference support mismatch
    #    fun-1 sec-5 sub-1: Check if mismatch supported
    #    fun-1 sec-5 sub-2: Move onto the next base
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-5 Sub-1: Check if mismatch supported
    ***************************************************************************/

    for(unsigned long uLngBase = 0;
        uLngBase < scoreStruct->cigarEntryULng;
        uLngBase++
    ) { /*for bases that have mismatches*/
        seqQInt = *scoreStruct->qScoreCStr - Q_ADJUST_CHAR; /*Get Q-score*/
        refQInt = **refQCStr - Q_ADJUST_CHAR; /*Get Q-score*/

        if(
           seqQInt < minReadStats->minQChar ||
           refQInt < minReadStats->minQChar
        ) /*If missmatch not support by reference & sequence*/
            samStruct->numIgnoreMisULng++;          /*Ignore, < min Q-score*/
        else
            samStruct->numMisULng++;                /*Keep, > min Q-score*/

        /*Get the stats (only for sequence median & mean Q)*/
        scoreStruct->seqQHistULng[seqQInt]++;       /*For aligned median*/
        scoreStruct->totalQScoreULng += seqQInt;    /*For aligned mean*/
        samStruct->readAligLenULng++;               /*incurment aligned length*/

        /*******************************************************************
        # Fun-1 Sec-5 Sub-2: Move onto the next base
        *******************************************************************/

        if(samStruct->flagUInt & 16)
        { /*If reference is reverse complement of read*/
            (*refSeqCStr)--;
            (*refQCStr)--;
        } /*If reference is reverse complement of read*/

        else
        { /*Else refernce and read are in same direction*/
            (*refSeqCStr)++;
            (*refQCStr)++;
        } /*Else refernce and read are in same direction*/

        scoreStruct->sequenceCStr++;       /*Move to next base in sequence*/
        scoreStruct->qScoreCStr++;         /*Move to next base in sequence*/
    } /*for bases that have mismatches*/

    return;
} /*refCheckMismatches*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept matches, total matches
##############################################################################*/
void refCheckMatches(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep mismatch*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
) /*Counts total matches, checks if matches are kept, & incurments pointers*/
{ /*refCheckMatches*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC:
    #    fun-2 sec-1: variable declerations
    #    fun-2 sec-2: Check if reference & sequence is missing the Q-score line
    #    fun-2 sec-3: Check if only the read has a Q-score line
    #    fun-2 sec-4: Check if only the reference has a Q-score line
    #    fun-2 sec-5: Check if both sequence & reference support match
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int
        seqQInt = 0,
        refQInt = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Check if has a Q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(
        *scoreStruct->qScoreCStr == '*' &&
        *(scoreStruct->qScoreCStr++) == '\t' &&
        **refQCStr == 0
    ) { /*If sequence & reference do not have a Q-score line*/

       /*Move to the bases after matches*/
       scoreStruct->sequenceCStr += scoreStruct->cigarEntryULng;

       /*Incurement the reference sequence (check if reverse complement (16))*/
       if(samStruct->flagUInt & 16)
           refSeqCStr -= scoreStruct->cigarEntryULng;
       else
           refSeqCStr += scoreStruct->cigarEntryULng;

       samStruct->numMatchULng++;
       samStruct->numKeptMatchULng++;
       samStruct->readAligLenULng++;         /*Add base to aligned length*/

       return;
    } /*If sequence & reference do not have a Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Check if only the read has a Q-score line
    #    fun-2 sec-3 sub-1: Check matches in the sequence (incurment pointer)
    #    fun-2 sec-3 sub-2: Incurment reference sequence pointer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-2 Sec-3 Sub-1: Check matches in the sequence (incurment pointer)
    ***************************************************************************/

    if(*refQCStr == 0)
    { /*if the reference does not have a Q-score line*/
        for(unsigned long uLngBase = 0;
            uLngBase < scoreStruct->cigarEntryULng;
            uLngBase++
        ) { /*for bases that have matches*/
            seqQInt = *scoreStruct->qScoreCStr - Q_ADJUST_CHAR; /*Get Q-score*/
            samStruct->numMatchULng++;

            if(seqQInt > minReadStats->minQChar)
                samStruct->numKeptMatchULng++;
    
            /*Get the stats*/
            scoreStruct->seqQHistULng[seqQInt]++;   /*For aligned median*/
            scoreStruct->totalQScoreULng += seqQInt;/*For aligned mean*/
            samStruct->readAligLenULng++;           /*incurment aligned length*/
    
            /*Move to the next base in the sequence*/
            scoreStruct->qScoreCStr++;
            scoreStruct->sequenceCStr++;

            /*******************************************************************
            # Fun-2 Sec-3 Sub-2: Incurment reference sequence pointer
            *******************************************************************/
    
            if(samStruct->flagUInt & 16) /*Reference is reverse complement*/
                refSeqCStr -= scoreStruct->cigarEntryULng;
            else
                refSeqCStr += scoreStruct->cigarEntryULng;
        } /*for bases that have matches*/

        return;
    } /*if the reference does not have a Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-4: Check if only the reference Q-score line
    #    fun-2 sec-4 sub-1: Check if should keep matches
    #    fun-2 sec-4 sub-2: Move onto the next base
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-2 Sec-4 Sub-1: Check if should keep matches
    ***************************************************************************/

    if(*scoreStruct->qScoreCStr == '*' && *(scoreStruct->qScoreCStr++) == '\t')
    { /*if read does not have a Q-score line*/
        for(unsigned long uLngBase = 0;
            uLngBase < scoreStruct->cigarEntryULng;
            uLngBase++
        ) { /*for bases that have matches*/
            refQInt = **refQCStr - Q_ADJUST_CHAR; /*Get Q-score*/
            samStruct->numMatchULng++;

            if(refQInt > minReadStats->minQChar)
                samStruct->numKeptMatchULng++;      /*Match worth keeping*/
    
            /*Get the stats*/
            scoreStruct->seqQHistULng[refQInt]++;   /*For aligned median*/
            scoreStruct->totalQScoreULng += refQInt;/*For aligned mean*/
            samStruct->readAligLenULng++;           /*incurment aligned length*/

            /*******************************************************************
            # Fun-2 Sec-4 Sub-2: Move onto the next base
            *******************************************************************/

            if(samStruct->flagUInt & 16)
            { /*If reference is reverse complement of read*/
                (*refSeqCStr)--;
                (*refQCStr)--;
            } /*If reference is reverse complement of read*/
    
            else
            { /*Else refernce and read are in same direction*/
                (*refSeqCStr)++;
                (*refQCStr)++;
            } /*Else refernce and read are in same direction*/
    
            scoreStruct->sequenceCStr++;       /*Move to next base in sequence*/
        } /*for bases that have matches*/

        return;
    } /*if read does not have a Q-score line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-5: Check if both sequence & reference support matches
    #    fun-2 sec-5 sub-1: Check if matches are supported
    #    fun-2 sec-5 sub-2: Move onto the next base
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-2 Sec-5 Sub-1: Check if matches are supported
    ***************************************************************************/

    for(unsigned long uLngBase = 0;
        uLngBase < scoreStruct->cigarEntryULng;
        uLngBase++
    ) { /*for bases that have matches*/
        seqQInt = *scoreStruct->qScoreCStr - Q_ADJUST_CHAR; /*Get Q-score*/
        refQInt = **refQCStr - Q_ADJUST_CHAR; /*Get Q-score*/
        samStruct->numMatchULng++;

        if(
           seqQInt >= minReadStats->minQChar &&
           refQInt >= minReadStats->minQChar
        ) /*If missmatch not support by reference & sequence*/
            samStruct->numKeptMatchULng++;

        /*Get the stats (only for sequence median & mean Q)*/
        scoreStruct->seqQHistULng[seqQInt]++;       /*For aligned median*/
        scoreStruct->totalQScoreULng += seqQInt;    /*For aligned mean*/
        samStruct->readAligLenULng++;               /*incurment aligned length*/

        /*******************************************************************
        # Fun-1 Sec-5 Sub-2: Move onto the next base
        *******************************************************************/

        if(samStruct->flagUInt & 16)
        { /*If reference is reverse complement of read*/
            (*refSeqCStr)--;
            (*refQCStr)--;
        } /*If reference is reverse complement of read*/

        else
        { /*Else refernce and read are in same direction*/
            (*refSeqCStr)++;
            (*refQCStr)++;
        } /*Else refernce and read are in same direction*/

        scoreStruct->sequenceCStr++;       /*Move to next base in sequence*/
        scoreStruct->qScoreCStr++;         /*Move to next base in sequence*/
    } /*for bases that have mismatches*/

    return;
} /*refCheckMatches*/

/*##############################################################################
# Output:
#    modifies: scoreReadStrut Q-score histogram, running Q-score total
#    incurments: sequence and Q-scores lines for read & reference
#    modifies: samStruct aligned length, kept mismatches, ignored insertions
# Note:
#    This function acts a little oddly, in that it only keeps a deletion if 
#      the reference has to little support. This is tricky, since I have no
#      Q-score for the read here
##############################################################################*/
void refCheckDeletions(
    struct scoreReadStruct *scoreStruct, /*Holds Q-score histogram & total*/
    char **refSeqCStr,                   /*Points to references sequence*/
    char **refQCStr,                     /*Points to references Q-score line*/
    struct minStats *minReadStats,       /*Holds min Q-score to keep insertion*/
    struct samEntry *samStruct           /*Holds reads sequence, Q, & cigar*/
) /*Uses refernce to checks if deletion is supported*/
{ /*checkInsertions*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #     fun-3 sec-1: variable declerations
    #     fun-3 sec-2: Find if is a valid indel or not
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpBaseCStr = 0;
    int refQInt = 0;
    unsigned long lenHomoULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Find if is a valid deletion or not
    #    fun-3 sec-2 sub-1: Check Q-scores eliminate insertion
    #    fun-3 sec-2 sub-2: Find the homopolymer length
    #    fun-3 sec-2 sub-3: Check if A insertion is > max A homopolymer length
    #    fun-3 sec-2 sub-4: Check if T insertion is > max T homopolymer length
    #    fun-3 sec-2 sub-5: Check if G insertion is > max G homopolymer length
    #    fun-3 sec-2 sub-6: Check if C insertion is > max C homopolymer length
    #    fun-3 sec-2 sub-7: Increment pointers and record Q-scores
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-3 Sec-2 Sub-1: Check Q-scores eliminate insertion
    ***************************************************************************/

    for(int intIndel = 0;
        intIndel < (*scoreStruct).cigarEntryULng;
        intIndel++
       )
    { /*For all deletions*/
        lenHomoULng = 1;             /*Have at least one base (the insertion)*/
        refQInt = **refQCStr - Q_ADJUST_CHAR;  /*get Q-score*/

        if(
           **refQCStr != 0 &&
           refQInt >= minReadStats->minQChar
          ) { /*if reference has strong support for its insertion*/

            samStruct->numIgnoreDelULng++;       /*Insertion > min Q*/

            if(samStruct->flagUInt & 16)
            { /*If reference is reverse complement of read*/
                (*refSeqCStr)--;
                (*refQCStr)--;
            } /*If reference is reverse complement of read*/
    
            else
            { /*Else refernce and read are in same direction*/
                (*refSeqCStr)++;
                (*refQCStr)++;
            } /*Else refernce and read are in same direction*/

            continue;
        } /*if reference has strong support for its insertion*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-2: Find the homopolymer length
        ***********************************************************************/

        /*for non-reverse complement, this would be the base behind insert*/
        tmpBaseCStr = *refSeqCStr - 1;

        while(
            checkIfBasesMatch(
                *scoreStruct->sequenceCStr,
                *tmpBaseCStr
            ) > 0
        ) { /*While there are bases in the homopolymer*/
            lenHomoULng++;                       /*Count bases in homopolymer*/
            tmpBaseCStr--;                       /*Move to next base*/
        } /*While there are bases in the homopolymer*/

        /*for non-reverse complement, this would be the base after insert*/
        tmpBaseCStr = *refSeqCStr + 1;

        while(
            checkIfBasesMatch(
                *scoreStruct->sequenceCStr,
                *tmpBaseCStr
            ) > 0
        ) { /*While there are bases in the homopolymer*/
            lenHomoULng++;                       /*Count bases in homopolymer*/
            tmpBaseCStr++;                       /*Move to next base*/
        } /*While there are bases in the homopolymer*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-3: Check if A insertion is > max A homopolymer length
        ***********************************************************************/

        if(
           checkIfBasesMatch(
               **refSeqCStr,
               'A'
           ) == 1
        ) { /*If was an Hompolymer of A's*/
           if(samStruct->flagUInt & 16)
           { /*If reference is reverse complement*/
               if(lenHomoULng > minReadStats->maxDelTHomoULng)
                    samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /*If reference is reverse complement*/

           else
           { /* Else reference is not reverse complement*/
               if(lenHomoULng > minReadStats->maxDelAHomoULng)
                   samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /* Else reference is not reverse complement*/
        } /*If was an Hompolymer of A's*/
    
        /***********************************************************************
        # Fun-3 Sec-2 Sub-4: Check if T insertion is > max T homopolymer length
        ***********************************************************************/

        else if(
           checkIfBasesMatch(
               **refSeqCStr,
               'T'
           ) == 1
        ) { /*If was an Hompolymer of T's*/
           if(samStruct->flagUInt & 16)
           { /*If reference is reverse complement*/
               if(lenHomoULng > minReadStats->maxDelAHomoULng)
                    samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /*If reference is reverse complement*/

           else
           { /* Else reference is not reverse complement*/
               if(lenHomoULng > minReadStats->maxDelTHomoULng)
                   samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /* Else reference is not reverse complement*/
        } /*If was an Hompolymer of T's*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-5: Check if G insertion is > max G homopolymer length
        ***********************************************************************/

        else if(
           checkIfBasesMatch(
               **refSeqCStr,
               'G'
           ) == 1
        ) { /*If was an Hompolymer of G's*/
           if(samStruct->flagUInt & 16)
           { /*If reference is reverse complement*/
               if(lenHomoULng > minReadStats->maxDelCHomoULng)
                    samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /*If reference is reverse complement*/

           else
           { /* Else reference is not reverse complement*/
               if(lenHomoULng > minReadStats->maxDelGHomoULng)
                   samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /* Else reference is not reverse complement*/
        } /*If was an Hompolymer of G's*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-6: Check if C insertion is > max C homopolymer length
        ***********************************************************************/

        else if(
           checkIfBasesMatch(
               **refSeqCStr,
               'C'
           ) == 1
        ) { /*If was an Hompolymer of C's*/
           if(samStruct->flagUInt & 16)
           { /*If reference is reverse complement*/
               if(lenHomoULng > minReadStats->maxDelGHomoULng)
                    samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /*If reference is reverse complement*/

           else
           { /* Else reference is not reverse complement*/
               if(lenHomoULng > minReadStats->maxDelCHomoULng)
                   samStruct->numIgnoreDelULng++;
                else
                    samStruct->numDelULng++;
           } /* Else reference is not reverse complement*/
        } /*If was an Hompolymer of C's*/

        /***********************************************************************
        # Fun-3 Sec-2 Sub-7: increment pointers and record Q-score
        ***********************************************************************/

        if(samStruct->flagUInt & 16)
        { /*If reference is reverse complement of read*/
            (*refSeqCStr)--;
            (*refQCStr)--;
        } /*If reference is reverse complement of read*/

        else
        { /*Else refernce and read are in same direction*/
            (*refSeqCStr)++;
            (*refQCStr)++;
        } /*Else refernce and read are in same direction*/
    } /*Loop through all insertions*/

    return;
} /*refCheckDeletions*/
