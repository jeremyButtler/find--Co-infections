/*######################################################################
# Use:
#   o Holds the Median Q-score functions for scoreReads & findCoInft
######################################################################*/

#include "FCIStatsFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: FCIStatsFun
'    fun-1 findQScores: find Q-score of read using q-score entry
'    fun-2 qHistToMed: Find median Q-score from a q-score histogram
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
    # Fun-1 Sec-1 Sub-1 TOC: findQScores
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
    # Fun-2 Sec-1 Sub-1 TOC: qHistToMedian
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
