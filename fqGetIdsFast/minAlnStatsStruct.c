/*######################################################################
# Use:
#   o Holds the minAlnStats structer and its supporting functions.
# Includes:
#   o defaultSettings.h
######################################################################*/

#include "minAlnStatsStruct.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: scoreReadsFun
'    fun-1 blankMinStats: set a minAlnStats structer to defaults
'    fun-2 blankMinStatsReadRead: For read/read settings
'    fun-3 blankMinStatsReadCon: For read/consensus settings
'    fun-4 blankMinStatsConCon: For consensus/consensus settings
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# Uses default read to reference mapping settings from defaultSettings.h
######################################################################*/
void blankMinStats(
    struct minAlnStats *minStats
) /*Sets minStats minimum requirements for sam alingemtns to defaults*/
{ /*blankMinStats*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1 Sub-1 TOC: blankMinStats
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    minStats->minMapqUInt = readRefMapq; /*min mapping quality*/
    minStats->minQChar = readRefMinBaseQ; /*default min Q-score*/
    minStats->minMedianQFlt = readRefMinMedQ;/*default median Q-score*/
    minStats->minMeanQFlt = readRefMinMeanQ; /*default mean Q-score*/
    minStats->minAlignedMedianQFlt = readRefMinAlnMedQ;/*med aligend Q*/
    minStats->minAlignedMeanQFlt = readRefMinAlnMeanQ;/*mean aligend Q*/
    minStats->maxReadLenULng = readRefMaxReadLen; /*max read length*/
    minStats->minReadLenULng = readRefMinReadLen; /*min read length*/
    
    for(uint8_t uCharCnt = 0; uCharCnt < 16; ++uCharCnt)
    { /*loop till have initialized the deltion and insertion arrays*/
        minStats->maxHomoInsAry[uCharCnt] = 0;
        minStats->maxHomoDelAry[uCharCnt] = 0;
    } /*loop till have initialized the deltion and insertion arrays*/

    minStats->maxHomoInsAry[0] = readRefMaxInsAHomo;
        /*Discard A insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[10] = readRefMaxInsTHomo;
        /*Discard T insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[1] = readRefMaxInsCHomo;
        /*Discard C insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[3] = readRefMaxInsGHomo;
        /*Discard G insertion if homopoymer in is larger than setting*/

    minStats->maxHomoDelAry[0] = readRefMaxDelAHomo;
        /*Discard A deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[10] = readRefMaxDelTHomo;
        /*Discard T deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[1] = readRefMaxDelCHomo;
        /*Discard C deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[3] = readRefMaxDelGHomo;
        /*Discard G deletion if homopoymer in is larger than setting*/

    /*These settings are for findCoInft*/
    minStats->minSNPsFlt = readRefMinPercSNPs;
    minStats->minDelsFlt = readRefMinPercDels;
    minStats->minInssFlt = readRefMinPercInss;
    minStats->minIndelsFlt = readRefMinPercIndels;
    minStats->minDiffFlt = readRefMinPercDiff;
    return;
} /*blankMinStats*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# Uses default read to read mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsReadRead(
    struct minAlnStats *minStats
) /*Sets minStats minimum requirements for sam alingemtns to defaults*/
{ /*blankMinStatsReadRead*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 Sub-1 TOC: blankMinStatsReadRead
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    minStats->minMapqUInt = readReadMapq;
        /*Mininum mapping quality to keep a read*/
    minStats->minQChar = readReadMinBaseQ;
         /*Minimum Q-score to keep a base or insertion*/
    minStats->minMedianQFlt = readReadMinMedQ;
        /*Minimum median Q-score to keep a read*/
    minStats->minMeanQFlt = readReadMinMeanQ;
        /*Minimum mean Q-score to keep a read*/
    minStats->minAlignedMedianQFlt = readReadMinAlnMedQ;
        /*Minimum median aligend Q-score to keep a read*/
    minStats->minAlignedMeanQFlt = readReadMinAlnMeanQ;
        /*Minimum mean aligend Q-score to keep a read*/
    minStats->maxReadLenULng = readReadMaxReadLen; /*max read length*/
    minStats->minReadLenULng = readReadMinReadLen; /*min read length*/
    
    for(uint8_t uCharCnt = 0; uCharCnt < 16; ++uCharCnt)
    { /*loop till have initialized the deltion and insertion arrays*/
        minStats->maxHomoInsAry[uCharCnt] = 0;
        minStats->maxHomoDelAry[uCharCnt] = 0;
    } /*loop till have initialized the deltion and insertion arrays*/

    minStats->maxHomoInsAry[0] = readReadMaxInsAHomo;
        /*Discard A insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[10] = readReadMaxInsTHomo;
        /*Discard T insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[1] = readReadMaxInsCHomo;
        /*Discard C insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[3] = readReadMaxInsGHomo;
        /*Discard G insertion if homopoymer in is larger than setting*/

    minStats->maxHomoDelAry[0] = readReadMaxDelAHomo;
        /*Discard A deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[10] = readReadMaxDelTHomo;
        /*Discard T deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[1] = readReadMaxDelCHomo;
        /*Discard C deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[3] = readReadMaxDelGHomo;
        /*Discard G deletion if homopoymer in is larger than setting*/

    /*These settings are for findCoInft*/
    minStats->minSNPsFlt = readReadMinPercSNPs;
    minStats->minDelsFlt = readReadMinPercDels;
    minStats->minInssFlt = readReadMinPercInss;
    minStats->minIndelsFlt = readReadMinPercIndels;
    minStats->minDiffFlt = readReadMinPercDiff;
    return;
} /*blankMinStatsReadRead*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# Uses default read to consensus mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsReadCon(
    struct minAlnStats *minStats
) /*Sets minStats minimum requirements for sam alingemtns to defaults*/
{ /*blankMinStatsReadCon*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 Sub-1 TOC: blankMinStatsReadCon
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    minStats->minMapqUInt = readConMapq;
        /*Mininum mapping quality to keep a read*/
    minStats->minQChar = readConMinBaseQ;
         /*Minimum Q-score to keep a base or insertion*/
    minStats->minMedianQFlt = readConMinMedQ;
        /*Minimum median Q-score to keep a read*/
    minStats->minMeanQFlt = readConMinMeanQ;
        /*Minimum mean Q-score to keep a read*/
    minStats->minAlignedMedianQFlt = readConMinAlnMedQ;
        /*Minimum median aligend Q-score to keep a read*/
    minStats->minAlignedMeanQFlt = readConMinAlnMeanQ;
        /*Minimum mean aligend Q-score to keep a read*/
    minStats->maxReadLenULng = readConMaxReadLen; /*max read length*/
    minStats->minReadLenULng = readConMinReadLen; /*min read length*/
    
    for(uint8_t uCharCnt = 0; uCharCnt < 16; ++uCharCnt)
    { /*loop till have initialized the deltion and insertion arrays*/
        minStats->maxHomoInsAry[uCharCnt] = 0;
        minStats->maxHomoDelAry[uCharCnt] = 0;
    } /*loop till have initialized the deltion and insertion arrays*/

    minStats->maxHomoInsAry[0] = readConMaxInsAHomo;
        /*Discard A insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[10] = readConMaxInsTHomo;
        /*Discard T insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[1] = readConMaxInsCHomo;
        /*Discard C insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[3] = readConMaxInsGHomo;
        /*Discard G insertion if homopoymer in is larger than setting*/

    minStats->maxHomoDelAry[0] = readConMaxDelAHomo;
        /*Discard A deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[10] = readConMaxDelTHomo;
        /*Discard T deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[1] = readConMaxDelCHomo;
        /*Discard C deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[3] = readConMaxDelGHomo;
        /*Discard G deletion if homopoymer in is larger than setting*/

    /*These settings are for findCoInft*/
    minStats->minSNPsFlt = readConMinPercSNPs;
    minStats->minDelsFlt = readConMinPercDels;
    minStats->minInssFlt = readConMinPercInss;
    minStats->minIndelsFlt = readConMinPercIndels;
    minStats->minDiffFlt = readConMinPercDiff;
    return;
} /*blankMinStatsReadCon*/

/*######################################################################
# Output:
#    modifes minStats to have default entries
# default consensus to consensus mapping settings from defaultSettings.h
######################################################################*/
void blankMinStatsConCon(
    struct minAlnStats *minStats
) /*Sets minStats minimum requirements for sam alingemtns to defaults*/
{ /*blankMinStatsConCon*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 Sub-1 TOC: blankMinStatsConCon
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    minStats->minMapqUInt = conConMapq;
        /*Mininum mapping quality to keep a con*/
    minStats->minQChar = conConMinBaseQ;
         /*Minimum Q-score to keep a base or insertion*/
    minStats->minMedianQFlt = conConMinMedQ;
        /*Minimum median Q-score to keep a con*/
    minStats->minMeanQFlt = conConMinMeanQ;
        /*Minimum mean Q-score to keep a con*/
    minStats->minAlignedMedianQFlt = conConMinAlnMedQ;
        /*Minimum median aligend Q-score to keep a con*/
    minStats->minAlignedMeanQFlt = conConMinAlnMeanQ;
        /*Minimum mean aligend Q-score to keep a con*/
    minStats->maxReadLenULng = conConMaxReadLen; /*max con length*/
    minStats->minReadLenULng = conConMinReadLen; /*min con length*/
    
    for(uint8_t uCharCnt = 0; uCharCnt < 16; ++uCharCnt)
    { /*loop till have initialized the deltion and insertion arrays*/
        minStats->maxHomoInsAry[uCharCnt] = 0;
        minStats->maxHomoDelAry[uCharCnt] = 0;
    } /*loop till have initialized the deltion and insertion arrays*/

    minStats->maxHomoInsAry[0] = conConMaxInsAHomo;
        /*Discard A insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[10] = conConMaxInsTHomo;
        /*Discard T insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[1] = conConMaxInsCHomo;
        /*Discard C insertion if homopoymer in is larger than setting*/
    minStats->maxHomoInsAry[3] = conConMaxInsGHomo;
        /*Discard G insertion if homopoymer in is larger than setting*/

    minStats->maxHomoDelAry[0] = conConMaxDelAHomo;
        /*Discard A deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[10] = conConMaxDelTHomo;
        /*Discard T deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[1] = conConMaxDelCHomo;
        /*Discard C deletion if homopoymer in is larger than setting*/
    minStats->maxHomoDelAry[3] = conConMaxDelGHomo;
        /*Discard G deletion if homopoymer in is larger than setting*/

    /*These settings are for findCoInft*/
    minStats->minSNPsFlt = conConMinPercSNPs;
    minStats->minDelsFlt = conConMinPercDels;
    minStats->minInssFlt = conConMinPercInss;
    minStats->minIndelsFlt = conConMinPercIndels;
    minStats->minDiffFlt = conConMinPercDiff;
    return;
} /*blankMinStatsConCon*/


