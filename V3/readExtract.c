/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF TOC: Start of Functions
'   fun-1 extractBestRead:
'     - Extract the best read from a bin
'   fun-2 fqOneIdExtract:
'     - Extract one read from a fastq file
'   fun-3 findBestXReads:
'     - Extract the top x (user specified) reads that mapped to the
'       selected best read.
'   fun-4 fqGetBestReadByMedQ:
'     - Extracts read with best medain Q-score from file.
'     - It also considers length if integer Q-scores are the same.
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "readExtract.h"

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Prints:                                                         |
|        - The read with the highest medain Q-score to the fastq file  |
|          name stored in binIn->bestReadCStr                          |
|    - Modifies:                                                       |
|        - Fastq file binIn-fqPathCStr to not have the best read       |
|        - Stats file binIn-statsPathCStr to not have the best read    |
|    - Returns:                                                        |
|        - 1: If sucessfull                                            |
|        - 2: For blank structer                                       |
|        - 4: For no fastq file                                        |
|        - 8: For no stats file                                        |
|        - 16: For error when extracting stats                         |
|        - 32: If could not open a temporary fastq file                |
|        - 128: If could not copy the fastq reads to their bins        |
\---------------------------------------------------------------------*/
uint8_t extractBestRead(
    struct readBin *binIn        /*Bin to extract best read from*/
) /*Splits bin into read with highest Q-score & other all other reads*/
{ /*extractBestRead*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: extractBestRead
    '    fun-1 sec-1: variable declerations                            \
    '    fun-1 sec-2: Check if Bin exists & set up best read file name /
    '    fun-1 sec-3: Check if fastq & stats file can be opened        \
    '    fun-1 sec-4: Open temporary stat file & write header          /
    '    fun-1 sec-5: Copy old stat file, except for best read         \
    '    fun-1 sec-6: Make temporay stat file the new bin stat file    /
    '    fun-1 sec-7: Extract the best read                            \
    '    fun-1 sec-8: Clean up, close files & make tmp file bin fq file/
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: variable declerations                               v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpFqCStr ="tmp-202301110827-reads-09876543211234567890.fastq",
        *tmpStatCStr ="tmp-202301110827-stats-09876543211234567890.tsv",
        *tmpCStr = 0;

    int8_t ignoreC = 0;

    uint8_t
        errUC = 0,     /*Holds error messages*/
        onHeaderBool = 1; /*Tells if need to ignore the first line*/

    struct readStat
        tmpRead,
        bestRead;            /*Holds best read to extract*/

    FILE 
        *inFILE = 0,
        *outFILE = 0,      /*Holds the read to polish with*/
        *otherOutFILE = 0; /*Holds every read except the polish read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-1 Sec-2: Check if Bin exists & set up best read file name    v
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(binIn == 0)
        return 2;    /*The structure has nothing*/

    strcpy(binIn->bestReadCStr, binIn->fqPathCStr); /*Copy the file name*/

    tmpCStr = binIn->bestReadCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    tmpCStr -= 6; /*Get on '.' part of ".fastq"*/

    strcpy(
        tmpCStr,
        "--best-read.fastq"
    ); /*Add in the best read ending*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-1 Sec-3: Check if can open fastq file & stats file           v
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    inFILE = fopen(binIn->statPathCStr, "r"); /*Open the stats file*/

    if(inFILE == 0)
        return 4;    /*The fastq file has nothing*/

    outFILE = fopen(binIn->fqPathCStr, "r"); /*Open the fastq file*/

    if(outFILE == 0)
    { /*If can not open the fastq file*/
        fclose(outFILE);
        return 8;
    } /*If can not open the fastq file*/

    fclose(outFILE); /*Was just a quick test to see if could open*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-4: Open temporary stat file & write header             v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the temporary stat file to write non-best read to*/
    outFILE = fopen(tmpStatCStr, "w");
    printStatHeader(outFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-5: Copy old stat file, except for best read            v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /*Read in the frist line (assume is best read*/
    errUC = readStatsFileLine(inFILE, &onHeaderBool, &bestRead);
        /*onHeaderBool tells if have read in the header*/

    if(!(errUC & 1))
    { /*If had a file error*/
        fclose(inFILE);
        return 16;     /*File error, return 16 so user knows*/
    } /*If had a file error*/

    /*Read in next line, so have something to compare*/
    errUC = readStatsFileLine(inFILE, &onHeaderBool, &tmpRead);

    if(!(errUC & 1))
    { /*If had a file error*/
        fclose(inFILE);
        return 16;     /*File error, return 16 so user knows*/
    } /*If had a file error*/
   
    while(errUC & 1)
    { /*While not at end of file or no problems*/
        ignoreC =
            (int8_t) (
              ((short) tmpRead.mapqUChar - (short) bestRead.mapqUChar)
              >> (sizeof(short) << 3)
            );
            /*tmpRead - bestRead is negative if bestRead > tmpRead
                - converting to hsort, so can have negatives
              (tmp - best) >> bitsInSht: Keeps only the negative bit
            */

        ignoreC |=
            (int8_t) (
                (
                  (int32_t) tmpRead.medianQFlt -
                  (int32_t) bestRead.medianQFlt
                )
              >> ((sizeof(int32_t) << 3) - 1) /*One off max size*/
            ) & ignoreC;
            /* tmp - best: is negative if best > tmp
               (tmp - best) >> 32: 1 if best > tmp; 0 if tmp>best
                   - change to sizeof(.* to deal with complier warnings
               & ingnoreUC: sets to 0 if mapq was best
            */

        ignoreC |=
            (int8_t) (
                (
                  (int32_t) tmpRead.readLenUInt -
                  (int32_t) bestRead.readLenUInt
                )
              >> ((sizeof(int32_t) << 3) - 1) /*One of max size*/
            ) & ignoreC;
            /* (int64_t) tmp - (int64_t) best: is negative if best > tmp
               (tmp - best) >> 64: 1 if negative, 0 if positive
               & ingnoreUC: sets to 0 if mapq was best
            */

        /*At this point ignoreC is 1 (discard) read or 0 (keep)*/
        if(ignoreC & 1)
            printReadStat(&tmpRead, outFILE);
        else
        { /*else have a new best mapq*/
            printReadStat(&bestRead, outFILE);
            cpReadStat(&bestRead, &tmpRead); /*copy the read*/
        } /*else have a new best mapq*/

        errUC =
            readStatsFileLine(
                inFILE,     /*File with line to grab*/
                &onHeaderBool, /*Tells if their is a header line*/
                &tmpRead   /*Will hold the stats from the stats file*/
        ); /*Read in the next line*/
    } /*While not at end of file or no problems*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-6: Make the temporay stat file the new bin stat file   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*Close the open files (no longer need to read from)*/
     fclose(inFILE);
     fclose(outFILE);

    if(errUC != 2) 
        return 16; /*File error*/

     remove(binIn->statPathCStr); /*Remove the old stats file*/

     /*Make temporary stats file without best read the stats file*/
     rename(tmpStatCStr, binIn->statPathCStr);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-7: Extract the best read                               v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the fastq file (already check if could open)*/
    inFILE = fopen(binIn->fqPathCStr, "r");

    if(inFILE == 0)
        return 32;

    /*Open the temporary file to hold the best read*/
    outFILE = fopen(binIn->bestReadCStr, "w");

    if(outFILE == 0)
    { /*If  could not open the output file*/
        fclose(inFILE);
        return 32;
    } /*If  could not open the output file*/

    /*File to hold all reads except the best read*/
    otherOutFILE = fopen(tmpFqCStr, "w");

    if(otherOutFILE == 0)
    { /*If could not open the other out file*/
        fclose(inFILE);
        fclose(outFILE);
        return 32;
    } /*If could not open the other out file*/

    errUC = 
        fqOneIdExtract(
            bestRead.queryIdCStr,
            inFILE,     /*fastq file to extract read from*/
            outFILE,
            otherOutFILE
    ); /*See if I can extract the read*/

    /*No longer need files open*/
    fclose(inFILE);
    fclose(outFILE);
    fclose(otherOutFILE);

    if(!(errUC & 1)) /*If could not extract the best read*/
        return 128;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-8: Clean up, close files & make tmp file bin fq file   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    --binIn->numReadsULng; /*Account for the extracted best read*/

    remove(binIn->fqPathCStr); /*Remove the old file*/

    /*Make the tempoaray fastq with all but teh best read to
      the bin fastq file*/
    rename(tmpFqCStr, binIn->fqPathCStr);

    return 1; /*Success*/
} /*extractBestRead*/

/*----------------------------------------------------------------------
# Output:
#    Prints: fastq entry to outFILE if finds
#    Returns:
#        - 1: if found and printed the id
#        - 2: If could not find the id (no printing)
#        - 4: If the fqFILE does not exist
#        - 8: If the outFILE does not exist
#        - 16: If an incomplete entry (EOF, but missing lines)
#        - 32: If the idCStr (id looking for) entry is incomplete
----------------------------------------------------------------------*/
uint8_t fqOneIdExtract(
    char *idCStr,  /*'\0' terminated read id to extract*/
    FILE *fqFILE,     /*fastq file to extract read from*/
    FILE *keptFILE,   /*File with the target read*/
    FILE *outFILE     /*fastq file to write read to*/
) /*Extracts one read id from a fastq file*/
{ /*fqOneIdExtract*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: fqOneIdExtract
    '    fun-2 sec-1: variable declerations
    '    fun-2 sec-2: Check if valid file and set up for read in
    '    fun-2 sec-3: Read in the header and compare read ids
    '    fun-2 sec-4: If not a match, move past sequence & spacer line
    '    fun-2 sec-5: If not a match, move past q-score lines
    '    fun-2 sec-6: Is match, print out header
    '    fun-2 sec-7: Is match, print sequence and spacer line
    '    fun-2 sec-8: Is match, print Q-score line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint16_t
        numSeqLineUSht = 0, /*Number of sequence lines in fastq file*/
        lenBuffUSht = 1024;

    char
        *tmpCStr = 0,
        *idIterCStr = 0,
        lineCStr[lenBuffUSht];

    FILE
        *tmpFILE = outFILE; /*File to write read to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Check if valid file and set up for read in
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(fqFILE == 0)
        return 4;
    if(outFILE == 0)
        return 8;
    if(keptFILE == 0)
        return 16;

    /*Make sure have marks set to determine if read in full line*/
    lineCStr[lenBuffUSht - 1] = '\0';
    lineCStr[lenBuffUSht - 2] = '\0';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Read in the header and compare read ids
    #    fun-2 sec-3 sub-1: Read in header and compare ids
    #    fun-2 sec-3 sub-2: Finsh reading in header if not a match
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-2 Sec-3 Sub-1: Read in header and compare ids
    *******************************************************************/

    while(fgets(lineCStr, lenBuffUSht, fqFILE))
    { /*While have a line to read in*/
        tmpCStr = lineCStr;
        idIterCStr = idCStr;

        /*Make sure off @ symbols for header*/
        if(*tmpCStr == '@')
            ++tmpCStr;

        if(*idIterCStr == '@')
            ++idIterCStr;

        while(*tmpCStr == *idIterCStr)
        { /*While the two strings are the same*/
            ++tmpCStr;
            ++idIterCStr;

            if(*idIterCStr == '\0') /*idIterCStr always ends with null*/
                break; /*Done with comparision*/
        } /*While the two strings are the same*/

        if(*idIterCStr == '\0')
            tmpFILE = keptFILE;
        else
            tmpFILE = outFILE;

        /*Print out the read in part of the header*/
        fprintf(tmpFILE, "%s", lineCStr);

        /***************************************************************
        # Fun-2 Sec-3 Sub-2: Finsh reading in header if not a match
        ***************************************************************/

        /*The header could have tabs, so can not do < 16*/
        /*I can not just seq till +\n, because Illumina fastq's break 
          into two lines, So I need to find the number of sequence lines
        */
        while(
            !(lineCStr[lenBuffUSht - 2] == '\0' ||
            lineCStr[lenBuffUSht - 2] == '\n')
        ) { /*While not on the next line*/
            lineCStr[lenBuffUSht - 2] = '\0';

            /*Grab next part of header from file*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/

            /*Print out next part of the header*/
            fprintf(tmpFILE, "%s", lineCStr);
        } /*While not on the next line*/

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Fun-2 Sec-4: If not a match, move past sequence & spacer line
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        lineCStr[lenBuffUSht - 2] = '\0'; /*Set up the maker*/

        /*Grab the first part of the sequence*/
        tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);
        numSeqLineUSht = 0; /*reset to 0 (number sequence lines)*/

        do { /*While not on the spacer entry*/

            /*Print out the read in part of the header*/
            fprintf(tmpFILE, "%s", lineCStr);

            /*I can get away with < 16 here because I know their are 
              no tabs in the sequence entry*/
            if(lineCStr[lenBuffUSht - 2] < 16)
                ++numSeqLineUSht; /*Record the number of lines read in*/

            lineCStr[lenBuffUSht - 2] = '\0'; /*Reset EOL marker*/

            /*Get next part of entry before spacer*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/

        } while(lineCStr[0] != '+'); /*While not on spacer*/
        /*fastq sequence line ends with \n+\nq-score line*/

        fwrite("+\n", sizeof(uint8_t), 2, tmpFILE);

        while(lineCStr[lenBuffUSht - 2] > 16)
        { /*While have extra heade entries*/
            lineCStr[lenBuffUSht - 2] = '\0'; /*reset EOL marker*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/
        } /*While have extra heade entries*/
        /*Just to remove extra user formating on spacer*/

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Fun-2 Sec-5: If not a match, move past q-score lines
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        for(uint16_t uShtCnt = 0; uShtCnt < numSeqLineUSht; ++uShtCnt)
        { /*Loop through all lines in the q-socre entry*/

             /*the Q-score entry will not have tabs, so < 16 works*/
        
             do { /*While not on the next line*/
                lineCStr[lenBuffUSht - 2] = '\0'; /*Reset marker*/

                /*Get next part of Q-score entry*/
                tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

                if(tmpCStr == 0)
                     return 32;    /*Premature end of file*/

                /*Print out the next part of the q-score entry*/
                fprintf(tmpFILE, "%s", lineCStr);
            } while(lineCStr[lenBuffUSht - 2] > 16); /*On next line?*/
        } /*Loop through all lines in the q-socre entry*/
    } /*While have a line to read in*/ 

    return 1; /*Sucess*/
} /*fqOneIdExtract*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|     - 1: if succeeded or no reads mapped (numReadsKept=0 for no reads)
|     - 4: if could not read reference
|     - 8: if could not open the fastq file
|     - 16: if minimap2 errored out or returned nothing
| Note:
|     - Score: percMult * (keptSNPs + keptIns + keptDels) / read length
|     - minSimUSht ranges from 1 (0.01%) to precMult (100%)            
\---------------------------------------------------------------------*/
uint8_t findBestXReads(
    const uint64_t *numReadConsULng, /*# reads for bulding a consensus*/
    uint64_t *numReadsKeptULng,  /*Number of reads binned to con*/
    char *threadsCStr,           /*Number threads to use with minimap2*/
    const char *useMapqBl,       /*1: use mapping quality in selection*/
    struct minAlnStats *minStats,/*Min stats to cluster reads together*/
    struct samEntry *samST,  /*Struct to use for reading sam file*/
    struct samEntry *refStruct,  /*holds the reference (0 to ignore)*/
    struct readBin *binTree,     /*Bin working on*/
    char noRefBl,                /*Only use fastq file*/
    char makeNameBl
       /*1: Make a name for the file using the fastq file; 0 use 
         binTree->topReadsCStr ('\0' for stdout)*/
) /*Extract the top reads that mapped to the selected best read*/
{ /*findBestXReads*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   \\ Fun-3 TOC: findBestXReads                                       /
   //    fun-3 sec-1: variable declerations                           \
   \\    fun-3 sec-2: Check if can open files                         /
   //    fun-3 sec-3: Map reads to best read & select top x           \
   \\    fun-3 sec-4: Trim, score, & select best mapped reads         /
   //    fun-3 sec-5: Set up the best x read file name                \
   \\    fun-3 sec-6: Set up hash table to extract reads              /
   //    fun-3 sec-7: Extract reads with fastq greps hash extract     \
   \\    fun-3 sec-8: Clean up                                        /
   //                                                                  \
   \\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-1: variable declerations                              v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int32_t lenBuffUInt = 1 << 15; /*wil make a 65536 byte array*/

    uint8_t errUC = 0; /*Holds error output*/
    uint8_t oneUChar = 1;
    uint8_t zeroUChar = 0;
    unsigned char lenStackUC = 128;
    unsigned char lenBigNum = 0;

    char minimapCmdCStr[2048];
    char *tmpCStr = 0;
    char buffCStr[lenBuffUInt];  /*Buffer to extract reads with*/

    int32_t lenIdUInt = 100;/*number of characters allowed for read id*/

    uint16_t topScoreUSht = 128 + (2020 * *useMapqBl); /*Max score*/
    uint16_t lowScoreUS = topScoreUSht; /*Lowest scoring kept read*/
    uint16_t scoreUS = 0;    /*Score of a single read*/

    uint64_t hashSizeULng = 0;/*Size of hash table for read extraction*/
    unsigned long majicNumULng = 0;     /*Number to use with hashing*/
    uint8_t digPerKeyUChar = 0;     /*(int) log2(hashSizeULng)*/

    struct readInfo **hashTbl = 0; /*Hash table for extraction*/
    struct readInfo *tmpRead = 0;
    struct readInfo *readOn = 0;
    struct readInfo *swapRead = 0;
    struct readInfo *scoresAry[topScoreUSht];
        /*look up table for scores*/

    struct readNodeStack searchStack[lenStackUC];  /*Used for fqGetIds*/

    FILE *testFILE = 0;
    FILE *bestReadsFILE = 0;  /*Holds the read to polish with*/
    FILE *stdinFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Check if can open files and copy reference
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(binTree == 0)
        return 2;    /*The structure has nothing*/

    /*Open reference file for read scoring*/
    if(!(noRefBl & 1))
    { /*If doing a mapping*/
        testFILE = fopen(binTree->bestReadCStr,"r");
        if(testFILE == 0) return 4;

        if(refStruct != 0)
        { /*If using the reference for deletion scoring*/
            blankSamEntry(refStruct);

            /*Read in the reference sequence*/
            errUC = readRefFqSeq(testFILE, refStruct, 0);

            if(!(errUC & 1))
            { /*If failed to read in the reference*/
                fclose(testFILE); /*No longer need open*/
                return 4;
            } /*If failed to read in the reference*/
        } /*If using the reference for deletion scoring*/

        fclose(testFILE); /*No longer need open*/
        testFILE = 0;
    } /*If doing a mapping*/

    /*See if the fastq file can be opened*/
    testFILE = fopen(binTree->fqPathCStr, "r");

    if(testFILE == 0)
        return 8;     /*Can not open the fastq file*/

    fclose(testFILE); /*Do not need open right know*/
    testFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-3: Map reads to best read & select top x              v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *numReadsKeptULng = 0; /*Make sure start at 0 reads*/

    if(!(noRefBl & 1))
    { /*If using a reference*/
        /*Build the command to run minimap2*/
        tmpCStr = cStrCpInvsDelm(minimapCmdCStr, minimap2CMD);
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
        tmpCStr =
           cpParmAndArg(
               tmpCStr,
               binTree->bestReadCStr,
               binTree->fqPathCStr
        ); /*Finsh off the minimap2 command*/
        
        /*Run minimap2 results*/
        stdinFILE = popen(minimapCmdCStr, "r");

        blankSamEntry(samST); /*Make sure start with blank*/

        /*Read in a single sam file line to check if valid (header)*/
        errUC = readSamLine(samST, stdinFILE);

        if(!(errUC & 1))
        { /*If an error occured*/
            pclose(stdinFILE);
            return 16;
        } /*If an error occured*/
    } /*If using a reference*/

    else
    { /*Else just using the fastq file*/
        stdinFILE = fopen(binTree->fqPathCStr, "r");
        errUC = readRefFqSeq(stdinFILE, samST, 0);
    } /*Else just using the fastq file*/

    /*Initalize my stack*/
    searchStack[0].readNode = 0;
    searchStack[lenStackUC - 1].readNode = 0;

    for(unsigned int IUElm = 0; IUElm < topScoreUSht; ++IUElm)
        scoresAry[IUElm] = 0; /*Intalize the scoring array*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-4: Trim, score, & select best mapped reads            v
    ^    fun-3 sec-4 sub-1: Trim read & remove if no sequence         v
    ^    fun-3 sec-4 sub-2: Score read                                v
    ^    fun-3 sec-4 sub-3: Check if Score meets requirements         v
    ^    fun-3 sec-4 sub-4: Check if at max number of reads to keep   v
    ^    fun-3 sec-4 sub-5: Check if read beats lowest kept read scorev
    ^    fun-3 sec-4 sub-6: Check if decided to keep the read         v
    ^    fun-3 sec-4 sub-7: Add kept read to the score list           v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    | Fun-3 Sec-4 Sub-1: Trim read & remove read if no sequence        |
    \******************************************************************/

    while(errUC & 1)
    { /*While their is a samfile entry to read in*/
        if(!(noRefBl & 1))
        { /*If using minimap2 input*/
            if(*samST->samEntryCStr == '@')
            { /*If was a header*/
                blankSamEntry(samST); /*Make sure start with blank*/
                errUC = readSamLine(samST,stdinFILE);
                continue; /*Is a header line, move to next line in file*/
            } /*If was a header*/

            if(samST->flagUSht & (2048 | 256 | 4))
            { /*If was a suplemental, secondary, or unmapped alignment*/
                blankSamEntry(samST); /*Make sure start with blank*/
                errUC = readSamLine(samST,stdinFILE);
                continue; /*Is a header line, move to next line*/
            } /*If was a suplemental, secondary, or unmapped alignment*/

            /*Convert & print out sam file entry*/
            errUC = trimSamEntry(samST);

            if(errUC >> 2)
            { /*If entry did not have a sequence, discard*/
                blankSamEntry(samST); /*Make sure start with blank*/
                errUC = readSamLine(samST, stdinFILE);
                continue;
            } /*If entry did not have a sequence, discard*/
        } /*If using minimap2 input*/

        /**************************************************************\
        * Fun-3 Sec-4 Sub-2: Score read
        \**************************************************************/

        findQScores(samST); /*Find the Q-scores*/

        if(!(noRefBl & 1))
        { /*If readning from minimap2*/
            scoreAln(
                minStats,
                samST,
                refStruct,   /*Reference struct to score deletions with*/ 
                &oneUChar,   /*Mapped read has Q-score*/
                &oneUChar    /*Make sure set to one if using reference*/
            ); /*Score the alignment*/

            /**************************************************************\
            * Fun-3 Sec-4 Sub-3: Check if Score meets requirements
            \**************************************************************/

            if(samST->mapqUChar < minStats->minMapqUInt ||
               !(checkIfKeepRead(minStats, samST) & 1) ||
               samST->readAligLenUInt < minStats->minReadLenULng
            ) { /*If the read Is to different from the reference*/
                /*Move to the next entry*/
                blankSamEntry(samST); /*Make sure start with blank*/
                errUC = readSamLine(samST,stdinFILE);
                continue;
            } /*If the read Is to different from the reference*/
        } /*If readning from minimap2*/

        else
        { /*Else just reading a fastq file*/
            if(samST->medianQFlt < minStats->minMedianQFlt ||
               samST->meanQFlt < minStats->minMeanQFlt ||
               samST->readLenUInt < minStats->minReadLenULng ||
               (
                   samST->readLenUInt > minStats->maxReadLenULng &&
                   minStats->maxReadLenULng != 0
               )
            ) { /*If discarding the read*/
                blankSamEntry(samST); /*Make sure start with blank*/
                errUC = readRefFqSeq(stdinFILE, samST, 0);
                continue;
            } /*If discarding the read*/
        } /*Else just reading a fastq file*/

        /**************************************************************\
        * Fun-3 Sec-4 Sub-4: Check if at max number of reads to keep
        \**************************************************************/

        tmpRead = 0;

        /*Discard decimal part of the Q-score*/
        scoreUS =
           (((uint16_t) samST->mapqUChar << 2)
              & (uint16_t) *useMapqBl
           ) + (uint16_t) samST->medianQFlt;
        /*Logic:
            mapping quality goes to 0 (ignored) if useMapqBl is 0.
            mapqUChar is 8 bits (1111 1111), mapqUChar << 2 max is 2027.
            The median Q-score is 7 bits (11 1111), with a max of 127).
            mapq << 2:  11 1111 1100
            median Q:       111 1111
                        || |||| ||||
                              16 4 1
                             32 8 2
                            64
            A median Q-score >= 4 will outweigh mapqs < 2
            A median Q-score >= 8 will outwight a mapq < 4
            A median Q-score >= 16 will outwight a mapq < 8
            A median Q-score >= 32 will outwight a mapq < 16
            A median Q-score >= 64 will outwight a mapq < 32
            A mapq >= 32 will outweigh all median Q-scores*/

        if(*numReadsKeptULng < *numReadConsULng)
        { /*If still accepting new reads*/
            readOn = makeBlankReadInfoStruct();
            readOn->idBigNum =
               buffToBigNum(samST->queryCStr, &tmpCStr, &lenBigNum);
               /*read id, pointer to id end, Number ints in big number*/

            /*Check if new lowest score*/
            if(scoreUS < lowScoreUS) lowScoreUS = scoreUS;

            if(scoresAry[scoreUS] == 0) scoresAry[scoreUS] = readOn;

            else
            { /*Else I already have reads with the same score*/
                tmpRead = scoresAry[scoreUS];
                scoresAry[scoreUS] = readOn;
                readOn->rightChild = tmpRead;
            } /*Else I already have reads with the same score*/

            ++(*numReadsKeptULng);
        } /*If still accepting new reads*/

        /**************************************************************\
        | Fun-3 Sec-4 Sub-5: Check if read beets lowest kept read score|
        \**************************************************************/

        else if(scoreUS > lowScoreUS)
        { /*Else if only keeping better reads*/

            tmpRead = scoresAry[lowScoreUS];

            /*Check if still have other reads at the low score*/
            if(tmpRead->rightChild != 0)
                scoresAry[lowScoreUS] = tmpRead->rightChild;
                     /*Mark remove the first low scoring read*/
            else
            { /*else need to find the next lowest score*/
                scoresAry[lowScoreUS] = 0;
                ++lowScoreUS; /*Move to the next score*/

                while(scoresAry[lowScoreUS] == 0)
                    ++lowScoreUS;  /*Move to read with higher score*/
            } /*else need to find the next lowest score*/

            /*Check if need to build the stack*/
            if(scoresAry[scoreUS] != 0)
                tmpRead->rightChild = scoresAry[scoreUS];
            else
                tmpRead->rightChild = 0;

            strToBackwardsBigNum(
                tmpRead->idBigNum,
                samST->queryCStr,
                &lenIdUInt
            ); /*Convert query id to a big number*/

            scoresAry[scoreUS] = tmpRead;
        } /*Else if only keeping better reads*/

        /*Move to the next read*/
        if(!(noRefBl & 1)) errUC = readSamLine(samST,stdinFILE);
        else errUC = readRefFqSeq(stdinFILE, samST, 0);
    } /*While their is a samfile entry to read in*/

    /*check which file close method I need to use*/
    if(!(errUC & 1)) pclose(stdinFILE);
    else             fclose(stdinFILE);

    stdinFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-5: Set up the best x read file name                   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Free the big numbers stored in the readInfo structs*/
    if(*numReadsKeptULng == 0)
        return 1;

    if(makeNameBl & 1)
    { /*If I need to make a name for the output file*/
        strcpy(binTree->topReadsCStr, binTree->fqPathCStr);
        tmpCStr = binTree->topReadsCStr;

        while(*tmpCStr != '\0')
            ++tmpCStr;

        tmpCStr -= 6; /*Get to . in .fastq*/
        strcpy(tmpCStr, "--top-reads.fastq"); /*Copy top read ending*/
    } /*If I need to make a name for the output file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-6: Set up hash table to extract reads                 v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readOn = scoresAry[lowScoreUS];

    for(uint32_t iRead = lowScoreUS + 1; iRead < topScoreUSht; ++iRead)
    { /*For all kept reads, set built the list for the hash table*/
       tmpRead = scoresAry[iRead];

       if(tmpRead != 0)
       { /*If I do not have a blank read*/
           swapRead = tmpRead;

           /*Move to the end of the list*/
           while(tmpRead->rightChild != 0)
               tmpRead = tmpRead->rightChild;

           tmpRead->rightChild = readOn;
           readOn = swapRead;
       } /*If I do not have a blank read*/
    } /*For all kept reads, set built the list for the hash table*/

    hashTbl = 
        readListToHash(
            readOn,
            numReadsKeptULng,
            searchStack,        /*Used for searching the hash table*/
            &hashSizeULng,      /*Will hold Size of hash table*/
            &digPerKeyUChar,    /*Power of two hash size is at*/
            &majicNumULng       /*Holds majick number for kunths hash*/
    ); /*Build the hash table*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-7: Extract reads with fastq greps hash extract        v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open fastq file and file to store the best reads in*/
    testFILE = fopen(binTree->fqPathCStr, "r");

    if(binTree->topReadsCStr[0] != '\0')
        bestReadsFILE = fopen(binTree->topReadsCStr, "w");
    else
        bestReadsFILE = stdout;

    tmpRead = 0; /*So extract reads knows not doing AVL tree search*/
 
    extractReads(
        testFILE,       /*File with reads to extract*/
        bestReadsFILE,  /*File to write reads to*/
        buffCStr,
        lenBuffUInt,
        majicNumULng,   /*Holds majick number for kunths hash*/
        digPerKeyUChar, /*Power of two hash size is at*/
        &zeroUChar,      /*Print the matches*/
        tmpRead,
        hashTbl
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-8: Clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    freeHashTbl(&hashTbl, &hashSizeULng, searchStack);

    fclose(testFILE);
    fclose(bestReadsFILE);

    return 1;
} /*findBestXReads*/

/*----------------------------------------------------------------------
| Output:
|    Files Created:
|      o A fastq file with only the best read
|    Files Modified:
|      o Removes the best read from the orignal fastq file
|    Returns:
|      o 1: if found and printed the id
|      o 4: If the fqFILE does not exist
|      o 8: If the outFILE does not exist
|      o 64: For memory allocation error
|        - In this case it does not change the old fastq file
----------------------------------------------------------------------*/
uint8_t fqGetBestReadByMedQ(
    struct readBin *clustOn,    /*Has fastq file to extract from*/
    struct samEntry *samST, /*Sam struct to use for extraction*/
    struct samEntry *bestRead   /*Holds best read till end*/
) /*Extracts read with best medain Q-score from file.
    It also considers length if integer Q-scores are the same.*/
{ /*fqGetBestReadByMedQ*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: fqOneIdExtract
    '   fun-4 sec-1: variable declerations
    '   fun-4 sec-2: Check if valid files were input
    '   fun-4 sec-3: Find best read & print out ignored reads
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t errUC = 0;
    uint8_t flagUC = 0;
    char *tmpFqCStr = "20232102-tmp10435y205y2034914701712097134.fastq";
    char *tmpCStr = 0;
    struct samEntry *swapStruct = 0;

    FILE *fqFILE = 0;
    FILE *tmpFqFILE = 0;
    FILE *bestReadFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-2: Check if valid file and set up for read in
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fqFILE = fopen(clustOn->fqPathCStr, "r");

    if(fqFILE == 0)
        return 4;

    tmpFqFILE = fopen(tmpFqCStr, "w");

    if(tmpFqFILE == 0)
    { /*If I could not open the temporary file*/
        fclose(fqFILE);
        return 4;
    } /*If I could not open the temporary file*/

    if(clustOn->bestReadCStr[0] == '\0')
    { /*If need to make a best read file*/
        tmpCStr =
            cStrCpInvsDelm(clustOn->bestReadCStr, clustOn->fqPathCStr);
        tmpCStr -= 6; /*get to "." in ".fastq" ending*/
        cStrCpInvsDelm(tmpCStr, "--best-read.fastq");
    } /*If need to make a best read file*/

    bestReadFILE = fopen(clustOn->bestReadCStr, "w");

    if(bestReadFILE == 0)
    { /*If I could not open the best read file for writing*/
        fclose(fqFILE);
        fclose(tmpFqFILE);
        return 4;
    } /*If I could not open the best read file for writing*/

    fclose(bestReadFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-3: Find best read & print out ignored reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errUC = 1; /*So I fire the loop*/
    blankSamEntry(samST);
    blankSamEntry(bestRead);

    while(errUC & 1)
    { /*While have a line to read in*/
        blankSamEntry(samST);
        errUC = readRefFqSeq(fqFILE, samST, 1);
            /*1 means to ignore new lines. In this case I am not worried
              about new lines, since all reads will have the samen number
              of new lines in their sequences. So, this will not affect
              the scoring to much, but also increase the speed of the
              step*/

        if(!(errUC & 1) || errUC == 0)
            continue; /*No more reads in file or other error*/

        findQScores(samST);

        flagUC =
            (int16_t) bestRead->medianQFlt <
            (int16_t) samST->medianQFlt;
            /*Set to 1 if samST->medianQFlt is better*/

        flagUC |=
           (
               (bestRead->readLenUInt < samST->readLenUInt) & 
               ((int16_t) bestRead->medianQFlt ==
               (int16_t) samST->medianQFlt)
           ); /*Set to 1 if flagUC is 1 or if new read length is better
                when median Q-scores are equal*/

        switch(flagUC)
        { /*switch, check if have a better read*/
            case 1: /*Have a better read & need to swap structures*/
                swapStruct = samST;
                samST = bestRead;
                bestRead = swapStruct;
        } /*switch, check if have a better read*/
            
        if(samST->seqCStr != 0) /*If have have something to print*/
            samToFq(samST, tmpFqFILE);
    } /*While have a line to read in*/ 

    fclose(fqFILE);
    fclose(tmpFqFILE);

    if(errUC == 64)
        return 64;

    bestReadFILE = fopen(clustOn->bestReadCStr, "w");
    samToFq(bestRead, bestReadFILE);
    fclose(bestReadFILE);

    /*Remove old fastq file and replace with fastq without best read*/
    remove(clustOn->fqPathCStr);
    rename(tmpFqCStr, clustOn->fqPathCStr);
    --clustOn->numReadsULng; /*Account for the removed read*/

    return 1; /*Sucess*/
} /*fqGetBestReadByMedQ*/
