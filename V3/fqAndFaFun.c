/*######################################################################
# Use:
#   o Holds functions for manipulating fastq and fasta files.
######################################################################*/

#include "fqAndFaFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' findCoInftMisc SOF:
'   fun-1 fqToFa:
'     o Converts a fastq file to a fasta file
'   fun-2 readRefFqSeq:
'     o Read in a single reference sequence from a fastq file
'   fun-3 addLineToBuff:
'     o Add characters from file to buffer, if needed resize.
'     o This will only read in till the end of the line
'   fun-4 getNumReadsInFq:
'     o Finds the number of reads in a fastq file
'   fun-5 copyFile:
'     o Makes a copy of a file (I am using this to copying fastq files)
'   fun-6 filterReads:
'     o Filters reads in a fastq file by length and mean/median Q-score
'   fun-7 moveToNextFastqEntry:
'     o Move to next entry in buffer holding data from a fastq file
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   o Creates:
|      - Fasta file with reads from the fastq file
|   o Returns:
|      - 1: Success
|      - 2: No fastq file (or invalid fastq file)
|      - 4: No fasta file
|      - 64: memory allocation error
\---------------------------------------------------------------------*/
unsigned char fqToFa(
    char *fqToCnvtCStr,         /*File name of fastq to convert*/
    char *outFaCStr,            /*File name of new fasta file*/
    struct samEntry *samStruct  /*Holds fastq file entries*/
) /*Converts a fastq file to a fasta file*/
{ /*fqToFa*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Sec-1 Sub-1: fqToFq
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *tmpCStr = 0;
    unsigned long lenFqUL = 0;
    FILE *fqFILE = 0;
    FILE *faFILE = 0;

    fqFILE = fopen(fqToCnvtCStr, "r");

    if(fqFILE == 0)
        return 2;

    fseek(fqFILE, 0, SEEK_END); /*Go to end of file*/
    lenFqUL = ftell(fqFILE);  /*Get offset (gives file size)*/
    fseek(fqFILE, 0, SEEK_SET); /*Go back to start of file*/

    faFILE = fopen(outFaCStr, "w");

    if(faFILE == 0)
    { /*If could not open the fasta file*/
        fclose(fqFILE);
        return 4;
    } /*If could not open the fasta file*/
    
    while(ftell(fqFILE) < lenFqUL)
    { /*While have more file to read in*/
        if(!(readRefFqSeq(fqFILE, samStruct, 0) & 1))
            return 64; /*Report the memory error (may not be)*/

        *samStruct->samEntryCStr = '>'; /*replace @ with >*/
        tmpCStr = samStruct->qCStr;

        while(*tmpCStr != '+')
            --tmpCStr;

        *tmpCStr = '\0'; /*Make the sequence into a c-string*/
        fprintf(faFILE, "%s", samStruct->samEntryCStr);
    } /*While have more file to read in*/

    fclose(fqFILE);
    fclose(faFILE);
    return 1;
} /*fqToFa*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: refStruct to hold the read in fastq entry & sets its
|              pointers
|    Returns:
|        - 0: if EOF
|        - 1: if succeded
|        - 2: If file was not a fastq file
|        - 130: If file had an invalide entry
|            - This error uses to flags, first it uses 2 to specify that
|              it is not a fastq file (invalid file).
|            - 2nd it uses 128 to specifty that it is not an blank file
|        - 64: If malloc failed to find memory
\---------------------------------------------------------------------*/
uint8_t readRefFqSeq(
    FILE *refFILE,      /*Pointer to fastq file to grab reference from*/
    struct samEntry *refStruct,/*Sam entry struct to hold reference*/
    char quickRunBl
      /*1: do not remove new lines, 0: remove any extra new lines*/
) /*Gets the frist reads sequence & q-score line from a fastq file*/
{ /*readRefFqSeq*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: readRefFqSeq
    '    fun-2 sec-1: Variable declarations
    '    fun-2 sec-2: Check if need to allocate memory for buffer
    '    fun-2 sec-3: Read in the first data
    '    fun-2 sec-4: If not at file end, see if have the full entry
    '    fun-2 sec-5: Read till end of file, check if valid fastq entry
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-2 Sec-1: Variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    unsigned char numLinesUC = 0; /*How many lines in sequence entry*/
    unsigned char errUC = 0;
    char *fqIterCStr = 0;         /*Marks spot working on in fastq*/
    char *oldIterCStr = 0;        /*For reallocs*/
    uint16_t extraBuffUS = 1024;
    unsigned long bytesInBuffUL = 0;
    unsigned long tmpBuffUL = 0;
    unsigned long headBuffUL = 0; /*To quickly find start of sequence*/
       /*How much to increase buffer for each new read*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-2 Sec-2: Read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(refFILE == 0)
        return 2;    /*No file provided*/

    errUC =
        addLineToBuff(
            &refStruct->samEntryCStr, /*Buffer to hold header*/
            &refStruct->lenBuffULng,  /*Size of samEntryCStr*/
            &bytesInBuffUL,           /*Number of bytes buffer holds*/
            extraBuffUS,          /*Amount to resize buffer by if full*/
            refFILE                   /*Fastq file to get header from*/
    ); /*Get the header*/

    headBuffUL = bytesInBuffUL; /*Get past header*/

    if(!(errUC & 1))
        return errUC;   /*have EOF (0) or memory allocation error (64)*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-2 Sec-3: Read in the sequence & spacer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*need to set this up so the loop does not error out*/
    oldIterCStr = refStruct->samEntryCStr + bytesInBuffUL;

    while(*oldIterCStr != '+')
    { /*While I have not reached the spacer entry*/
        tmpBuffUL = bytesInBuffUL; /*So can get to this position later*/

        errUC =
            addLineToBuff(
                &refStruct->samEntryCStr, /*Buffer to hold header*/
                &refStruct->lenBuffULng,  /*Size of samEntryCStr*/
                &bytesInBuffUL,         /*Number of bytes buffer holds*/
                extraBuffUS,      /*Amount to resize buffer by if full*/
                refFILE              /*Fastq file to get sequence from*/
        ); /*Get the header*/

        if(errUC & 64) return errUC;   /*memory allocation error (64)*/
        if(errUC == 0) return 2 + 128; /*Invalid fastq entry*/

        /*Get on first character in the new buffer*/
        oldIterCStr = refStruct->samEntryCStr + tmpBuffUL;
        ++numLinesUC; /*Count number of new lines in sequence entry*/
    } /*While I have not reached the spacer entry*/

    --numLinesUC; /*Account for the overcounting*/

    /*Check if user wants new line smoothing*/
    if(quickRunBl & 1)
        refStruct->readLenUInt = tmpBuffUL - headBuffUL - 1;
        /*-1 to account for +, head buff to addjust for header*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-2 Sec-4: Read in the Q-score entry
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    while(numLinesUC > 0)
    { /*While I have not reached the spacer entry*/
        errUC =
            addLineToBuff(
                &refStruct->samEntryCStr, /*Buffer to hold header*/
                &refStruct->lenBuffULng,  /*Size of samEntryCStr*/
                &bytesInBuffUL,         /*Number of bytes buffer holds*/
                extraBuffUS,      /*Amount to resize buffer by if full*/
                refFILE              /*Fastq file to get sequence from*/
        ); /*Get the header*/

        if(errUC & 64) return errUC;   /*memory allocation error (64)*/
        if(errUC == 0) return 2 + 128; /*Invalid fastq entry*/

        --numLinesUC; /*Count number of new lines in sequence entry*/
    } /*While I have not reached the spacer entry*/

    /*Check if user wanted to remove extra new lines*/
    if(quickRunBl & 1)
    { /*If not removing new lines, set sequence and q-score pointers*/
        refStruct->refCStr = refStruct->samEntryCStr + 1;
        refStruct->queryCStr = refStruct->samEntryCStr + 1;
        refStruct->seqCStr = refStruct->samEntryCStr + headBuffUL;
        refStruct->qCStr = refStruct->seqCStr + refStruct->readLenUInt;

        /*I need to get off the spacer entry after the Q-score entry*/
        ++refStruct->qCStr; /*get off the new line*/
        while(*refStruct->qCStr != '\n')
            ++refStruct->qCStr; /*get off the spacer entry*/
        ++refStruct->qCStr; /*get off the new line after spacer*/
    
        return errUC; /*Is 1 for another entry or 0 for EOF*/
    } /*If not removing new lines, set sequence and q-score pointers*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun-2 Sec-5: Find the read length and remove extra newlines
    ^   o fun-2 sec-5 sub-1: Set header pointers and get off header
    ^   o fun-2 sec-5 sub-2: Get seq length, remove '\n''s, set pointer
    ^   o fun-2 sec-5 sub-3: Set Q-score pointer &remove extra new lines
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /******************************************************************\
    * Fun-2 Sec-5 Sub-1: Set header pointers and get off header
    \******************************************************************/

    fqIterCStr = refStruct->samEntryCStr;
    refStruct->refCStr = fqIterCStr + 1;
    refStruct->queryCStr = fqIterCStr + 1;

    /*Go to start of the sequence entry*/
    fqIterCStr = refStruct->samEntryCStr + headBuffUL;
    oldIterCStr = fqIterCStr;

    /******************************************************************\
    * Fun-2 Sec-5 Sub-2: Get sequence length, remove '\n''s, set pointer
    \******************************************************************/

    refStruct->seqCStr = fqIterCStr;

    while(*fqIterCStr != '+')
    { /*While still on the sequence entry*/
        if(*fqIterCStr > 32)
        { /*If not on a new line entry*/
            *oldIterCStr = *fqIterCStr;
            ++refStruct->readLenUInt;
            ++oldIterCStr;
        } /*If not on a new line entry*/

        ++fqIterCStr;
    } /*While still on the sequence entry*/

    /*Copy over space, but make sure no funny stuff from user in it*/
    *oldIterCStr = '\n';
    ++oldIterCStr;
    *oldIterCStr = '+';
    ++oldIterCStr;
    *oldIterCStr = '\n';
    ++oldIterCStr;

    ++fqIterCStr; /*Get off the new line*/

    /******************************************************************\
    * Fun-2 Sec-5 Sub-3: Set Q-score pointer & remove extra new lines
    \******************************************************************/

    refStruct->qCStr = oldIterCStr;

    /*Get off the spacer entry*/
    while(*fqIterCStr != '\n')
        ++fqIterCStr;

    ++fqIterCStr; /*Get off the new line*/

    /*Gov over the Q-score entry and remove new lines*/
    while(*fqIterCStr != '\0')
    { /*While still on the sequence entry*/
        if(*fqIterCStr > 32)
        { /*If not on a new line entry*/
            *oldIterCStr = *fqIterCStr;
            ++oldIterCStr;
        } /*If not on a new line entry*/

        ++fqIterCStr;
    } /*While still on the sequence entry*/

    return errUC; /*Is 1 for another entry or 0 for EOF*/
} /*readRefFqSeq*/

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
    unsigned long *lenBuffUL, /*Size of the buffer*/
    unsigned long *curBuffUL, /*Length buffer with valid data*/
    unsigned long resBuffUL,  /*Amount to resize buffer by if full*/
    FILE *inFILE              /*File to grab data from*/
) /*Add characters from file to buffer, if needed resize.
    This will only read in till the end of the line*/
{ /*addToBuff*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: addToBuff
    '   o fun-3 sec-1: variable declerations
    '   o fun-3 sec-2: Check if need to resize the buffer
    '   o fun-3 sec-3: Read in the next line in the buffer
    '   o fun-3 sec-4: If at end of file, update read in lengths
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    unsigned long spareBuffUL = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-2: Check if need to resize the buffer
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
    ^ Fun-3 Sec-3: Read in the next line in the buffer
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
    ^ Fun-3 Sec-4: If at end of file, update read in lengths
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
} /*addToBuff*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 0 if the fastq file could not be opened
|     o The number of reads in the fastq file
| Note:
|   - fqFILE will be set back to its starting position at the end
\---------------------------------------------------------------------*/
unsigned long getNumReadsInFq(
    char *fqFileCStr /*Path to fastq file to get number of reads in*/
) /*Find the number of reads in a fastq file*/
{ /*getnUmReadsInFq*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: Sec-1 Sub-1: getNumReadsInFq
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *tmpCS = 0;
    char buffCS[1024];
    uint64_t tmpUL = 1024;
    unsigned long numReadsUL = 0;
    FILE *fqFILE = fopen(fqFileCStr, "r");

    if(fqFileCStr == 0)
        return 0; /*No fastq file to work with*/

    tmpCS = buffCS;
    buffCS[0] = '\0';

     while(moveToNextFastqEntry(buffCS, &tmpCS,1024, &tmpUL,fqFILE) & 2)
         ++numReadsUL;

     return numReadsUL;
} /*getnUmReadsInFq*/

/*---------------------------------------------------------------------\
| Output:
|   o Writes a copy of the original file (orgFqCStr):
|   o Returns:
|     - 1 if no problems happened
|     - 2 if could not open orgFqCStr
|     - 4 if could not open newFqCStr
\---------------------------------------------------------------------*/
unsigned char copyFile(
    char *orgFqCStr,  /*Path to fastq file to copy (orignal)*/
    char *newFqCStr   /*Path to the new duplicate Fastq file*/
) /*Copies a fastq file to a new file*/
{ /*copyFile*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-5 TOC: Sec-1 Sub-1: copyFqFile
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    unsigned short lenBuffUS = 1024;
    char buffCStr[lenBuffUS];
    unsigned long bytesReadInUL = 0;
    FILE *fqFILE = 0;
    FILE *cpFqFILE = 0;

    /*Make sure can open and write both files*/
    fqFILE = fopen(orgFqCStr, "r");

    if(fqFILE == 0)
        return 2;

    cpFqFILE = fopen(newFqCStr, "w");

    if(cpFqFILE == 0)
    { /*If could not open the copy file*/
        fclose(fqFILE);
        return 4;
    } /*If could not open the copy file*/

    /*Copy the input file*/
    bytesReadInUL = fread(buffCStr, sizeof(char), lenBuffUS, fqFILE);

    while(bytesReadInUL == lenBuffUS)
    { /*While I have not reached the end of the file*/
        fwrite(buffCStr, sizeof(char), bytesReadInUL, cpFqFILE);
        bytesReadInUL = fread(buffCStr,sizeof(char), lenBuffUS, fqFILE);
    } /*While I have not reached the end of the file*/

    if(bytesReadInUL > 0)
        fwrite(buffCStr, sizeof(char), bytesReadInUL, cpFqFILE);

    fclose(fqFILE);
    fclose(cpFqFILE);
    return 1;
} /*copyFile*/

/*---------------------------------------------------------------------\
| Output:
|    Creates File: from outFqPathCStr with filtered reads
|    Returns:
|        - 1: if succeded
|        - 2: If file was not a fastq file
|        - 130: If file had an invalide entry
|            - This error uses to flags, first it uses 2 to specify that
|              it is not a fastq file (invalid file).
|            - 2nd it uses 128 to specifty that it is not an blank file
|        - 64: If malloc failed to find memory
\---------------------------------------------------------------------*/
unsigned char filterReads(
    char *fqCStr,/*Fastq file with reads to filter (null for stdin)*/
    char *outCStr,/*Name of fastq to write reads to (null for stdout)*/
    struct samEntry *samST, /*For reading in lines for the fastq file*/
    struct minAlnStats *minStats
        /*Has min/max lengths & mean/median Q-score*/
) /*Filters reads in a fastq file by length and mean/median Q-score*/
{ /*filterReads*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-6 TOC: Sec-1 Sub-1: filterReads
    '   o fun-6 sec-1: variable declerations
    '   o fun-6 sec-2: Check if input files are valid
    '   o fun-6 sec-3: Filter reads
    '   o fun-6 sec-4: Clean up
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned char errUC = 0;
    FILE *fqFILE = 0;
    FILE *outFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-2: Check if input files are valid
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(fqCStr == 0)
        fqFILE = stdin;
    else
        fqFILE = fopen(fqCStr, "r");

    if(fqFILE == 0)
        return 2;

    if(*outCStr == 0)
        outFILE = stdout;
    else
        outFILE = fopen(outCStr, "w");

    if(outFILE == 0)
        return 4;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-3: Filter reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Read in the first line*/
    blankSamEntry(samST);
    errUC = readRefFqSeq(fqFILE, samST, 0);
       /*Want to remove extra new lines for Q-score calculations*/

    while(errUC & 1)
    { /*While have entries to check*/

        findQScores(samST);

        /*Check if keeping the read*/
        if(samST->medianQFlt < minStats->minMedianQFlt)
            goto readNextLine;

        if(samST->meanQFlt < minStats->minMeanQFlt)
            goto readNextLine;

        if(samST->readLenUInt < minStats->minReadLenULng)
            goto readNextLine;

        if(minStats->maxReadLenULng == 0)
            samToFq(samST, outFILE); /*Save the read*/

        else if(samST->readLenUInt <= minStats->maxReadLenULng)
            samToFq(samST, outFILE); /*Save the read*/

        readNextLine:
            blankSamEntry(samST);
            errUC = readRefFqSeq(fqFILE, samST, 0);
    } /*While have entries to check*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-4: Clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fclose(fqFILE);
    fclose(outFILE);

    if(errUC != 0)
       return errUC; /*If errored out*/

    return 1;
} /*filterReads*/

/*##############################################################################
# Output:
#    Modifies: bufferCStr to have the next buffer if empty
#    Modifies: incurments pointInBufferCStr to start of next read
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
uint8_t moveToNextFastqEntry(
    char *bufferCStr,  /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr, /*position working on in buffer*/
    uint32_t buffSizeInt,        /*Size of buffer to work on*/
    uint64_t *lenInputULng,      /*Length of input from fread*/
    FILE *fastqFile              /*Fastq file to get data from*/
) /*Moves to next fastq read, without printing out*/
{ /*moveToNextFastqEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 TOC:
    #     fun-7 sec-1: Variable declerations
    #     fun-7 sec-2: Move past header
    #     fun-7 sec-3: Find number of newlines in sequence & move to spacer
    #     fun-7 sec-4: Move past spacer
    #     fun-7 Sec-5: Move past Q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t
        numSeqlinesULng = 0;  /*Record the number of new lines in the sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-2: Move past header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on header, move to sequence line*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on header, move to sequence line*/

    ++(*pointInBufferCStr); /*get off the new line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-3: Find number new lines in seqence line & move to spacer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(**pointInBufferCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**pointInBufferCStr == '\n')
            numSeqlinesULng++;    /*Record number new lines, for q-score entry*/

        (*pointInBufferCStr)++;
    } /*While on sequence line, count number new lines & move to spacer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-4: Move past spacer entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on the spacer entry*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on the spacer entry*/

    ++(*pointInBufferCStr); /*get off the new line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-5: Move past q-score entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputULng < buffSizeInt)
            { /*If at end of the file*/
                if(numSeqlinesULng > 1)
                    return 0;         /*Missing q-score lines for entry*/
                else
                    return 4;         /*At end of file & printed last line*/
            } /*If at end of the file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**pointInBufferCStr == '\n')
            numSeqlinesULng--;    /*Record number new lines, for q-score entry*/

        (*pointInBufferCStr)++;
    } /*While have q-score entry lines to print out*/

    return 2; /*Copied name sucessfully*/
} /*moveToNextFastqEntry*/
