/*##############################################################################
# Name: fastqGrepFastqFun
# Use:
#    Functions to read input from fastq, find fastq header, & print out entries
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepFastqFun TOC:
#   fun-1: parseFastqHeader: Sets pointer to start of fastq header
#   fun-2: printFastqEntry: Print out an entry in a fastq file
#   fun-3: moveToNextFastqEntry: Move to next fastq entry without printing
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fqGetIdsFqFun.h"

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdInt to hold the length of the read id
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
uint8_t parseFastqHeader(
    char *bufferCStr,     /*buffer to hold fread input (can have data)*/
    char **startNameCStr, /*Will hold the start of the name*/
    char **endNameCStr,   /*start of read name, will point to end*/
    uint64_t *lenInputULng,  /*Length of input from fread*/
    uint32_t buffSizeInt,    /*Size of buffer to work on*/
    int32_t *lenIdInt,     /*Length of the read id*/
    struct bigNum *idBigNum,  /*Will hold big number found*/
    FILE *fastqFile          /*Fastq file to get data from*/
) /*Reads input from file & sets pointer to start of read name*/
{ /*parseFastqHeader*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: Move to start of read id
    #    fun-1 sec-2: loop through buffer & check if need to get input from file
    #    fun-1 sec-3: Convert c-string read id to bigNum read id (U long array)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t charBit = 0;
    unsigned long *elmOnPtrULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Move to start of read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *lenIdInt = 0; /*Make sure starts at 0*/
    idBigNum->lenUsedElmChar = 0;

    if(**endNameCStr == '\n')
        (*endNameCStr)++;

    if(**endNameCStr == '\0')
    { /*If at the end of the buffer, but not at start of read*/
        if(*lenInputULng < buffSizeInt)
          return 4;                    /*Done with file*/

        *lenInputULng = fread(
                           bufferCStr,
                           sizeof(char),
                           buffSizeInt,
                           fastqFile
        ); /*Read in more of the file*/

      *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
      *endNameCStr = bufferCStr;
    } /*If at the end of the buffer, but not at start of read*/

    *startNameCStr = *endNameCStr;    /*know at the start of the read name*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Convert c-string read id to bigNum read id (U long array)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
      
     do { /*While still on the read name part of header*/

        /*Graph unsigned long element working on*/
        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        (idBigNum->lenUsedElmChar)++; /*Track number ULngs acctualy used*/
        *elmOnPtrULng = 0;
        charBit = 0;

        while(charBit < (sizeof(unsigned long) << 3))
        { /*while empty bits in the current big number unsigned long element*/
            if(hexTblCharAry[(unsigned char) **endNameCStr] & 64)
                return 2; /*If have finshed converting the hex string*/

            if(!(hexTblCharAry[(unsigned char) **endNameCStr] & 32))
            { /*If is a hex character*/
                (*elmOnPtrULng) +=
                    (hexTblCharAry[(unsigned char) **endNameCStr] << charBit);
                charBit += 4;
            } /*If is a hex character*/

            (*lenIdInt)++;
            (*endNameCStr)++;

            if(**endNameCStr == '\0')
            { /*If ran out of buffer & need to read in more of the file*/
                if(*lenInputULng < buffSizeInt)
                    return 0;     /*At end of file, but no sequence*/

                /*Put file pointer back to the start of read name*/
                fseek(fastqFile, *lenIdInt * -1, SEEK_CUR);

                *lenInputULng = fread(bufferCStr,
                                     sizeof(char),
                                     buffSizeInt,
                                     fastqFile
                ); /*Read in more of the file*/

                *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
                *startNameCStr = bufferCStr;
                *endNameCStr = bufferCStr + *lenIdInt;
                continue;
            } /*If ran out of buffer & need to read more of the file*/
        } /*while empty bits in the current big number unsigned long element*/
    } while((unsigned char) **endNameCStr > 32);
    /*While still on the read name part of header*/

    return 2; /*Copied name sucessfully*/
} /*parseFastqHeader*/

/*##############################################################################
# Output:
#    stdout: prints id & the remaing three lines in fastq for read
#    Modifies: bufferCStr to have the next buffer if empty
#    Modifies: incurments pointInBufferCStr to start of next read
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
uint8_t printFastqEntry(
    char *bufferCStr,  /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr, /*locatoin working on in buffer*/
    char **readStartCStr,     /*start of read name*/
    uint32_t buffSizeInt,        /*Size of buffer to work on*/
    uint64_t *lenInputULng,      /*Length of input from fread*/
    FILE *outFILE,               /*File to print reads to*/
    FILE *fastqFile              /*Fastq file to get data from*/
) /*Reads input from file & marks end of read name*/
{ /*printFastqEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC:
    #     fun-2 sec-1: Variable declerations
    #     fun-2 sec-2: Print out the header
    #     fun-2 sec-3: Print out & find number new lines in seqence line
    #     fun-2 sec-4: Print out the spacer entry
    #     fun-2 Sec-5: Print out the q-scores entry (same # lines as sequence)
    #     fun-2 sec-6: Print out remaning parts of entry in buffer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t
        numSeqlinesULng = 0;  /*Record the number of new lines in the sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Print out the header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on header, move to sequence line*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            /*print out old buffer*/
            fprintf(outFILE, "%s", *readStartCStr);

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on header, move to sequence line*/

    ++(*pointInBufferCStr); /*get off the new line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Print out & find number new lines in seqence line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(**pointInBufferCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            /*print out old buffer*/
            fprintf(outFILE, "%s", *readStartCStr);

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**pointInBufferCStr == '\n')
            numSeqlinesULng++;    /*Record number new lines, for q-score entry*/

        (*pointInBufferCStr)++;
    } /*While on sequence line, count number new lines & move to spacer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-4: Print out the spacer entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on the spacer entry*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            /*print out old buffer*/
            fprintf(outFILE, "%s", *readStartCStr);

            if(*lenInputULng < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputULng = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputULng) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on the spacer entry*/

    ++(*pointInBufferCStr); /*get off the new line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-5: Print out the q-scores entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            /*print out old buffer*/
            fprintf(outFILE, "%s", *readStartCStr);

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
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**pointInBufferCStr == '\n')
            numSeqlinesULng--;    /*Record number new lines, for q-score entry*/

        (*pointInBufferCStr)++;
    } /*While have q-score entry lines to print out*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-6: Print out remaning parts of entry in buffer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    (*pointInBufferCStr)--;      /*Get back on the new line*/
    **pointInBufferCStr = '\0'; /*Turn '\n' into '\0', so print stops at line*/

    /*Finish printing out the read entry*/
    fprintf(outFILE, "%s\n", *readStartCStr);

    **pointInBufferCStr = '\n'; /*Convert back to \n so no longer c-string*/
    
    return 2; /*Copied name sucessfully*/
} /*printFastqEntry*/

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
    # Fun-3 TOC:
    #     fun-3 sec-1: Variable declerations
    #     fun-3 sec-2: Move past header
    #     fun-3 sec-3: Find number of newlines in sequence & move to spacer
    #     fun-3 sec-4: Move past spacer
    #     fun-3 Sec-5: Move past Q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t
        numSeqlinesULng = 0;  /*Record the number of new lines in the sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Move past header
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
    # Fun-3 Sec-3: Find number new lines in seqence line & move to spacer
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
    # Fun-3 Sec-4: Move past spacer entry
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
    # Fun-3 Sec-5: Move past q-score entry (same number lines as sequence)
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
