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

#include "fastqGrepFastqFun.h"

/*Make look up table to look if character is valid hex character
    64 is invisivle character
    32 is non-hex character (printable)
*/
char hexTblCharAry[] =
    {
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,/*0-32 (invisible)*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, /*non-#*/

     0, 1, 2, 3, 4, 5, 6, 7, 8, 9, /*48-57 Numbers*/

     32, 32, 32, 32, 32, 32, 32, /*Between numbers & uppcase letters*/

     10, 11, 12, 13, 14, 15,     /*A-F*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32,  /*G-` (71 to 96)*/

     10, 11, 12, 13, 14, 15,     /*a-f (97 to 102)*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32 /*g to ascii limit (103 to 256)*/
    };

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdULng to hold the length of the read id
#    Returns:
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char parseFastqHeader(
    char *bufferCStr,        /*buffer to hold fread input (can have data)*/
    char **startNameCStr,    /*Will hold the start of the name*/
    char **endNameCStr,      /*Points to start of read name, will point to end*/
    int *lenInputInt,        /*Length of input from fread*/
    int buffSizeInt,         /*Size of buffer to work on*/
    unsigned long *lenIdULng,/*Lengtho of the read id*/
    struct bigNum *idBigNum,  /*Will hold big number found*/
    FILE *fastqFile          /*Fastq file to get data from*/
) /*Reads input from file & sets pointer to start of read name*/
{ /*parseFastqHeader*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: loop through buffer & check if need to get input from file
    #    fun-1 sec-3: Copy header over, recored length, & incurment pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char charBit = 0;
    unsigned long *elmOnPtrULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: loop through buffer & check if need to get input from file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *lenIdULng = 0; /*Make sure starts at 0*/
    idBigNum->lenUsedElmChar = 0;

    if(**endNameCStr == '\n')
        (*endNameCStr)++;

    if(**endNameCStr == '\0')
    { /*If at the end of the buffer, but not at start of read*/
        if(*lenInputInt < buffSizeInt)
          return 4;                    /*Done with file*/

        *lenInputInt = fread(
                           bufferCStr,
                           sizeof(char),
                           buffSizeInt,
                           fastqFile
        ); /*Read in more of the file*/

      *endNameCStr = bufferCStr;
    } /*If at the end of the buffer, but not at start of read*/

    *startNameCStr = *endNameCStr;    /*know at the start of the read name*/
      
     do { /*While still on the read name part of header*/

        /*Graph unsigned long element working on*/
        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        *elmOnPtrULng = 0;
        charBit = 0;

        while(charBit < (sizeof(unsigned long) << 3))
        { /*while empty bits in the current big number unsigned long element*/
            if(hexTblCharAry[**endNameCStr] & 64)
                break; /*If have finshed converting the hex string*/

            else if(!(hexTblCharAry[**endNameCStr] & 32)) /*Was hex char*/
            { /*Else if is a hex character*/
                *elmOnPtrULng =
                    *elmOnPtrULng +
                    (hexTblCharAry[**endNameCStr] << charBit);
                charBit += 4;   /*make sureo only recored conversion*/
            } /*Else if is a hex character*/

            (*lenIdULng)++;
            (*endNameCStr)++;

            if(**endNameCStr == '\0')
            { /*If ran out of buffer & need to read in more of the file*/
                if(*lenInputInt < buffSizeInt)
                    return 0;     /*At end of file, but no sequence*/

                /*Put file pointer back to the start of read name*/
                fseek(fastqFile, *lenIdULng * -1, SEEK_CUR);

                *lenInputInt = fread(bufferCStr,
                                     sizeof(char),
                                     buffSizeInt,
                                     fastqFile
                ); /*Read in more of the file*/

                *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
                *startNameCStr = bufferCStr;
                *endNameCStr = bufferCStr + *lenIdULng;
                continue;
            } /*If ran out of buffer & need to read more of the file*/
        } /*while empty bits in the current big number unsigned long element*/

        (idBigNum->lenUsedElmChar)++; /*Track number ULngs acctualy used*/
    } while(**endNameCStr > 32); /*While still on the read name part of header*/

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
char printFastqEntry(
    char *bufferCStr,            /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr,     /*Points to locatoin working on in buffer*/
    char **readStartCStr,         /*Points to start of read name*/
    int buffSizeInt,              /*Size of buffer to work on*/
    int *lenInputInt,             /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
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

    unsigned long
        numSeqlinesULng = 0;  /*Record the number of new lines in the sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Print out the header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on header, move to sequence line*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            printf("%s", *readStartCStr); /*print out old buffer*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on header, move to sequence line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Print out & find number new lines in seqence line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(**pointInBufferCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            printf("%s", *readStartCStr); /*print out old buffer*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
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

            printf("%s", *readStartCStr); /*print out old buffer*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            *readStartCStr = *pointInBufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on the spacer entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-5: Print out the q-scores entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            printf("%s", *readStartCStr); /*print out old buffer*/

            if(*lenInputInt < buffSizeInt)
            { /*If at end of the file*/
                if(numSeqlinesULng > 1)
                    return 0;         /*Missing q-score lines for entry*/
                else
                    return 4;         /*At end of file & printed last line*/
            } /*If at end of the file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
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
    printf("%s\n", *readStartCStr); /*Print out the read name*/
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
char moveToNextFastqEntry(
    char *bufferCStr,            /*buffer to hold fread input (can have data)*/
    char **pointInBufferCStr,     /*Points to locatoin working on in buffer*/
    int buffSizeInt,              /*Size of buffer to work on*/
    int *lenInputInt,             /*Length of input from fread*/
    FILE *fastqFile               /*Fastq file to get data from*/
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

    unsigned long
        numSeqlinesULng = 0;  /*Record the number of new lines in the sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Move past header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**pointInBufferCStr != '\n')
    { /*While on header, move to sequence line*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on header, move to sequence line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Find number new lines in seqence line & move to spacer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(**pointInBufferCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
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

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*pointInBufferCStr)++;
    } /*While on the spacer entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-5: Move past q-score entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(**pointInBufferCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
            { /*If at end of the file*/
                if(numSeqlinesULng > 1)
                    return 0;         /*Missing q-score lines for entry*/
                else
                    return 4;         /*At end of file & printed last line*/
            } /*If at end of the file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *pointInBufferCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**pointInBufferCStr == '\n')
            numSeqlinesULng--;    /*Record number new lines, for q-score entry*/

        (*pointInBufferCStr)++;
    } /*While have q-score entry lines to print out*/

    return 2; /*Copied name sucessfully*/
} /*moveToNextFastqEntry*/
