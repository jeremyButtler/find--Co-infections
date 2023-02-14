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

/*##############################################################################
# Output:
#    Modifies: startNameCStr to point to the start of the read name
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
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
    FILE *fastqFile          /*Fastq file to get data from*/
) /*Reads input from file & sets pointer to start of read name*/
{ /*parseFastqHeader*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: loop through buffer & check if need to get input from file
    #    fun-1 sec-3: Copy header over, recored length, & incurment pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: loop through buffer & check if need to get input from file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long charFromStartULng = 0;

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
      
    while(**endNameCStr != ' ' &&          /*End of read name*/
          **endNameCStr != '\n'            /*End of line*/
    ) { /*While still on the read name part of header*/

        if(**endNameCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            if(*lenInputInt < buffSizeInt)
                return 0;                    /*At end of file, but no sequence*/

            /*Go to the start of the name*/
            fseek(
                fastqFile,
                charFromStartULng * -1,
                SEEK_CUR
            ); /*Put file pointer back to the start of read name*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *startNameCStr = bufferCStr;
            *endNameCStr = bufferCStr;
            continue;
        } /*If ran out of buffer & need to read more of the file*/

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Fun-1 Sec-3: Copy header over, recored length, & incurment pointers
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        charFromStartULng++;
        (*endNameCStr)++;
    } /*While still on the read name part of header*/

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
