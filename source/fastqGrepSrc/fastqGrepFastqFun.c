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
#    Modifies: idBigNum to holder head as hex big number (set to 0 for 
#    Modifies: lenBigNumAryUChar when string needs more longs, set to 0 when
#              memory allocatoin failed
#    Returns:
#        8: memory allocation failed
#        4: If the end of the file
#        2: if nothing went wrong
#        0: If ran out of file
##############################################################################*/
char parseFastqHeader(
    char *buffCStr,              /*buffer to hold fread input (can have data)*/
    char **startNameCStr,        /*Start of read id*/
    int *lenInInt,               /*Length of input from fread*/
    int lenBuffInt,             /*Size of buffer to work on*/
    struct bigNum *idBigNum,   /*Big number struct to hold read id*/
    unsigned char *lenBigNumAryUChar,  /*Length of U longs in big number*/
    FILE *fqFILE              /*Fastq file to get data from*/
) /*Reads input from file & sets pointer to start of read name*/
{ /*parseFastqHeader*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: Make sure the bigNum struct has been initalized
    #    fun-1 sec-3: Find an move past @ makring header start
    #    fun-1 sec-4: loop through buffer & check if need to get input from file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *seqIterCStr = *startNameCStr;

    unsigned long
        lenIdLng = 0,   /*Number of offset fseek by, if need to regrab*/
        *elmOnPtrULng = 0; /*For number conversion*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Make sure the bigNum struct has been initalized
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make sure have an array large enough to store number*/
    if(idBigNum->lenAllElmChar < *lenBigNumAryUChar)
    { /*If I need to make the array biger*/
        if(idBigNum->bigNumAryULng == 0)
            idBigNum->bigNumAryULng =
                malloc(sizeof(unsigned long) * (*lenBigNumAryUChar));
        else
            idBigNum->bigNumAryULng =
                realloc(
                    idBigNum->bigNumAryULng,
                    sizeof(unsigned long) * (*lenBigNumAryUChar)
            ); /*Need to reallocate memory*/

        if(idBigNum->bigNumAryULng == 0)
        { /*If memory reallocation failed*/
            idBigNum->lenUsedElmChar = 0;
            fprintf(
                stderr,
               "Memory allocation failed: Fun-1 fastqFastqRun.c line 60 to 63\n"
            ); /*Print error message to user*/

            return 8;
        } /*If memory reallocation failed*/

        idBigNum->lenAllElmChar = *lenBigNumAryUChar;
    } /*If I need to make the array biger*/

    else
        *lenBigNumAryUChar = idBigNum->lenAllElmChar; /*more longs in stuct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Find an move past @ makring header start
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*seqIterCStr != '@')
    { /*While not at header*/

        if(*seqIterCStr == '\0')
        { /*If at the end of the buffer, but not at start of read*/
            if(*lenInInt < lenBuffInt)
                return 4;                    /*Done with file*/

            *lenInInt=fread(buffCStr, sizeof(char), lenBuffInt, fqFILE);

            seqIterCStr = buffCStr;
            *startNameCStr = buffCStr;
        } /*If at the end of the buffer, but not at start of read*/

        else
            seqIterCStr++; /*Move to the header*/
    } /*While not at header*/

    seqIterCStr++; /*Move off the header symbol*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: loop through id & convert to big number
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    idBigNum->lenUsedElmChar = 0; /*Make sure 0, since copying new number*/

    do
    { /*While still on the read name part of header*/

        if(idBigNum->lenUsedElmChar >= idBigNum->lenAllElmChar)
        { /*If need to resize the unsigned long array*/
            idBigNum->bigNumAryULng =
                realloc(
                    idBigNum->bigNumAryULng,
                    sizeof(unsigned long) * (*lenBigNumAryUChar)
            ); /*Need to make the unsigned long array biger*/

            if(idBigNum->bigNumAryULng == 0)
            { /*If memory reallocation failed*/
                idBigNum->lenUsedElmChar = 0;
                fprintf(
                    stderr,
                   "Memory allocation failed: Fun-1 fastqGrepFastqFun.c 121\n"
                ); /*Print error message to user*/
    
                return 8;
            } /*If memory reallocation failed*/

            idBigNum->lenAllElmChar++;
            (*lenBigNumAryUChar)++;
        } /*If need to resize the unsigned long array*/

        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        *elmOnPtrULng = 0;

        for(
            unsigned char charBit = 0;       /*Using ULng to make easiy*/
            charBit < (sizeof(unsigned long) << 3); /*While bits to fill in*/
            charBit += 4                     /*Bits used per hex character*/
        ) { /*For empty bits in the current big number unsigned long element*/
            if(*seqIterCStr < 33)
                break; /*If have finshed converting the hex string*/

            if(*seqIterCStr > 47 && *seqIterCStr < 71) /*0-9 or A-F, (0-15)*/
                *elmOnPtrULng = *elmOnPtrULng + (((*seqIterCStr)-48) << charBit);

            else if(*seqIterCStr > 96 && *seqIterCStr < 103) /*a-f, (10-15)*/
                *elmOnPtrULng = *elmOnPtrULng + (((*seqIterCStr)-87) << charBit);
            else
                charBit -= 4;   /*make sureo only recored conversion*/

            seqIterCStr++;  /*move to next character in id*/
            lenIdLng++; /*count number of characters read in*/

            if(*seqIterCStr == '\0')
            { /*If ran out of buffer & need to read in more of the file*/
                if(*lenInInt < lenBuffInt)
                    return 0; /*At end of file, but no sequence*/

                lenIdLng++; /*Is one off the @ symbol*/

                fseek(fqFILE, lenIdLng * -1 /*start of header*/, SEEK_CUR);
                *lenInInt = fread(buffCStr, sizeof(char), lenBuffInt, fqFILE);

                *(buffCStr + *lenInInt) = '\0';/*make sure a c-string*/
                seqIterCStr = buffCStr + lenIdLng; /*move back to position at*/
                *startNameCStr = buffCStr;  /*Start of read id*/
            } /*If ran out of buffer & need to read more of the file*/
        } /*For empty bits in the current big number unsigned long element*/

        (idBigNum->lenUsedElmChar)++;
    } while(*seqIterCStr > 32); /*32 is space, < 33 '\n', '\t', '\r', & \0*/
       /*Ensured that the loop never starts at \0*/

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

    char *seqIterCStr = *readStartCStr;
    unsigned long numSeqlinesULng = 0; /*Holds number of new lines in sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Print out the header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*seqIterCStr != '@')
    { /*While not at start of fastq header*/
        if(*seqIterCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            seqIterCStr = bufferCStr;
            *readStartCStr = bufferCStr;
        } /*If ran out of buffer & need to read in more of the file*/

        else
            seqIterCStr++; /*Move off any new lines*/
    } /*While not at start of fastq header*/
 
    while(*seqIterCStr != '\n')
    { /*While on header, move to sequence line*/
        if(*seqIterCStr == '\0')
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
            seqIterCStr = bufferCStr;
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        seqIterCStr++; /*Move to next character in header*/
    } /*While on header, move to sequence line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Print out & find number new lines in seqence line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*seqIterCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(*seqIterCStr == '\0')
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
            seqIterCStr = bufferCStr;
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(*seqIterCStr == '\n')
            numSeqlinesULng++;    /*Record number new lines, for q-score entry*/

        seqIterCStr++;
    } /*While on sequence line, count number new lines & move to spacer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-4: Print out the spacer entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(*seqIterCStr != '\n')
    { /*While on the spacer entry*/

        if(*seqIterCStr == '\0')
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
            seqIterCStr = bufferCStr;
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        seqIterCStr++;
    } /*While on the spacer entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-5: Print out the q-scores entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(*seqIterCStr == '\0')
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
            seqIterCStr = bufferCStr;
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(*seqIterCStr == '\n')
            numSeqlinesULng--;    /*Record number new lines, for q-score entry*/

        seqIterCStr++;
    } /*While have q-score entry lines to print out*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-6: Print out remaning parts of entry in buffer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    seqIterCStr--;               /*Get back on the new line*/
    *seqIterCStr = '\0';         /*Turn '\n' into '\0', so print stops at line*/
    printf("%s\n", *readStartCStr); /*Print out the read name*/
    *readStartCStr = seqIterCStr + 1; /*Move to header in next entry*/
    
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
    char **readStartCStr,        /*Points to locatoin working on in buffer*/
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

    unsigned long numSeqlinesULng = 0; /*Holds number of new lines in sequence*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Move past header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**readStartCStr != '@')
    { /*While not at start of fastq header*/
        if(**readStartCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/
            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *readStartCStr = bufferCStr;
        } /*If ran out of buffer & need to read in more of the file*/

        else
            (*readStartCStr)++; /*Move off any new lines*/
    } /*While not at start of fastq header*/


    while(**readStartCStr != '\n')
    { /*While on header, move to sequence line*/

        if(**readStartCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*readStartCStr)++;
    } /*While on header, move to sequence line*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Find number new lines in seqence line & move to spacer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(**readStartCStr != '+')
    { /*While on sequence line, count number new lines & move to spacer*/

        if(**readStartCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**readStartCStr == '\n')
            numSeqlinesULng++;    /*Record number new lines, for q-score entry*/

        (*readStartCStr)++;
    } /*While on sequence line, count number new lines & move to spacer*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-4: Move past spacer entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
    while(**readStartCStr != '\n')
    { /*While on the spacer entry*/

        if(**readStartCStr == '\0')
        { /*If ran out of buffer & need to read in more of the file*/

            if(*lenInputInt < buffSizeInt)
                return 0;         /*Is not a complete fastq file*/

            *lenInputInt = fread(bufferCStr,
                                 sizeof(char),
                                 buffSizeInt,
                                 fastqFile
            ); /*Read in more of the file*/

            *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        (*readStartCStr)++;
    } /*While on the spacer entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-5: Move past q-score entry (same number lines as sequence)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(numSeqlinesULng > 0)
    { /*While have q-score entry lines to print out*/

        if(**readStartCStr == '\0')
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
            *readStartCStr = bufferCStr;
            continue;                           /*so can check if '\n'*/
        } /*If ran out of buffer & need to read more of the file*/

        if(**readStartCStr == '\n')
            numSeqlinesULng--;    /*Record number new lines, for q-score entry*/

        (*readStartCStr)++;
    } /*While have q-score entry lines to print out*/

    return 2; /*Copied name sucessfully*/
} /*moveToNextFastqEntry*/
