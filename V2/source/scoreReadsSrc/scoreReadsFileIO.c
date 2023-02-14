/*##############################################################################
# TOC:
#   void openFile: Check if is a file, open file, and return max line length
#       fun-1 sec-1: variable declarations
#       fun-1 sec-2: check if input file can be opened
#       fun-1 sec-3: find the file size
#       fun-1 sec-4: Find the max line length
#   void stdoutReadScore: Print out the stats for the samfile entry
#       fun-2 Sec-1: Print the samfile stats to stdout
#    fun-3: readFirstSeqInFastq: read in fastq 1st sequence & Q-score lines
##############################################################################*/

#include "scoreReadsFileIO.h"

/*##############################################################################
# Name: openFile
# Use: Opens an input sam file and finds the longest line
# Input: 
#    samPathCStr: Path/Name of sam file to open (c-string)
#    maxLineLenULng: Length of the longest line in the same file (unsigned long)
#    buffSizeInt: Size of buffer used to scan through the file (int)
# Output:
#    Returns: pointer to the file (returns 0 if nothing if file or not valid)
#    Modifies: maxLineLenULng to hold the length of the longest line
##############################################################################*/
FILE * openFile(char *fileInCStr,
                unsigned long *maxLineLenULng,
                int buffSizeInt)
{ /*openFile*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-1 Sec-1: Variable declarations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Set up buffers to read each line*/
    char buffCStr[buffSizeInt + 1],  /*Buffer to read in part of the file*/
         *tmpCStr = buffCStr,        /*Pointer to advance through the buffer*/
         tabCntChar = 0;             /*Number of tabs in a line*/

    /*Set up file length variables*/
    unsigned long lenFileULng = 0,   /*Length of the file in characters*/
                  posInFileULng = 0, /*position in the file currently at*/
                  lenLineULng = 0;   /*Length of the line on in the file*/
    FILE *samFile = 0;               /*points to the sam file*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-1 Sec-2: Check if the input can be opened
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    samFile = fopen(fileInCStr, "r"); /*open the sam file*/

    if(samFile == 0)
        return 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-1 Sec-3: Find the file size
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    fseek(samFile, 0 , SEEK_END); /*find the end of the file*/
    lenFileULng = ftell(samFile); /*get the length of the file*/
    fseek(samFile, 0, SEEK_SET); /*go back to the start of the file*/

    if(lenFileULng < 2)
    { /*If there was nothing in the file, quite*/
        fclose(samFile); /*close the file*/
        return 0;
    } /*If there was nothing in the file, quite*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-1 Sec-4: Find the max line length
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    buffCStr[buffSizeInt] = '\0'; /*make a c-string*/

    while(posInFileULng < lenFileULng)
    { /*Loop till at the end of the file*/
        /*Check if at the end of the file*/
        if(posInFileULng + buffSizeInt > lenFileULng)
            buffSizeInt = lenFileULng - posInFileULng;

        /*read in the file into the buffer*/
        fread(buffCStr, sizeof(char), buffSizeInt, samFile); /*read file*/
            /*buffer to read into, size of single element (char), buffer size,
              file to to read from*/
        posInFileULng = posInFileULng + buffSizeInt;
        tmpCStr = buffCStr; /*reset pointer to start of buffer*/

        /*find the length of the new line*/
        while(*tmpCStr != '\0')
        { /*Loop through buffer until new line or null found*/
            lenLineULng++; /*Add in the new character (since not '\0')*/

            if(*tmpCStr == '\n')
            { /*If on a new line*/
                if(lenLineULng > *maxLineLenULng) /*check if longest line*/
                    *maxLineLenULng = lenLineULng;
                lenLineULng = 0; /*Reset the current line length*/
                tabCntChar = 0; /*Reset for the next line*/
            } /*If on a new line*/

            if(*tmpCStr == '\t')
                tabCntChar++;

            tmpCStr++; /*Move to the next character in buffer*/
        } /*Loop through buffer until new line or null found*/
    } /*Loop till at the end of the file*/

    fseek(samFile, 0, SEEK_SET); /*go back to the start of the file*/
    return samFile;
} /*openFile*/

/*##############################################################################
# Name: stdoutReadScore
# Use: Outputs the read score line to stdout
# Input: 
#    samEntryStruct: Structer with stats to print out
#    printHeaderChar: 1 = print header line, 0 do not print header
# Output:
#     stdout: line with stats from samEntryStruct
#     Modifies: printHeaderChar to be 0 if set to 1
##############################################################################*/
void stdoutReadScore(samEntry *samEntryStruct, char *printHeaderChar)
{ /*outputReadScore*/
    char *tmpQueryCStr = (*samEntryStruct).queryCStr,
         *tmpRefCStr = (*samEntryStruct).refCStr;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-2 Sec-1: Print the samfile stats to stdout
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(*printHeaderChar == 1)
    { /*If wanted to print the header*/
        printf("Read\tRef\tMAPQ\treadLength\talignedLength\tmatches");
        printf("\tkeptMatches\tmismatches"); 
        printf("\tinsertions\tdeletions\tmedianQ\tmeanQ\talignedMedianQ"); 
        printf("\talignedMeanQ\tignoredMismatches\tignoredIndels"); 
        printf("\tignoredDeletions\n");
        *printHeaderChar = 0;
    } /*If wanted to print the header*/

    while(*tmpQueryCStr != '\t')
        tmpQueryCStr++; /*move to end of query name*/
    *tmpQueryCStr = '\0'; /*turn into c-string*/

    while(*tmpRefCStr != '\t')
        tmpRefCStr++; /*move to end of reference name*/
    *tmpRefCStr = '\0'; /*turn into c-string*/

    /*Print out the entry stats*/
  printf(
   "%s\t%s\t%u\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%f\t%f\t%f\t%f\t%lu\t%lu\t%lu\n",
       samEntryStruct->queryCStr, 
       samEntryStruct->refCStr,
       samEntryStruct->mapqUInt,
       samEntryStruct->readLenULng,
       samEntryStruct->readAligLenULng,
       samEntryStruct->numMatchULng,
       samEntryStruct->numKeptMatchULng,
       samEntryStruct->numMisULng,
       samEntryStruct->numInsULng,
       samEntryStruct->numDelULng,
       samEntryStruct->medianQDbl,
       samEntryStruct->meanQDbl,
       samEntryStruct->medianAligQDbl,
       samEntryStruct->meanAligQDbl,
       samEntryStruct->numIgnoreMisULng,
       samEntryStruct->numIgnoreInsULng,
       samEntryStruct->numIgnoreDelULng
    ); /*printf: print out stats*/

    *tmpQueryCStr = '\t'; /*turn back into single string*/
    *tmpRefCStr = '\t'; /*turn back into singel string*/
} /*outputReadScore*/

/*
 Output:
    Allocates: seqCStr to have the first sequence from the fastq file
        - Set to 0 & freed if no sequence line
    Allocates: qCStr to have the first sequences Q-score line from the fastq
        - Set to 0 & freed if no q-score line
    Returns: 1: if no malloc errors
             0: If malloc failed to find memory
*/
char readFirstSeqInFastq(
    char *fileInCStr, /*C-string with path & name of file to open*/
    char **seqCStr,   /*points to c-string with sequence, c-string is on heap*/
    char **qCStr,     /*points to c-string with q-score line, c-string on heap*/
    int buffSizeInt   /*Size of buffer to read file in with*/
) /*Gets the frist reads sequence & q-score line from a fastq file*/
{ /*openFile*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: Variable declarations
    #    fun-3 sec-2: Check if the input can be opened
    #    fun-3 sec-3: Find end of header (1st) line (start of sequence line)
    #    fun-3 sec-4: If read file in one go, find start & length of sequence
    #    fun-3 sec-5: Check if there is a sequence line in file
    #    fun-3 sec-6: Get sequence length when sequence not in buffer
    #    fun-3 sec-7: Get sequence length when file read in buffer
    #    fun-3 sec-8: allocate memory for the sequence & Q-score line
    #    fun-3 sec-9: Copy sequence when entire file is not in buffer
    #    fun-3 sec-10: Copy sequence when entire file is in buffer
    #    fun-3 sec-11: Get line between sequence and Q-score (+\n)
    #    fun-3 sec-12: Check if line between seqence and Q-score is present
    #    fun-3 sec-13: Read in Q-score when entire file is in buffer
    #    fun-3 sec-14: Clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-1: Variable declarations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Set up buffers to read each line*/
    char 
        buffCStr[buffSizeInt + 1], /*Buffer to read in part of the file*/
        *tmpCStr = buffCStr,       /*Pointer to advance through the buffer*/
        *tmpFastqCStr = 0,         /*Points to sequence or Q-score buffers*/
        enitreFileBool = 0;        /*1: read in entire file*/

    unsigned long 
        lenSeqULng = 0,            /*Length of sequence, same as q-score line*/
        seqStartULng = 0;          /*Marks start of sequence line in file*/

    FILE 
        *fastqFile = 0;            /*points to the fastq file*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-2: Check if the input can be opened
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    fastqFile = fopen(
                    fileInCStr,    /*File to open*/
                    "r"            /*Only read the file (as characters)*/
    ); /*open the file*/

    if(fastqFile == 0)
    { /*If there was nothing in the file, quite*/
        fclose(fastqFile);         /*close the file*/
        *seqCStr = 0;
        *qCStr = 0;
        fprintf(
            stderr,
            "%s is not a file\n",
            fileInCStr
        ); /*tell user not a valid file*/
        return 1;
    } /*If there was nothing in the file, quite*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-3: Find end of header (1st) line (start of sequence line)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    buffCStr[buffSizeInt] = '\0'; /*make a c-string*/

    while(                           /*Find the sequence line in the file*/
        fread(
            buffCStr,                /*buffer to put file input into*/
            sizeof(char),            /*size of each element in the file*/
            buffSizeInt,             /*size of my buffer*/
            fastqFile                /*file to read*/
        ) /*Read file in*/
        == buffSizeInt
    ) { /*While I have not found the sequence line*/
        tmpCStr = buffCStr;
        seqStartULng = 0;

        while(
            *tmpCStr != '\0' &&     /*Null marks end of the buffer*/
            *tmpCStr != '\n'        /*\n marks start of the sequence line*/
        ) { /*While not at end of first line*/
            tmpCStr++;
            seqStartULng++;
        } /*While not at end of first line*/

        if(*tmpCStr == '\n')
        { /*If on a new line*/
            seqStartULng++;         /*Move of the new line char*/
            break;                  /*At the start of the sequence line*/
        } /*If on a new line*/

        else
            lenSeqULng = ftell(fastqFile); /*reusing variable*/
    } /*While I have not found the sequence line*/

    lenSeqULng = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-4: If read file in one go, find start & length of sequence
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(seqStartULng == 0)
    { /*If read in entire file*/
        enitreFileBool = 1;
        tmpCStr = buffCStr;

        while(
            *tmpCStr != '\0' &&              /*End of buffer*/
            *tmpCStr != '\n'                 /*End of header line*/
        ) /*while not at the start of the sequence*/
            tmpCStr++;

    } /*If read in entire file*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-5: Check if there is a sequence line in file
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(*tmpCStr == '\0')
    { /* If is not a fastq file*/
        *seqCStr = 0;
        *qCStr = 0;
        fclose(fastqFile);
        fprintf(
            stderr,
            "%s is not a fastq file\n",
            fileInCStr
        ); /*tell user not a valid file*/

        return 1;
    } /* If is not a fastq file*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-6: Get sequence length when sequence not in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(enitreFileBool == 0)
    { /*If did not read in the entire file in one go*/
        fseek(
            fastqFile,    /*File to find postion in*/
            seqStartULng, /*My offset to first base of the sequence*/
            SEEK_SET /*File start. SEEK_SET, SEEK_END, & SEEK_CUR only options*/
        ); /*Find the start of sequence in the file*/
    
        while(                           /*Find the sequence line in the file*/
            fread(
                buffCStr,                /*buffer to put file input into*/
                sizeof(char),            /*size of each element in the file*/
                buffSizeInt,             /*size of my buffer*/
                fastqFile                /*file to read*/
            ) /*Read file in*/
            == buffSizeInt
        ) { /*While I have not found the sequence line*/
            tmpCStr = buffCStr;
    
            while(
                *tmpCStr != '\0' &&      /*Null marks end of the buffer*/
                *tmpCStr != '\n'         /*\n marks start of the sequence line*/
            ) { /*While not at end of first line*/
                tmpCStr++;
                lenSeqULng++;            /*Count each base in the sequence*/
            } /*While not at end of first line*/
    
            if(*tmpCStr == '\n')
                break;           /*At end of the sequence line*/
        } /*While I have not found the sequence line*/
    } /*If did not read in the entire file in one go*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-7: Get sequence length when file read in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    else
    { /*Else the entire sequence is in the buffer*/
        tmpCStr++;                       /*move onto sequence line*/

         while(
             *tmpCStr != '\0' &&         /*End of buffer*/
             *tmpCStr != '\n'            /*End of header line*/
         ) { /*While not at end of first line*/
            tmpCStr++;
            lenSeqULng++;                /*Get the sequence length*/
         } /*While not at end of first line*/

         /*Go back to the start of the sequence*/
         tmpCStr--;

         while(*tmpCStr != '\n')
             tmpCStr--;
         tmpCStr++;                      /*Get of \n*/
    } /*Else the entire sequence is in the buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-8: allocate memory for the sequence & Q-score line
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    *seqCStr = malloc(sizeof(char) * (lenSeqULng + 1));
    if(seqCStr == 0)
    { /*If malloc could not allocate memory for the sequence*/
        *seqCStr = 0;
        *qCStr = 0;
        fclose(fastqFile);
        fprintf(
            stderr,
            "Malloc failed fun-3: readFirstSeqInFastq: scoreReadsFileIO.c:350\n"
        ); /*Print error to user*/
        return 0;
    } /*If malloc could not allocate memory for the sequence*/

    *qCStr = malloc(sizeof(char) * (lenSeqULng + 1));
    if(qCStr == 0)
    { /*If malloc could not allocate memory for the sequence*/
        free(*seqCStr);
        *seqCStr = 0;
        *qCStr = 0;
        fclose(fastqFile);
        fprintf(
            stderr,
            "Malloc failed fun-3: readFirstSeqInFastq: scoreReadsFileIO.c:367\n"
        ); /*Print error to user*/
        return 0;
    } /*If malloc could not allocate memory for the sequence*/

    *(*seqCStr + lenSeqULng) = '\0';
    *(*qCStr + lenSeqULng) = '\0';

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-9: Copy sequence when entire file is not in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(enitreFileBool == 0)
    { /*If did not read in the entire file in one go*/
        fseek(
            fastqFile,    /*File to find postion in*/
            seqStartULng, /*My offset to first base of the sequence*/
            SEEK_SET /*File start. SEEK_SET, SEEK_END, & SEEK_CUR only options*/
        ); /*Find the start of sequence in the file*/
    
        fread(
            *seqCStr,                /*Copy sequence to */
            sizeof(char),           /*size of each element in the file*/
            lenSeqULng,             /*size of my buffer*/
            fastqFile                /*file to read*/
        ); /*Read file in*/
    } /*If did not read in the entire file in one go*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-10: Copy sequence when entire file is in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    else
    { /*Else entire file is in buffer*/

        tmpFastqCStr = *seqCStr;

        while(
            *tmpCStr != '\n' &&
            *tmpCStr != '\0'
        ) { /*While there are bases to be read in*/
            *tmpFastqCStr = *tmpCStr;
            tmpCStr++;
            tmpFastqCStr++;
        } /*While there are bases to be read in*/

        *tmpFastqCStr = '\0';        /*Make sure c-string includes no junk*/
    } /*Else entire file is in buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-11: Get line between sequence and Q-score (+\n)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(enitreFileBool == 0)
    { /*if the entire file is not in the buffer*/
        if(
            fread(
                buffCStr,                /*Copy sequence to */
                sizeof(char),            /*size of each element in the file*/
                3,                       /*size of my buffer*/
                fastqFile                /*file to read*/
            ) /*Read file in*/
            != 3
         ) /*If nothing else in file*/
             buffCStr[0] = '\0';         /*Marks there is no Q-score line*/

        tmpCStr = buffCStr;
    } /*if the entire file is not in the buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-12: Check if line between seqence and Q-score is present
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(
        *tmpCStr != '\n' ||
        *(tmpCStr + 1) != '+' ||
        *(tmpCStr + 2) != '\n'
    )
    { /*If has no Q-score line*/
        free(*qCStr);
        qCStr = 0;
        fclose(fastqFile);
        return 1;
    } /*If has no Q-score line*/

    if(enitreFileBool == 1) tmpCStr += 3;      /*Get to Q-score line*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-11: Read in Q-score when entire file not in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(enitreFileBool == 0)
    { /*If did not read in the entire file in one go*/
        fread(
            *qCStr,                /*Copy sequence to */
            sizeof(char),           /*size of each element in the file*/
            lenSeqULng,             /*size of my buffer*/
            fastqFile                /*file to read*/
        ); /*Read file in*/
    } /*If did not read in the entire file in one go*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-13: Read in Q-score when entire file is in buffer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    else
    { /*Else need to copy Q-score line from buffer*/
        tmpFastqCStr = *qCStr;

        while(
            *tmpCStr != '\n' &&
            *tmpCStr != '\0'
        ) { /*While there are bases to be read in*/
            *tmpFastqCStr = *tmpCStr;
            tmpCStr++;
            tmpFastqCStr++;
        } /*While there are bases to be read in*/

        *tmpFastqCStr = '\0';        /*Make sure c-string includes no junk*/
    } /*Else need to copy Q-score line from buffer*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Fun-3 Sec-14: Clean up
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    fclose(fastqFile);
    return 1;
} /*openFile*/
