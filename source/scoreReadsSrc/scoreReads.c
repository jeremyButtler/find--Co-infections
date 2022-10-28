/*##############################################################################
# Name: scoreReads.c
# Use: Scores and prints out scores for reads in a samfile
# Input:
#    -f:
#        - Samfile to score
#        - Required if -stdin not used
#        - Needs the --eqx cigar (minimap2 --eqx)
#    -stdin:
#        - Take piped input from the command line (-f command ignored)
#        - Requires: len-max-seq to set the input buffer
#    -ref:
#        - Fastq with reference reads were mapped against
#        - Note: Only the top sequence is used
#        - uses reference in checking matches, mismatches, & indels
#        - Default: Do not use reference
#    -ref-del:
#        - Fastq with reference reads were mapped against
#        - Note: Only the top sequence is used 
#        - uses reference to validate deletions
#        - Default: Do not use reference
#    -len-seq:
#        - The length of the longest sequence in your fastq/sam file
#        - Default: 100000000 (not used when -f input)
#    -min-q: Min Q-score to replace low quality bases with (0 to 97)
#        Default 13
#    -min-map-q: Min mapping quality Q-score needed to keep read
#        Defualt: 20
#    -min-read-length: Min aligned read length to keep read 
#        Default: 600
#    -max-read-length: Max read length to keep a read (0 is any length)
#        Default: 1000
#    -min-mean-q: Min mean Q-score needed to keep a read
#        Default: 13
#    -min-median-q: Min median Q-score needed to keep a read
#        Default: 13
#    -min-aligned-mean-q: Min mean Q-score of the alinged sections of a read
#                         needed to keep a read
#        Default: 13
#    -min-aligned-median-q: Min median Q-score of the alinged sections of a read
#                         needed to keep a read
#        Default: 13
#
#    -ins-A: Max A homopolymer length to keep a insertion
#        Default: 2
#    -ins-T: Max T homopolymer length to keep a insertion
#        Default: 2
#    -ins-G: Max G homopolymer length to keep a insertion
#        Default: 1
#    -ins-C: Max C homopolymer length to keep a insertion
#        Default: 1
#    -del-A: Max A homopolymer length to keep a deletion
#        Default: 0
#    -del-T: Max T homopolymer length to keep a deletion
#        Default: 0
#    -del-G: Max G homopolymer length to keep a deletion
#        Default: 0
#    -del-C: Max C homopolymer length to keep a deletion
#        Default: 0
# Output:
#    stdout: Line with the read, query, and scores
# NOTE: ADD STDIN INPUT OPITON LATER
##############################################################################*/

/*##############################################################################
# TOC:
#   fun-1: main: Main function to glue everything together
#   char checkInput: Checks the user input [returns: 0 if invalid]
##############################################################################*/

/* Line to build program
egcc -g -Wall \
	*.c \
	-o scoreReads
*/
#include <string.h> /*strcmp function*/

#include "scoreReadsScoringWrappers.h"
#include "scoreReadsConversion.h" /*conversion functions*/

/*##############################################################################
# Output: Modifies: Each input variable to hold user input
##############################################################################*/
char checkInput(
    int *lenArgsInt,               /*Number arguments & paramteters user input*/
    char *argsCStr[],              /*Array with user arguments & parameters*/
    char **samPathCStr,            /*file name & path to same file to score*/
    char **refPathCStr,            /*file name & path to mapping reference*/
    char *refForDelChar,           /*1: Only use reference for deletions*/
    struct minStats *minReadStats, /*holds min thresholds user provides*/
    unsigned long *longestSeqULng, /*Holds user provided longest seq length*/
    char *stdinChar                /*Set 1: if input comes from stdin
                                     Set 0: If input comes from file*/
); /*Checks & extracts user input*/

int main(int lenArgsInt, char *argsPtrCStr[])
{ /*main function*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1: main
    #    fun-1 sec-1: Variable declarations
    #    fun-1 sec-2: Read in and check user input
    #    fun-1 sec-3: Call the read scoring functions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declarations
    #    fun-1 sec-1 sub-1: normal variable declerations
    #    fun-1 sec-1 sub-2: help message
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-1 Sub-1: normal variable declerations
    ***************************************************************************/

    char
        *samPathCStr = 0,            /*Path to sam file to work on*/
        *refPathCStr = 0,            /*Fastq file with mapping reference*/
        refForDelChar = 0,          /*1: only checking deltions with reference*/
        stdinChar = 0;               /*Char makring if input is from stdin*/
    /*max line size and buffer size limits*/

    int
        buffSizeInt = 10000;          /*Max buffer size for openFile*/

    unsigned long
        maxLineLenULng = 100000000; /*Max line length in same file*/

    FILE
        *samFile = 0; /*Points to file to get data from*/

    struct minStats
        minReadStats;

    struct samEntry
        refEntry;                   /*Holds reference stats*/

    /***************************************************************************
    # Fun-1 Sec-1 Sub-2: help message
    ***************************************************************************/

    char
        *helpMesgCStr = "\
            \n Command: scoreReads -f reads.sam [options ...]\
            \n Use:\
            \n   Finds the aligned length, Q-scores, number of mismatches, &\
            \n     number of indels for all reads in the input sam file.\
            \n   Mismatches & indels are only kept if meet input thresholds\
            \n Output:\
            \n   - stdout: Stats printed in tsv format\
            \n Input:\
            \n  -stdin:\
            \n    - Take input from command line instd of file     [Not used]\
            \n  -len-seq:\
            \n    - Length of longest sequence in reads.sam        [100000000]\
            \n  -ref:\
            \n    - Fastq with reference used in mapping           [None]\
            \n    - Keeps mismatch/indel if ref & read support\
            \n    - Note: Only top reference is used\
            \n  -ref-del:\
            \n    - Fastq with reference used in mapping           [None]\
            \n    - Keeps deletion only if reference supports\
            \n    - Note: Only the top sequence is used\
            \n  -min-q:\
            \n    - Min Q-score to keep snp or indel               [13]\
            \n  -min-map-q:\
            \n    - Min mapping quality needed for read            [20]\
            \n  -min-read-length:\
            \n    - Min aligned read length                        [600]\
            \n  -max-read-length:\
            \n    - Max read length (0 = no max)                   [1000]\
            \n  -min-mean-q:\
            \n    - Min mean Q-score to keep read                  [13]\
            \n  -min-median-q:\
            \n    - Min median Q-score to keep read                [13]\
            \n  -min-aligned-mean-q:\
            \n    - Min mean Q of alinged read                     [1]\
            \n  -min-aligned-median-q:\
            \n    - Min read med aligned Q                         [13]\
            \n  -ins-A:\
            \n    - Max A homopolymer length to keep insert        [2]\
            \n  -ins-T:\
            \n    - Max T homopolymer length to keep insert        [2]\
            \n  -ins-G:\
            \n    - Max G homopolymer length to keep insert        [1]\
            \n  -ins-C:\
            \n    - Max C homopolymer length to keep insert        [1]\
            \n  -del-A:\
            \n    - Max A homopolymer length to keep del           [0]\
            \n  -del-T:\
            \n    - Max T homopolymer length to keep del           [0]\
            \n  -del-G:\
            \n    - Max G homopolymer length to keep del           [0]\
            \n  -del-C:\
            \n    - Max C homopolymer length to keep del           [0]\
            \n";

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Read in and check user input
    #    fun-1 sec-2 sub-1: get user input & check if need to print help message
    #    fun-1 sec-2 sub-2: Check if stdin buffer length input is valid
    #    fun-1 sec-2 sub-3: Check if min Q-score is valid
    #    fun-1 sec-2 sub-4: Check if can open file (not using stdin)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-1: get user input & check if need to print help message
    ***************************************************************************/

    blankMinStats(&minReadStats); /*Intalize my min stats to default*/
    blankSamFileEntry(&refEntry); /*Intalize my refence entry struct*/

    /*Check the user input*/
    if(
        checkInput(
            &lenArgsInt,
            argsPtrCStr,
            &samPathCStr,
            &refPathCStr,
            &refForDelChar,
            &minReadStats,
            &maxLineLenULng,
            &stdinChar
        ) /*Check & get user input*/
        < 1
    ) { /*if the user input an invalid input*/
        fprintf(
            stderr,
            "%s\n",
            helpMesgCStr
        );
        exit(-1);
    } /*if the user input an invalid input*/

    /***************************************************************************
    # Fun-1 Sec-2 Sub-2: Check if user provided buffer length for stdin is valid
    ***************************************************************************/

    if(
        stdinChar == 1 &&
        maxLineLenULng < 100
    ) { /*If user is not providing max line length for stdin*/
        fprintf(
            stderr,
            "%s\nFor -stdin -len-seq must be greater than 100:\n",
            helpMesgCStr
        ); /*print out low sequence length*/

        exit(-1);
    } /*If user is not providing max line length for stdin*/

    /**************************************************************************
    # Fun-1 Sec-2 Sub-3: Check if min Q-score is valid
    ***************************************************************************/

    if(minReadStats.minQChar < 0)
    { /*if the user input an invalid q-score*/
        fprintf(
            stderr,
            "Q-score must be greater than -1, but %i input\n",
            minReadStats.minQChar
        ); /*Print out invalid minimum Q-score*/
        return 1;
    } /*if the user input an invalid q-score*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Open sam file for reading
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(stdinChar != 1)
    { /*If using a file for input*/
        samFile = openFile(
                      samPathCStr,
                      &maxLineLenULng,
                      buffSizeInt
        ); /*Open same file*/
    
        if(samFile == 0)
        { /*If was unable to open the sam file*/
            fprintf(
                stderr,
                "%s\nCan not open file or file has nothing (%s)\n",
                helpMesgCStr,
                samPathCStr
            ); /*If file was not valid*/

            exit(-1);
        } /*If was unable to open the sam file*/
    } /*If using a file for input*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Open reference fastq for reading
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refPathCStr != 0)
    { /*If user provded the reference the reads mapped to*/
        if(
            readFirstSeqInFastq(
                refPathCStr,       /*Path & name of file to open*/
                &refEntry.seqCStr, /*Will point to the sequence line*/
                &refEntry.qCStr,   /*Will point to the Q-score line*/
                1000               /*Size of buffer (how many times read file)*/
            ) /*Read in the reference sequence and Q-score lines*/
            == 0
        ) /*Open same file*/
        { /*If malloc failed to create memory*/
            fprintf(
                stderr,
                "Malloc failed to allocate memory: main: scoreReads.c:305"
            );
            if(samFile != 0)
                fclose(samFile);   /*Need to release the file*/
            exit(-1);
        } /*If malloc failed to create memory*/
    
        else if(refEntry.seqCStr == 0)
        { /*If was unable to open the sam file*/
            fprintf(
                stderr,
                "%s\nCan not open file or file has nothing (%s)\n",
                helpMesgCStr,
                samPathCStr
            ); /*If file was not valid*/
            if(samFile != 0)
                fclose(samFile);   /*Need to release the file*/
            exit(-1);
        } /*If was unable to open the sam file*/
    } /*If user provded the reference the reads mapped to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-5: Call stdin read scoring functions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(stdinChar == 1)
    { /*If taking input from stdin*/
        maxLineLenULng = maxLineLenULng * 3; /*Worst case cigar = seq = q line*/

        if(refPathCStr == 0)       /*Check if a mapping refeence was provided*/
        { /*If no mapping reference provided*/
            scoreReadsStdin(
                maxLineLenULng,
                &minReadStats
            ); /*Score each alignment using only the cigar*/
        } /*If no mapping reference provided*/

        else
        { /*Else a mapping reference was provided*/
            refScoreReadsStdin(
                maxLineLenULng,
                &minReadStats,
                refForDelChar,
                &refEntry
            ); /*Score each alignment using the provided reference & cigar*/
        } /*Else a mapping reference was provided*/
    } /*If taking input from stdin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-6: Call file read scoring functions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    else
    { /*Else if getting input from a file*/
        if(refPathCStr == 0)
        { /*If no mapping reference was provided*/
            scoreReadsFile(
               samFile,
               maxLineLenULng,
               &minReadStats
            ); /*Sore reads using input file & cigar*/
        } /*If no mapping reference was provided*/

        else
        { /*Else an mapping reference was provided*/
            refScoreReadsFile(
                samFile,
                maxLineLenULng,
                &minReadStats,
                refForDelChar,
                &refEntry
            ); /*Score reads in file using cigar & reference*/
        } /*Else an mapping reference was provided*/
    } /*Else if getting input from a file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-7: Clean up
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refEntry.seqCStr != 0)
        free(refEntry.seqCStr);

    if(refEntry.qCStr != 0)
        free(refEntry.qCStr);

    exit(0);
} /*main function*/

/*##############################################################################
# Output: Modifies: Each input variable to hold user input
##############################################################################*/
char checkInput(
    int *lenArgsInt,               /*Number arguments & paramteters user input*/
    char *argsCStr[],              /*Array with user arguments & parameters*/
    char **samPathCStr,            /*file name & path to same file to score*/
    char **refPathCStr,            /*file name & path to mapping reference*/
    char *refForDelChar,           /*1: Only use reference for deletions*/
    struct minStats *minReadStats, /*holds min thresholds user provides*/
    unsigned long *longestSeqULng, /*Holds user provided longest seq length*/
    char *stdinChar                /*Set 1: if input comes from stdin
                                     Set 0: If input comes from file*/
) /*Checks & extracts user input*/
{ /*checkInput*/
    char *tmpCStr = 0, *singleArgCStr = 0;

    if(*lenArgsInt < 2)
        return 0; /*no arguments input*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
    { /*loop through all user input arguments*/           /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1);          /*supplied argument*/
        tmpCStr = *(argsCStr + intArg);                   /*Paramter*/

        if(strcmp(tmpCStr, "-f") == 0)
            *samPathCStr = singleArgCStr;
        else if(strcmp(tmpCStr, "-stdin") == 0)
        { /*If if taking input from stdin*/
            *stdinChar = 1;
            intArg--; /*Account for incurment at end of loop*/
        } /*If taking input from stdin*/
        else if(strcmp(tmpCStr, "-ref") == 0)
            *refPathCStr = singleArgCStr;
        else if(strcmp(tmpCStr, "-ref-del") == 0)
        { /*Else if have refence, but only want to check for deletions*/
            *refPathCStr = singleArgCStr;
            *refForDelChar = 1;
        } /*Else if have refence, but only want to check for deletions*/

        else if(strcmp(tmpCStr, "-len-seq") == 0)
            *longestSeqULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-min-q") == 0)
            (*minReadStats).minQChar = cStrToChar(singleArgCStr, '\0');
        else if(strcmp(tmpCStr, "-min-map-q") == 0)
            cStrToUInt(singleArgCStr, &(*minReadStats).minMapqUInt, '\0');

        /*Read lengths*/
        else if(strcmp(tmpCStr, "-min-read-length") == 0)
            (*minReadStats).minReadLenULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-max-read-length") == 0)
            (*minReadStats).maxReadLenULng = strtoul(singleArgCStr, NULL, 10);

        /*Get the min read Q-score values*/
        else if(strcmp(tmpCStr, "-min-mean-q") == 0)
            sscanf(singleArgCStr, "%lf", &(*minReadStats).minMeanQDbl);
        else if(strcmp(tmpCStr, "-min-median-q") == 0)
            sscanf(singleArgCStr, "%lf", &(*minReadStats).minMedianQDbl);
        else if(strcmp(tmpCStr, "-min-aligned-mean-q") == 0)
            sscanf(singleArgCStr, "%lf", &(*minReadStats).minAlignedMeanQDbl);
        else if(strcmp(tmpCStr, "-min-aligned-median-q") == 0)
            sscanf(singleArgCStr, "%lf", &(*minReadStats).minAlignedMedianQDbl);

        /*Max Insertion lengths*/
        else if(strcmp(tmpCStr, "-ins-A") == 0)
            (*minReadStats).maxInsAHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-ins-T") == 0)
            (*minReadStats).maxInsTHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-ins-G") == 0)
            (*minReadStats).maxInsGHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-ins-C") == 0)
            (*minReadStats).maxInsCHomoULng = strtoul(singleArgCStr, NULL, 10);

        /*Max deletion length*/
        else if(strcmp(tmpCStr, "-del-A") == 0)
            (*minReadStats).maxDelAHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-del-T") == 0)
            (*minReadStats).maxDelTHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-del-G") == 0)
            (*minReadStats).maxDelGHomoULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-del-C") == 0)
            (*minReadStats).maxDelCHomoULng = strtoul(singleArgCStr, NULL, 10);

        else
            return 0;
        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    tmpCStr = *samPathCStr;

    if(*stdinChar != 1)
    { /*If taking input from a file*/
        while(*tmpCStr != '\0') /*find end of c-string*/
            tmpCStr++;

        if(*(tmpCStr - 1)!='m' || *(tmpCStr - 2)!='a' || *(tmpCStr - 3)!='s')
        { /*If file is not a sam file*/
            fprintf(
                stderr,
                "%s is not a sam file\n",
                *samPathCStr
            ); /*Tell user input was not a sam file*/
            return 0; /*If file does not end in sam*/
        } /*If file is not a sam file*/
    } /*If taking input from a file*/

    return 1; /*input is valid*/
} /*checkInput*/
