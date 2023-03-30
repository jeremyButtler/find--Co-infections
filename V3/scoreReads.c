/*######################################################################
# Name: scoreReads.c
# Use: Scores and prints out scores for reads in a samfile
# Input:
#    -file:
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
#    -min-aligned-mean-q: Min mean Q-score of the alinged sections of a
#                         read needed to keep a read
#        Default: 13
#    -min-aligned-median-q: Min median Q-score of the alinged sections
#                           of a read needed to keep a read
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
# Includes:
#    - "scoreReadsFun.h"
#        - "minAlnStatsStructs.h"
#        - "fqAndFqFun.h"
#        - "FCIStatsFun.h"
#        o "samEntryStruct.h"
#            - <stdlib.h>
#            - "cStrToNumberFun.h"
#                - <sdtint.h>
#            - "printError.h"
#                - <stdio.h>
######################################################################*/

/*######################################################################
# TOC:
#   main: Main function to glue everything together
#   fun-2 checkInput: Checks the user input [returns: 0 if invalid]
######################################################################*/

#include <string.h> /*strcmp function*/
#include "scoreReadsFun.h"/*Structs & functions specific to scoreReads*/

/*######################################################################
# Output: Modifies: Each input variable to hold user input
######################################################################*/
char checkInput(
    int *lenArgsInt,          /*size of argsCStr*/
    char *argsCStr[],         /*Array with user arguments & parameters*/
    char **samPathCStr,          /*file path to sam file to score*/
    char **refPathCStr,          /*file path to mapping reference*/
    uint8_t *refForDelUC,      /*Set to 1 if: use ref for dels only*/
    struct minAlnStats *minStats, /*min thresholds user provides*/
    char *stdinChar                /*Set 1: if input comes from stdin
                                     Set 0: If input comes from file*/
); /*Checks & extracts user input*/

int main(int lenArgsInt, char *argsPtrCStr[])
{ /*main function*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main TOC:
    #    main sec-1: Variable declarations
    #    main sec-2: Read in and check user input
    #    main sec-3: Call the read scoring functions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-1: Variable declarations
    #    main sec-1 sub-1: normal variable declerations
    #    main sec-1 sub-2: help message
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Main Sec-1 Sub-1: normal variable declerations
    *******************************************************************/

    char *samPathCStr = 0;   /*Path to sam file to work on*/
    char *refPathCStr = 0;   /*Fastq file with mapping reference*/
    char stdinChar = 0;      /*Char makring if input is from stdin*/
    uint8_t refForDelUC = 0;  /*1: checking deltions with reference*/

    FILE *samFILE = 0;       /*Points to file to get data from*/
    FILE *refFILE = 0;       /*reference file to use in comparison*/
    FILE *outFILE = stdout;  /*File to put output in*/

    struct minAlnStats minStats;    /*Holds user input*/

    /*******************************************************************
    # Main Sec-1 Sub-2: help message
    *******************************************************************/

    char
        *helpMesgCStr = "\
            \n Command: scoreReads -file reads.sam [options ...]\
            \n Use:\
            \n   Finds the aligned length, Q-scores, number of\
            \n     mismatches, & number of indels for all reads in the\
            \n     input sam file. Mismatches & indels are only kept if\
            \n     they meet the min input thresholds\
            \n Output:\
            \n   - stdout: Stats printed in tsv format\
            \n Input:\
            \n  -file:\
            \n    - Take input from a input file            [Default]\
            \n  -stdin:\
            \n    - Take input from command line            [Use file]\
            \n  -ref:\
            \n    - Fastq with reference used in mapping    [None]\
            \n    - Keeps a mismatch/indel if ref & read support\
            \n    - Note: Only the first reference is used\
            \n    - Note: Deletions are treated like insertions\
            \n  -ref-del:\
            \n    - Like -ref, but only uses reference for  [None]\
            \n      for evaluating deletions\
            \n  -min-q:\
            \n    - Min Q-score to keep SNP or indel        [13]\
            \n  -min-map-q:\
            \n    - Min mapping quality needed for read     [20]\
            \n  -min-read-length:\
            \n    - Min aligned read length                 [600]\
            \n  -max-read-length:\
            \n    - Max read length (0 = no max)            [1000]\
            \n  -min-mean-q:\
            \n    - Min mean Q-score to keep read           [13]\
            \n  -min-median-q:\
            \n    - Min median Q-score to keep read         [13]\
            \n  -min-aligned-mean-q:\
            \n    - Min mean Q of alinged read              [1]\
            \n  -min-aligned-median-q:\
            \n    - Min read med aligned Q                  [13]\
            \n  -ins-A:\
            \n    - Max A homopolymer length to keep insert [2]\
            \n  -ins-T:\
            \n    - Max T homopolymer length to keep insert [2]\
            \n  -ins-G:\
            \n    - Max G homopolymer length to keep insert [1]\
            \n  -ins-C:\
            \n    - Max C homopolymer length to keep insert [1]\
            \n  -del-A:\
            \n    - Max A homopolymer length to keep del    [0]\
            \n  -del-T:\
            \n    - Max T homopolymer length to keep del    [0]\
            \n  -del-G:\
            \n    - Max G homopolymer length to keep del    [0]\
            \n  -del-C:\
            \n    - Max C homopolymer length to keep del    [0]\
            \n";

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-2: Read in and check user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankMinStats(&minStats); /*Intalize my min stats to default*/

    if(
        checkInput(
            &lenArgsInt,
            argsPtrCStr,
            &samPathCStr,
            &refPathCStr,
            &refForDelUC,
            &minStats,
            &stdinChar
        ) == 0 /*Check if input is valid (also get input)*/
    ) { /*if the user input an invalid input*/
        if(
            strcmp(samPathCStr, "-h") ||
            strcmp(samPathCStr, "--h") ||
            strcmp(samPathCStr, "-help") ||
            strcmp(samPathCStr, "--help")
        ) { /*If user requested the help message*/
            fprintf(stdout, "%s\n", helpMesgCStr);
            exit(0);
        } /*If user requested the help message*/

        fprintf(stderr, "%s\n", helpMesgCStr);
        exit(-1);
    } /*if the user input an invalid input*/

    if(minStats.minQChar < 0)
    { /*if the user input an invalid q-score*/
        fprintf(stderr, "Min q-score (%i) < 0\n",minStats.minQChar);
        return 1;
    } /*if the user input an invalid q-score*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Open sam file for reading
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(stdinChar != 1)
    { /*If using a file for input*/
        samFILE = fopen(samPathCStr, "r");
    
        if(samFILE == 0)
        { /*If was unable to open the sam file*/
            fprintf(stderr, "Can not open file (%s)\n", samPathCStr);
            exit(-1);
        } /*If was unable to open the sam file*/
    } /*If using a file for input*/

    else
        samFILE = stdin;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-4: Open reference fastq for reading
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refPathCStr != 0)
    { /*If user provded the reference the reads mapped to*/
        refFILE = fopen(refPathCStr, "r");
    
        if(refFILE == 0)
        { /*If was unable to open the reference file*/
            fprintf(stderr, "Can not open file (%s)\n", refPathCStr);
            exit(-1);
        } /*If was unable to open the reference file*/
    } /*If user provded the reference the reads mapped to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-5: Call stdin read scoring functions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    scoreReads(&minStats, &refForDelUC, samFILE, refFILE, outFILE);

    fclose(samFILE);

    if(outFILE != 0)
        fclose(outFILE);

    if(refFILE != 0)
        fclose(refFILE);

    exit(0);
} /*main function*/

/*######################################################################
# Output: Modifies: Each input variable to hold user input
######################################################################*/
char checkInput(
    int *lenArgsInt,          /*size of argsCStr*/
    char *argsCStr[],         /*Array with user arguments & parameters*/
    char **samPathCStr,          /*file path to sam file to score*/
    char **refPathCStr,          /*file path to mapping reference*/
    uint8_t *refForDelUC,        /*Set to 1 if: use ref for dels only*/
    struct minAlnStats *minStats, /*min thresholds user provides*/
    char *stdinChar                /*Set 1: if input comes from stdin
                                     Set 0: If input comes from file*/
) /*Checks & extracts user input*/
{ /*checkInput*/
    char *tmpCStr = 0, *singleArgCStr = 0;

    if(*lenArgsInt < 2)
        return 0; /*no arguments input*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg += 2)
    { /*loop through all user input arguments*/    /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1);   /*supplied argument*/
        tmpCStr = *(argsCStr + intArg);            /*Paramter*/

        if(strcmp(tmpCStr, "-file") == 0)
            *samPathCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-stdin") == 0)
        { /*If if taking input from stdin*/
            *stdinChar = 1;
            --intArg; /*Account for +2 incurment*/
        } /*If taking input from stdin*/

        else if(strcmp(tmpCStr, "-ref") == 0)
            *refPathCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-ref-del") == 0)
        { /*Else if have refence, but only want to check for deletions*/
            *refPathCStr = singleArgCStr;
            *refForDelUC = 1;
        } /*Else if have refence, but only want to check for deletions*/

        else if(strcmp(tmpCStr, "-min-q") == 0)
            cStrToUChar(singleArgCStr, &(minStats->minQChar));
        else if(strcmp(tmpCStr, "-min-map-q") == 0)
            cStrToUInt(singleArgCStr, &(minStats->minMapqUInt));

        /*Read lengths*/
        else if(strcmp(tmpCStr, "-min-read-length") == 0)
            minStats->minReadLenULng = strtoul(singleArgCStr,NULL,10);
        else if(strcmp(tmpCStr, "-max-read-length") == 0)
            minStats->maxReadLenULng = strtoul(singleArgCStr, NULL, 10);

        /*Get the min read Q-score values*/
        else if(strcmp(tmpCStr, "-min-mean-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMeanQFlt);
        else if(strcmp(tmpCStr, "-min-median-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMedianQFlt);
        else if(strcmp(tmpCStr, "-min-aligned-mean-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minAlignedMeanQFlt);
        else if(strcmp(tmpCStr, "-min-aligned-median-q") == 0)
            sscanf(singleArgCStr,"%f",&minStats->minAlignedMedianQFlt);

        /*Max Insertion lengths*/
        else if(strcmp(tmpCStr, "-ins-A") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[0]);
        else if(strcmp(tmpCStr, "-ins-T") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[10]);
        else if(strcmp(tmpCStr, "-ins-G") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[1]);
        else if(strcmp(tmpCStr, "-ins-C") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[3]);

        /*Max Delertion lengths*/
        else if(strcmp(tmpCStr, "-del-A") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[0]);
        else if(strcmp(tmpCStr, "-del-T") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[10]);
        else if(strcmp(tmpCStr, "-del-G") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[1]);
        else if(strcmp(tmpCStr, "-del-C") == 0)
            cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[3]);

        else
        { /*Else invalid input*/
            *samPathCStr = tmpCStr;
            return 0;
        } /*Else invalid input*/
    } /*loop through all user input arguments*/

    tmpCStr = *samPathCStr;

    if(*stdinChar != 1)
    { /*If taking input from a file*/
        while(*tmpCStr != '\0') /*find end of c-string*/
            tmpCStr++;

        if(
            *(tmpCStr - 3)!='s' ||
            *(tmpCStr - 2)!='a' ||
            *(tmpCStr - 1)!='m')
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
