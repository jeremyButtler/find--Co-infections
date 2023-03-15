/*######################################################################
# Name: filterReads
# Use:
#   o Selecs the top reads in a file
# Input:
#    -fastq file.fastq:                                      [Required]
#      o Fastq file to extract reads from
#    -ref                                                    [None]
#      o Fasta file with reference to map reads to.
#      o Default is to use length and median Q-scores.
#    -out:                                                   [stdout]
#      o output file to print filtered reads to
#    -threads:                                               [3]
#      o Number of threads to use with minimap2.
#      o Only applied when -ref is selected.
#    -num-reads:                                             [300]
#      o Maximum number of reads to extract
#    -min-mapq                                               [20]
#      o min mapping quality to keep a read (-ref only)
#    -min-median-q:                                          [10]
#      o Remove reads with median Q-scores under
#        this Q-score
#      o Not applied when -ref is selected.
#    -min-mean-q:                                            [10]
#      o Remove reads with mean Q-scores under
#        this Q-score
#      o Not applied when -ref is selected.
#    -max-length:                                            [1000]
#      o Remove reads over this length (0 = all lengths)
#      o Not applied when -ref is selected.
#    -min-length:                                            [500]
#      o Remove reads under this length
#      o Is alinged length instead of read length for -ref
#    -min-q:                                                 [10]
#      o Minimum quality score to count an SNP or inseertion
#      o Only applied when -ref is selected.
#    -max-snps:                                              [0.07=7%]
#      o Maximum percent difference in SNPs to keep a read
#      o Only applied when -ref is selected.
#    -max-diff:                                              [1=100%]
#      o Maximum percent difference to keep a read
#      o Only applied when -ref is selected.
#    -max-inss:                                              [1=100%]
#      o Maximum percent difference in insertions to keep a read
#      o Only applied when -ref is selected.
#    -max-dels:                                              [1=100%]
#      o Maximum percent difference in deletions to keep a read
#      o Only applied when -ref is selected.
#    -max-indels:                                            [1=100%]
#      o Maximum percent difference in indels to keep a read
#      o Only applied when -ref is selected.
#    -max-ins-a-homo:                                        [1]
#      o Maximum hompolymer size an A insertion can be part of
#        to keep the insertion.
#      o Only applied when -ref is selected.
#    -max-ins-t-homo:                                        [1]
#      o Maximum hompolymer size an T insertion can be part of
#        to keep the insertion.
#      o Only applied when -ref is selected.
#    -max-ins-g-homo:                                        [1]
#      o Maximum hompolymer size an G insertion can be part of
#        to keep the insertion.
#      o Only applied when -ref is selected.
#    -max-ins-c-homo:                                        [1]
#      o Maximum hompolymer size an C insertion can be part of
#      o Only applied when -ref is selected.
#    -max-del-a-homo:                                        [1]
#      o Maximum hompolymer size an A deletion can be part of
#        to keep the delertion.
#      o Only applied when -ref is selected.
#    -max-del-t-homo:                                        [1]
#      o Maximum hompolymer size an T deletion can be part of
#        to keep the delertion.
#      o Only applied when -ref is selected.
#    -max-del-g-homo:                                        [1]
#      o Maximum hompolymer size an G deletion can be part of
#        to keep the delertion.
#      o Only applied when -ref is selected.
#    -max-del-c-homo:                                        [1]
#      o Maximum hompolymer size an C deletion can be part of
#      o Only applied when -ref is selected.
#    -v:
#      o Print version & exit
# Output:
#    stdout: prints out the top x reads reads to stdout
# Non c-standard includes:
#   - "readExtrac.h"
#   o "defaultSettings.h"
#   o "cStrFun.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
#   o "FCIStatsFun.h"
#   o "minAlnStats.h"
#   o "samEntryStruct.h"
#   o "trimSam.h"
#   o "fqAndFaFun.h"
#   o "scoreReadsFun.h"
#   o "findCoInftBinTree.h"
#   o "findCoInftChecks.h"
#   o "fqGetIdsSearchFq"
#   o "fqGetIdsFqFun.h"
#   o "fqGetIdsStructs.h"
#   o "fqGetIdsHash.h"
#   o "fqGetIdsAVLTree.h"
# C standard Includes (all though non c-standard includes):
#   o <string.h>
#   o <stdlib.h>
#   o <sdtint.h>
#   o <stdio.h>
######################################################################*/

#include "readExtract.h" /*Holds functions to do read extraction*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start Of Program
'  o Main: function that runs everything
'  o fun-1: checkInput: check and process the user input (TO BE WRITTEN)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: 0 if no errors, pointer to argumet errored on for errors
|    Modifies: Every input varible to hold user input
\---------------------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,        /*Number of arugments user input*/
    char *argsCStr[],       /*Argumenas & parameters input*/
    char **fqFileCStr,      /*Will hold path to reads fastq file*/
    char **refFaFileCStr,   /*Reference file to use in selecting reads*/
    char **outFileCStr,     /*Will hold path of output file*/
    char *threadsCStr,      /*Number of threads to use with minimap2*/
    uint64_t *numReadsToExtUL, /*Number of reads to extract*/
    struct minAlnStats *minStats /*Holds mininum stats to keep a read*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-1 TOC: Sec-1 Sub-1: checkInput
   |  - Checks user input & puts input into variables for later use
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int main(int lenArgsInt, char *argsCStr[])
{ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Main TOC: main function
  '  o main sec-1: variable declerations
  '  o main sec-2: Check user input
  '  o main sec-3: Check if can open the primer fasta file
  '  o main sec-4: Check if can open fastq file
  '  o main sec-5: Check if can open output file
  '  o main sec-6: Run function to trim reads
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *fqFileCStr = 0;     /*file to open*/
    char *refFaFileCStr = 0;  /*File to hold input reference*/
    char *outFileCStr = 0;    /*file to write to*/
    char *inputChar = 0;      /*Holds arguemnt that had input error*/
    char threadsCStr[64];     /*Number of threads for minimap2*/
    char oneC = 1;            /*For passing 1 to functions*/
    char noRefBl = 0;           /*1: using reference; 0 I am not*/
    uint64_t numReadsToExtractUL = 300;
    uint64_t numReadsExtractedUL = 0;

    unsigned char errUC = 0; /*For error messages*/

    FILE *testFILE = 0;       /*Input file*/

    struct minAlnStats minStats;
    struct samEntry samST;       /*For file reading*/
    struct readBin binTree;      /*Holds fastq file*/

    char *helpMesgCStr = "\
        \n filterReads -fastq rads.fastq [options ...]\
        \n Use: Uses primer mappings to trim reads in a fastq file.\
        \n Input:\
        \n   -fastq file.fastq:                              [Required]\
        \n     o Fastq file to extract reads from.\
        \n   -ref                                            [None]\
        \n     o Fasta file with reference to map reads to.\
        \n     o Default is to use length and median Q-scores.\
        \n   -out:                                           [stdout]\
        \n     o output file to print filtered reads to\
        \n   -threads:                                       [3]\
        \n     o Number of threads to use with minimap2.\
        \n     o Only applied when -ref is selected.\
        \n   -num-reads:                                     [300]\
        \n     o Maximum number of reads to extract\
        \n   -min-mapq                                       [20]\
        \n     o min mapping quality to keep a read.\
        \n     o Only applied when -ref is selected.\
        \n   -min-median-q:                                  [10]\
        \n     o Remove reads with median Q-scores under\
        \n       this Q-score\
        \n     o Not applied when -ref is selected.\
        \n   -min-mean-q:                                    [10]\
        \n     o Remove reads with mean Q-scores under\
        \n       this Q-score.\
        \n     o Not applied when -ref is selected.\
        \n   -max-length:                                    [1000]\
        \n     o Remove reads over this length (0=all lengths)\
        \n     o Not applied when -ref is selected.\
        \n   -min-length:                                    [500]\
        \n     o Remove reads under this length.\
        \n     o Alinged length is used instead of read\
        \n       length for -ref.\
        \n   -min-q:                                         [10]\
        \n     o Minimum quality score to keep an SNP or\
        \n       insertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-snps:                                      [0.07=7%]\
        \n     o Maximum percent difference in SNPs\
        \n       between read and reference to keep read.\
        \n     o Only applied when -ref is selected.\
        \n   -max-diff:                                      [1=100%]\
        \n     o Maximum percent difference between read\
        \n       reference to keep read.\
        \n     o Only applied when -ref is selected.\
        \n   -max-inss:                                      [1=100%]\
        \n     o Maximum percent difference in insertions\
        \n       between read and reference to keep read.\
        \n     o Only applied when -ref is selected.\
        \n   -max-dels:                                      [1=100%]\
        \n     o Maximum percent difference in deletions\
        \n       between read and reference to keep read.\
        \n     o Only applied when -ref is selected.\
        \n   -max-indels:                                    [1=100%]\
        \n     o Maximum percent difference in indels\
        \n       between read and reference to keep read.\
        \n     o Only applied when -ref is selected.\
        \n   -max-ins-a-homo:                                [1]\
        \n     o Maximum hompolymer size an A insertion can\
        \n       be part of to keep the insertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-ins-t-homo:                                [1]\
        \n     o Maximum hompolymer size an T insertion can\
        \n       be part of to keep the insertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-ins-g-homo:                                [1]\
        \n     o Maximum hompolymer size an G insertion can\
        \n       be part of to keep the insertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-ins-c-homo:                                [1]\
        \n     o Maximum hompolymer size an C insertion can\
        \n       be part of to keep the insertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-del-a-homo:                                [1]\
        \n     o Maximum hompolymer size an A deletion can\
        \n       be part of to keep the delertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-del-t-homo:                                [1]\
        \n     o Maximum hompolymer size an T deletion can\
        \n       be part of to keep the delertion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-del-g-homo:                                [1]\
        \n     o Maximum hompolymer size an G deletion can\
        \n        be part of to keep the deletion.\
        \n     o Only applied when -ref is selected.\
        \n   -max-del-c-homo:                                [1]\
        \n     o Maximum hompolymer size an C deletion can\
        \n       be part of to keep the deletion.\
        \n     o Only applied when -ref is selected.\
        \n   -v:\
        \n     o Print version & exit\
        \n Output:\
        \n   stdout: prints out the kept reads to stdout\
         "; /*Help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Check user input
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankMinStats(&minStats); /*Set up default values (read/ref map)*/
    initSamEntry(&samST); /*Blank the sam entry struct (avoids errors)*/
    strcpy(threadsCStr, "3");

    inputChar =
        checkInput(
            &lenArgsInt,
            argsCStr,
            &fqFileCStr,
            &refFaFileCStr,
            &outFileCStr,
            threadsCStr,
            &numReadsToExtractUL,
            &minStats
    ); /*Get the user input*/

    if(inputChar != 0)
    { /*If I have an error or a non-run request*/
        if(strcmp(inputChar, "-h") == 0 ||
           strcmp(inputChar, "--h") == 0 ||
           strcmp(inputChar, "-help") == 0 ||
           strcmp(inputChar, "--help") == 0 ||
           strcmp(inputChar, "help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpMesgCStr);
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputChar, "-V") == 0 ||
           strcmp(inputChar, "-v") == 0 ||
           strcmp(inputChar, "--V") == 0 ||
           strcmp(inputChar, "--v") == 0 ||
           strcmp(inputChar, "--version") == 0 ||
           strcmp(inputChar, "--Version") == 0 ||
           strcmp(inputChar, "-version") == 0 ||
           strcmp(inputChar, "-Version") == 0
        ) { /*if the user wanted the version number*/
            fprintf(
                stdout,
                "extractTopReads from findCoInft version: %.8f\n",
                defVersion
            ); /*Print out the closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        else if(inputChar != 0)
        { /*If user had invalid input*/
            fprintf(
                stderr,
                "%s\n%s is invalid\n",
                helpMesgCStr,
                inputChar
            ); /*Print out the problem*/
            exit(1); /*Let user know their was an error*/
        } /*If user had invalid input*/
     } /*If I have an error or a non-run request*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-3: Check if can open fastq file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(fqFileCStr == 0)
    { /*If no file input*/
       fprintf(stderr, "No fastq file was input with -fastq\n");
       exit(-1);
    } /*If no file input*/

    testFILE = fopen(fqFileCStr, "r");

    if(testFILE == 0)
    { /*If no file was oppened*/
        fprintf(
            stderr,
            "Fastq file (-fastq %s) could not be opened\n",
            fqFileCStr
        ); /*If the filter file was invalid, let the user know*/

        exit(-1);
    } /*If no file was oppened*/

    fclose(testFILE);
    testFILE = 0;
    strcpy(binTree.fqPathCStr, fqFileCStr);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4: Check if can open reference file (if one was input)
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refFaFileCStr != 0)
    { /*If the user provide a reference file*/
        testFILE = fopen(refFaFileCStr, "r");

        if(testFILE == 0)
        { /*If no file was oppened*/
            fprintf(
                stderr,
                "Reference file (-ref %s) could not be opened\n",
                refFaFileCStr
            ); /*If the filter file was invalid, let the user know*/

            exit(-1);
        } /*If no file was oppened*/

        fclose(testFILE);
        testFILE = 0;
        strcpy(binTree.bestReadCStr, refFaFileCStr);
    } /*If the user provide a reference file*/
 
     else noRefBl = 1; /*Mark that I am not using the reference*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Check if can open output file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(outFileCStr != 0)
    { /*If taking input from file*/
        strcpy(binTree.topReadsCStr, outFileCStr);

        testFILE = fopen(outFileCStr, "w");    /*Re-using the out file*/
    
        if(testFILE == 0)
        { /*If no file was oppened*/
            fprintf(
                stderr,
                "Could not open output file (-out %s)\n",
                outFileCStr
            ); /*If the filter file was invalid, let the user know*/

            exit(-1);
        } /*If no file was oppened*/

        fclose(testFILE);
        testFILE = 0;
    } /*If taking input from file*/

    else binTree.topReadsCStr[0] = '\0'; /*Flag for stdout*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Run function to trim reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errUC = 
        findBestXReads(
            &numReadsToExtractUL,
            &numReadsExtractedUL,
            threadsCStr,
            &oneC,
            &minStats,
            &samST,
            0,  /*Do not use reference in scoring (needs to be fastq)*/
            &binTree,     /*Bin working on*/
            noRefBl,      /*Only use fastq file*/
            0             /*Use input name*/
    ); /*Extract the top reads that mapped to the selected best read*/

    freeStackSamEntry(&samST); /*Free the buffer in samST*/

    if(errUC & 16)
    { /*If had a memory allocation error*/
        fprintf(stderr, "Something happened (minimap2 errored out?)\n");
        exit(-1);
    } /*If had a memory allocation error*/

    fprintf(stderr, "Extracted %ju reads\n", numReadsExtractedUL);
    exit(0);
} /*main function*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: 0 if no errors, pointer to argumet errored on for errors
|    Modifies: Every input varible to hold user input
\---------------------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,             /*Number of arugments user input*/
    char *argsCStr[],            /*Argumenas & parameters input*/
    char **fqFileCStr,           /*Will hold path to reads fastq file*/
    char **refFaFileCStr,   /*Reference file to use in selecting reads*/
    char **outFileCStr,          /*Will hold path of output file*/
    char *threadsCStr,      /*Number of threads to use with minimap2*/
    uint64_t *numReadsToExtUL, /*Number of reads to extract*/
    struct minAlnStats *minStats /*Holds mininum stats to keep a read*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-1 TOC: Sec-1 Sub-1: checkInput
   |  - Checks user input & puts input into variables for later use
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *parmCStr = 0;
    char *singleArgCStr = 0;

    if(*lenArgsInt < 2)
        return parmCStr; /*no arguments input*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
    { /*loop through all user input arguments*/  /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1); /*supplied argument*/
        parmCStr = *(argsCStr + intArg);          /*Paramter*/

        if(strcmp(parmCStr, "-fastq") == 0)
            *fqFileCStr = singleArgCStr;

        else if(strcmp(parmCStr, "-ref") == 0)
            *refFaFileCStr = singleArgCStr;

        else if(strcmp(parmCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        else if(strcmp(parmCStr, "-threads") == 0)
            strcpy(threadsCStr, singleArgCStr);

        else if(strcmp(parmCStr, "-num-reads") == 0)
            *numReadsToExtUL = strtoul(singleArgCStr, &parmCStr, 10);

        /*General filtering settings*/
        else if(strcmp(parmCStr, "-min-median-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMedianQFlt);

        else if(strcmp(parmCStr, "-min-mean-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMeanQFlt);
            
       else if(strcmp(parmCStr, "-min-length") == 0)
           minStats->minReadLenULng=strtoul(singleArgCStr,&parmCStr,10);

       else if(strcmp(parmCStr, "-max-length") == 0)
           minStats->maxReadLenULng=strtoul(singleArgCStr,&parmCStr,10);

       else if(strcmp(parmCStr, "-min-q") == 0)
           cStrToUChar(singleArgCStr, &minStats->minQChar);

       else if(strcmp(parmCStr, "-min-mapq") == 0)
           cStrToUInt(singleArgCStr, &minStats->minMapqUInt);

       /*Percent difference settings*/
       else if(strcmp(parmCStr, "-min-snps") == 0)
           sscanf(singleArgCStr, "%f", &minStats->minSNPsFlt);

       else if(strcmp(parmCStr, "-min-diff") == 0)
           sscanf(singleArgCStr, "%f", &minStats->minDiffFlt);

       else if(strcmp(parmCStr, "-min-inss") == 0)
           sscanf(singleArgCStr, "%f", &minStats->minInssFlt);

       else if(strcmp(parmCStr, "-min-dels") == 0)
           sscanf(singleArgCStr, "%f", &minStats->minDelsFlt);

       else if(strcmp(parmCStr, "-min-indels") == 0)
           sscanf(singleArgCStr, "%f", &minStats->minIndelsFlt);

       /*Homopolymer settings (insertions)*/
       else if(strcmp(parmCStr, "-max-ins-a-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[0]);
       else if(strcmp(parmCStr, "-max-ins-t-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[10]);
       else if(strcmp(parmCStr, "-max-ins-c-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[1]);
       else if(strcmp(parmCStr, "-max-ins-g-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoInsAry[3]);

       /*Homopolymer settings (deletions)*/
       else if(strcmp(parmCStr, "-max-del-a-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[0]);
       else if(strcmp(parmCStr, "-max-del-t-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[10]);
       else if(strcmp(parmCStr, "-max-del-c-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[1]);
       else if(strcmp(parmCStr, "-max-del-g-homo") == 0)
           cStrToUInt(singleArgCStr, &minStats->maxHomoDelAry[3]);

       else
           return parmCStr;

        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
