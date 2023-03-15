/*######################################################################
# Name: trimPrimers
# Use:
#   o Uses primer coordinatest to trims reads in a fastq file
# Input:
#    -primers primers.fasta:                                 [Required]
#      o Fasta file with primers to map & trime reads with
#    -fastq file.fastq:                                      [Required]
#      o Fastq file to filter reads from
#    -no-hash:
#      o Do search with only the AVL tree (no hashing)       [Hashing]
#      o Takes lonber, but uses slightly (~10%) less memory. 
#    -v:
#      o Print version & exit
# Output:
#    stdout: prints out the trimmed reads to stdout
# Includes:
#   - "trimPrimersSeach.h"
#   o "trimPrimersHash.h"
#   o "fqGetIdsHash.h"
#   o "trimPrimersAVLTree.h"
#   o "defaultSettings.h"
#   o "cStrFun.h"
#   o "fqGetIdsFqFun.h"
#   o "trimPrimersStructs.h"
#   o "fqGetIdsStructs.h"
#   o "cStrToNumberFun.h"
#   o "fqAndFaFun.h"
#   o "FCIStatsFun.h"     (fqAndFqFun.h)
#   o "minAlnStats.h"     (fqAndFaFun.h)
#   o "samEntryStruct.h"  (fqAndFaFun.h->FCIstatsFun.h)
#   o "printError.h"      (fqAndFaFun.h)
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#include "trimPrimersSearch.h" /*Holds functions to do read extraction*/

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
    int *lenArgsInt,      /*Number of arugments user input*/
    char *argsCStr[],     /*Argumenas & parameters input*/
    char **primFileCStr,  /*Will hold path to primer fasta file*/
    char **pafFileCStr,   /*Holds paf file with read/primer mappings*/
    char *stdinPafBl,     /*Sets to 1 if taking paf file from stdin*/
    char **fqFileCStr, /*Will hold path to reads fastq file*/
    char **outFileCStr,   /*Will hold path of output file*/
    char *hashBl,    /*Set to 0 if user wants tree search*/
    char *threadsCStr /*Number of threads to use with minimap2*/
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

    char *primFileCStr = 0;/*file to open*/
    char *pafFileCStr = 0; /*Points to paf file to use for mappings*/
    char stdinPafBl = 0;   /*1: take paf file from stdin; 0 do not*/
    char *fqFileCStr = 0;  /*file to open*/
    char *outFileCStr = 0; /*file to write to*/
    char hashBl = 1;       /*Holds if user wanted hashing [1: use hash]*/
    char *inputChar = 0;   /*Holds arguemnt that had input error*/
    char threadsCStr[128]; /*Holds number of threads to use*/

    unsigned char errUC = 0; /*For error messages*/

    FILE *testFILE = 0;       /*Input file*/

    char *helpMesgCStr = "\
        \n trimPrimers -primers primers.fasta -fastq reads.fastq [...]\
        \n Use: Uses primer mappings to trim reads in a fastq file.\
        \n Input:\
        \n   -primers primers.fasta:                         [Required]\
        \n     o Fasta file with primers to map & trime reads with\
        \n   -paf                                            [None]\
        \n     o Paf file with reads mapped to primers. This\
        \n       replaces -primers and skips minimap2.\
        \n   -stdin-paf:                                     [No]\
        \n     o Take the paf file (-paf) from stdin.\
        \n   -fastq file.fastq:                              [Required]\
        \n     o Fastq file to filter reads from\
        \n   -no-hash:                                       [Hashing]\
        \n     o Do search with only the AVL tree (no hashing)\
        \n     o Takes lonber, but uses slightly (~10%) less memory.\
        \n   -v:\
        \n     o Print version & exit\
        \n   -out:                                           [stdout]\
        \n     o Name of file to output reads to\
        \n Output:\
        \n   stdout: prints out the kept reads to stdout\
         "; /*Help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Check user input
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    strcpy(threadsCStr, defThreads);

    inputChar = checkInput(&lenArgsInt,
                           argsCStr,
                           &primFileCStr,
                           &pafFileCStr,
                           &stdinPafBl,
                           &fqFileCStr,
                           &outFileCStr,
                           &hashBl,
                           threadsCStr
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
                "trimPrimers from findCoInft version: %.8f\n",
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
    ^ Main Sec-3: Check if can open the filter file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(!(stdinPafBl & 1) && pafFileCStr == 0)
    { /*If have to do mapping with minimap2*/
        if(primFileCStr != 0) testFILE = fopen(primFileCStr, "r");

        else
        { /*Else no file input*/
            fprintf(stderr, "No primer file was input with -primers\n");
            exit(-1);
        } /*Else no file input*/

        if(testFILE == 0)
        { /*If no file was oppened*/
            fprintf(
                 stderr,
                "Input primer fasta file (%s) could not be opened\n",
                primFileCStr
            ); /*If the filter file was invalid, let the user know*/

            exit(-1);
        } /*If no file was oppened*/

        fclose(testFILE);
        testFILE = 0;
    } /*If have to do mapping with minimap2*/

    else
    { /*Else user is providing the mappig files*/
        if(pafFileCStr != 0)
        { /*If user provided a paf file*/
            testFILE = fopen(pafFileCStr, "r");

            if(testFILE == 0)
            { /*If an invalid file was provided*/
                fprintf(
                    stderr,
                    "paf file (-paf %s) does not exist\n",
                     pafFileCStr
                ); /*Let user know the paf file did not exist*/

                exit(-1);
            } /*If an invalid file was provided*/

            fclose(testFILE);
        } /*If user provided a paf file*/

        testFILE = 0;
    } /*Else user is providing the mappig files*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4: Check if can open fastq file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(fqFileCStr != 0) testFILE = fopen(fqFileCStr, "r");

    else
    { /*Else no file input*/
        fprintf(stderr, "No fastq file was provided with -fastq\n");
        exit(-1);
    } /*Else no file input*/

    if(testFILE == 0)
    { /*If no file was oppened*/
        fprintf(
            stderr,
            "Input fastq file (%s) could not be opened\n",
            fqFileCStr
        ); /*If the filter file was invalid, let the user know*/

        exit(-1);
    } /*If no file was oppened*/

    fclose(testFILE);
    testFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Check if can open output file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(outFileCStr != 0)
    { /*If taking input from file*/
        testFILE = fopen(outFileCStr, "w");    /*Re-using the out file*/
    
        if(testFILE == 0)
        { /*If no file was oppened*/
            fprintf(
                stderr,
                "Could not make the output file (%s)\n",
                outFileCStr
            ); /*If the filter file was invalid, let the user know*/

            exit(-1);
        } /*If no file was oppened*/

        fclose(testFILE);
        testFILE = 0;
    } /*If taking input from file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Run function to trim reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Trim the reads using the input primers*/
    errUC =
      trimPrimers(
          primFileCStr,
          pafFileCStr,
          stdinPafBl,
          fqFileCStr,
          outFileCStr,
          threadsCStr,
          hashBl
    );

    if(errUC & 64)
    { /*If had a memory allocation error*/
        fprintf(stderr, "Memory error (ran out of memory)\n");
        exit(-1);
    } /*If had a memory allocation error*/

    if(errUC & 32)
    { /*If had a memory allocation error*/
        fprintf(stderr, "%s is not a fastq file\n", fqFileCStr);
        exit(-1);
    } /*If had a memory allocation error*/

    exit(0);
} /*main function*/

/*---------------------------------------------------------------------\
| Output:
|    Returns: 0 if no errors, pointer to argumet errored on for errors
|    Modifies: Every input varible to hold user input
\---------------------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,      /*Number of arugments user input*/
    char *argsCStr[],     /*Argumenas & parameters input*/
    char **primFileCStr,  /*Will hold path to primer fasta file*/
    char **pafFileCStr,   /*Holds paf file with read/primer mappings*/
    char *stdinPafBl,     /*Sets to 1 if taking paf file from stdin*/
    char **fqFileCStr, /*Will hold path to reads fastq file*/
    char **outFileCStr,   /*Will hold path of output file*/
    char *hashBl,    /*Set to 0 if user wants tree search*/
    char *threadsCStr /*Number of threads to use with minimap2*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   | Fun-1 TOC: Sec-1 Sub-1: checkInput
   |  - Checks user input & puts input into variables for later use
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *tmpCStr = 0, *singleArgCStr = 0;

    if(*lenArgsInt < 2)
        return tmpCStr; /*no arguments input*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
    { /*loop through all user input arguments*/  /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1); /*supplied argument*/
        tmpCStr = *(argsCStr + intArg);          /*Paramter*/

        if(strcmp(tmpCStr, "-primers") == 0)
            *primFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-paf") == 0)
            *pafFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-stdin-paf") == 0)
        { /*Else if taking the paf file from stdin*/
            *stdinPafBl = 1;
            --intArg;
        } /*Else if taking the paf file from stdin*/

        else if(strcmp(tmpCStr, "-fastq") == 0)
            *fqFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-threads") == 0)
            strcpy(threadsCStr, singleArgCStr);

        else if(strcmp(tmpCStr, "-no-hash") == 0)
        { /*If user wants to do tree search instead*/
            *hashBl = 0;
            --intArg;         /*Account for incurment at end of loop*/
        } /*If user wants to do tree search instead*/

        else
            return tmpCStr;

        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
