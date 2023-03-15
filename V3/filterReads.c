/*######################################################################
# Name: filterReads
# Use:
#   o filters reads using length, median Q-score, and mean Q-score
# Input:
#    -fastq file.fastq:                                      [Required]
#      o Fastq file to filter reads from
#    -stdin:                                                 [No]
#      o Take input fastq file from stdin
#    -min-median-q:                                          [10]
#      o Remove reads with median Q-scores under
#        this Q-score
#    -min-mean-q:                                            [10]
#      o Remove reads with mean Q-scores under
#        this Q-score
#    -max-length:                                            [1000]
#      o Remove reads over this length (0 = all lengths)
#    -min-length:                                            [500]
#      o Remove reads under this length
#    -out:                                                   [stdout]
#      o output file to print filtered reads to
#    -v:
#      o Print version & exit
# Output:
#    stdout: prints out the kept reads to stdout
# Non c-standard includes:
#   - "fqAndFaFun.h"
#   o "defaultSettings.h"
#   o "minAlnStatsStruct.h"
#   o "FCIStatsFun.h"
#   o "samEntryStruct.h"
#   o "cStrToNumberFun.h"
#   o "printError.h"
# C standard Includes (all though non c-standard includes):
#   - <string.h>
#   o <stdlib.h>
#   o <sdtint.h>
######################################################################*/

#include <string.h>
#include "fqAndFaFun.h" /*Holds functions to do read extraction*/

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
    int *lenArgsInt,             /*Number of arugments user input*/
    char *argsCStr[],            /*Argumenas & parameters input*/
    char **fqFileCStr,           /*Will hold path to reads fastq file*/
    char *stdinBl,          /*Sets to 1 if taking paf file from stdin*/
    char **outFileCStr,          /*Will hold path of output file*/
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

    char stdinBl = 0;   /*1: take paf file from stdin; 0 do not*/
    char *fqFileCStr = 0;  /*file to open*/
    char *outFileCStr = 0; /*file to write to*/
    char *inputChar = 0;   /*Holds arguemnt that had input error*/

    unsigned char errUC = 0; /*For error messages*/

    FILE *testFILE = 0;       /*Input file*/

    struct minAlnStats minStats;
    struct samEntry samST;       /*For file reading*/

    char *helpMesgCStr = "\
        \n filterReads -fastq rads.fastq [options ...]\
        \n Use: Uses primer mappings to trim reads in a fastq file.\
        \n Input:\
        \n   -fastq file.fastq:                              [Required]\
        \n     o Fastq file to filter reads from\
        \n   -stdin:                                         [No]\
        \n     o Take input fastq file from stdin.\
        \n   -min-median-q:                                  [10]\
        \n     o Remove reads with median Q-scores under\
        \n       this Q-score.\
        \n   -min-mean-q:                                    [10]\
        \n     o Remove reads with mean Q-scores under\
        \n       this Q-score.\
        \n   -max-length:                                    [1000]\
        \n     o Remove reads over this length (0 = all lengths).\
        \n   -min-length:                                    [500]\
        \n     o Remove reads under this length.\
        \n   -out:                                           [stdout]\
        \n     o Name of file to output reads to\
        \n   -v:\
        \n     o Print version & exit\
        \n Output:\
        \n   stdout: prints out the kept reads to stdout\
         "; /*Help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Check user input
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankMinStats(&minStats); /*Set up default values*/
    initSamEntry(&samST); /*Blank the sam entry struct (avoids errors)*/

    inputChar =
        checkInput(
            &lenArgsInt,
            argsCStr,
            &fqFileCStr,
            &stdinBl,
            &outFileCStr,
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
                "filterReads from findCoInft version: %.8f\n",
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

    if(!(stdinBl & 1))
    { /*If I need to check the fastq file*/
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
    } /*If I need to check the fastq file*/

    else
        fqFileCStr = 0; /*So stdin input fired in filterReads()*/

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
                "Could not open output file (-out %s)\n",
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

    /*Filter reads*/
    errUC = filterReads(fqFileCStr, outFileCStr, &samST, &minStats);

    freeStackSamEntry(&samST); /*Free the buffer in samST*/

    if(errUC & 64)
    { /*If had a memory allocation error*/
        fprintf(stderr, "Memory error (ran out of memory)\n");
        exit(-1);
    } /*If had a memory allocation error*/

    if(errUC & 2)
    { /*If had a memory allocation error*/
        fprintf(stderr, "-fastq %s is not a fastq file\n", fqFileCStr);
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
    int *lenArgsInt,             /*Number of arugments user input*/
    char *argsCStr[],            /*Argumenas & parameters input*/
    char **fqFileCStr,           /*Will hold path to reads fastq file*/
    char *stdinBl,          /*Sets to 1 if taking paf file from stdin*/
    char **outFileCStr,          /*Will hold path of output file*/
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

        if(strcmp(parmCStr, "-stdin") == 0)
        { /*Else if taking the paf file from stdin*/
            *stdinBl = 1;
            --intArg;
        } /*Else if taking the paf file from stdin*/

        else if(strcmp(parmCStr, "-fastq") == 0)
            *fqFileCStr = singleArgCStr;

        else if(strcmp(parmCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        else if(strcmp(parmCStr, "-min-median-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMedianQFlt);

        else if(strcmp(parmCStr, "-min-mean-q") == 0)
            sscanf(singleArgCStr, "%f", &minStats->minMeanQFlt);
            
       else if(strcmp(parmCStr, "-min-length") == 0)
           minStats->minReadLenULng=strtoul(singleArgCStr,&parmCStr,10);

         else if(strcmp(parmCStr, "-max-length") == 0)
            minStats->maxReadLenULng =
                strtoul(singleArgCStr, &parmCStr, 10);

        else
            return parmCStr;

        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
