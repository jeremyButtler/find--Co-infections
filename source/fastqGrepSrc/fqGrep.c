/*##############################################################################
# Name: clustGraph.c
# Use: Graps selected reads from a fastq file
# Input:
#    -f file.txt:
#      - File with reads to copy from a fastq file           [Required]
#      - Each read name should be on a separate line
#    -stdin-fastq:
#      - Take fastq input from stdin                         [Default: not set]
#      - Do not use with stdin-filt (will break)
#    -stdin-filt:
#      - Take filter file from stdin                         [Default: not set]
#      - Do not use with stdin-fastq (will break)
#    -fastq file.fastq:
#      - Fastq file to filter reads from                     [Required]
#      -no-hash:
#          - Use a tree search instead of hashing.           [Default: hashing]
#              - Default search is hash combined with tree. 
#          - Adds more time, but uses slightly less memory. 
#    -v:
#      - Print reads not provided by -f                      [Default: not set]
#    -V:
#      - Print version & exit
# Output:
#    stdout: prints out the kept reads to stdout
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    fun-1: main: function that runs everything
#    fun-2: checkInput: check and process the user input (TO BE WRITTEN)
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fqGrepSearchFq.h" /*Holds functions to do read extraction*/
    /*
      Includes:
          - "fastqGrepFastqFun.h"    # Functoions to read fastq file
              - <stdio.h>            # FILE, fread, printf, & fprintf 
          - "fastqGrepHash.h"        # functions to build my hash table
              - "fastqGrepAVLTree.h" # functions to build AVL tree
                  - <string.h>       # strcmp()
                  - "fastqGrepStructs.h" # Structers used in this code
                      - <stdlib.h>  # Malloc & free
                      - <stdio.h>   # fprintf
    */

/*##############################################################################
# Output:
#    Returns: 0 if no errors, pointer to argumet errored on for errors
#    Modifies: filtFileCStr, fastqFileCStr, stdinChar, & flipChar to hold or
#              point to user input
##############################################################################*/
char * checkInput(int *lenArgsInt,        /*Number of arugments user input*/
                  char *argsCStr[],       /*Argumenas & parameters input*/
                  char **filtFileCStr,    /*Holds path/name of filter file*/
                  char **fastqFileCStr,   /*Holds path/name of fastq file*/
                  char *stdinFastqChar,   /*1: stdin input for fastq*/
                  char *stdinFiltChar,    /*1: stdin input for filter*/
                  char *useHashChar,      /*Set to 0 if user wants tree search*/
                  char *flipChar          /*1: print reads not in filtFileCStr
                                            0: print reads in filtFileCStr*/
); /*Checks user input & puts input into variables for later use*/

int main(int lenArgsInt, char *argsCStr[])
{ /*main function*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main TOC: main function
    #    main sec-1: variable declerations
    #    main sec-2: Check user input and index file
    #    main sec-3: Check if can open the filter file
    #    main sec-4: Check if can open fastq file
    #    main sec-5: Run function to extract reads by id and exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long
        buffSizeULng = 1 << 19;  /*16^4-1, buffer size for file reading*/

    char
        *filtFileCStr = 0,    /*file to open*/
        *fastqFileCStr = 0,   /*file to open*/
        stdinFastqChar = 0,   /*If 1 taking input from stdin*/
        stdinFiltChar = 0,    /*If 1 taking input from stdin*/
        printReverseChar = 0, /*Print sequences not in filter file*/
        useHashChar = 1,      /*Holds if user wanted hashing [1: use hash]*/
        *inputChar = 0;       /*Holds arguemnt that had input error*/

    unsigned char
        sizeReadStackUChar = 200;

    FILE
        *fastqFile = 0,       /*Input file*/
        *filtFile = 0;        /*File to hold the filter*/

    char
        *helpMesgCStr = "grepFastq -f filter.txt -fastq file.fastq [options...]\
            \n Use: grabs selected reads from a fastq file\
            \n Input:\
            \n    -f file.txt:\
            \n      - File with read names to copy from fastq file   [Required]\
            \n      - Each read name should be on a separate line\
            \n    -stdin-fastq:\
            \n      - Take fastq input from stdin            [Default: not set]\
            \n    -stdin-filt:\
            \n      - Take filter file from stdin            [Default: not set]\
            \n    -fastq file.fastq:\
            \n      - Fastq file to filter reads from        [Required]\
            \n    -no-hash:\
            \n      - Use a tree search instead of hashing.  [Default: hashing]\
            \n          - Default search is hash combined with tree.\
            \n      - Adds more time, but uses slightly less memory.\
            \n    -v:\
            \n      - Print reads not provided by -f         [Default: not set]\
            \n    -V:\
            \n      - Print version and exit\
            \n Output:\
            \n   stdout: prints out the kept reads to stdout\
         "; /*Help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-2: Check user input and index file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    inputChar = checkInput(&lenArgsInt,
                           argsCStr,
                           &filtFileCStr,
                           &fastqFileCStr,
                           &stdinFastqChar,
                           &stdinFiltChar,
                           &useHashChar,
                           &printReverseChar
    ); /*Get the user input*/

    if(
       inputChar != 0 &&
       strcmp(inputChar, "-h") == 0
    ) { /*If user wanted the help message*/
        fprintf(
            stdout,       /*Print to stdout so user can redirect easily*/
            "%s\n",       /*Print out help messge & new line*/
            helpMesgCStr  /*Holds help message*/
        ); /*Print out help message*/

        exit(0);
    } /*If user wanted the help message*/

    else if(
        inputChar != 0 &&
        strcmp(
            inputChar, /*Paremeter checkInput exited on*/
            "-V"       /*Input for version number*/
        ) == 0         /*String was a match*/
    ) { /*Else if the user wanted the version number*/
        fprintf(
            stdout,                  /*stdout so user can pipe & grab easily*/
            "fastqGrep version: 1\n" /*Version number*/
        ); /*Print out the version number*/
        exit(0);
    } /*Else if the user wanted the version number*/

    else if(inputChar != 0)
    { /*If user had invalid input*/
        fprintf(stderr, "%s\n%s is invalid\n", helpMesgCStr, inputChar);
        exit(-1);
    } /*If user had invalid input*/

    if(stdinFastqChar == 1 && stdinFiltChar == 1)
    { /*If user takeing all input from stdin*/
        fprintf(stderr, "both fastq file & filter file are from stdin\n");
        exit(-1);
    } /*If user takeing all input from stdin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Check if can open the filter file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if filter input comming from stdin or file*/
    if(stdinFiltChar == 0)
        filtFile = fopen(filtFileCStr, "r");
    else
        filtFile = stdin;

    if(filtFile == 0)
    { /*If no file was oppened*/
        fprintf(
            stderr,
            "Input filter file (%s) could not be opened\n",
            filtFileCStr
        ); /*If the filter file was invalid, let the user know*/

        exit(-1);
    } /*If no file was oppened*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-4: Check if can open fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(stdinFastqChar == 0)
    { /*If taking input from file*/
        fastqFile = fopen(fastqFileCStr, "r");      /*Re-using the fastq file*/
    
        if(fastqFile == 0)
        { /*If no file was oppened*/
            fprintf(
                stderr,
                "Input fastq file (%s) could not be opened\n",
                fastqFileCStr
            ); /*If the filter file was invalid, let the user know*/

            exit(-1);
        } /*If no file was oppened*/
    } /*If taking input from file*/

    else
        fastqFile = stdin;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-5: Run function to extract reads by id and exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(
        fastqExtract(
            filtFile,           /*File with read targets to keep or ignore*/
            fastqFile,          /*Fastq File with reads to extract*/
            sizeReadStackUChar, /*Number of elements to use in stack*/
            buffSizeULng,       /*Size of buffer to read input with*/
            useHashChar,        /*Tells if doing search with hashing*/
            printReverseChar /*Tells if extracting filter ids (1) or other (0)*/
        ) == 0
    ) { /*If fastq file was not a valid fastq file, issue a warning*/
        fprintf(
            stderr,
            "fastq file (%s) supplied with -fastq is not a fastq file\n",
            fastqFileCStr
        ); /*Warn user fastq is not formatted corectly*/
    } /*If fastq file was not a valid fastq file, issue a warning*/

    exit(0);
} /*main function*/

/*##############################################################################
# Output:
#    Returns: 0 if no errors, pointer to argumet errored on for errors
#    Modifies: filtFileCStr, fastqFileCStr, stdinChar, & flipChar to hold or
#              point to user input.
#    Sets: useHashChar to 0 if user did not want to using hashing
##############################################################################*/
char * checkInput(int *lenArgsInt,        /*Number of arugments user input*/
                  char *argsCStr[],       /*Argumenas & parameters input*/
                  char **filtFileCStr,    /*Holds path/name of filter file*/
                  char **fastqFileCStr,   /*Holds path/name of fastq file*/
                  char *stdinFastqChar,   /*1: stdin input for fastq*/
                  char *stdinFiltChar,    /*1: stdin input for filter*/
                  char *useHashChar,      /*Set to 0 if user wants tree search*/
                  char *flipChar          /*1: print reads not in filtFileCStr
                                            0: print reads in filtFileCStr*/
) /*Checks user input & puts input into variables for later use*/
{ /*checkInput*/
    char *tmpCStr = 0, *singleArgCStr = 0;

    if(*lenArgsInt < 2)
        return tmpCStr; /*no arguments input*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
    { /*loop through all user input arguments*/           /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1);          /*supplied argument*/
        tmpCStr = *(argsCStr + intArg);                   /*Paramter*/

        if(strcmp(tmpCStr, "-f") == 0)
            *filtFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-fastq") == 0)
            *fastqFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-no-hash") == 0)
        { /*If user wants to do tree search instead*/
            *useHashChar = 0;
            intArg--;                  /*Account for incurment at end of loop*/
        } /*If user wants to do tree search instead*/

        else if(strcmp(tmpCStr, "-stdin-fastq") == 0)
        { /*If if taking input from stdin*/
            *stdinFastqChar = 1;
            intArg--;                  /*Account for incurment at end of loop*/
        } /*If taking input from stdin*/

        else if(strcmp(tmpCStr, "-stdin-filt") == 0)
        { /*If if taking input from stdin*/
            *stdinFiltChar = 1;
            intArg--;                  /*Account for incurment at end of loop*/
        } /*If if taking input from stdin*/

        else if(strcmp(tmpCStr, "-v") == 0)
        { /*If if taking input from stdin*/
            *flipChar = 1;
            intArg--;                  /*Account for incurment at end of loop*/
        } /*If taking input from stdin*/

        else if(strcmp(tmpCStr, "-V") == 0)
            return tmpCStr; /*Main function handles version number*/

        else if(strcmp(tmpCStr, "-h") == 0)
            return tmpCStr; /*Main function handles help message*/

        else
            return tmpCStr;
        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    if(*fastqFileCStr == 0)
        return 0;            /*Program will catch latter*/

    tmpCStr = *fastqFileCStr;

    while(*tmpCStr != '\0')
        tmpCStr++;            /*find the end of the file name*/

    if(
        *(tmpCStr - 1) != 'q' ||
        *(tmpCStr - 2) != 't' ||
        *(tmpCStr - 3) != 's' ||
        *(tmpCStr - 4) != 'a' ||
        *(tmpCStr - 5) != 'f' ||
        *(tmpCStr - 6) != '.'
    ) /*If not a fastq file*/
    { /*If input is not a fastq file*/
        fprintf(
            stderr,
            "%s is not a fastq file (should end in .fastq)\n",
            *fastqFileCStr
        ); /*Warn user that provided fastq file is not a fastq file*/

        return tmpCStr;
    } /*If input is not a fastq file*/

    return 0; /*input is valid*/
} /*checkInput*/
