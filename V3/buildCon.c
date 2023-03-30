/*######################################################################
# Use:
#   o buildCon is the driver file for the buildCon program.
#   o It will buid a consensus from a fastq file with or without a 
#     reference.
# Requires:
#   o cStrFun.c/h
#   o cStrToNumberFun.c/h
#   o defaultSettings.h
#   o printError.c/h
#   o samEntryStruct.c/h
#   o trimSam.c/h
#   o fqGetIdsStructs.c/h
#   o fqGetIdsAVLTree.c/h
#   o fqGetIdsHash.c/h
#   o fqGetIdsFqFun.c/h
#   o fqGetIdsSearchFq.c/h
#   buildConFun.c/h
# C libaries:
#   o string.h
#   o stdlib.h
#   o stdio.h
#   o stdint.h
######################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' HOP: Head of program
'   header:
'     o has header files and includes
'   main:
'     o The main function
'   fun-1 getUserInput:
'     o processes usser input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
^ Header:
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "buildConFun.h"

/*---------------------------------------------------------------------\
| Output:
|    - Modifies input variables to hold user input
|    - Returns:
|        - 0: if nothing went wrong
|        - pionter to invalid paramter
|            - if no agruments input, returns ponter to argsCStr
\---------------------------------------------------------------------*/
char * getUserInput(
    int32_t lenArgsInt,
    char *argsCStr[],
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refPathCStr, /*Holds path to references*/
    char **statsPathCStr, /*Holds path to scoreReads output to use*/
    char *prefixCStr,   /*Prefix to name everything*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    struct conBuildStruct *conSet, 
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *minReadConStats  /*Read cluster socring set*/
); /*Reads in user input*/

int main(
    int32_t lenArgsInt, /*Number of parameters & arguments user input*/
    char *argsCStr[]  /*List of parameters & arguements user input*/
) /*Function to run everything*/
{ /*main*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Main TOC:
    '   main sec-1: Variable declerations
    '   main sec-2: Set up default settings in structures
    '   main sec-3: Get user input & set up commands
    '   main sec-4: Check user input
    '   main sec-5: Build the consensus
    '   main sec-6: Clean up, Print out version numbers, and exit
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-1: Variable declerations
    ^    main sec-1 sub-1: General variables
    ^    main sec-1 sub-2: Help message
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-1 Sub-1: General variables
    \******************************************************************/

    /*User Input or related variables (non file names)*/
    char *fqPathCStr = 0;
    char *refPathCStr = 0;
    char *statsPathCStr = 0;      /*Stats file to use in read seletion*/
    char prefixCStr[64] = defPrefix;
    char tmpPathCStr[128];               /*copy of Fastq file working on*/
    char threadsCStr[16] = defThreads; /*# of threads for system calls*/

    /*Variables for the majority consensus step*/
    struct conBuildStruct conSetting;

    /*C-strings that hold commands*/
    char tmpCmdCStr[1024];      /*Holds a quick system command*/
    char minimap2VersionCStr[256];
    char raconVersionCStr[256];
    char medakaVersionCStr[2048]; /*For gpu error messages*/
    char *tmpCStr = 0;         /*For string manipulation*/

    /*Miscalanious variables*/
    unsigned char errUC = 0;         /*Holds error messages*/

    /*Holds thresholds to keep alignments, errors, & matches*/
    struct samEntry samStruct; /*For reading files*/
    struct samEntry refStruct; /*For reading files*/

    struct minAlnStats readToReadMinStats; /*for read to consensus map*/
    struct minAlnStats minReadConStats; /*for read to consensus map*/

    struct readBin fastqStruct;         /*Holds my fastq file*/

    FILE *stdinFILE = 0;
    FILE *cpFILE = 0;  /*For making a copy of the input fastq file*/
        /*This avoids modifications of the input fastq file*/

    /******************************************************************\
    * Main Sec-1 Sub-2: Help message
    \******************************************************************/

    char *helpCStr = "\
            \n buildCon -fastq reads.fastq [Options ...]\
            \n Use: Builds a consensus from a fastq file.\
            \n    -fastq:\
            \n        - Fastq file with reads to search      [Required]\
            \n    -ref:\
            \n        - Reference used to build consensus   [Best read]\
            \n    -prefix:                                  [out]\
            \n        - Prefix to name output file\
            \n    -stats:                                    [Not used]\
            \n        - tsv file output by score reads to use\
            \n          to select the best read by mapping\
            \n          quality.\
            \n        - Default is to use the read with the\
            \n          best median Q-score.\
            \n    -threads:\
            \n        - Number of threads to use             [3]\
            \n    -min-reads-per-bin:\
            \n        - Min number of reads needed to keep   [100]\
            \n          a bin or a cluster\
            \n    -max-reads-per-con:\
            \n        - Max numver of reads to use in        [300]\
            \n          a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
            \n    -min-con-length:                               [500]\
            \n       - Discard consensuses that are under the\
            \n         input length.\
            \n       - If you lower this your should also lower\
            \n         -min-read-read-map-length &\
            \n         -min-read-con-map-length\
            \n Additional Help messages:\
            \n    -h-build-consensus:\
            \n        - Print paramaters for building the consensus\
            \n    -h-read-con-map:\
            \n        - Print out the parameters for the read to\
            \n          consensus mapping step.\
            \n    -h-read-read-map:\
            \n        - Print out the parameters for the read to read\
            \n          mapping step.\
            \n Output:\
            \n    - File:\
            \n      o Contents: Consensus built from the fasta file\
            \n      o Name: fastqFileName--cluster-0--consensus.fasta\
            \n      o Stdout: Minimap2, Racon, and Medaka versions.\
            \n Requires:\
            \n    - Minimap2\
            \n Optional dependencies:\
            \n    - Racon\
            \n        o Required if using racon in consensus building.\
            \n    - Medaka\
            \n        o Required if using Medaka in consensus building.\
            \n        o Can be installed at ~/medaka through the python\
            \n          virtual enviroment (git hub install) or by\
            \n           miniconda\
        "; /*Help message*/

    char *conBuildHelpCStr = "\
            \n buildCon can use several different consensuses\
            \n   building methods to build a consensus. These methods\
            \n   can be run separately or can be combined together.\
            \n   When run together the majority consensus step builds\
            \n   the consensus using the selected read, while the Racon\
            \n   and Medaka steps are used to polish the consensus\
            \n   built with the majority consensus step. Medaka is used\
            \n   used to polish the consensus build by Racon when Racon\
            \n   is used to build a consensus.\
            \n\
            \n Note: The majority consensus step is the fastest\
            \n   consensus building method. It handles SNPs and matches\
            \n   well, provided their is enough read depth, but is weak\
            \n   for indels.\
            \n    -min-con-length:                               [500]\
            \n       - Discard consensuses that are under the\
            \n         input length.\
            \n       - If you lower this your should also lower\
            \n         -min-read-read-map-length &\
            \n         -min-read-con-map-length\
            \n    -max-reads-per-con:\
            \n        - Max number of reads to use in        [300]\
            \n          a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
            \n    -disable-majority-consensus: [Use majority consensus]\
            \n        - Build a consensus using a simple\
            \n          majority consensus. This consensus\
            \n          will be polished with Racon or\
            \n          Medaka if Racon and Medaka set.\
            \n    -maj-con-min-bases                         [0.35=35%]\
            \n        - When building the majority consesus\
            \n          make a deletion in positions that\
            \n          have less than x\% of bases (35%).\
            \n    -maj-con-min-base-q:                       [7]\
            \n        - Minimum q-score to keep a base when\
            \n          building a majority consensus.\
            \n    -maj-con-min-ins                           [0.3=30%]\
            \n        - When building the majority consesus\
            \n          ingore insertions that have support\
            \n          from less than x\% of reads (30%).\
            \n    -maj-con-min-ins-q:                        [5]\
            \n        - Minimum q-score to keep a insertion\
            \n          when building a majority consensus.\
            \n    -enable-racon:                                [No]\
            \n        - Do not use racon to polish the\
            \n          consensus.\
            \n    -rounds-racon:\
            \n        - Number of rounds to polish a         [4]\
            \n          consensus with racon\
            \n    -enable-medaka:                               [No]\
            \n        - Do not use medaka to polish the\
            \n          consensus.\
            \n    -model:\
            \n        - Model to use with medaka   [r941_min_high_g351]\
            \n          (calling medaka_consensus)\
       "; /*consensus building parameters*/

    char
        *readReadMapHelpCStr ="\
            \n The read mapping step uses a higher quality (by median\
            \n    Q-score) read in the bin to identify other high\
            \n    qaulity reads (by median Q-score) that likely came\
            \n    from the same virus. This is done by mapping reads\
            \n    to the selected read and removing all reads that are\
            \n    beneath the quality thresholds. Only the top X\
            \n    (default 300) are kept to build a consensus with.\
            \n\
            \n     -min-read-read-map-length:                     [500]\
            \n        - Minimum aligned read length needed\
            \n          to keep a read to read mapping.\
            \n        - Values less than this will not be\
            \n          extracted.\
            \n        - Change this will also require changing\
            \n          -min-read-con-map-length.\
            \n     -read-read-snps:                      [0.021 = 2.1%]\
            \n        - Minimum percentage of snps needed to\
            \n          discard a read during the read to\
            \n          read clustering step.\
            \n     -read-read-diff:                       [0.02 = 2%]\
            \n        - Minimum percent difference needed to\
            \n          discard a read during the read to\
            \n          read clustering step.\
            \n     -read-read-dels:                       [1 = 100%]\
            \n        - Minimum percentage of deletions\
            \n          needed to discard a read during\
            \n          the read to read clustering step.\
            \n     -read-read-inss:                       [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          needed to discard a read during\
            \n          the read to read clustering step.\
            \n     -read-read-indels:                     [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          and deletions needed to discard a\
            \n          read during the read to read\
            \n          clustering step.\
            \n\
            \n\
            \n The scoring settings for read mapping only counts SNPs\
            \n     and indels in reads that are at or above the user\
            \n     input thresholds. These counts are used to find the\
            \n     percentages in the comparison step.\
            \n\
            \n     -read-read-min-base-q:                 [10]\
            \n        - Minimum Q-score needed to keep an\
            \n          SNP or insertion when comparing\
            \n          reads.\
            \n     -read-read-min-mapq:                   [20]\
            \n        - Minimum mapping quality needed to\
            \n          keep a read when comparing reads.\
            \n\
            \n     -read-read-max-a-ins-homo:                [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-t-ins-homo:                [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-c-ins-homo:                [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-g-ins-homo:                [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            \n     -read-read-max-a-del-homo:                [0]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-t-del-homo:                [0]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-c-del-homo:                [0]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-g-del-homo:                [0]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            "; /*Read mapping help message*/


    char
        *readConMapHelpCStr = "\
            \n These paramaters are for controlling which reads are\
            \n    kept when mapped to the best read or the consensus.\
            \n\
            \n     -min-read-con-map-length:                   [500]\
            \n        - Minimum aligned read length needed to\
            \n          keep a read to consensus mapping.\
            \n        - Reads that have an alinged length less\
            \n          than this will not be extracted.\
            \n        - Change this will also require changing\
            \n          -min-read-read-map-length.\
            \n     -read-con-snps:                       [0.015 = 1.5%]\
            \n        - Minimum percentage of snps needed to\
            \n          discard a read during the read to\
            \n          consensus clustering step.\
            \n     -read-con-diff:                       [0.02 = 2%]\
            \n        - Minimum percent difference needed to\
            \n          discard a read during the read to\
            \n          consensus clustering step.\
            \n     -read-con-dels:                       [1 = 100%]\
            \n        - Minimum percentage of deletions\
            \n          needed to discard a read during\
            \n          the read to consensus clustering.\
            \n     -read-con-inss:                       [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          needed to discard a read during\
            \n          the read to consensus clustering.\
            \n     -read-con-indels:                     [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          and deletions needed to discard a\
            \n          read during the read to consensus\
            \n          clustering step.\
            \n\
            \n\
            \n The scoring settings for only counts SNPs and indels\
            \n     in reads that are at or above the input thresholds.\
            \n\
            \n     -read-con-min-base-q:                 [10]\
            \n        - Minimum Q-score needed to keep an\
            \n          SNP or insertion when clustering\
            \n          reads.\
            \n     -read-con-min-mapq:                   [20]\
            \n        - Minimum mapping quality needed to\
            \n          keep a read when clustering reads.\
            \n\
            \n     -read-con-max-a-ins-homo:             [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-t-ins-homo:             [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-c-ins-homo:             [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-g-ins-homo:             [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            \n     -read-con-max-a-del-homo:             [0]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-t-del-homo:             [0]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-c-del-homo:             [0]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-g-del-homo:             [0]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            "; /*Clustering help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Set up default settings in structures
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Set default consensus step settings*/
    initConBuildStruct(&conSetting);
    initSamEntry(&samStruct);
    initSamEntry(&refStruct);

    /*Set up filter and scoring settings*/
    blankMinStats(&readToReadMinStats);
    blankMinStatsReadCon(&minReadConStats);

    /*Make sure the struct holding the fastq file does not have noise*/
    fastqStruct.refIdCStr[0] = '\0';
    fastqStruct.statPathCStr[0] = '\0';
    fastqStruct.bestReadCStr[0] = '\0';
    fastqStruct.topReadsCStr[0] = '\0';
    fastqStruct.consensusCStr[0] = '\0';
    fastqStruct.numReadsULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Get user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr =
        getUserInput(
            lenArgsInt,
            argsCStr,
            &fqPathCStr,
            &refPathCStr,
            &statsPathCStr,
            prefixCStr,
            threadsCStr,
            &conSetting,
            &readToReadMinStats, /*to keep a read/read map*/
            &minReadConStats   /*To keep read/consensus map*/
    ); /*Get user input*/

    if(tmpCStr != 0)
    { /*If have an error*/

        if(tmpCStr == 0)
        { /*If no user input was supplied*/
            fprintf(
                stderr,
                "%s\n\nNo user input was supplied\n",
                helpCStr
            );

            exit(1);
        } /*If no user input was supplied*/

       if(
            strcmp(tmpCStr, "-v") == 0 ||
            strcmp(tmpCStr, "-V") == 0 ||
            strcmp(tmpCStr, "-version") == 0
        ) { /*If the user is requesting the version number*/
            fprintf(stdout, "buildCon built with findCoInft version:");
            fprintf(stdout, " %.8f\n", defVersion);
            exit(0);
        } /*If the user is requesting the version number*/

        if(strcmp(tmpCStr, "-h-build-consensus") == 0)
        { /*If the user wants the consensus building options*/
            fprintf(stdout, "%s", conBuildHelpCStr);
            exit(0);
        } /*If the user wants the consensus building options*/

        if(strcmp(tmpCStr, "-h-read-read-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readReadMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/
        
        if(strcmp(tmpCStr, "-h-read-con-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readConMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/

        if(
            strcmp(tmpCStr, "-h") == 0 ||
            strcmp(tmpCStr, "-help") == 0 ||
            strcmp(tmpCStr, "--h") == 0 ||
            strcmp(tmpCStr, "--help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);/*Print out help message*/
            exit(0);
        } /*If user wanted the help message*/

        fprintf(
            stderr,
            "%s\n\n%s is an invalid parameter\n",
            helpCStr,   /*Print out the help message*/
            tmpCStr     /*Print out the error*/
        ); /*Let user know about the invalid parameter*/

        exit(1);
    } /*If have an error*/

    if(!(
        conSetting.majConSet.useMajConBl |
        conSetting.raconSet.useRaconBl |
        conSetting.medakaSet.useMedakaBl
    )) { /*If the user said to ingore all consensus building steps*/
        printf("Current settings have turned off all consensus");
        printf(" building methods.\nSelect a consensus step by");
        printf(" removing: -enable-medaka or -enable-racon\n");
    } /*If the user said to ingore all consensus building steps*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4: Get minimap2, racon, & medaka versions
    ^   main sec-4 sub-1: Check if fastq file exists
    ^   main sec-4 sub-1: Find minimap2 version & check if exists
    ^   main sec-4 sub-2: If using Racon, find version & check if exists
    ^   main sec-4 sub-3: If using Medaka, find version
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-4 Sub-1: Check if the fastq file exists
    \******************************************************************/

    stdinFILE = fopen(fqPathCStr, "r");

    if(stdinFILE == 0)
    { /*If fastq file was invalid*/
        fprintf(
            stdout,
            "Could not open provided fastq file (%s)\n",
             fqPathCStr
         ); /*Let user know about the issue*/

         exit(1);
    } /*If fastq file was invalid*/

    fclose(stdinFILE);

    if(statsPathCStr != 0)
    { /*If the user provided a stats file for read selection*/
        stdinFILE = fopen(statsPathCStr, "r");
        if(stdinFILE == 0)
        { /*If I could not open the provided file*/
            fprintf(
                stdout,
                "Could not open stats file (%s) from scoreReads\n",
                 statsPathCStr
             ); /*Let user know about the issue*/

             exit(1);
        } /*If I could not open the provided file*/

        fclose(stdinFILE);
    } /*If the user provided a stats file for read selection*/


    /*Check if the reference file exists*/
    if(refPathCStr != 0)
    { /*If the user provided a reference file*/
        stdinFILE = fopen(refPathCStr, "r");

        if(stdinFILE == 0)
        { /*If fastq file was invalid*/
            fprintf(
                stdout,
                "Could not open provided reference file (%s)\n",
                 refPathCStr
             ); /*Let user know about the issue*/

             exit(1);
        } /*If fastq file was invalid*/

        fclose(stdinFILE);
    } /*If the user provided a reference file*/

    /******************************************************************\
    * Main Sec-4 Sub-1: Find minimap2 version & check if exists
    \******************************************************************/

    /*Set up minimap 2 check*/
    tmpCStr = cpParmAndArg(tmpCmdCStr, "minimap2", "--version");

    stdinFILE = popen(tmpCmdCStr, "r");
    fgets(minimap2VersionCStr, 256, stdinFILE);
    fclose(stdinFILE);

    if(minimap2VersionCStr[0] == '\0')
    { /*If could not find minimap2*/
        fprintf(stdout, "Minimap2 could not be found\n");
        exit(1);
    } /*If could not find minimap2*/

    /******************************************************************\
    * Main Sec-4 Sub-2: If using Racon, find version & check if exists
    \******************************************************************/

    if(conSetting.raconSet.useRaconBl & 1)
    { /*If using racon, get the version used*/
        /*Set up racon check*/
        tmpCStr = cpParmAndArg(tmpCmdCStr, "racon", "--version");

        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(raconVersionCStr, 256, stdinFILE);
        fclose(stdinFILE);

        if(raconVersionCStr[0] == '\0')
        { /*If racon does not exist*/
            fprintf(stderr, "Racon could not be found\n");
            exit(1);
        } /*If racon does not exist*/
    } /*If using racon, get the version used*/

    /******************************************************************\
    * Main Sec-4 Sub-3: If using Medaka, find version & check if exists
    \******************************************************************/

    if(conSetting.medakaSet.useMedakaBl & 1)
    { /*If using medaka, check version*/
        /*Set up non-miniconda command*/
        tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medakaCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, medakaCMDEnd, "");

        stdinFILE = popen(tmpCmdCStr, "r");
        fgets(medakaVersionCStr, 2048, stdinFILE);
        fclose(stdinFILE);

        if(medakaVersionCStr[0] == '\0')
        { /*If python virtual enviorment medaka does not exist*/
            /*Set up miniconda medaka enviroment command*/
            tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medCondCMD);
            tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
            tmpCStr = cStrCpInvsDelm(tmpCStr, medCondCMDEnd);

            stdinFILE = popen(tmpCmdCStr, "r");
            fgets(medakaVersionCStr, 2048, stdinFILE);
            fclose(stdinFILE);

            if(medakaVersionCStr[0] == '\0')
            { /*If medaka could not be found*/
                fprintf(stderr, "Medaka could not be found\n");
                exit(1);
            } /*If medaka could not be found*/

            /*Mark that I am using medaka from miniconda*/
            conSetting.medakaSet.condaBl = 1;
        } /*If python virtual enviorment medaka does not exist*/
    } /*If using medaka, check version*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Build the consensus
    ^   o main sec-5 sub-1: Make copy of fastq file so orignal is safe
    ^   o main sec-5 sub-2: Make copy of stats file so orignal is safe
    ^   o main sec-5 sub-3: Build the consensus
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-5 Sub-1: Make copy of fastq file so orignal is safe
    \******************************************************************/

    tmpCStr = cStrCpInvsDelm(tmpPathCStr, prefixCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".fastq");
    copyFile(fqPathCStr, tmpPathCStr);
    strcpy(fastqStruct.fqPathCStr, tmpPathCStr);

    /******************************************************************\
    * Main Sec-5 Sub-2: Make copy of stats file so orignal is safe
    \******************************************************************/

    if(statsPathCStr != 0)
    { /*If I need to copy over the stats file*/
        tmpCStr = cStrCpInvsDelm(tmpPathCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--stats.tsv");
        stdinFILE = fopen(statsPathCStr, "r");
        cpFILE = fopen(tmpPathCStr, "w");

        while(fgets(tmpCmdCStr, 1024, stdinFILE))
            fprintf(cpFILE, "%s", tmpCmdCStr);

        fclose(stdinFILE);
        fclose(cpFILE);
        strcpy(fastqStruct.statPathCStr, tmpPathCStr);

        conSetting.useStatBl = 1;
    } /*If I need to copy over the stats file*/

    /******************************************************************\
    * Main Sec-5 Sub-2: Build the consensus
    \******************************************************************/

    errUC =
        buildCon(
            &fastqStruct,     /*Fastq file with reads*/
            refPathCStr,      /*Fasta file with reference, 0 to ignore*/
            threadsCStr,      /*Number threads to use with system calls*/
            &conSetting,
            &samStruct,
            &refStruct,
            &readToReadMinStats,  /*Min stats to keep mapped reads*/
            &minReadConStats    /*Min stats to keep mapped reads*/
    ); /*Build the consensus if possible*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Clean up, Print out version numbers, and exit
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*remove(fastqStruct.bestReadCStr); *//*Make sure no extra files*/
        /*buildCon function already did this*/
    remove(fastqStruct.fqPathCStr);
    remove(fastqStruct.bestReadCStr);

    if(statsPathCStr != 0)
        remove(fastqStruct.statPathCStr); /*Make sure no extra files*/

    freeStackSamEntry(&samStruct);
    freeStackSamEntry(&refStruct);

    if(errUC & 64)
    { /*If had a memory allocation error*/
        fprintf(stdout, "Memory allocation error\n");
        exit(1);
    } /*If had a memory allocation error*/

    fprintf(stdout, "Minimap2 version: %s\n", minimap2VersionCStr);

    if(conSetting.raconSet.useRaconBl & 1)
        fprintf(stdout, "Racon version: %s\n", raconVersionCStr);
    if(conSetting.medakaSet.useMedakaBl & 1)
        fprintf(stdout, "Medaka version: %s\n", medakaVersionCStr);

    if(errUC & 16)
    { /*If had a memory allocation error*/
        fprintf(stdout, "Unable to build a consensus\n");
        exit(0);
    } /*If had a memory allocation error*/

    exit(0);
} /*main*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies input variables to hold user input
|    - Returns:
|        - 0: if nothing went wrong
|        - pionter to invalid paramter
|            - if no agruments input, returns ponter to argsCStr
\---------------------------------------------------------------------*/
char * getUserInput(
    int32_t lenArgsInt,
    char *argsCStr[],
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refPathCStr, /*Holds path to references*/
    char **statsPathCStr, /*Holds path to scoreReads output to use*/
    char *prefixCStr,   /*Prefix to name everything*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    struct conBuildStruct *conSet, 
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *minReadConStats  /*Read cluster socring set*/
) /*Reads in user input*/
{ /*getUserInput*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: getUserInput
    '    fun-1 sec-1: variable declerations
    '    fun-1 sec-2: Get user input
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   char
       *tmpCStr = 0,
       *inputCStr = 0, /*Points to user input part of line*/
       *parmCStr = 0;   /*Points to argument part of parameter*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Get user input
    ^    fun-1 sec-2 sub-1: General input
    ^    fun-1 sec-2 sub-2: Percent difference settings
    ^    fun-1 sec-2 sub-3: scoreReads read to reference settings
    ^    fun-1 sec-2 sub-4: scoreReads read to read settings
    ^    fun-1 sec-2 sub-5: scoreReads read to consensus settings
    ^    fun-1 sec-2 sub-6: scoreReads consensus to consensus settings
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-1 Sec-2 Sub-1: General input
    \******************************************************************/

    if(lenArgsInt < 2)
        return 0;       /*Nothing input*/

    for(int32_t intArg = 1; intArg < lenArgsInt; intArg += 2)
    { /*Loop through all user arguments (for)*/
        parmCStr = *(argsCStr + intArg); /*Get the parameter used*/
        inputCStr = *(argsCStr + intArg + 1); /*setting*/

        if(strcmp(parmCStr, "-fastq") == 0)
            *fqPathCStr = inputCStr;  /*Fastq file to check*/

        else if(strcmp(parmCStr, "-stats") == 0)
            *statsPathCStr = inputCStr; /*Stats file with reads to use*/

        else if(strcmp(parmCStr, "-ref") == 0)
            *refPathCStr = inputCStr;  /*references*/

        else if(strcmp(parmCStr, "-prefix") == 0)
            strcpy(prefixCStr, inputCStr); /*Prefix to name files with*/

        else if(strcmp(parmCStr, "-threads") == 0)
            strcpy(threadsCStr, inputCStr); /*Number of threads*/

        else if(strcmp(parmCStr, "-model") == 0)
            strcpy(conSet->medakaSet.modelCStr, inputCStr);

        else if(strcmp(parmCStr, "-min-reads-per-bin") == 0)
            conSet->minReadsToBuildConUL=strtoul(inputCStr,&tmpCStr,10);

        else if(strcmp(parmCStr, "-max-reads-per-con") == 0)
            conSet->maxReadsToBuildConUL=strtoul(inputCStr,&tmpCStr,10);

        else if(strcmp(parmCStr, "-rounds-racon") == 0)
            cStrToUChar(inputCStr, &conSet->raconSet.rndsRaconUC);
 
        else if(strcmp(parmCStr, "-extra-consensus-steps") == 0)
            cStrToUInt(inputCStr, &conSet->numRndsToPolishUI);

        else if(strcmp(parmCStr, "-min-con-length") == 0)
            conSet->minConLenUI = strtoul(inputCStr, &tmpCStr, 10);

        else if(strcmp(parmCStr, "-disable-majority-consensus") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            conSet->majConSet.useMajConBl = 0;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/

        else if(strcmp(parmCStr, "-enable-racon") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            conSet->raconSet.useRaconBl = 1;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/

        else if(strcmp(parmCStr, "-enable-medaka") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            conSet->medakaSet.useMedakaBl = 1;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/
            
        else if(strcmp(parmCStr, "-maj-con-min-bases") == 0)
          sscanf(inputCStr,"%f",&conSet->majConSet.minReadsPercBaseFlt);

        else if(strcmp(parmCStr, "-maj-con-min-ins") == 0)
           sscanf(inputCStr,"%f",&conSet->majConSet.minReadsPercInsFlt);

        else if(strcmp(parmCStr, "-maj-con-min-base-q") == 0)
            cStrToUChar(inputCStr, &conSet->majConSet.minBaseQUC);

        else if(strcmp(parmCStr, "-maj-con-min-ins-q") == 0)
            cStrToUChar(inputCStr, &conSet->majConSet.minInsQUC);
            
        /**************************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-min-read-con-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            minReadConStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);
        } /*Else if the user provided a minimum read length*/

        else if(strcmp(parmCStr, "-read-con-snps") == 0)
            sscanf(inputCStr, "%f", &minReadConStats->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-con-diff") == 0)
            sscanf(inputCStr, "%f", &minReadConStats->minDiffFlt);
        else if(strcmp(parmCStr, "-read-con-dels") == 0)
            sscanf(inputCStr, "%f", &minReadConStats->minDelsFlt);
        else if(strcmp(parmCStr, "-read-con-inss") == 0)
            sscanf(inputCStr, "%f", &minReadConStats->minInssFlt);
        else if(strcmp(parmCStr, "-read-con-indels") == 0)
            sscanf(inputCStr, "%f", &minReadConStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to consensus settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-con-min-base-q") == 0)
            cStrToUChar(inputCStr, &minReadConStats->minQChar);

        else if(strcmp(parmCStr, "-read-con-min-mapq") == 0)
            cStrToUInt(inputCStr, &minReadConStats->minMapqUInt);

         else if(strcmp(parmCStr, "-read-con-max-a-ins-homo") == 0)
           cStrToUInt(inputCStr, &minReadConStats->maxHomoInsAry[0]);

         else if(strcmp(parmCStr, "-read-con-max-t-ins-homo") == 0)
           cStrToUInt(inputCStr,&minReadConStats->maxHomoInsAry[10]);

         else if(strcmp(parmCStr, "-read-con-max-c-ins-homo") == 0)
            cStrToUInt(inputCStr,&minReadConStats->maxHomoInsAry[1]);

         else if(strcmp(parmCStr, "-read-con-max-g-ins-homo") == 0)
            cStrToUInt(inputCStr,&minReadConStats->maxHomoInsAry[3]);

         else if(strcmp(parmCStr, "-read-con-max-a-del-homo") == 0)
           cStrToUInt(inputCStr, &minReadConStats->maxHomoDelAry[0]);

         else if(strcmp(parmCStr, "-read-con-max-t-del-homo") == 0)
           cStrToUInt(inputCStr,&minReadConStats->maxHomoDelAry[10]);

         else if(strcmp(parmCStr, "-read-con-max-c-del-homo") == 0)
            cStrToUInt(inputCStr,&minReadConStats->maxHomoDelAry[1]);

         else if(strcmp(parmCStr, "-read-con-max-g-del-homo") == 0)
            cStrToUInt(inputCStr,&minReadConStats->maxHomoDelAry[3]);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-min-read-read-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            readToReadMinStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);
        } /*Else if the user provided a minimum read length*/

        else if(strcmp(parmCStr, "-read-read-snps") == 0)
            sscanf(inputCStr, "%f", &readToReadMinStats->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-read-diff") == 0)
            sscanf(inputCStr, "%f", &readToReadMinStats->minDiffFlt);
        else if(strcmp(parmCStr, "-read-read-dels") == 0)
            sscanf(inputCStr, "%f", &readToReadMinStats->minDelsFlt);
        else if(strcmp(parmCStr, "-read-read-inss") == 0)
            sscanf(inputCStr, "%f", &readToReadMinStats->minInssFlt);
        else if(strcmp(parmCStr, "-read-read-indels") == 0)
            sscanf(inputCStr, "%f", &readToReadMinStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to readsensus settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-read-min-base-q") == 0)
            cStrToUChar(inputCStr, &readToReadMinStats->minQChar);

        else if(strcmp(parmCStr, "-read-read-min-mapq") == 0)
            cStrToUInt(inputCStr, &readToReadMinStats->minMapqUInt);

         else if(strcmp(parmCStr, "-read-read-max-a-ins-homo") == 0)
           cStrToUInt(inputCStr, &readToReadMinStats->maxHomoInsAry[0]);

         else if(strcmp(parmCStr, "-read-read-max-t-ins-homo") == 0)
           cStrToUInt(inputCStr,&readToReadMinStats->maxHomoInsAry[10]);

         else if(strcmp(parmCStr, "-read-read-max-c-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToReadMinStats->maxHomoInsAry[1]);

         else if(strcmp(parmCStr, "-read-read-max-g-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToReadMinStats->maxHomoInsAry[3]);

         else if(strcmp(parmCStr, "-read-read-max-a-del-homo") == 0)
           cStrToUInt(inputCStr, &readToReadMinStats->maxHomoDelAry[0]);

         else if(strcmp(parmCStr, "-read-read-max-t-del-homo") == 0)
           cStrToUInt(inputCStr,&readToReadMinStats->maxHomoDelAry[10]);

         else if(strcmp(parmCStr, "-read-read-max-c-del-homo") == 0)
            cStrToUInt(inputCStr,&readToReadMinStats->maxHomoDelAry[1]);

         else if(strcmp(parmCStr, "-read-read-max-g-del-homo") == 0)
            cStrToUInt(inputCStr,&readToReadMinStats->maxHomoDelAry[3]);

        else
            return parmCStr;
    } /*Loop through all user arguments (for)*/

    return 0;
} /*getUserInput*/
