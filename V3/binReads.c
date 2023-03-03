/*######################################################################
# Use:
#   o Bins reads to a set of input references
# Requires:
#   o cStrFun.c/h
#   o cStrToNumberFun.c/h
#   o defaultSettings.h
#   o printError.c/h
#   o samEntryStruct.c/h
#   o trimSam.c/h
#   o scoreReadsFun.c/h
#   o findCoInftBinTree.c/h
#   o findCoInftChecks.c/h
#   o binReadFun.c/h
# C libaries:
#   o string.h
#   o stdlib.h
#   o stdio.h
#   o stdint.h
######################################################################*/

#include "binReadsFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' HOP:
'   o header: Has all the includes and other misc stuff
'   o main: Main function to drive everything
'   o getUserInput: Process the user supplied input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
    char *prefCStr,                      /*Holds user supplied prefix*/
    char **fqPathCStr,                   /*Holds path to fastq file*/
    char **refsPathCStr,                 /*Holds path to references*/
    char *threadsCStr,                   /*Number threads for minimap2*/
    char *rmSupAlnBl,       /*Remove reads with supplemenat alignments*/
    char *trimBl,                        /*1 trim reads, 0 do not*/
    unsigned long *minReadsPerBinUL,       /*Min # reads to keep a bin*/
    struct minAlnStats *readToRefMinStats  /*Binning scoring settings*/
); /*Reads in user input*/

int main(
    int32_t lenArgsInt, /*Number of parameters & arguments user input*/
    char *argsCStr[]  /*List of parameters & arguements user input*/
) /*Function to run everything*/
{ /*main*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Main TOC:
    '   o main sec-1: Variable declerations
    '   o main sec-2: Set structer default settings & initalize structs
    '   o main sec-3: Get user input
    '   o main sec-4: Check the version of minimap2
    '   o main sec-5: Bin reads
    '   o main sec-6: Remove low read count bins
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-1: Variable declerations
    ^   o main sec-1 sub-1: Program variables
    ^   o main sec-1 sub-2: Help messages
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-1 Sub-1: Program variables
    \******************************************************************/

    char *fqPathCStr = 0;        /*Fastq file to bin*/
    char *refsPathCStr = 0;      /*References to bin with*/
    char trimBl = 0;             /*1: trim reads, 0: do not*/
    char prefixCStr[128] = defPrefix;  /*Prefix to name the bins*/
    char threadsCStr[16] = defThreads; /*Number of threads to use*/
    char rmSupAlnBl = rmReadsWithSupAln;
       /*if rmSupAlnBl = 1, Remove reads with supplementary alignments*/

    char minimap2CmdCStr[1024]; /*To check if minimap2 exists*/
    char minimap2VersionCStr[1024]; /*holds minimap2 version*/
    char *inputErrC = 0;
    char *tmpCStr = 0;
    char readCntNameCStr[256]; /*Name of the read count file*/
    unsigned char errUC = 0;

    unsigned long minReadsPerBinUL = minReadsPerBin;

    FILE *readCountFILE = 0;

    struct readBin *binTree = 0; /*Holds output bins, will be freeded*/

    struct samEntry samStruct;
    struct samEntry oldStruct;

    struct minAlnStats minStats;

    /******************************************************************\
    * Main Sec-1 Sub-2: Help messages
    \******************************************************************/

    char *helpCStr = "\
            \n findCoInfc -fastq reads.fastq -ref refs.fasta [Options]\
            \n    -fastq:\
            \n        - Fastq file with reads to search      [Required]\
            \n    -ref:\
            \n        -Fasta file with references for bining [Required]\
            \n    -prefix:\
            \n        - Prefix to add to file names          [Out]\
            \n    -threads:\
            \n        - Number of threads to use             [3]\
            \n    -min-reads-per-bin:\
            \n        - Min number of reads needed to keep   [100]\
            \n          a bin or a cluster\
            \n    -read-ref-min-read-length:                 [600]\
            \n       - Discard reads with read lengths under\
            \n         this when binning.\
            \n       - Applied after the trimming step.\
            \n    -read-ref-max-read-length:                 [1000]\
            \n       - Discard reads with read lengths over\
            \n         this when binning.\
            \n       - Applied after the trimming step.\
            \n    -trim:                                     [No]\
            \n       - Trim reads to reference.\
            \n    -rm-sup-reads                 [No]\
            \n       - Removes any read that has a\
            \n         supplemental alignment. These might\
            \n         be chimeras or could just be repeated\
            \n         genes.\
            \n       - If this setting is off, then supplemental\
            \n         alignments are just ignored.\
            \n       - This setting requires running minimap2 with\
            \n         just one thread for binning.\
            \n Additional Help messages:\
            \n    -h-bin:\
            \n        - Print out the parameters for the binning step.\
            \n Output:\
            \n    - fastq files: With the reads for each co-infection\
            \n    - prefix--read-counts.tsv:\
            \n        o File with number of reads per kept cluster\
            \n        o Also has the read counts for each discard bin\
            \n          or cluster. This makes this file a good for\
            \n          checking this programs progress. It will often\
            \n          be printing bining read counts while\
            \n          clustering.\
            \n Requires:\
            \n    - Minimap2\
        "; /*main help message*/

    char
         *binHelpCStr = "\
            \n The binning step compares a read to the best reference\
            \n    its mapped to. If the read meets the qualitity\
            \n    thresholds then it is assigned to the references bin.\
            \n    The read is discarded if it does not meet the quality\
            \n    thresholds.\
            \n\
            \n -read-ref-snps:                           [0.02 = 2%]\
            \n    - Minimum percentage of snps needed to\
            \n      discard a read during the read to\
            \n      reference binning step.\
            \n -read-ref-diff:                           [0.03 = 3%]\
            \n    - Minimum percent difference needed to\
            \n      discard a read during the read to\
            \n      reference binning step.\
            \n -read-ref-dels:                           [1 = 100%]\
            \n    - Minimum percentage of deletions\
            \n      needed to discard a read during\
            \n      the read to reference binning step.\
            \n -read-ref-inss:                           [1 = 100%]\
            \n    - Minimum percentage of insertions\
            \n      needed to discard a read during\
            \n      the read to reference binning step.\
            \n -read-ref-indels:                         [1 = 100%]\
            \n    - Minimum percentage of insertions\
            \n      and deletions needed to discard a\
            \n      read during the read to reference\
            \n      binning step.\
            \n\
            \n -read-ref-min-mapq:                        [20]\
            \n    - Minimum mapping quality needed to\
            \n      keep a read when binning.\
            \n -read-ref-min-read-length:                 [600]\
            \n    - Discard reads with read lengths under\
            \n      this when binning.\
            \n    - Applied after the trimming step.\
            \n -read-ref-max-read-length:                 [1000]\
            \n    - Discard reads with read lengths over\
            \n      this when binning.\
            \n    - Applied after the trimming step.\
            \n\
            \n -read-ref-min-median-q:                    [13]\
            \n    - Minimum read median quality score\
            \n      needed to keep a read when binning.\
            \n -read-ref-min-mean-q:                      [13]\
            \n    - Minimum read mean quality score\
            \n      needed to keep a read when binning.\
            \n -read-ref-min-aligned-median-q:            [13]\
            \n    - Minimum read median quality score\
            \n      of the aligned region of a read\
            \n      needed to keep a read when binning.\
            \n -read-ref-min-aligned-mean-q:              [13]\
            \n    - Minimum read mean quality score\
            \n      of the aligned region of a read\
            \n      needed to keep a read when binning.\
            \n\
            \n\
            \n The scoring settings for binning only counts SNPs and\
            \n     or indels that are at or above the user input\
            \n     thresholds. These counts are used to find the\
            \n     percentages in the comparison step.\
            \n\
            \n -read-ref-min-base-q:                      [13]\
            \n    - Minimum Q-score needed to keep an\
            \n      SNP or insertion when binning.\
            \n\
            \n -read-ref-max-a-ins-homo:                  [1]\
            \n    - Maximum A homopolymer size to keep\
            \n      an insertion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-t-ins-homo:                  [1]\
            \n    - Maximum T homopolymer size to keep\
            \n      an insertion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-c-ins-homo:                  [1]\
            \n    - Maximum C homopolymer size to keep\
            \n      an insertion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-g-ins-homo:                  [1]\
            \n    - Maximum G homopolymer size to keep\
            \n      an insertion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n\
            \n -read-ref-max-a-del-homo:                  [0]\
            \n    - Maximum A homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-t-del-homo:                  [0]\
            \n    - Maximum T homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-c-del-homo:                  [0]\
            \n    - Maximum C homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-g-del-homo:                  [0]\
            \n    - Maximum G homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            "; /*binning parameters help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Set default settings for structers / initalize structs
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankMinStats(&minStats);
    initSamEntry(&samStruct); /*Remove old stats in sam file*/
    initSamEntry(&oldStruct); /*Remove old stats in sam file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-3: Get user input
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    inputErrC =
        getUserInput(
            lenArgsInt,
            argsCStr,
            prefixCStr,           /*Holds user supplied prefix*/
            &fqPathCStr,          /*Holds path to fastq file*/
            &refsPathCStr,        /*Holds path to references*/
            threadsCStr,          /*Number threads for minimap2*/
            &rmSupAlnBl,    /*Remove reads with supplemenat alignments*/
            &trimBl,              /*1 trim reads, 0 do not*/
            &minReadsPerBinUL,    /*Min # reads to keep a bin*/
            &minStats             /*scoreReads variables*/
    ); /*Get user input*/

    if(inputErrC != 0)
    { /*If had a problematic parameter*/
        if(inputErrC == 0)
        { /*If no user input was supplied*/
            fprintf(stderr, "%s\n\nNo input was supplied\n", helpCStr);
            exit(1);
        } /*If no user input was supplied*/

        if(
            strcmp(inputErrC, "-v") == 0 ||
            strcmp(inputErrC, "-V") == 0 ||
            strcmp(inputErrC, "-version") == 0
        ) { /*If the user is requesting the version number*/
            fprintf(
                stdout,
                "binReads built from findCoInft version: %f\n",
                defVersion
            );
            exit(0);
        } /*If the user is requesting the version number*/

        if(strcmp(inputErrC, "-h-bin") == 0)
        { /*If user wants to know about the binning parameters*/
            fprintf(stdout, "%s", binHelpCStr);
            exit(0);
        } /*If user wants to know about the binning parameters*/

        if(
            strcmp(inputErrC, "-h") == 0 ||
            strcmp(inputErrC, "-help") == 0 ||
            strcmp(inputErrC, "--h") == 0 ||
            strcmp(inputErrC, "--help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);/*Print out help message*/
            exit(0);
        } /*If user wanted the help message*/

        fprintf(
            stderr,
            "%s\n\n %s is an invalid parameter\n",
            helpCStr,   /*Print out the help message*/
            inputErrC     /*Print out the error*/
        ); /*Let user know about the invalid parameter*/

        exit(1);
    } /*If have an error*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4: Check if minimap2 exists
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = cpParmAndArg(minimap2CmdCStr, "minimap2", "--version");
    readCountFILE = popen(minimap2CmdCStr, "r");

    minimap2VersionCStr[0] = '\0';
    fgets(minimap2VersionCStr, 1024, readCountFILE);
    pclose(readCountFILE);

    if(minimap2VersionCStr[0] == '\0')
    { /*If could not find minimap2*/
        fprintf(stdout, "Minimap2 could not be found\n");
        exit(1);
    } /*If could not find minimap2*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Bin reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    binTree =
        binReads(
            fqPathCStr,        /*Fastq file to bin*/
            refsPathCStr,      /*References to bin with*/
            prefixCStr,
            threadsCStr,       /*Number threads to use with minimap2*/
            rmSupAlnBl,   /*Remove reads with supplementary alignments*/
            trimBl,       /*1: trim reads, 0: do not*/
            &samStruct, /*Holds minimap2 output*/
            &oldStruct, /*Holds previous line of minimap2 output*/
            &minStats,
            &errUC     /*Reports any errors*/
    ); /*Bin reads using the provided references*/

    /*No longer need the samEntry structures*/
    freeStackSamEntry(&samStruct);
    freeStackSamEntry(&oldStruct);

    if(binTree == 0)
    { /*If had an error*/
        if(errUC & 64)
            fprintf(stdout, "Memory error: not enough memory\n");
        else if(errUC & 2)
            fprintf(
                stdout,
                "Fastq file (%s) could not be opened\n",
                fqPathCStr
            ); /*Let user know fastq file is invalid*/
        else
            fprintf(stdout, "Could not create fastq or stats files\n");

        exit(1);
    } /*If had an error*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Remove low read count bins
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the read count file for writing*/
    tmpCStr = cStrCpInvsDelm(readCntNameCStr, prefixCStr);
    strcpy(tmpCStr, "--read-counts.tsv");
    readCountFILE = fopen(readCntNameCStr, "w");

    if(readCountFILE == 0)
    { /*If could not open the file*/
        fprintf(stdout, "Could not make file to right bin counts to\n");
        exit(1);
    } /*If could not open the file*/

    cnvtBinTreeToList(&binTree); /*Conver to linked list*/

    while(binTree != 0)
    { /*While I have bins to free*/
        fprintf(
            readCountFILE,
            "%s\t%lu",
            binTree->refIdCStr,
            binTree->numReadsULng
        ); /*Recording the number of reads in the bin*/

        if(binTree->numReadsULng >= minReadsPerBinUL)
        { /*If this bin is good, make sure files doe not get deleted*/
            binTree->fqPathCStr[0] = '\0';
            binTree->statPathCStr[0] = '\0';
            fprintf(readCountFILE, "\tkept\n");
        } /*If this bin is good, make sure files doe not get deleted*/

        else
            fprintf(readCountFILE, "\tremoved\n");

        fflush(readCountFILE); /*make sure io printed out*/
        rmBinFromList(&binTree); /*Also removes all files*/
    } /*While I have bins to free*/

    /*Let user know the minimap2 version used*/
    fprintf(stdout, "Minimap2 version: %s\n", minimap2VersionCStr);

    fclose(readCountFILE);
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
    char *prefCStr,                      /*Holds user supplied prefix*/
    char **fqPathCStr,                   /*Holds path to fastq file*/
    char **refsPathCStr,                 /*Holds path to references*/
    char *threadsCStr,                   /*Number threads for minimap2*/
    char *rmSupAlnBl,       /*Remove reads with supplemenat alignments*/
    char *trimBl,                        /*1 trim reads, 0 do not*/
    unsigned long *minReadsPerBinUL,     /*Min # reads to keep a bin*/
    struct minAlnStats *readToRefMinStats  /*Binning scoring settings*/
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
   char *tmpCStr = 0;
   char *inputCStr = 0; /*Points to user input part of line*/
   char *parmCStr = 0;   /*Points to argument part of parameter*/

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

        else if(strcmp(parmCStr, "-ref") == 0)
            *refsPathCStr = inputCStr;  /*references*/

        else if(strcmp(parmCStr, "-prefix") == 0)
            strcpy(prefCStr, inputCStr);     /*Have prefix to use*/

        else if(strcmp(parmCStr, "-threads") == 0)
            strcpy(threadsCStr, inputCStr); /*Number of threads*/

        else if(strcmp(parmCStr, "-min-reads-per-bin") == 0)
            *minReadsPerBinUL = strtoul(inputCStr, &tmpCStr, 10);

        else if(strcmp(parmCStr, "-rm-sup-reads") == 0)
        { /*Else if removing supplemntal alignments*/
            *rmSupAlnBl = !(*rmSupAlnBl);
            --intArg;
        } /*Else if removing supplemntal alignments*/

        else if(strcmp(parmCStr, "-trim") == 0)
        { /*Else if trimming the reads*/
            *trimBl = !(*trimBl);
            --intArg;
        } /*Else if trimming the reads*/

        /**************************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-ref-snps") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-ref-diff") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minDiffFlt);
        else if(strcmp(parmCStr, "-read-ref-dels") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minDelsFlt);
        else if(strcmp(parmCStr, "-read-ref-inss") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minInssFlt);
        else if(strcmp(parmCStr, "-read-ref-indels") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-3: scoreReads read to reference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-ref-min-base-q") == 0)
            cStrToUChar(inputCStr, &readToRefMinStats->minQChar);

        else if(strcmp(parmCStr, "-read-ref-min-mapq") == 0)
            cStrToUInt(inputCStr, &readToRefMinStats->minMapqUInt);

        else if(strcmp(parmCStr, "-read-ref-min-median-q") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minMedianQFlt);

        else if(strcmp(parmCStr, "-read-ref-min-mean-q") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minMeanQFlt);

        else if(strcmp(parmCStr, "-read-ref-min-aligned-median-q") == 0)
            sscanf(
                inputCStr,
                "%f",
                &readToRefMinStats->minAlignedMedianQFlt
            );

        else if(strcmp(parmCStr, "-read-ref-min-aligned-mean-q") == 0)
            sscanf(
                inputCStr,
                "%f",
                &readToRefMinStats->minAlignedMeanQFlt
            );

        else if(strcmp(parmCStr, "-read-ref-min-read-length") == 0)
            readToRefMinStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);

         else if(strcmp(parmCStr, "-read-ref-max-read-length") == 0)
            readToRefMinStats->maxReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);

         else if(strcmp(parmCStr, "-read-ref-max-a-ins-homo") == 0)
            cStrToUInt(inputCStr, &readToRefMinStats->maxHomoInsAry[0]);

         else if(strcmp(parmCStr, "-read-ref-max-t-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoInsAry[10]);

         else if(strcmp(parmCStr, "-read-ref-max-c-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoInsAry[1]);

         else if(strcmp(parmCStr, "-read-ref-max-g-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoInsAry[3]);

         else if(strcmp(parmCStr, "-read-ref-max-a-del-homo") == 0)
            cStrToUInt(inputCStr, &readToRefMinStats->maxHomoDelAry[0]);

         else if(strcmp(parmCStr, "-read-ref-max-t-del-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoDelAry[10]);

         else if(strcmp(parmCStr, "-read-ref-max-c-del-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoDelAry[1]);

         else if(strcmp(parmCStr, "-read-ref-max-g-del-homo") == 0)
            cStrToUInt(inputCStr,&readToRefMinStats->maxHomoDelAry[3]);

        else
            return parmCStr;
    } /*Loop through all user arguments (for)*/

    return 0;
} /*getUserInput*/
