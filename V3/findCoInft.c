/*#####################################################################\
# Name: findCoInfect
# Use:
#    - Detects co-infections in nanopore sequenced reads using Minimap2
# Optinal Dependencies (depends on user input settings):
#    - Minimap2 (Required for -enable-medaka)
#    - Racon    (Required for -enable-racon)
# Internal Requirments:
#    - Every *.c file that has a *.h file in this directory
#      o buildCon.c, binReads.c, scoreReads.c, trimSamFile.c,
#        fqGetIds.c, and fqGetIdsThread.c are not needed, but do allow
#        for building some internal commands as stand alone programs.
# Libraries used (in various files):
#     - <stdlib.h>
#     - <sdtint.h>
#     - <stdio.h>
#     - <string.h>
\#####################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start of program
'    header:
'      - Header section for functions, structures, & libraries
'    main:
'      - main function that runs everything
'    fun-1 getUserInput:
'      - function to get user input with
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header: Included files
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*My includes*/
#include "buildConFun.h" /*Various dependencies through readExtract.h*/
#include "binReadsFun.h" /*Functions for binning reads*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies input variables to hold user input
|    - Returns:
|        - 1: if nothing went wrong
|        - 2: if no arguments were supplied
\---------------------------------------------------------------------*/
char * getUserInput(
    int32_t lenArgsInt,
    char *argsCStr[],
    char *prefCStr,  /*Holds user supplied prefix*/
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refsPathCStr, /*Holds path to references*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    char *rmSupAlnBl,  /*1: Remove reads with supplemental alignments*/
    char *skipBinBl,   /*1: Skip the binning step, 0 do not*/
    char *skipClustBl, /*1: Skip the clusterin step, 0 do not*/
    double *minReadsDbl,
    struct conBuildStruct *conSet,   /*Settings for consensus building*/
    struct minAlnStats *readToRefMinStats, /*Binning scoring settings*/
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *readToConMinStats, /*Read cluster socring set*/
    struct minAlnStats *conToConMinStats   /*Consensus comparison set*/
); /*Reads in user input*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start of functions
'   o main: driver function to run find co-infections
'   o getUserInput: Process user input provided by command line
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
    '   main sec-5: Put non-required input into log for user
    '   main sec-6: Find initial bins with references
    '   main sec-7: Non-reference based binning (clustering?)
    '   main sec-8: Compare all consensus to remove false positives
    '   main sec-9: Use longest read to build consensus (TO DO)
    '   main sec-10: Clean up and exit
    '
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
    char rmSupAlnBl = rmReadsWithSupAln;
    char skipBinBl = defSkipBinBl;     /*Skip binning step?*/
    char skipClustBl = defSkipClustBl; /*Skip clustering step?*/
    char prefCStr[100];        /*Holds the user prefix*/
    char threadsCStr[7];      /*Number of threads for minimap2 & racon*/
    double minReadsDbl = defMinPercReads;
        /*What percentage of reads should the final cluster have*/

    char *fqPathCStr = 0;      /*Holds fastq file to process*/
    char *refsPathCStr = 0;    /*Holds references for binning*/
    char logFileCStr[256];    /*Holds the name of the log file*/
    char readCntFileCStr[256]; /*Holds Number of reads per bin/cluster*/

    /*C-strings that hold commands*/
    char tmpCmdCStr[1024];      /*Holds a quick system command*/
    char *tmpCStr = 0;          /*For string manipulation*/

    /*Miscalanious variables*/
    uint8_t errUC = 0;         /*Holds error messages*/

    char *inutErrCStr = 0; /*holds user input error*/

    unsigned long totalKeptReadsUL = 0;

    /*FILES opened*/
    FILE *logFILE = 0;      /*Holds the log*/
    FILE *stdinFILE = stdin;/*Points to piped input from minimap2*/
    FILE *statFILE = 0;     /*File to output stats from score reads to*/
    /*FILE *tmpFILE = 0;*/
    /*FILE *fqBinFILE = 0;*/  /*Points to file adding binned reads to*/

    /*My structures*/
    struct samEntry samStruct;    /*Holds sam file entry from minimap2*/
    struct samEntry refStruct;    /*Holds the reference sequence*/

    /*Holds thresholds to keep alignments, errors, & matches*/
    struct minAlnStats readToRefMinStats;
        /*map read to reference thresholds*/
    struct minAlnStats readToReadMinStats;
        /*read to read mapping thresolds*/
    struct minAlnStats readToConMinStats;
        /*map read to reference thresholds*/
    struct minAlnStats conToConMinStats; 
        /*Consensus to consensus mapping thresholds*/
 
    struct conBuildStruct conSet; /*Settings for building a consensus*/

    struct readBin *binTree = 0;   /*Tree of read bins*/
    struct readBin *bestBin = 0;   /*Bin with most similar consensus*/
    struct readBin *clustOn = 0;   /*Bin in list/tree am working on*/
    struct readBin *lastClust = 0; /*Last cluster worked on for a bin*/
    struct readBin *lastBin = 0;   /*For keeping list in order*/
    struct readBin *tmpBin = 0;    /*Pionts to a readBin to work on*/
       
    /******************************************************************\
    * Main Sec-1 Sub-2: Help message
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
            \n    -min-per-reads:                          [0.003=0.3%]\
            \n        - Minimum percentage of reads expected\
            \n          to keep a cluster.\
            \n        - Is the percentage of all reads that\
            \n          were assigned to a cluster.\
            \n    -min-reads-per-bin:\
            \n        - Min number of reads needed to keep   [100]\
            \n          a bin or a cluster\
            \n    -max-reads-per-con:\
            \n        - Max number of reads to use in        [300]\
            \n          a consensus.\
            \n    -min-con-length:                               [500]\
            \n        - Minimum length to keep a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
            \n    -min-read-length:                          [600]\
            \n       - Discard reads or consensuses built by\
            \n         the majority consensus step that are\
            \n         under the input length.\
            \n       - If you lower this your should also lower\
            \n         -min-read-read-length &\
            \n         -min-read-con-length\
            \n    -max-read-length:                          [1000]\
            \n       - Discard reads with read lengths over\
            \n         input setting (0 to ignore)\
            \n       - When binning is done, this is applied\
            \n         after the trimming step.\
            \n    -pick-read-with-med-q                      [No]\
            \n       - Picks the read to start consensus\
            \n         building by median Q-score.\
            \n       - Default is by mapping quality to the\
            \n         binned reference.\
            \n       - This setting requires more time to run.\
            \n    -min-median-q:                             [10]\
            \n       - Minimum read median quality score\
            \n         needed to keep a read when binning.\
            \n    -min-mean-q:                               [10]\
            \n       - Minimum read mean quality score\
            \n         needed to keep a read when binning.\
            \n Additional Help messages:\
            \n    -h-build-consensus:\
            \n        - Print paramaters for building the consensus\
            \n    -h-bin:\
            \n        - Print out the parameters for the binning step.\
            \n    -h-read-map:\
            \n        - Print out the parameters for the read mapping\
            \n          step.\
            \n    -h-clust:\
            \n        - Print out the parameters for the clustering\
            \n          step.\
            \n    -h-con:\
            \n        - Print out the parameters for the consensus\
            \n          comparison steps.\
            \n Output:\
            \n    - fastq files: With the reads for each co-infection\
            \n    - fasta files: With consensuses for each co-infection\
            \n    - prefix--log.txt: Log file\
            \n    - prefix--read-counts.tsv:\
            \n        o File with number of reads per kept cluster\
            \n        o Also has the read counts for each discard bin\
            \n          or cluster. This makes this file a good for\
            \n          checking this programs progress. It will often\
            \n          be printing bining read counts while\
            \n          clustering.\
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
            \n Find co-infections can use several different consensuses\
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
            \n        - Minimum length to keep a consensus.\
            \n    -disable-majority-consensus:                   [Yes]\
            \n        - Build a consensus using a simple\
            \n          majority consensus. This consensus\
            \n          will be polished with Racon or\
            \n          Medaka if Racon and Medaka set.\
            \n    -max-reads-per-con:\
            \n        - Max number of reads to use in        [300]\
            \n          a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
            \n    -maj-con-min-bases                         [0.4=40%]\
            \n        - When building the majority consesus\
            \n          make a deletion in positions that\
            \n          have less than x\% of bases (40%).\
            \n    -maj-con-min-base-q:                       [10]\
            \n        - Minimum q-score to keep a base when\
            \n          building a majority consensus.\
            \n    -maj-con-min-ins                           [0.3=30%]\
            \n        - When building the majority consesus\
            \n          ingore insertions that have support\
            \n          from less than x\% of reads (40%).\
            \n    -maj-con-min-ins-q:                        [5]\
            \n        - Minimum q-score to keep a insertion\
            \n          when building a majority consensus.\
            \n    -enable-racon:                                [No]\
            \n        - Do not use racon to polish the\
            \n          consensus.\
            \n    -rounds-racon:\
            \n        - Number of rounds to polish a         [4]\
            \n          consensus with racon\
       "; /*consensus building parameters*/

    char
         *binHelpCStr = "\
            \n The binning step compares a read to the best reference\
            \n    its mapped to. If the read meets the qaulity\
            \n    thresholds then it is assigned to the references bin.\
            \n    The read is discarded if it does not meet the quality\
            \n    thresholds.\
            \n\
            \n -skip-bin:                                [No]\
            \n    - Skip the binning step.\
            \n -rm-sup-reads                             [No]\
            \n    - Removes any read that has a\
            \n      supplemental alignment. These might\
            \n      be chimeras or could just be repeated\
            \n      genes.\
            \n    - If this setting is off, then supplemental\
            \n      alignments are just ignored.\
            \n    - This setting requires running minimap2 with\
            \n      just one thread for binning.\
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
            \n\
            \n -read-ref-min-aligned-median-q:            [10]\
            \n    - Minimum read median quality score\
            \n      of the aligned region of a read\
            \n      needed to keep a read when binning.\
            \n -read-ref-min-aligned-mean-q:              [10]\
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
            \n -read-ref-min-base-q:                      [10]\
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
            \n -read-ref-max-a-del-homo:                  [1]\
            \n    - Maximum A homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-t-del-homo:                  [1]\
            \n    - Maximum T homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-c-del-homo:                  [1]\
            \n    - Maximum C homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            \n -read-ref-max-g-del-homo:                  [1]\
            \n    - Maximum G homopolymer size to keep\
            \n      an deletion (1 = no hompolymer,\
            \n      0 = always discard).\
            "; /*binning parameters help message*/

    char
        *readMapHelpCStr ="\
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
            \n     -read-read-snps:                       [0.015 = 1.5%]\
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
            \n     -read-read-min-base-q:                 [13]\
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
            \n     -read-read-max-a-del-homo:                [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-t-del-homo:                [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-c-del-homo:                [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-read-max-g-del-homo:                [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            "; /*Read mapping help message*/


    char
        *clustHelpCStr = "\
            \n The clustering step uses the consensus built after the\
            \n    read mapping step to identify reads that came from\
            \n    the same strain. Reads are first mapped to the\
            \n    consensus and any read that meets the user thresholds\
            \n    is considered part of the cluster. The discarded\
            \n    reads are saved for another round of clustering.\
            \n\
            \n     -skip-clust:                                [No]\
            \n        - Skip the clusterin step.\
            \n     -min-read-con-map-length:                   [500]\
            \n        - Minimum aligned read length needed to\
            \n          keep a read to consensus mapping.\
            \n        - Reads that have an alinged length less\
            \n          than this will not be extracted.\
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
            \n The scoring settings for clustering only counts SNPs\
            \n     and indels in reads that are at or above the user\
            \n     input thresholds. These counts are used to find if\
            \n     the read is a good quality read (belongs to cluster)\
            \n     or a low quality read (discard).\
            \n\
            \n     -read-con-min-base-q:                 [13]\
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
            \n     -read-con-max-a-del-homo:             [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-t-del-homo:             [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-c-del-homo:             [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -read-con-max-g-del-homo:             [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            "; /*Clustering help message*/

    char
        *conCompHelpCStr = "\
            \n The consensus comparison steps are used to identify\
            \n    potential false positive clusters. Consensus are\
            \n    mapped to each other and any pair of consensus that\
            \n    are to similar (meets user ciriteria) are considered\
            \n    the same. The clusters of consensus considered the\
            \n    same are merged into one cluster.\
            \n\
            \n     -con-con-snps:                        [0.015 = 1.5%]\
            \n        - Minimum percentage of snps needed to\
            \n          merge clusters during the consensus\
            \n          to consensus comparison steps.\
            \n     -con-con-diff:                        [0.023 = 2.3%]\
            \n        - Minimum percent difference needed to\
            \n          merge clusters during the consensus\
            \n          to consensus comparison steps.\
            \n     -con-con-dels:                        [1 = 100%]\
            \n        - Minimum percentage of deletions\
            \n          needed to merge consensuses during\
            \n          the consensus to consensus\
            \n          comparison steps.\
            \n     -con-con-inss:                        [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          needed to merge consensuses during\
            \n          the consensus to consensus\
            \n          comparison steps.\
            \n     -con-con-indels:                      [1 = 100%]\
            \n        - Minimum percentage of insertions\
            \n          and deletions needed to merge \
            \n          consensuses during the consensus to\
            \n          consensus comparison steps.\
            \n\
            \n\
            \n The scoring settings for the consensus comparison steps\
            \n     only count SNPs and indels in consensuses that are\
            \n     at or above the user input thresholds. These counts\
            \n     are used to find if the two clusters should be\
            \n     merged (consensus very similar).\
            \n\
            \n     -con-con-max-a-ins-homo:              [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-t-ins-homo:              [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-c-ins-homo:              [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-g-ins-homo:              [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an insertion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n\
            \n     -con-con-max-a-del-homo:              [1]\
            \n        - Maximum A homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-t-del-homo:              [1]\
            \n        - Maximum T homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-c-del-homo:              [1]\
            \n        - Maximum C homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            \n     -con-con-max-g-del-homo:              [1]\
            \n        - Maximum G homopolymer size to keep\
            \n          an deletion (1 = no hompolymer,\
            \n          0 = always discard).\
            "; /*difference help settings*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-2: Set up default settings in structures
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankMinStats(&readToRefMinStats);
    blankMinStatsReadRead(&readToReadMinStats);
    blankMinStatsReadCon(&readToConMinStats);
    blankMinStatsConCon(&conToConMinStats);

    initConBuildStruct(&conSet);          /*default consensus settings*/
    conSet.useStatBl = 1;/*Select read with stats file from scoreReads*/

    initSamEntry(&samStruct);
    initSamEntry(&refStruct);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Get user input and set up commands 
    #    sec-3 sub-1: Get user input & check for errors
    #    sec-3 sub-2: Add prefix to file names
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Main Sec-3 Sub-1: Get user input & check for errors
    *******************************************************************/

    /*Set up default values*/
    strcpy(threadsCStr, defThreads); /*Setup default number of threads*/
    strcpy(prefCStr, defPrefix);

    inutErrCStr =
        getUserInput(
            lenArgsInt,
            argsCStr,
            prefCStr,
            &fqPathCStr,
            &refsPathCStr,
            threadsCStr,
            &rmSupAlnBl,
            &skipBinBl,   /*1: Skip the binning step, 0 do not*/
            &skipClustBl, /*1: Skip the clusterin step, 0 do not*/
            &minReadsDbl,
            &conSet,                /*consensus building settings*/
            &readToRefMinStats,
            &readToReadMinStats,
            &readToConMinStats,
            &conToConMinStats
    ); /*Get user input*/

    if(inutErrCStr != 0)
    { /*If have an error*/
        if(inutErrCStr == 0)
        { /*If no user input was supplied*/
            fprintf(
                stderr,
                "%s\n\nNo user input was supplied\n",
                helpCStr
            );

            exit(1);
        } /*If no user input was supplied*/

        if(
            strcmp(inutErrCStr, "-v") == 0 ||
            strcmp(inutErrCStr, "-V") == 0 ||
            strcmp(inutErrCStr, "-version") == 0
        ) { /*If the user is requesting the version number*/
            fprintf(stdout, "%.8f\n", defVersion);
            exit(0);
        } /*If the user is requesting the version number*/

        if(strcmp(inutErrCStr, "-h-build-consensus") == 0)
        { /*If the user wants the consensus building options*/
            fprintf(stdout, "%s", conBuildHelpCStr);
            exit(0);
        } /*If the user wants the consensus building options*/
        
        if(strcmp(inutErrCStr, "-h-bin") == 0)
        { /*If user wants to know about the binning parameters*/
            fprintf(stdout, "%s", binHelpCStr);
            exit(0);
        } /*If user wants to know about the binning parameters*/

        if(strcmp(inutErrCStr, "-h-read-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/

        if(strcmp(inutErrCStr, "-h-clust") == 0)
        { /*If the user wants to know about the clustering parameters*/
            fprintf(stdout, "%s", clustHelpCStr);
            exit(0);
        } /*If the user wants to know about the clustering parameters*/

        if(strcmp(inutErrCStr, "-h-con") == 0)
        { /*If user wants the consensus comparison parameters*/
            fprintf(stdout, "%s", conCompHelpCStr);
            exit(0);
        } /*If user wants the consensus comparison parameters*/

        if(
            strcmp(inutErrCStr, "-h") == 0 ||
            strcmp(inutErrCStr, "-help") == 0 ||
            strcmp(inutErrCStr, "--h") == 0 ||
            strcmp(inutErrCStr, "--help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);/*Print out help message*/
            exit(0);
        } /*If user wanted the help message*/

        fprintf(
            stderr,
            "%s\n\n %s is an invalid parameter\n",
            helpCStr,   /*Print out the help message*/
            inutErrCStr     /*Print out the error*/
        ); /*Let user know about the invalid parameter*/

        exit(1);
    } /*If have an error*/

    if(!(conSet.majConSet.useMajConBl |
         conSet.raconSet.useRaconBl |
         conSet.medakaSet.useMedakaBl
    )) { /*If the user said to ingore all consensus building steps*/
        printf("Current settings have turned off all consensus");
        printf(" building methods.\nSelect a consensus step by");
        printf(" removing: -enable-medaka or -enable-racon\n");

        exit(1);
    } /*If the user said to ingore all consensus building steps*/

    if((skipBinBl & skipClustBl) & 1)
    { /*If skipping both the binning and clustering step*/
        printf("Binning and clustering steps has been turned of, so");
        printf(" nothing can be done. Please enable the binning step");
        printf(" or clustering step\n");
        exit(1);
    } /*If skipping both the binning and clustering step*/

    /*******************************************************************
    # Main Sec-3 Sub-2: Add prefix to file names
    *******************************************************************/

    tmpCStr = cStrCpInvsDelm(logFileCStr, prefCStr);
    strcpy(tmpCStr , "--log.txt"); /*Finsh log name*/

    tmpCStr = cStrCpInvsDelm(readCntFileCStr, prefCStr);
    strcpy(tmpCStr, "--read-counts.tsv");

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-4: Open files and check user input
    ^    main sec-4 sub-1: Check if can open the log file
    ^    main sec-4 sub-2: Check if minimap2 exists
    ^    main sec-4 sub-3: Check if racon exists
    ^    main sec-4 sub-4: Check if medaka exists
    ^    main sec-4 sub-5: Check if the fastq file of reads exists
    ^    main sec-4 sub-6: Check if the reference fasta file extists
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Main Sec-4 Sub-1: Check if can open the log file
    *******************************************************************/

    logFILE = fopen(logFileCStr, "a"); /*Test if can open the log file*/

    if(logFILE == 0)
    { /*If could not set up the log file*/
        fprintf(
            stderr,
            "Could not open the log file (%s)\n",
            logFileCStr
        ); /*Let user know that the log file was not opened*/

        exit(1);
    } /*If could not set up the log file*/

    /*******************************************************************
    # Main Sec-4 Sub-2: Check if minimap2 exists
    *******************************************************************/

    /*Set up minimap 2 check*/
    tmpCStr = cpParmAndArg(tmpCmdCStr, "minimap2", "--version");
    tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);

    fprintf(logFILE, "minimap2 version: ");
    fclose(logFILE); /*Closing to avoid system appending to open file*/
    logFILE = 0;

    if(system(tmpCmdCStr) != 0)
    { /*If minimap2 does not exist*/
        fprintf(stderr, "Minimap2 could not be found\n");
        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Minimap2 could not be found\n");
        fclose(logFILE);

        exit(1);
    } /*If minimap2 does not exist*/

    /*******************************************************************
    # Main Sec-4 Sub-3: Check if racon exists
    *******************************************************************/

    if(conSet.raconSet.useRaconBl & 1)
    { /*If using racon, get the version used*/
        /*Set up racon check*/
        tmpCStr = cpParmAndArg(tmpCmdCStr, "racon", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);

        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Racon version: ");
        fclose(logFILE); /*Closing to avoid system appending to file*/
        logFILE = 0;

        if(system(tmpCmdCStr) != 0)
        { /*If racon does not exist*/
            fprintf(stderr, "Racon could not be found\n");
            logFILE = fopen(logFileCStr, "a");
            fprintf(logFILE, "Racon could not be found\n");
            fclose(logFILE);

            exit(1);
        } /*If racon does not exist*/
    } /*If using racon, get the version used*/

    else
    { /*Else if not using racon*/
        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Not using Racon\n");
        fclose(logFILE);
        logFILE = 0;
    } /*Else if not using racon*/

    /******************************************************************\
    * Main Sec-4 Sub-4: Check if medaka exists
    \******************************************************************/

    if(conSet.medakaSet.useMedakaBl & 1)
    { /*If using medaka, check version*/
        /*Set up non-miniconda command*/
        tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medakaCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);
        tmpCStr = cpParmAndArg(tmpCStr, medakaCMDEnd, "");

        logFILE = fopen( logFileCStr, "a");
        fprintf(logFILE, "Medaka version: ");
        fclose(logFILE);
        logFILE = 0;

        if(system(tmpCmdCStr) != 0)
        { /*If python virtual enviorment medaka does not exist*/
            /*Set up miniconda medaka enviroment command*/
            tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medCondCMD);
            tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
            tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);
            tmpCStr = cStrCpInvsDelm(tmpCStr, medCondCMDEnd);
           
            if(system(tmpCmdCStr) != 0)
            { /*If medaka could not be found*/
                fprintf(stderr, "Medaka could not be found\n");
                logFILE = fopen(logFileCStr, "a");
                fprintf(logFILE, "Racon could not be found\n");
                fclose(logFILE);

                exit(1);
            } /*If medaka could not be found*/

            logFILE = fopen(logFileCStr, "a");

            fprintf(
                logFILE,
                "    - Using medaka installed by miniconda\n"
            );
            conSet.medakaSet.condaBl = 1; /*Using medaka from conda*/

            fclose(logFILE);
            logFILE = 0;
        } /*If python virtual enviorment medaka does not exist*/

        else
        { /*Else if found the python virtual enviorment medaka*/
            logFILE = fopen(logFileCStr, "a");

            conSet.medakaSet.condaBl = 0; /*Installed by python env*/
            fprintf(logFILE, "    - Using medaka installed by python\n");

            fclose(logFILE);
            logFILE = 0;
        } /*Else if found the python virtual enviorment medaka*/
    } /*If using medaka, check version*/

    else
    { /*Else if not using medaka*/
        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Not using Medaka\n");
        fclose(logFILE);
        logFILE = 0;
    } /*Else if not using medaka*/

    /*******************************************************************
    # Main Sec-4 Sub-5: Check if the fastq file of reads extists
    *******************************************************************/

    logFILE = fopen(logFileCStr, "a"); /*So can print out the error*/

    fprintf(logFILE, "findCoInft version: %f\n", defVersion);
    fprintf(logFILE, "\n\nfindCoInft settings:\n");
    fprintf(logFILE, "  findCoInft \\\n");

    stdinFILE = fopen(fqPathCStr, "r"); /*Check if can open fastq file*/

    if(stdinFILE == 0)
    { /*If no fastq file was provided*/
        fprintf(
            stderr,
            "The provided fastq file (%s) does not exist\n",
            fqPathCStr
        ); /*Let user know the fastq file does not exist*/

        fprintf(
            logFILE,
            "File provided by -fastq (%s) does not exist\n",
            fqPathCStr
        ); /*Print error to log*/

        fclose(logFILE);
        exit(1);
    } /*If no fastq file was provided*/

    fclose(stdinFILE); /*No longer need open*/

    /*Print out file used for user*/
    fprintf(logFILE, "    -fastq %s \\\n", fqPathCStr);

    /*******************************************************************
    # Main Sec-4 Sub-5: Check if the reference fasta file extists
    *******************************************************************/

    if(!(skipBinBl & 1))
    { /*If binning reads*/
        stdinFILE = fopen(refsPathCStr, "r"); /*Open reference file*/

        if(stdinFILE == 0)
        { /*If no fastq file was provided*/
            fprintf(
                stderr,
                "The provided reference file (%s) does not exist\n",
                refsPathCStr
            ); /*Let user know the fastq file does not exist*/

            fprintf(
                logFILE,
                "File provided by -ref (%s) does not exist\n",
                refsPathCStr
            ); /*Print error to log*/

            fclose(logFILE);
            exit(1);
        } /*If no fastq file was provided*/

        fclose(stdinFILE); /*No longer need open*/

        /*Print out file used for user*/
        fprintf(logFILE, "    -ref %s \\\n", refsPathCStr);
    } /*If binning reads*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-5: Put non-required input into log for user
    ^     main sec-5 sub-1: Print out general settings
    ^     main sec-5 sub-2: Print out Comparision settings for binning
    ^     main sec-5 sub-3: Print Comparision settings for read mapping
    ^     main sec-5 sub-4: Print comparision setting, consensus compare
    ^     main sec-5 sub-5: Print our scoring settings for binning
    ^     main sec-5 sub-6: Print out scoring settings for read mapping
    ^     main sec-5 sub-7: Print out scoring settings for clustering
    ^     main sec-5 sub-8: Print scoring settings for consensus compare
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-5 Sub-1: Print out general settings
    \******************************************************************/

    /*Print out the number of threads input*/
    fprintf(logFILE, "    -threads %s \\\n", threadsCStr);

    /*Print out the general settings*/
    fprintf(logFILE, "    -prefix %s \\\n", prefCStr);

    fprintf(
        logFILE,
        "    -min-reads-per-bin %ju \\\n",
        conSet.minReadsToBuildConUL
    );

    fprintf(
        logFILE,
        "    -max-reads-per-con %ju \\\n",
        conSet.maxReadsToBuildConUL
    );

    fprintf(
        logFILE,
        "    -extra-consensus-steps %u \\\n",
        conSet.numRndsToPolishUI
    );

    fprintf(logFILE, "    -min-perc-reads %f \\\n", minReadsDbl);

    if(skipBinBl & 1)
        fprintf(logFILE, "    -skip-bin \\\n");

    if(skipClustBl & 1)
        fprintf(logFILE, "    -skip-clust \\\n");

    if(rmSupAlnBl & 1)
        fprintf(logFILE, "    -rm-sup-reads \\\n");

    if(!(conSet.useStatBl & 1) || skipBinBl & 1)
    { /*If using the median Q-score*/
        fprintf(logFILE, "    -pick-read-with-med-q \\\n");

        /*If user is skipping binning, make sure not using stat file*/
        if(conSet.useStatBl & 1)
            conSet.useStatBl = 0; 
    } /*If using the median Q-score*/

    if(conSet.majConSet.useMajConBl & 1)
    { /*If using the majority consensus  step*/
        fprintf(logFILE, "    -disable-majority-consensus \\\n");

        fprintf(
            logFILE,
            "    -maj-con-min-bases %f \\\n",
            conSet.majConSet.minReadsPercBaseFlt
        ); /*min percentage of read support to keep a match/SNP*/

        fprintf(
            logFILE,
            "    -maj-con-min-base-q %u \\\n",
            conSet.majConSet.minBaseQUC
        );

        fprintf(
            logFILE,
            "    -maj-con-min-ins %f \\\n",
            conSet.majConSet.minReadsPercInsFlt
        ); /*Min % of read support needed to keep an insertion*/

        fprintf(
            logFILE,
            "    -maj-con-min-ins-q %u \\\n",
            conSet.majConSet.minInsQUC
        );
    } /*If using the majority consensus  step*/

    if(conSet.raconSet.useRaconBl & 1)
    { /*If using racon*/
        fprintf(logFILE, "    -enable-racon \\\n");

        fprintf(
           logFILE,
           "    -rounds-racon %u \\\n",
           conSet.raconSet.rndsRaconUC
        );
    } /*If using racon*/

    if(conSet.medakaSet.useMedakaBl & 1)
    { /*If using medaka*/
        fprintf(logFILE, "    -enable-medaka \\\n");
       fprintf(logFILE,"    -model %s \\\n",conSet.medakaSet.modelCStr);
    } /*If using medaka*/

    /******************************************************************\
    * Main Sec-5 Sub-2: Print out Comparision settings for binning
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-ref-snps %f \\\n",
        
        readToRefMinStats.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-diff %f \\\n",
        
        readToRefMinStats.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-dels %f \\\n",
        
        readToRefMinStats.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-inss %f \\\n",
        
        readToRefMinStats.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-indels %f \\\n",
        
        readToRefMinStats.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-3: Print out Comparision settings for read mapping
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-read-snps %f \\\n",
        
        readToReadMinStats.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-read-diff %f \\\n",
        
        readToReadMinStats.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-read-dels %f \\\n",
        
        readToReadMinStats.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-read-inss %f \\\n",
        
        readToReadMinStats.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-read-indels %f \\\n",
        
        readToReadMinStats.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-3: Print out Comparision settings for clustering
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-con-snps %f \\\n",
        
        readToConMinStats.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-con-diff %f \\\n",
        
        readToConMinStats.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-con-dels %f \\\n",
        
        readToConMinStats.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-con-inss %f \\\n",
        
        readToConMinStats.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-con-indels %f \\\n",
        
        readToConMinStats.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-4: Print Comparision settings for consensus compare
    \******************************************************************/

    fprintf(
        logFILE,
        "    -con-con-snps %f \\\n",
        
        conToConMinStats.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -con-con-diff %f \\\n",
        
        conToConMinStats.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -con-con-dels %f \\\n",
        
        conToConMinStats.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -con-con-inss %f \\\n",
        
        conToConMinStats.minInssFlt
    );

    fprintf(
        logFILE,
        "    -con-con-indels %f \\\n",
        
        conToConMinStats.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-5: Print our scoring settings for binning
    \******************************************************************/

     fprintf(
         logFILE,
         "    -read-ref-min-base-q %u \\\n",
         readToRefMinStats.minQChar
     );

     fprintf(
         logFILE,
         "    -read-ref-min-mapq %u \\\n",
         readToRefMinStats.minMapqUInt
     );

     fprintf(
         logFILE,
         "    -min-median-q %f \\\n",
         readToRefMinStats.minMedianQFlt
     );

     fprintf(
         logFILE,
         "    -min-mean-q %f \\\n",
         readToRefMinStats.minMeanQFlt
     );

     fprintf(
         logFILE,
         "    -read-ref-min-aligned-median-q %f \\\n",
         readToRefMinStats.minAlignedMedianQFlt
     );

     fprintf(
         logFILE,
         "    -read-ref-min-aligned-mean-q %f \\\n",
         readToRefMinStats.minAlignedMedianQFlt
     );

     fprintf(
         logFILE,
         "    -min-read-length %u \\\n",
         readToRefMinStats.minReadLenULng
     );

     fprintf(
         logFILE,
         "    -max-read-length %u \\\n",
         readToRefMinStats.maxReadLenULng
     );

     fprintf(
         logFILE,
         "    -read-ref-max-a-ins-homo %u \\\n",
         readToRefMinStats.maxHomoInsAry[0]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-t-ins-homo %u \\\n",
         readToRefMinStats.maxHomoInsAry[10]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-c-ins-homo %u \\\n",
         readToRefMinStats.maxHomoInsAry[1]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-g-ins-homo %u \\\n",
         readToRefMinStats.maxHomoInsAry[3]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-a-del-homo %u \\\n",
         readToRefMinStats.maxHomoDelAry[0]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-t-del-homo %u \\\n",
         readToRefMinStats.maxHomoDelAry[10]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-c-del-homo %u \\\n",
         readToRefMinStats.maxHomoDelAry[1]
     );

     fprintf(
         logFILE,
         "    -read-ref-max-g-del-homo %u \\\n",
         readToRefMinStats.maxHomoDelAry[3]
     );

    /******************************************************************\
    * Main Sec-5 Sub-6: Print out scoring settings for read mapping
    \******************************************************************/

     fprintf(
         logFILE,
         "    -read-read-min-base-q %u \\\n",
         readToReadMinStats.minQChar
     );

     fprintf(
         logFILE,
         "    -read-read-min-mapqq %u \\\n",
         readToReadMinStats.minMapqUInt
     );

     fprintf(
         logFILE,
         "    -read-read-max-a-ins-homo %u \\\n",
         readToReadMinStats.maxHomoInsAry[0]
     );

     fprintf(
         logFILE,
         "    -read-read-max-t-ins-homo %u \\\n",
         readToReadMinStats.maxHomoInsAry[10]
     );

     fprintf(
         logFILE,
         "    -read-read-max-c-ins-homo %u \\\n",
         readToReadMinStats.maxHomoInsAry[1]
     );

     fprintf(
         logFILE,
         "    -read-read-max-g-ins-homo %u \\\n",
         readToReadMinStats.maxHomoInsAry[3]
     );

     fprintf(
         logFILE,
         "    -read-read-max-a-del-homo %u \\\n",
         readToReadMinStats.maxHomoDelAry[0]
     );

     fprintf(
         logFILE,
         "    -read-read-max-t-del-homo %u \\\n",
         readToReadMinStats.maxHomoDelAry[10]
     );

     fprintf(
         logFILE,
         "    -read-read-max-c-del-homo %u \\\n",
         readToReadMinStats.maxHomoDelAry[1]
     );

     fprintf(
         logFILE,
         "    -read-read-max-g-del-homo %u \\\n",
         readToReadMinStats.maxHomoDelAry[3]
     );

    /******************************************************************\
    * Main Sec-5 Sub-7: Print out scoring settings for clustering
    \******************************************************************/

     fprintf(
         logFILE,
         "    -read-con-min-base-q %u \\\n",
         readToConMinStats.minQChar
     );

     fprintf(
         logFILE,
         "    -read-con-min-mapqq %u \\\n",
         readToConMinStats.minMapqUInt
     );

     fprintf(
         logFILE,
         "    -read-con-max-a-ins-homo %u \\\n",
         readToConMinStats.maxHomoInsAry[0]
     );

     fprintf(
         logFILE,
         "    -read-con-max-t-ins-homo %u \\\n",
         readToConMinStats.maxHomoInsAry[10]
     );

     fprintf(
         logFILE,
         "    -read-con-max-c-ins-homo %u \\\n",
         readToConMinStats.maxHomoInsAry[1]
     );

     fprintf(
         logFILE,
         "    -read-con-max-g-ins-homo %u \\\n",
         readToConMinStats.maxHomoInsAry[3]
     );

     fprintf(
         logFILE,
         "    -read-con-max-a-del-homo %u \\\n",
         readToConMinStats.maxHomoDelAry[0]
     );

     fprintf(
         logFILE,
         "    -read-con-max-t-del-homo %u \\\n",
         readToConMinStats.maxHomoDelAry[10]
     );

     fprintf(
         logFILE,
         "    -read-con-max-c-del-homo %u \\\n",
         readToConMinStats.maxHomoDelAry[1]
     );

     fprintf(
         logFILE,
         "    -read-con-max-g-del-homo %u \\\n",
         readToConMinStats.maxHomoDelAry[3]
     );

    /******************************************************************\
    * Main Sec-5 Sub-8: Print out scoring settings for consensus compare
    \******************************************************************/

     fprintf(
         logFILE,
         "    -con-con-max-a-ins-homo %u \\\n",
         conToConMinStats.maxHomoInsAry[0]
     );

     fprintf(
         logFILE,
         "    -con-con-max-t-ins-homo %u \\\n",
         conToConMinStats.maxHomoInsAry[10]
     );

     fprintf(
         logFILE,
         "    -con-con-max-c-ins-homo %u \\\n",
         conToConMinStats.maxHomoInsAry[1]
     );

     fprintf(
         logFILE,
         "    -con-con-max-g-ins-homo %u \\\n",
         conToConMinStats.maxHomoInsAry[3]
     );

     fprintf(
         logFILE,
         "    -con-con-max-a-del-homo %u \\\n",
         conToConMinStats.maxHomoDelAry[0]
     );

     fprintf(
         logFILE,
         "    -con-con-max-t-del-homo %u \\\n",
         conToConMinStats.maxHomoDelAry[10]
     );

     fprintf(
         logFILE,
         "    -con-con-max-c-del-homo %u \\\n",
         conToConMinStats.maxHomoDelAry[1]
     );

     fprintf(
         logFILE,
         "    -con-con-max-g-del-homo %u;\n",
         conToConMinStats.maxHomoDelAry[3]
     );

    fclose(logFILE); /*Flush output to the log file*/
    logFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Find initial bins with references
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(!(skipBinBl & 1))
    { /*If binning reads*/
        binTree =
            binReads(
                fqPathCStr,        /*Fastq file to bin*/
                refsPathCStr,      /*References to bin with*/
                prefCStr,          /*prefix to name all bins with*/
                threadsCStr,      /*Number of threads to use with minimap2*/
                rmSupAlnBl,         /*Remove supplementary alignments*/
                1,                  /*1: trim reads, 0: do not*/
                &samStruct,
                &refStruct,
                &readToRefMinStats,
                &errUC              /*Reports any errors*/
        );

        if(binTree == 0)
        { /*If had an error*/
            logFILE = fopen(logFileCStr, "a");

            if(errUC & 64)
            { /*If the binning step errored out*/
                fprintf(stdout, "Memory error: not enough memory\n");
                fprintf(logFILE, "Binning step ran out of memory\n");
            } /*If the binning step errored out*/

            else if(errUC & 2)
            { /*If the fastq file could not be opened*/
                fprintf(
                    stdout,
                    "Fastq file (%s) could not be opened\n",
                    fqPathCStr
                ); /*Let user know fastq file is invalid*/
                fprintf(
                    logFILE,
                    "Binning: Fastq file (%s) could not be opened\n",
                    fqPathCStr
                ); /*Let user know fastq file is invalid*/
            } /*If the fastq file could not be opened*/

            else
            { /*Else could not create a file*/
                fprintf(stdout, "Could not create fastq/stats files\n");
                fprintf(logFILE, "Binning: Could not create files\n");
            } /*Else could not create a file*/

            freeStackSamEntry(&samStruct);
            freeStackSamEntry(&refStruct);

            if(statFILE != 0)
                fclose(statFILE);
            if(logFILE != 0)
                fclose(logFILE);

            exit(1);
        } /*If had an error*/
    } /*If binning reads*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Main Sec-7: Non-reference based binning steps
    ^    main sec-7 sub-1: Set up for clustering or consnsus buidling
    ^    main sec-7 sub-2: Print read counts for bins & decide if keep
    ^    main sec-7 sub-3: Set up readBin to hold a clusters files
    ^    main sec-7 sub-4: Build the consensus
    ^    main sec-7 sub-5: Copy best read back into to its orignal bin
    ^    main sec-7 sub-6: Bin reads to the consensus
    ^    main sec-7 sub-7: Compare new consensus to old consensuses
    ^    main sec-7 sub-8: Update list of clusters in bin
    ^    main sec-7 sub-9: Add clusters to bin & move to next bin
    ^        - Before going in the bin list is an balanced tree
    ^        - after this the bins are in a list, with the left pointer
    ^          pointing towards the bins & the right pointer pointing
    ^          towards clusters (unless user did not want clustering)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-7 Sub-1: Set up for clustering or consensus building
    \******************************************************************/

    if(skipBinBl & 1)
    { /*If skipping the binning step, then need to create a single bin*/
        binTree = malloc(sizeof(struct readBin));

        /*Blank the variables in binTree*/
        blankReadBin(binTree);
        strcpy(binTree->refIdCStr, "clust");
        binTree->balUChar = 1; /*To mark keeping*/

        /*Make a copy of the fastq file to protect the original*/
        tmpCStr = cStrCpInvsDelm(binTree->fqPathCStr, prefCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--clust.fastq");

        filterReads(
            fqPathCStr,
            binTree->fqPathCStr,
            &samStruct,
            &readToRefMinStats
        ); /*Remove low quality reads*/

        /*Get the number of reads in the copied fastq file*/
        binTree->numReadsULng = getNumReadsInFq(binTree->fqPathCStr);
    } /*If skipping the binning step, then need to create a single bin*/

    else
        cnvtBinTreeToList(&binTree);
        /*Convert our bin tree to a list (No longer need AVL tree)*/

    clustOn = binTree;
    lastBin = binTree; /*so I can reset pointers when removing bin*/
    tmpBin = 0;

    /*Blank the rad count file*/
    statFILE = fopen(readCntFileCStr, "w");/*File to recored counts to*/
    fprintf(statFILE, "\nBins\t*\t0\t*\t*\n");

    /******************************************************************\
    * Main Sec-7 Sub-2: Print out read counts for bins & decide if keep
    \******************************************************************/

    while(clustOn != 0)
    { /*While have reads to cluster*/
        conSet.clustUC = 0;
        lastClust = clustOn; /*Head of cluster list*/

        fprintf(
            statFILE,
            "%s\t%lu",
            clustOn->refIdCStr,
            clustOn->numReadsULng
        ); /*Recording the number of reads in the bin*/

        if(clustOn->numReadsULng < conSet.minReadsToBuildConUL)
        { /*if have to remove a bin*/
            if(binTree == clustOn)
            { /*If removing bin removes the head of the list*/
                rmBinFromList(&binTree); /*Remove head bin from list*/
                clustOn = binTree;
                lastBin = binTree;
            } /*If removing bin removes the head of the list*/

            else
            { /*Else needo to remove the cluster from the list*/
                lastBin->leftChild = clustOn->leftChild;
                rmBinFromList(&clustOn); /*Remove the bin from list*/
            } /*Else needo to remove the cluster from the list*/
           
            /*Let user know why bin was not kept*/
            fprintf(statFILE, "\tremoved\tfirst-binning\n");
            fflush(statFILE); /*Make sure io printed out*/
            continue;
        } /*if have to remove a bin*/

        fprintf(statFILE, "\tkept\tfirst-binning\n");
        fflush(statFILE); /*make sure io printed out*/

        /**************************************************************\
        * Main Sec-7 Sub-3: Set up readBin to hold a clusters files
        \**************************************************************/
 
        while(clustOn->numReadsULng >= conSet.minReadsToBuildConUL)
        { /*While have reads to bin*/
            errUC = 0; /*reset*/

            if(tmpBin == 0)
                tmpBin = malloc(sizeof(struct readBin));

            if(tmpBin != 0)
            { /*If I need to blank the bin*/
                blankReadBin(tmpBin);
                tmpBin->balUChar = 1; /*To mark keeping*/
            } /*If I need to blank the bin*/

            else
            { /*Else memory allocation error*/
                logFILE = fopen(logFileCStr, "a");

                fprintf(
                    stderr,
                    "Memory error in non-ref binning (main sec-7)\n"
                ); /*Let user know about memory issue*/

                fprintf(
                    logFILE,
                    "Memory error in non-ref binning (main sec-7)\n"
                );

                fclose(logFILE);
                freeStackSamEntry(&samStruct);
                freeStackSamEntry(&refStruct);

                freeBinTree(&binTree);

                exit(1); 
            } /*Else memory allocation error*/

            /**********************************************************\
            * Main Sec-7 Sub-4: Buld the consensus
            \**********************************************************/

             errUC  = 
                 buildCon(
                     clustOn,
                     0,        /*Path to fasta file with reference*/
                     threadsCStr, /*# threads to use with system calls*/
                     &conSet,     /*settings for building a consensus*/
                     &samStruct,  /*Will hold sam file data*/
                     &refStruct,  /*For read median Q extraction*/
                     &readToReadMinStats,
                     &readToConMinStats
            ); /*Builds a consensus using fastq file & best read*/

            if(errUC & 16)
                continue;
                /*Unable to build consensus, let loop terminate*/

            if(errUC & 64)
            { /*If had a memory allocation error*/
                logFILE = fopen(logFileCStr, "a");

                fprintf(
                    stderr,
                    "Memory error in cluster step (main sec-7 sub-2)\n"
                ); /*Let user know about memory issue*/

                fprintf(
                    logFILE,
                    "Memory error in cluster step (main sec-7 sub-2)\n"
                );

                fclose(logFILE);
                freeStackSamEntry(&samStruct);
                freeStackSamEntry(&refStruct);

                freeBinTree(&binTree);
                exit(1); 
            } /*If had a memory allocation error*/

            /**********************************************************\
            * Main Sec-7 Sub-5: Copy best read back to its oringal bin
            * buildCon already did this
            \**********************************************************/

            /*++clustOn->numReadsULng;*/ /*Account for the best read*/
            /*tmpFILE = fopen(clustOn->bestReadCStr, "r");
            fqBinFILE = fopen(clustOn->fqPathCStr, "a");

            while(fgets(
                    refStruct.samEntryCStr,
                    refStruct.lenBuffULng,
                    tmpFILE
                )*/ /*Read in file line*/
            /*) fprintf(tmpFILE, "%s", refStruct.samEntryCStr);
    
            fclose(tmpFILE);
            fclose(fqBinFILE);

            remove(clustOn->bestReadCStr);*/
                 /*Remove best read consensus*/

            /**********************************************************\
            * Main Sec-7 Sub-6: Bin reads to the consensus
            \**********************************************************/

            /*Copy the consensus name to the clusters bin*/
            strcpy(tmpBin->consensusCStr, clustOn->consensusCStr);

            if(skipClustBl & 1)
            { /*If not clustering, move onto the next bin*/
                lastClust->rightChild = tmpBin;
                /*So can merge bins durning consensus comparisions*/
                tmpBin->numReadsULng = lastClust->numReadsULng;
                totalKeptReadsUL += tmpBin->numReadsULng;

                strcpy(tmpBin->fqPathCStr, lastClust->fqPathCStr);
                lastClust->fqPathCStr[0] = '\0';/*avoid deleting atEnd*/
                tmpBin = 0;
                break;        /*If not clusterin, move to next bin*/
            } /*If not clustering, move onto the next bin*/

            binReadToCon(
                &conSet.clustUC,    /*Cluster on*/
                clustOn,            /*Bin working on*/
                tmpBin,             /*Bin to hold the cluster*/
                &samStruct,         /*To hold temporary input*/
                &readToConMinStats, /*Settings to keep read to con*/
                threadsCStr         /*# threads to use with Minimap2*/
            ); /*Find reads that mapp to the consensus*/

            /*Find how many reads were kept in clustering*/
            totalKeptReadsUL += tmpBin->numReadsULng;

            /**********************************************************\
            * Main Sec-7 Sub-7: Compare new consensus to old consensuses
            \**********************************************************/

            *clustOn->consensusCStr = '\0';

            bestBin =
                cmpCons(
                    tmpBin,       /*Bin with consensus to compare*/
                    clustOn,     /*Other clusters to comapre to*/
                    &samStruct,  /*Struct to hold input from minimap2*/
                    &refStruct,  /*Struct to hold input from minimap2*/
                    &conToConMinStats, /*Cons to consensus thresholds*/
                    threadsCStr  /*Number threads to use with Minimap2*/
            ); /*Compares a consenses to other consensuses*/

            /**********************************************************\
            * Main Sec-7 Sub-8: Update list of clusters in bin
            \**********************************************************/

            if(bestBin != 0)
            { /*If the consensuses are to similar (the same?)*/
                mergeBins(bestBin, tmpBin);
                continue; /*This cluster is not worth keeping*/
            } /*If the consensuses are to similar (the same?)*/

            /*else add this bin to the end of the list*/
            lastClust->rightChild = tmpBin;
            lastClust = tmpBin;
            tmpBin = 0;
            ++conSet.clustUC;
        } /*While have reads to bin*/

        /**************************************************************\
        * Main Sec-7 Sub-9: Add clusters to bin & move to next bin
        \**************************************************************/

        lastBin = clustOn; /*For reording the list*/
        clustOn = clustOn->leftChild;
    } /*While have reads to cluster*/

    if(tmpBin != 0)
        freeReadBin(&tmpBin); /*Make sure no loose ends*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Main Sec-8: Compare all consensus to remove false positives
    ^    main sec-8 sub-1: Remove empty bins (no clusters) from list
    ^    - Also remove uneeded files
    ^    main sec-8 sub-2: Find most similar consensus to the cluster
    ^    main sec-8 sub-3: If consensuses are to similar, merge clusters
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-8 Sub-1: Remove empty bins (no clusters) from list
    *    - Also remove uneeded files
    \******************************************************************/

    fprintf(statFILE, "\nBins-and-clusters\t*\t0\t*\t*\n");
    clustOn = binTree;
    lastBin = binTree;

    while(clustOn != 0)
    { /*While have bins to check*/
        if(clustOn->rightChild != 0)
        { /*If is a bin I want to keep*/
            tmpBin = clustOn->rightChild;

            while(tmpBin != 0)
            { /*While have clusters to print out*/
                fprintf(
                    statFILE,
                    "%s\t%lu\tkept\tclustering\n",
                    tmpBin->consensusCStr,
                    tmpBin->numReadsULng
                ); /*While have clusters to print out*/

                fflush(statFILE);          /*make sure io printed out*/
                tmpBin = tmpBin->rightChild; /*Move to next cluster*/
            } /*While have clusters to print out*/

            lastBin = clustOn;           /*Last bin was good*/
            clustOn = clustOn->leftChild;
            continue;
        } /*If is a bin I want to keep*/

        fprintf(
            statFILE,
            "%s\t%lu\tdiscared\tno-clusters\n",
            clustOn->refIdCStr,
            clustOn->numReadsULng
        ); /*Recored that the bin was discarded*/

        if(binTree == clustOn)
        { /*If removing bin removes the head of the list*/
            rmBinFromList(&binTree); /*Remove head bin from list*/
            clustOn = binTree;
            lastBin = binTree;
        } /*If removing bin removes the head of the list*/

        else
        { /*Else needo to remove the cluster from the list*/
            lastBin->leftChild = clustOn->leftChild;
            rmBinFromList(&clustOn); /*Remove the bin from list*/
        } /*Else needo to remove the cluster from the list*/
    } /*While have bins to check*/

    /******************************************************************\
    * Main Sec-8 Sub-2: Find the most similar consensus to the cluster
    \******************************************************************/

    /*Just rerun this step if just did clustering, should end pretty
      quickly*/
    clustOn = binTree;
    lastClust = 0;     /*marks when have to leave consensus*/

    while(clustOn != 0)
    { /*While have bins to compare*/
        tmpBin = clustOn->rightChild; /*First cluster in bin*/

        while(tmpBin != 0)
        { /*While have clusters to compare to other clusters*/
            if(tmpBin->balUChar < 0)
            { /*If have already merged this bin*/
                tmpBin = tmpBin->rightChild;
                continue; /*Already merged this bin*/
            } /*If have already merged this bin*/

            if(tmpBin->numReadsULng < conSet.minReadsToBuildConUL)
            { /*If not enough reads to keep*/
                binDeleteFiles(tmpBin); /*Remove its files*/
                tmpBin->balUChar = -1;
                tmpBin = tmpBin->rightChild;
                continue;
            } /*If not enough reads to keep*/

            if((double) tmpBin->numReadsULng / (double) totalKeptReadsUL
               < minReadsDbl
            ) { /*If discarding the bin*/
                binDeleteFiles(tmpBin); /*Remove its files*/
                tmpBin->balUChar = -1; /*mark for removal*/
                tmpBin = tmpBin->rightChild;
                continue;
            } /*If discarding the bin*/

            bestBin =
                cmpCons(
                    tmpBin,         /*Consensus to check*/
                    clustOn->leftChild, /*Other bins & their clusters*/
                    &samStruct,   /*Struct to hold input from minimap2*/
                    &refStruct,   /*Struct to hold input from minimap2*/
                    &conToConMinStats,  /*Cons to consensus thresholds*/
                    threadsCStr  /*Number threads to use with Minimap2*/
            ); /*Compares a consenses to other consensuses*/

            /**********************************************************\
            * Main Sec-8 Sub-3: If consensuses to similar, mergeClusters
            \**********************************************************/

            while(bestBin != 0)
            { /*While have clusters with highly similar consensuses*/
                if(tmpBin->numReadsULng >= bestBin->numReadsULng)
                { /*If the current cluster has more reads*/
                    mergeBins(tmpBin, bestBin);
                    bestBin->balUChar = -1;
                } /*If the current cluster has more reads*/

                else
                { /*else, the best bin has more reads*/
                    mergeBins(bestBin, tmpBin);
                    tmpBin->balUChar = -1;
                    break; /*Will hit the best bin later*/
                } /*else, the best bin has more reads*/

                /*Restart search (no idea about best bin cluster)*/
                bestBin =
                    cmpCons(
                        tmpBin,         /*Consensus to check*/
                        clustOn->leftChild, /*Other bins*/
                        &samStruct,         /*hold minimap2 output*/
                        &refStruct,         /*hold input from minimap2*/
                        &conToConMinStats,  /*min thresholds*/
                        threadsCStr     /*Number threads with Minimap2*/
                ); /*Compares a consenses to other consensuses*/
            } /*While have clusters with highly similar consensuses*/

           tmpBin = tmpBin->rightChild; /*Move to next cluster*/
                /*First node is old bin, so is not a cluster*/
        } /*While have clusters to compare to other clusters*/

        clustOn = clustOn->leftChild; /*Move to next set of clusters*/
    } /*While have bins to compare*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-9: Build longest read consensuses
    ^    main sec-9 sub-1: print out final read counts
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-9 Sub-1: print out final read counts
    \******************************************************************/

    /*Print out final read counts*/
    fprintf(statFILE, "\nFinal-clusters\t*\t0\t*\t*\n");
    clustOn = binTree;

    while(clustOn != 0)
    { /*While have bins to check*/
        tmpBin = clustOn->rightChild;

        while(tmpBin != 0)
        { /*While have cluster stats to print out*/
            if(tmpBin->balUChar > -1)
            { /*If kept the cluster*/
                fprintf(
                    statFILE,
                    "%s\t%lu\tkept\tfinal-check\n",
                    tmpBin->consensusCStr,
                    tmpBin->numReadsULng
                ); /*Recored that the bin was discarded*/
            } /*If kept the cluster*/

            tmpBin = tmpBin->rightChild;
        } /*While have cluster stats to print out*/

        clustOn = clustOn->leftChild;
    } /*While have bins to check*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-10: Clean up and exit
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(logFILE != 0)
        fclose(logFILE);
    if(statFILE != 0)
        fclose(statFILE);
    freeStackSamEntry(&samStruct);
    freeStackSamEntry(&refStruct);

    clustOn = binTree;

    while(clustOn != 0)
    { /*While have bins to free*/
        tmpBin = clustOn->rightChild;
        lastBin = clustOn;
        clustOn = clustOn->leftChild;
        lastBin->consensusCStr[0] ='\0'; /*So cluster consenus kept*/
        binDeleteFiles(lastBin);         /*Remove uneeded files*/
        freeReadBin(&lastBin);           /*Free the bin*/

        while(tmpBin != 0)
        { /*While have clusters to free*/
            lastClust = tmpBin;
            tmpBin = tmpBin->rightChild;
            freeReadBin(&lastClust); /*Free list of bins*/
        } /*While have clusters to free*/
    } /*While have bins to free*/

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
    char *prefCStr,  /*Holds user supplied prefix*/
    char **fqPathCStr, /*Holds path to fastq file*/
    char **refsPathCStr, /*Holds path to references*/
    char *threadsCStr, /*Number threads for minimap2 & racon*/
    char *rmSupAlnBl,  /*1: Remove reads with supplemental alignments*/
    char *skipBinBl,   /*1: Skip the binning step, 0 do not*/
    char *skipClustBl, /*1: Skip the clusterin step, 0 do not*/
    double *minReadsDbl,
    struct conBuildStruct *conSet,   /*Settings for consensus building*/
    struct minAlnStats *readToRefMinStats, /*Binning scoring settings*/
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *readToConMinStats, /*Read cluster socring set*/
    struct minAlnStats *conToConMinStats   /*Consensus comparison set*/
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

        else if(strcmp(parmCStr, "-ref") == 0)
            *refsPathCStr = inputCStr;  /*references*/

        else if(strcmp(parmCStr, "-prefix") == 0)
            strcpy(prefCStr, inputCStr);     /*Have prefix to use*/

        else if(strcmp(parmCStr, "-threads") == 0)
            strcpy(threadsCStr, inputCStr);

        else if(strcmp(parmCStr, "-skip-bin") == 0)
        { /*Else if skipping the binning step*/
            *skipBinBl = 1;
            --intArg; /*Account for this being a true or false*/
        } /*Else if skipping the binning step*/

        else if(strcmp(parmCStr, "-skip-clust") == 0)
        { /*Else if skipping the clustering step*/
            *skipClustBl = 1;
            --intArg; /*Account for this being a true or false*/
        } /*Else if skipping the clustering step*/

        else if(strcmp(parmCStr, "-min-perc-reads") == 0)
            sscanf(inputCStr, "%lf", minReadsDbl);

        else if(strcmp(parmCStr, "-pick-read-with-med-q") == 0)
        { /*If user wanted to use the median Q-score instead*/
            conSet->useStatBl = 0;
            --intArg;
        } /*If user wanted to use the median Q-score instead*/

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

        else if(strcmp(parmCStr, "-rm-sup-reads") == 0)
           *rmSupAlnBl = !(*rmSupAlnBl);
           
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

        else if(strcmp(parmCStr, "-read-con-snps") == 0)
            sscanf(inputCStr, "%f", &readToConMinStats->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-con-diff") == 0)
            sscanf(inputCStr, "%f", &readToConMinStats->minDiffFlt);
        else if(strcmp(parmCStr, "-read-con-dels") == 0)
            sscanf(inputCStr, "%f", &readToConMinStats->minDelsFlt);
        else if(strcmp(parmCStr, "-read-con-inss") == 0)
            sscanf(inputCStr, "%f", &readToConMinStats->minInssFlt);
        else if(strcmp(parmCStr, "-read-con-indels") == 0)
            sscanf(inputCStr, "%f", &readToConMinStats->minIndelsFlt);

        else if(strcmp(parmCStr, "-con-con-snps") == 0)
            sscanf(inputCStr, "%f", &conToConMinStats->minSNPsFlt);
        else if(strcmp(parmCStr, "-con-con-diff") == 0)
            sscanf(inputCStr, "%f", &conToConMinStats->minDiffFlt);
        else if(strcmp(parmCStr, "-con-con-dels") == 0)
            sscanf(inputCStr, "%f", &conToConMinStats->minDelsFlt);
        else if(strcmp(parmCStr, "-con-con-inss") == 0)
            sscanf(inputCStr, "%f", &conToConMinStats->minInssFlt);
        else if(strcmp(parmCStr, "-con-con-indels") == 0)
            sscanf(inputCStr, "%f", &conToConMinStats->minIndelsFlt);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-3: scoreReads read to reference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-ref-min-base-q") == 0)
            cStrToUChar(inputCStr, &readToRefMinStats->minQChar);

        else if(strcmp(parmCStr, "-read-ref-min-mapq") == 0)
            cStrToUInt(inputCStr, &readToRefMinStats->minMapqUInt);

        else if(strcmp(parmCStr, "-min-median-q") == 0)
            sscanf(inputCStr, "%f", &readToRefMinStats->minMedianQFlt);

        else if(strcmp(parmCStr, "-min-mean-q") == 0)
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

        else if(strcmp(parmCStr, "-min-read-length") == 0)
        { /*Else if the user provided a minimum read length*/
            readToRefMinStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);
            conSet->minMajConLenUI = readToRefMinStats->minReadLenULng;
        } /*Else if the user provided a minimum read length*/

         else if(strcmp(parmCStr, "-max-read-length") == 0)
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

        /**************************************************************\
        * Fun-1 Sec-2 Sub-4: scoreReads read to read settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-min-read-read-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            readToReadMinStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);
        } /*Else if the user provided a minimum read length*/

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

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to consensus settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-min-read-con-map-length") == 0)
        { /*Else if the user provided a minimum read length*/
            readToConMinStats->minReadLenULng =
                strtoul(inputCStr, &tmpCStr, 10);
        } /*Else if the user provided a minimum read length*/

        else if(strcmp(parmCStr, "-read-con-min-base-q") == 0)
            cStrToUChar(inputCStr, &readToConMinStats->minQChar);

        else if(strcmp(parmCStr, "-read-con-min-mapq") == 0)
            cStrToUInt(inputCStr, &readToConMinStats->minMapqUInt);

         else if(strcmp(parmCStr, "-read-con-max-a-ins-homo") == 0)
           cStrToUInt(inputCStr, &readToConMinStats->maxHomoInsAry[0]);

         else if(strcmp(parmCStr, "-read-con-max-t-ins-homo") == 0)
           cStrToUInt(inputCStr,&readToConMinStats->maxHomoInsAry[10]);

         else if(strcmp(parmCStr, "-read-con-max-c-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToConMinStats->maxHomoInsAry[1]);

         else if(strcmp(parmCStr, "-read-con-max-g-ins-homo") == 0)
            cStrToUInt(inputCStr,&readToConMinStats->maxHomoInsAry[3]);

         else if(strcmp(parmCStr, "-read-con-max-a-del-homo") == 0)
           cStrToUInt(inputCStr, &readToConMinStats->maxHomoDelAry[0]);

         else if(strcmp(parmCStr, "-read-con-max-t-del-homo") == 0)
           cStrToUInt(inputCStr,&readToConMinStats->maxHomoDelAry[10]);

         else if(strcmp(parmCStr, "-read-con-max-c-del-homo") == 0)
            cStrToUInt(inputCStr,&readToConMinStats->maxHomoDelAry[1]);

         else if(strcmp(parmCStr, "-read-con-max-g-del-homo") == 0)
            cStrToUInt(inputCStr,&readToConMinStats->maxHomoDelAry[3]);

        /**************************************************************\
        * Fun-1 Sec-2 Sub-6: scoreReads consensus to consensus settings
        \**************************************************************/

         else if(strcmp(parmCStr, "-con-con-max-a-ins-homo") == 0)
           cStrToUInt(inputCStr, &conToConMinStats->maxHomoInsAry[0]);

         else if(strcmp(parmCStr, "-con-con-max-t-ins-homo") == 0)
           cStrToUInt(inputCStr,&conToConMinStats->maxHomoInsAry[10]);

         else if(strcmp(parmCStr, "-con-con-max-c-ins-homo") == 0)
            cStrToUInt(inputCStr,&conToConMinStats->maxHomoInsAry[1]);

         else if(strcmp(parmCStr, "-con-con-max-g-ins-homo") == 0)
            cStrToUInt(inputCStr,&conToConMinStats->maxHomoInsAry[3]);

         else if(strcmp(parmCStr, "-con-con-max-a-del-homo") == 0)
           cStrToUInt(inputCStr, &conToConMinStats->maxHomoDelAry[0]);

         else if(strcmp(parmCStr, "-con-con-max-t-del-homo") == 0)
           cStrToUInt(inputCStr,&conToConMinStats->maxHomoDelAry[10]);

         else if(strcmp(parmCStr, "-con-con-max-c-del-homo") == 0)
            cStrToUInt(inputCStr,&conToConMinStats->maxHomoDelAry[1]);

         else if(strcmp(parmCStr, "-con-con-max-g-del-homo") == 0)
            cStrToUInt(inputCStr,&conToConMinStats->maxHomoDelAry[3]);

        else
            return parmCStr;
    } /*Loop through all user arguments (for)*/

    return 0;
} /*getUserInput*/
