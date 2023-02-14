/*#####################################################################\
# Name: findCoInfect
# Use:
#    - Detects co-infectoins in nanopore sequenced reads using Minimap2
#      & Racon
# External Requirments:
#    - Minimap2
#    - Racon
# Internal Requirments:
#    scoreReadsFun.c/h   # Has read scoring functions
#    trimSamFile.c/h     # Has functions to trim sam file alignments
#    samEntryStruct.c/h  # Has structure trimSamFile & scoreReadsFun use
#    printErrors.c/h     # For error messages
#    cStrToNumberFum.c/h # Number conversion functions
#    findCoInftBinTree.c/h # For building and AVL tree with the bins
#    fastqGrepAVLTree.c/h     # Needed for fastq id extraction
#    fastqGrepFastqFun.c/h    # Needed for fastq id extraction
#    fastqGrepHash.c/h        # Needed for fastq id extraction
#    fastqGrepSearchFastq.c/h # Needed for fastq id extractoin
#    fastqGrepStructs.c/h     # Needed for fastq id extraction
# Libraries used (in various files):
#     - <stdlib.h>
#     - <sdtint.h>
#     - <stdio.h>
#     - <string.h>
# Note:
#    - fastqGrep.c, trimSam.c, & scoreReads.c are here so that the user
#      can build these programs for personal use
#        - fastqGrep:
#            - Is for read id extraction from a fastq file. It
#              is faster than seqkit when extracting 30% to 50% of reads
#              from a file and only a few seconds slower for extracting
#              under 1% of reads.
#            - Requires
#                - fastqGrep.c
#                - fastqGrepAVLTree.c/h
#                - fastqGrepFastqFun.c/h
#                - fastqGrepHash.c/h
#                - fastqGrepSearchFastq.c/h
#                - fastqGrepStructs.c/h
#        - trimSam:
#            - Trims soft clipping of ends of sam alignments
#            - Requires
#                - trimSam.c
#                - trimSamFile.c/h
#                - samEntryStruct.c/h
#                - printErrors.c/h
#                - cStrToNumberFum.c/h
#        - scoreReads:
#            - Counts number of matches, SNPs, & indels in a read. It
#              also keeps track of the number matches, SNPs, & indels
#              above a certian q-score. For indels it will also ignore
#              indels in longer (user specified) homopolymers.
#            - Requires
#                - scoreReads.c
#                - scoreReads.c
#                - samEntryStruct.c/h
#                - printErrors.c/h
#                - cStrToNumberFum.c/h
\#####################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' TOC: Table of Contents                                               |
'    header toc:                                                       |
'        - Header section for functions, structures, & libraries       |
'    sof toc:                                                          |
'        - Start of functions                                          |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' Header TOC: Header section for functions and libaries
'    header sec-1: Libaries, default settings, & external program cmds
'    header sec-2: Structers
'    header sec-3: Functions
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
^ Header Sec-1: Libaries and definitions                               v
*    header sec-1 Sub-1: Included files
*    header sec-1 Sub-2: Version number & Default settings
*    header sec-1 Sub-3: System commands
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/**********************************************************************\
* Header Sec-1 Sub-1: Included files
\**********************************************************************/

/*My includes*/
#include "findCoInftBinTree.h" /*<string.h>, <stdio.h>, <stdlib.h>*/
#include "scoreReadsFun.h"
#include "trimSam.h"
#include "fqGetIdsSearchFq.h"     /*fqGrep files*/

/*Both scoreReadsFun.h and trimSamFile.h include:
    - samEntryStruct.h
        - printErrors.h
        - cStrToNumberFun.h
*/

/**********************************************************************\
* Header Sec-1 Sub-2: Version number & Default settings
\**********************************************************************/

#define bitsInFlt sizeof(float) << 3 /*Number bits in a float*/
#define bitsInSht sizeof(short) << 3 /*Number of bits in a short*/

#define defVersion 3.0       /*The Version number of this program*/
#define defPrefix "out"      /*Default prefix to use*/
#define defThreads "3"       /*Default number of threads to use*/
#define defReadersPerCon 300 /*max # reads used to build a consensus*/
#define minReadsPerBin 100   /*Min number reads to keep a bin*/
#define defNumPolish 2      /*Number of extra consensus building steps*/

#define defUseMajCon 1    /*1: Use majority in consnesus building*/
#define majConMinBaseQ 10 /*Min Q-score to keep base in majority con*/
#define percBasesPerPos 0.4 /*majority consensus (-maj-con-min-bases)*/
#define majConMinInsQ 5   /*Min Q-score to keep insertion in majority*/
#define percInsPerPos 0.3 /*% of supporting reads to keep insertion*/

#define defUseRaconCon 0  /*1: Racon in consensus building*/
#define defRoundsRacon 4    /*Default number of racon rounds/consensus*/

#define defUseMedakaCon 0 /*1: Use medaka in consensus building*/
#define defMedakaModel "r941_min_high_g351" /*Model to use with medaka*/


/*Read & consensus mapping filter settings*/
    /*As a rule I do not trust indels. SNPs are more likely in
      variants that are not quasi species
    */

/*Read mapped to refrence filters*/
#define readRefMinPercSNPs 0.07 /*Min % of SNPs to keep a read*/
#define readRefMinPercDiff 1    /*Min % difference to keep a read*/
#define readRefMinPercDels 1    /*Min % of deletions to keep a read*/
#define readRefMinPercInss 1    /*Min % of insertions to keep a read*/
#define readRefMinPercIndels 1  /*Min % of indels to keep a read*/

#define readRefMapq 20          /*Min mapq to keep read mapped ref*/

/*Read mapped to read filters*/
#define readReadMinPercSNPs 0.021 /*Min % of SNPs to keep a read*/
                                 /*0.015 works well for good databases,
                                   but 0.02 works better for sparse*/
#define readReadMinPercDiff 0.03 /*Min % difference to keep a read*/
#define readReadMinPercDels 1     /*Min % of deletions to keep a read*/
#define readReadMinPercInss 1     /*Min % of insertions to keep a read*/
#define readReadMinPercIndels 1   /*Min % of indels to keep a read*/

/*Read mapped to consensus filters*/
#define readConMinPercSNPs 0.015 /*Min % of SNPs to keep a read*/
#define readConMinPercDiff 0.02  /*Min % difference to keep a read*/
#define readConMinPercDels 1     /*Min % of deletions to keep a read*/
#define readConMinPercInss 1     /*Min % of insertions to keep a read*/
#define readConMinPercIndels 1   /*Min % of indels to keep a read*/

/*Consensus mapped to consensus filters*/
#define conConMinPercSNPs 0.01  /*Min % of SNPs to keep a consensus*/
#define conConMinPercDiff 0.03  /*Min % difference to keep a consensus*/
#define conConMinPercDels 1     /*% of deletions to keep a consensus*/
#define conConMinPercInss 1     /*% of insertions to keep consensus*/
#define conConMinPercIndels 1   /*Min % of indels to keep a consensus*/

/**********************************************************************\
* Header Sec-1 Sub-3: System commands
\**********************************************************************/

#define minimap2CMD "minimap2 --eqx --secondary=no -a -x map-ont"
#define raconCMD "racon -m 8 -x 6 -g -8 -w 500"

#define medakaCMD "/bin/bash -c \"source ~/medaka/venv/bin/activate;"
#define medakaCMDEnd "; deactivate\"" /*Telling to quite enviroment*/
#define medCondCMD "eval \"$(conda shell.bash hook)\"; conda activate medaka;"
#define medCondCMDEnd "; conda deactivate"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-2: Structers
^     header sec-2 struct-1: readStat   (Holds stats from scoreReads)
^     header sec-2 struct-2: readScore  (Holds query id and score)
^     header sec-2 struct-3: minDiff (min settings to keep reads/cons)
^     header sec-2 struct-4: baseStruct (error type, base, & counter)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*---------------------------------------------------------------------\
| Header Sec-2 Struct-1: readStat
| Use:
|    - Holds query id, reference id, & stats for a single read
\---------------------------------------------------------------------*/
typedef struct readStat
{ /*readStat*/
    char
        queryIdCStr[100],
        refIdCStr[100];

    unsigned char
        mapqUChar;

    float
       medianQFlt,     /*holds the median read Q-score*/
       medianAligQFlt, /*holds the aligned read median Q-score*/
       meanQFlt,       /*holds the mean read Q-score*/
       meanAligQFlt;   /*holds mean aligned read Q-score (no low Q)*/

    uint32_t
        readLenUInt,       /*Holds the read length of sam entry*/
        readAligLenUInt,   /*Aligned read length of sam entry*/
        numMatchUInt,      /*Holds number of matches*/
        numKeptMatchUInt,  /*Holds number matches with Q-score > min Q*/
        numKeptSNPUInt,    /*number of kept mismatches in sam entry*/
        numSNPUInt,        /*total number of mismatches in sam entry*/
        numKeptDelUInt,    /*number of kept deletions in sam entry*/
        numDelUInt,        /*total number of deletions in sam entry*/
        numKeptInsUInt,    /*number of kept insertions in sam entry*/
        numInsUInt;        /*total number of insertions in sam entry*/
}readStat;

/*---------------------------------------------------------------------\
| Header Sec-2 Struct-2: readScore
| Use:
|    - Holds the read id & score (as 32 bit integer) of read
\---------------------------------------------------------------------*/
typedef struct readScore
{ /*readScore*/
    char
        queryIdCStr[100];

    uint32_t
        readScoreUInt;

    struct readScore
        *nextRead;
}readScore;

/*---------------------------------------------------------------------\
| Header Sec-2 Struct-3: minDiff
| Use:
|    - Holds the min difference (SNPs, indels, or % sim) needed to 
|      consdier a read the same
\---------------------------------------------------------------------*/
typedef struct minDiff
{ /*readScore*/
    float
        minSNPsFlt,
        minDelsFlt,
        minInssFlt,
        minIndelsFlt,
        minDiffFlt;
}minDiff;

/*---------------------------------------------------------------------\
| Header Sec-2 Struct-4: baseStruct
| Use:
|    - Holds the base, error type (if insertion or deletion), & a
|      counter for number of reads with the same base.
\---------------------------------------------------------------------*/
typedef struct baseStruct
{ /*baseStruct*/
    char baseChar;     /*What is the base (0 for deletion)*/
    char errTypeChar;  /*0 = match, 1 = SNP, 2 = insertion*/
    
   unsigned long numBasesUL; /*Number of times this base is repeated*/

   struct baseStruct *altBase; /*Alternate options for the base*/
   struct baseStruct *nextBase; /*For linked lists*/
}baseStruct;

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-3: Function headers                                       v
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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
    unsigned int *numPolishUI, /*Number to times to rebuild consensus*/
    char *medakaModelCStr, /*Model to use with medaka*/
    char *useMajConBl,   /*Use majority consensus to build a consensus*/
    float *minBasePerFlt, /*Majority con min # reads to keep base*/
    float *minInsPerFlt, /*Majority con min # reads to keep insertion*/
    unsigned char *majConBaseQUC, /*Maj con min Q-score to keep base*/
    unsigned char *majConInsQUC, /*Maj con min Q-score to keep insert*/
    char *useRaconConBl, /*Use racon to build a consensus*/
    char *useMedakaConBl,/*Use medaka to build a consensus*/
    struct minDiff *minReadRefDiff, /*min difference for read to ref
                                      alignment to be different*/
    struct minDiff *minReadReadDiff, /*minReadRefDiff, for read read*/
    struct minDiff *minReadConDiff,/*minReadRefDiff, for read consenus*/
    struct minDiff *minConConDiff,/*minReadRef, for consensus consenus*/
    uint64_t *maxReadsPerConULng, /*Max reads for consensus building*/
    uint64_t *minReadsPerBinULng, /*min reads to keep a bin*/
    unsigned char *rndsRaconUC,            /*rounds to run racon*/
    struct minAlnStats *readToRefMinStats, /*Binning scoring settings*/
    struct minAlnStats *readToReadMinStats,/*Read pull scoring setting*/
    struct minAlnStats *readToConMinStats, /*Read cluster socring set*/
    struct minAlnStats *conToConMinStats   /*Consensus comparison set*/
); /*Reads in user input*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: minMap2CmdCStr to have the contents of refCStr and fqCStr
\---------------------------------------------------------------------*/
void addMinimap2Files(
    char *minMap2CmdCStr, /*Pointer to file portion of minimap2 cmd*/
    char *refFileCStr,    /*Has name of the file with references*/
    char *readFileCStr   /*Has name of fastq file with reads*/
); /*Adds the reference and read file names to the minimap2 command*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: raconCmdCStr to hold the file names input
\---------------------------------------------------------------------*/
void addRaconFiles(
    char *raconCmdCStr, /*Points to end of racon command*/
    char *refFileCStr,  /*Holds file with references*/
    char *samFileCStr,  /*Holds sam file to use with racon*/
    char *readsFileCStr, /*Holds reads to use*/
    char *outFileCStr    /*Holds file name to output consensus to*/
); /*Adds the input files to the racon command*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        -0: If does not meet min stats
|        -1: If meets min stats
\---------------------------------------------------------------------*/
uint8_t checkRead(
    struct minAlnStats *minStats, /*Structer holding min requirments*/
    struct samEntry *samStruct    /*Structer with the alignemtn stats*/
); /*Checks if the sam entry meets the min user requirements*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - readToBlank to have all stats set to 0 & c-strins to start 
|          with '\0'
\---------------------------------------------------------------------*/
void blankReadList(
    struct readStat *readToBlank
); /*Blanks a read list struct*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: newReadList to have the same values as oldReadList
\---------------------------------------------------------------------*/
void cpReadList(
    struct readStat *newReadList, /*Read to copy stats to*/
    struct readStat *oldReadList /*Read to copy stats from*/
); /*Copies stats from one read list to another*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Prints:                                                         |
|        - The read with the highest medain Q-score to the fastq file  |
|          name stored in binIn->bestReadCStr                          |
|    - Modifies:                                                       |
|        - Fastq file binIn-fqPathCStr to not have the best read       |
|        - Stats file binIn-statsPathCStr to not have the best read    |
|    - Returns:                                                        |
|        - 1: If sucessfull                                            |
|        - 2: For blank structer                                       |
|        - 4: For no fastq file                                        |
|        - 8: For no stats file                                        |
|        - 16: For error when extracting stats                         |
|        - 32: If could not open a temporary fastq file                |
|        - 64: If could not copy the fastq reads to their bins         |
\---------------------------------------------------------------------*/
uint8_t extractBestRead(
    struct readBin *binIn        /*Bin to extract best read from*/
); /*Splits bin into read with highest Q-score & other all other reads*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Modifies:                                                       |
|        - readStruct to have stats from the next line in statsFILE    |
\---------------------------------------------------------------------*/
uint8_t readStatsFileLine(
    FILE *statsFILE,            /*File with line to grab*/
    uint8_t *onHeaderBool,      /*1: skip one line, 0: grab first line*/
    struct readStat *readStruct /*Holds the stats from the stats file*/
); /*Reads single line from printSamStats function in samEntryStruct.c*/

/*---------------------------------------------------------------------\
| Output:
|    Prints: fastq entry to outFILE if finds
|    Returns:
|        - 1: if found and printed the id
|        - 2: If could not find the id (no printing)
|        - 4: If the fqFILE does not exist
|        - 8: If the outFILE does not exist
|        - 16: If the keptFILE (file for target read) does not exist
|        - 32: If their is an incomplete fastq entry
\---------------------------------------------------------------------*/
uint8_t fqOneIdExtract(
    char *idCStr, /*Read id to extract from the fastq file*/
    FILE *fqFILE,     /*fastq file to extract read from*/
    FILE *keptFILE,   /*File with the target read*/
    FILE *outFILE     /*fastq file to write read to*/
); /*Extracts one read id from a fastq file*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|   Prints: Prints out variables in readToPrint structer               |
\---------------------------------------------------------------------*/
void printReadStat(
    struct readStat *readToPrint, /*Read to print stats out for*/
    FILE *outFILE                  /*File to print read to*/
); /*Prints out the stats in a readStat structure to a file*/

/*---------------------------------------------------------------------\
| Output:
| Note:
|     - Score: percMult * (keptSNPs + keptIns + keptDels) / read length
|     - minSimUSht ranges from 1 (0.01%) to percMult (100%)
\---------------------------------------------------------------------*/
uint8_t findBestXReads(
    const uint64_t *numReadConsULng, /*# reads for bulding a consensus*/
    uint64_t *numReadsKeptULng,  /*Number of reads binned to con*/
    struct minDiff *maxDifference,
    char *threadsCStr,           /*Number threads to use with minimap2*/
    const char *useMapqBl,       /*1: use mapping quality in selection*/
    struct minAlnStats *minStats,/*Min stats to cluster reads together*/
    struct samEntry *samStruct,  /*Struct to use for reading sam file*/
    struct samEntry *refStruct,  /*holds the reference (0 to ignore)*/
    struct readBin *binTree      /*Bin working on*/
); /*Extract the top reads that mapped to the selected best read*/

/*---------------------------------------------------------------------\
| Output:
|    Uses: Racon to build a consensus (file name in bin->consensusCStr)
\---------------------------------------------------------------------*/
void buildConWithRacon(
    const uint8_t *clustUChar,         /*Cluster on*/
    const unsigned char *rndsRaconUC, /*Number of rounds to run racon*/
    const char *useFqConBl,   /*1 consensus is fastq, 0 use best read*/
    char *threadsCStr,             /*Number threads to use with racon*/
    struct readBin *binTree            /*Bin working on*/
); /*Builds a consensus using racon*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 1: if succeded
|    Modifies:
|        - fastq in binClust->fqPathCStr to be the fastq for the cluster
|        - fastq in binTree->fqPathCStr to not have clustered reads
\---------------------------------------------------------------------*/
uint8_t binReadToCon(
    struct minDiff *maxDifference,  /*Has min thresholds to keep reads*/
    const uint8_t *clustUChar,       /*Cluster on*/
    struct readBin *binTree,         /*Bin working on*/
    struct readBin *binClust,        /*Bin to assign reads & consensus*/
    struct samEntry *samStruct,      /*To hold temporary input*/
    struct minAlnStats *minStats,    /*Min stats needed to keep a read*/
    char *minimapCmdCStr,         /*Minimap2 command to run*/
    char *addFileToMiniCmdCStr    /*Ponts to position to add file
                                       names to in minimap2 command*/
); /*Maps reads to consensus and keeps reads that meet user criteria*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        - readBin Struct with the closest consensus
|    - Modifies:
|        - closesConUint to hold the score of the most similar consensus
| Note:
|    - clusters with ->balUChar < 0 will be ignored
\---------------------------------------------------------------------*/
struct readBin * cmpCons(
    struct minDiff *maxDifference,  /*Has min thresholds to keep reads*/
    struct readBin *conBin,       /*Bin with consensus to compare*/
    struct readBin *conBinTree,   /*Tree of consensus to compare to*/
    struct samEntry *samStruct,   /*Struct to hold input from minimap2*/
    struct samEntry *refStruct,   /*Struct to hold input from minimap2*/
    struct minAlnStats *minStats, /*Min stats needed to keep a error*/
    char *minimapCmdCStr,      /*Minimap2 command to run*/
    char *addFileToMiniCmdCStr /*Pionts to position to add file names
                                    in minimap2 command*/
); /*Compares a consenses to a another consensus. This will do a
    recursive call if conBinTree has children*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - binToKeep->fqPathCStr to hold binToMerge reads
|    Deletes:
|        - File named after binToMerge->fqPathCStr
\---------------------------------------------------------------------*/
void mergeBins(
    struct readBin *binToKeep,  /*Bin to merge into*/
    struct readBin *binToMerge /*Bin to merge into binToKeep*/
); /*Merge two bins togetehr*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 0: Read does not meet minimum thresholds in maxDifference
|        - 1: Keep the read (meets minimum thresholds)
\---------------------------------------------------------------------*/
uint8_t checkIfKeepRead(
    struct minDiff *maxDifference, /*Has thresholds to keep reads*/
    struct samEntry *samStruct     /*Has stats for read*/
); /*Checks if read does not meet one threshold in maxDifference*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies:
|        - newBin to have stats in samStruct
\---------------------------------------------------------------------*/
void samEntryToReadStat(
    struct readStat *newBin,   /*Read bin to hold stats from samStruct*/
    struct samEntry *samStruct /*copy stats from this struct*/
); /*Copies stats from a samEntry struct to a readStat struct*/

/*---------------------------------------------------------------------\
| Output:
|   Returns:
|     - 1: if both alignments are for the same query
|     - 0: if alignments are for different querys
\---------------------------------------------------------------------*/
uint8_t isSamAlnDup(
    struct samEntry *samAln,    /*sam alignemnt to check if duplicate*/
    struct samEntry *lastSamAln /*last sam alignment to compare to*/
); /*Checks if both sam entrys (alignments) have the same querys*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        binTree->consensusCStr to hold the polished consensus
|    Returns:
|        - 1 If built a consensus
|        - 2 If input file does not exists
|        - 4 IF consensus not built
\---------------------------------------------------------------------*/
unsigned char medakaPolish(
    const char *condaBl,    /*Tells if using miniconda medaka*/
    char *medakaModelCStr,  /*Model to use with medaka*/
    char *threadsCStr,      /*Number threads to use*/
    struct readBin *binTree /*bin with consensus & best reads*/
); /*Polish a consensus with medaka using the best reads*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold the copied C-string
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cStrCpInvsDelm(
    char *cpToCStr,  /*C-string to copy values to*/
    char *cpFromCStr /*C-string to copy*/
); /*Copy one c-string till an tab, newline, or '\0' (keeps spaces)*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold space, parameter, space, & argument
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cpParmAndArg(
    char *cpToCStr,   /*Holds copied parameter and argement*/
    char *cpParmCStr, /*Paramater to copy*/
    char *cpArgCStr   /*Argument to copy*/
); /*Copies adds a space and copies a paramater and a agrugment*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        o 1 If succeded
|        o 2 if could not open or read in the reference sequence
|        o 4 if could not open the bestReads file
|        o 16 if had to few sequences map to the selected read
|        o 32 if minimap2 crashed
|        o 64 for memory allocation errors
\---------------------------------------------------------------------*/

unsigned char simpleMajCon(
    const unsigned char *clustUChar, /*Number of cluster working on*/
    struct readBin *binStruct, /*Has best read & top reads files*/
    struct samEntry *samStruct,/*For reading in sam file entries*/
    unsigned char minBaseQUC,  /*Has min mapq to keep an base*/
    unsigned char minInsQUC,   /*Has min mapq to keep an insertion*/
    char *threadsCStr,         /*Number threads to use with minimap2*/
    const float *minBasePerFlt, /*Min percentage of bases needed to
                                       not be a deletion (counts all
                                       bases that pass the quality
                                       threshold)*/
    const float *minInsPercFlt  /*Min percentage of supporting reads
                                  needed to keep an insertion*/
); /*Builds a majority consensus from the best reads & top read*/

/*---------------------------------------------------------------------\
| Output:
|   o Creates:
|      - Fasta file with reads from the fastq file
|   o Returns:
|      - 1: Success
|      - 2: No fastq file (or invalid fastq file)
|      - 4: No fasta file
|      - 64: memory allocation error
\---------------------------------------------------------------------*/
unsigned char fqToFa(
    char *fqToCnvtCStr,         /*File name of fastq to convert*/
    char *outFaCStr,            /*File name of new fasta file*/
    struct samEntry *samStruct  /*Holds fastq file entries*/
); /*Converts a fastq file to a fasta file*/

/*---------------------------------------------------------------------\
| Output:
|    o Creates:
|      - Fasta file with the consensus
|    o Modifies:
|      - clustOn->consensusCStr to have file name of created consensus 
|    o Returns:
|      - 1 for success
|      - 2 or 4 for file errors
|      - 16 if had to few sequences map to the selected read
|      - 32 if minimap2 crashed
|      - 64 for memory allocation error
\---------------------------------------------------------------------*/
unsigned char buildCon(
    char *threadsCStr,       /*Number threads to use with system calls*/
    unsigned char *clustUC,  /*Cluter number*/
    char useMajConBl,        /*Use the majority (maj) consensus step*/
    float *minBasePerFlt,     /*For majority call del for bases < %*/
    float *minInsPercFlt,  /*for majority call min bases to keep ins*/
    unsigned char minBaseQUC, /*min Q-score to keep base (majority)*/
    unsigned char minInsQUC, /*Has min mapq to keep an insertion*/
    char useRaconConBl,      /*Use Racon in building a consnesus*/ 
    const unsigned char *rndsRaconUC, /*Number of rounds to run racon*/
    char useMedakaConBl,     /*Use the Medaka in building a consensus*/ 
    char *modelCStr,         /*Model to use with medaka*/
    const char *condaBl,     /*1: Use conda medaka install 0: python*/
    struct readBin *clustOn, /*Has fastq file to build the consensus*/
    struct samEntry *samStruct /*For reading in sequences or sam file*/
); /*Builds a consensus using the best read & top reads in clustOn*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies
|    - refStruct: To hold the sequence (no header)
|  - Returns
|    - 1 if succeeded
|    - 2 if file does not exist
|    - 4 invalid file
|    - 64 memory allocation error
| Note:
|  - Fasta file should only have one sequence
\---------------------------------------------------------------------*/
unsigned char readInConFa(
    char *conFaToReadCStr, /*Name of fasta file with the consensus*/
    struct samEntry *refStruct /*Sam struct to hold consensus*/
); /*Reads in reference sequence in fasta file*/

/*---------------------------------------------------------------------\
| Output:
|   o Returns:
|     - The next base in the list (if altBase != 0, returns altBase)
|   o Frees:
|     - The baseToFree from memory
|   o Modifies:
|     - If their is an alternative base (baseToFree->altBase != 0)
|         o lastBase->nextBase is set to baseToFree->altBase
|         o bastToFree is set to bastToFree->altBase
|     - If their is not alternative base (bastToFree->altBase == 0)
|         o lastBase->nextBase is set to baseToFree->nextBase
|         o bastToFree is set to 0
| Note:
|    o This function assumes that the next base pointer (nextBase) is
|      is alwasy 0 (not set) for teh alternate base pointer (altBase==0)
\---------------------------------------------------------------------*/
struct baseStruct * freeBaseStruct(
    struct baseStruct **baseToFree, /*Insertion list to free*/
    struct baseStruct *lastBase     /*Base to assign pointers to*/
); /*Frees an base from a linked list of bases*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF TOC: Start of Functions
'    main:
'        - main function to run everything
'    fun-1 getUserInput:
'        - Gets user input
'    fun-2 addMinimap2Files:
'        - Adds file names to the minimap2 command
'    fun-3 addRaconFiles:
'        - Adds file names to the racon command
'    fun-4 checkRead:
'        - Checks if should keep the sam alignment
'    fun-5 blankReadList:
'        - Blank a readList structer
'    fun-6 cpReadList:
'        - Copy a readList structer
'    fun-7 extractBestRead:
'        - Extract the best read from a bin
'    fun-8 readStatsFileLine:
'        - Read a single line from the stats file output by scoreReads
'    fun-9 fqOneIdExtract:
'        - Extract one read from a fastq file
'    fun-10 printReadStat:
'        - Prints out read id, ref id, & stats in a readStat structure
'    fun-11 findBestXReads:
'        - Extract the top x (user specified) reads that mapped to
'          the selected best read.
'    fun-12 buildConWithRacon:
'        - Buids a consensus using racon
'    fun-13 binReadToCon:
'        - Bins reads to a consensus to from a cluster
'    fun-14 cmpCons:
'        - Compares two consensus (does recursive call if tree input
'    fun-15 mergeBins:
'        - Merges reads in two bins & delets one of the bines fq file
'    fun-16 checkIfKeepRead:
'        - Checks if stats in reads samEntry struct meet the min
'          thresholds to keep a read
'    fun-17 samEntryToReadStat:
'        - Copies stats from samEntry to readStat struct (not used)
'    fun-18 isSamAlnDup:
'        - Checks if two sam alignmnents are aligned to the same query
'    fun-19 medakaPolish:
'        - Polish a consensus with medaka using the best reads
'    fun-20 cStrCpInvsDelm:
'        - Copy one c-string till an tab, newline, or '\0' (keeps ' ')
'    fun-21 cpParmAndArg:
'        - Copies adds a space and copies a paramater and a agrugment
'    fun-22 simpleMajCon:
'        - Builds a majority consensus (ignores insertions)
'    fun-23 fqToFa:
'        - Convert fastq file to fasta file (old fastq file kept)
'    fun-24 buildCon:
'        - Builds a consensus using: majority, racon, medaka functions
'    fun-26 freeBaseStruct:
'        - Frees a base struct from a linked list.
'        - The freeded base is set to the alternate base if their is an
'          alternate base or 0 if no alternative base is present.
'        - Assumes alternateive bases (altBase) have nextBase set to
'          0 (null).
'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
    char useMajConBl = defUseMajCon;     /*Use majority consensus step*/
    char useRaconConBl = defUseRaconCon;  /*Use racon step*/
    char useMedakaConBl = defUseMedakaCon;/*Use medaka step*/
    char prefCStr[100];        /*Holds the user prefix*/
    char threadsCStr[7];      /*Number of threads for minimap2 & racon*/
    char medakaModelCStr[128]; /*Model to use with medaka*/


    /*Number of times to build the extra conesnsus before clustering*/
    unsigned int numPolishUI = defNumPolish;

    /*Variables for the majority consensus step*/
    unsigned char minBaseQUC = majConMinBaseQ; /*min q to keep base*/
    unsigned char minInsQUC = majConMinInsQ; /*min q to keep insertion*/
    float minBasePerFlt = percBasesPerPos; /*% of reads to keep base*/
    float minInsPercFlt = percInsPerPos; /*% of reads to keep insert*/

    /*Programs finds these settings automaticly*/
    char condaBl = 0;            /*Tells if using conda medaka install*/
    unsigned char lenPrefUChar = 0;           /*Holds length of prefix*/
    unsigned char rndsRaconUC = defRoundsRacon; /*Number racon rounds*/

    uint64_t
        numReadConsULng = defReadersPerCon,
            /*max number of reads used to build a consensus*/
        numReadsKeptULng = 0,  /*Number of reads extracted*/
        minReadsPerBinULng = minReadsPerBin;
            /*Min number reads to keep a bin*/

    /*FILE names*/
    char
        *fqPathCStr = 0,      /*Holds fastq file to process*/
        *refsPathCStr = 0,    /*Holds references for binning*/
        binFileCStr[200],     /*Holds name of the assigned bin*/
        logFileCStr[200],    /*Holds the name of the log file*/
        statFileCStr[200],    /*Holds name of stat file for a bin*/
        readCntFileCStr[200]; /*Holds Number of reads per bin/cluster*/

    /*C-strings that hold commands*/
    char
        useMapqBl = 0,         /*1: Tells findBestXReads to use mapq*/
        tmpCmdCStr[1024],      /*Holds a quick system command*/
        minimap2CmdCStr[1024], /*Buffer for minimap2 command*/
        *tmpMiniMapCStr = 0, /*Points to file entry in minimap2CmdCStr*/
        raconCmdCStr[1024],    /*Holds command to run racon*/
        *tmpRaconCStr = 0, /*Points to file entry part of raconCmdCStr*/
        *tmpCStr = 0,         /*For string manipulation*/
        *cpTmpCStr = 0;       /*C-string to copy to*/


    /*Miscalanious variables*/
    uint8_t
        clustUChar = 0,       /*cluster on*/
        zeroUChar = 0,        /*Handy 0 to pass around*/
        errUChar = 0,         /*Holds error messages*/
        dupBool = 0,          /*Tells if was a duplicate*/
        printStatsHeadUChar = 0; /*Tells if to print header for stats*/

    char
        *inputErrCSStr = 0; /*holds user input error*/

    int64_t
        tmpULng = 0;

    /*FILES opened*/
    FILE
        *logFILE = 0,       /*Holds the log*/
        *stdinFILE = stdin, /*Points to piped input from minimap2*/
        *statFILE = 0,      /*File to output stats from score reads to*/
        *tmpFILE = 0,
        *fqBinFILE = 0;     /*Points to file adding binned reads to*/

    /*My structures*/
    struct samEntry
        *tmpSam = 0, /*For swapping newSam and oldSam pointers*/
        *newSam = 0, /*For duplicate detection*/
        *oldSam = 0, /*For duplicate detection*/
        samStruct, /*Holds sam file entry from minimap2*/
        refStruct, /*Holds the reference sequence*/
        *samZeroStruct = 0; /*Handy zero to pass around*/

    /*Holds thresholds to keep alignments, errors, & matches*/
    struct minAlnStats
        readToRefMinStats,     /*map read to reference thresholds*/
        readToReadMapMinStats, /*read to read mapping thresolds*/
        readToConMinStats,     /*map read to reference thresholds*/
        conToConMinStats;  /*Consensus to consensus mapping thresholds*/

    
    struct readBinStack binStack[200]; /*Stack for readBin avl tree*/

    struct readBin *binTree = 0;   /*Tree of read bins*/
    struct readBin *bestBin = 0;   /*Bin with most similar consensus*/
    struct readBin *clustOn = 0;   /*Bin in list/tree am working on*/
    struct readBin *lastClust = 0; /*Last cluster worked on for a bin*/
    struct readBin *lastBin = 0;   /*For keeping list in order*/
    struct readBin *tmpBin = 0;    /*Pionts to a readBin to work on*/
       
    struct minDiff minReadRefDiff;
    struct minDiff minReadReadDiff;
    struct minDiff minReadConDiff;
    struct minDiff minConConDiff;

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
            \n    -min-reads-per-bin:\
            \n        - Min number of reads needed to keep   [100]\
            \n          a bin or a cluster\
            \n    -max-reads-per-con:\
            \n        - Max number of reads to use in        [300]\
            \n          a consensus.\
            \n    -extra-consensus-steps:                    [2]\
            \n        - Number of times to rebuild the\
            \n          consensus using a new set of best\
            \n          reads.\
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
            \n    -enable-majority-consensus:                   [Yes]\
            \n        - Build a consensus using a simple\
            \n          majority consensus. This consensus\
            \n          will be polished with Racon or\
            \n          Medaka if Racon and Medaka set.\
            \n    -maj-con-min-bases                         [0.4=40%]\
            \n        - When building the majority consesus\
            \n          make a deletion in positions that\
            \n          have less than x\% of bases (40%).\
            \n    -maj-con-min-base-q:                       [10]\
            \n        - Minimum q-score to keep a base when\
            \n          building a majority consensus.\
            \n    -maj-con-min-ins                           [0.4=40%]\
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
            \n    -enable-medaka:                               [No]\
            \n        - Do not use medaka to polish the\
            \n          consensus.\
            \n    -model:\
            \n        - Model to use with medaka   [r941_min_high_g351]\
            \n          (calling medaka_consensus)\
       "; /*consensus building parameters*/



    char
         *binHelpCStr = "\
            \n The binning step compares a read to the best reference\
            \n    its mapped to. If the read meets the qaulity\
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
            \n    - Minimum read length to keep a read\
            \n      when binning (applies only to trimmed\
            \n      reads).\
            \n -read-ref-max-read-length:                 [1000]\
            \n    - Maximum read length to keep a read\
            \n      when binning (applies only to\
            \n      trimmed reads).\
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

    /*read to reference (inital binning) mapping*/
    minReadRefDiff.minSNPsFlt = readRefMinPercSNPs;
    minReadRefDiff.minDiffFlt = readRefMinPercDiff;
    minReadRefDiff.minDelsFlt = readRefMinPercDels;
    minReadRefDiff.minInssFlt = readRefMinPercInss;
    minReadRefDiff.minIndelsFlt = readRefMinPercIndels;

    /*Read to read mapping*/
    minReadReadDiff.minSNPsFlt = readReadMinPercSNPs;
    minReadReadDiff.minDiffFlt = readReadMinPercDiff;
    minReadReadDiff.minDelsFlt = readReadMinPercDels;
    minReadReadDiff.minInssFlt = readReadMinPercInss;
    minReadReadDiff.minIndelsFlt = readReadMinPercIndels;

    /*read to consensus mapping*/
    minReadConDiff.minSNPsFlt = readConMinPercSNPs;
    minReadConDiff.minDiffFlt = readConMinPercDiff;
    minReadConDiff.minDelsFlt = readConMinPercDels;
    minReadConDiff.minInssFlt = readConMinPercInss;
    minReadConDiff.minIndelsFlt = readConMinPercIndels;

    /*consensus to consensus mapping*/
    minConConDiff.minSNPsFlt = conConMinPercSNPs;
    minConConDiff.minDiffFlt = conConMinPercDiff;
    minConConDiff.minDelsFlt = conConMinPercDels;
    minConConDiff.minInssFlt = conConMinPercInss;
    minConConDiff.minIndelsFlt = conConMinPercIndels;

    initSamEntry(&samStruct);
    initSamEntry(&refStruct);

    /*Need to add indel settings in*/
    blankMinStats(&readToRefMinStats);
    blankMinStats(&readToReadMapMinStats);
    blankMinStats(&readToConMinStats);
    blankMinStats(&conToConMinStats);

    readToRefMinStats.minMapqUInt = readRefMapq;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Main Sec-3: Get user input and set up commands 
    #    sec-1: Get user input & check for errors
    #    sec-3: Set up minimap2 & racon commands
    #    sec-3: Add prefix to file names
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Main Sec-3 Sub-1: Get user input & check for errors
    *******************************************************************/

    /*Set up default values*/
    strcpy(threadsCStr, defThreads); /*Setup default number of threads*/
    strcpy(medakaModelCStr, defMedakaModel); /*default model*/
    strcpy(prefCStr, defPrefix);

    inputErrCSStr =
        getUserInput(
            lenArgsInt,
            argsCStr,
            prefCStr,
            &fqPathCStr,
            &refsPathCStr,
            threadsCStr,
            &numPolishUI, /*Number to times to rebuild consensus*/
            medakaModelCStr,
            &useMajConBl,
            &minBasePerFlt,
            &minInsPercFlt,
            &minBaseQUC,
            &minInsQUC,
            &useRaconConBl,
            &useMedakaConBl,
            &minReadRefDiff,
            &minReadReadDiff,
            &minReadConDiff,
            &minConConDiff,
            &numReadConsULng,
            &minReadsPerBinULng,
            &rndsRaconUC,
            &readToRefMinStats,
            &readToReadMapMinStats,
            &readToConMinStats,
            &conToConMinStats
    ); /*Get user input*/

    if(inputErrCSStr != 0)
    { /*If have an error*/

        if(inputErrCSStr == 0)
        { /*If no user input was supplied*/
            fprintf(
                stderr,
                "%s\n\nNo user input was supplied\n",
                helpCStr
            );

            exit(1);
        } /*If no user input was supplied*/

        if(
            strcmp(inputErrCSStr, "-v") == 0 ||
            strcmp(inputErrCSStr, "-version") == 0
        ) { /*If the user is requesting the version number*/
            fprintf(stdout, "%f\n", defVersion);
            exit(0);
        } /*If the user is requesting the version number*/

        if(strcmp(inputErrCSStr, "-h-build-consensus") == 0)
        { /*If the user wants the consensus building options*/
            fprintf(stdout, "%s", conBuildHelpCStr);
            exit(0);
        } /*If the user wants the consensus building options*/
        
        if(strcmp(inputErrCSStr, "-h-bin") == 0)
        { /*If user wants to know about the binning parameters*/
            fprintf(stdout, "%s", binHelpCStr);
            exit(0);
        } /*If user wants to know about the binning parameters*/

        if(strcmp(inputErrCSStr, "-h-read-map") == 0)
        { /*If user wants to know about the read mapping parameters*/
            fprintf(stdout, "%s", readMapHelpCStr);
            exit(0);
        } /*If user wants to know about the read mapping parameters*/

        if(strcmp(inputErrCSStr, "-h-clust") == 0)
        { /*If the user wants to know about the clustering parameters*/
            fprintf(stdout, "%s", clustHelpCStr);
            exit(0);
        } /*If the user wants to know about the clustering parameters*/

        if(strcmp(inputErrCSStr, "-h-con") == 0)
        { /*If user wants the consensus comparison parameters*/
            fprintf(stdout, "%s", conCompHelpCStr);
            exit(0);
        } /*If user wants the consensus comparison parameters*/

        if(
            strcmp(inputErrCSStr, "-h") == 0 ||
            strcmp(inputErrCSStr, "-help") == 0 ||
            strcmp(inputErrCSStr, "--h") == 0 ||
            strcmp(inputErrCSStr, "--help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);/*Print out help message*/
            exit(0);
        } /*If user wanted the help message*/

        fprintf(
            stderr,
            "%s\n\n %s is an invalid parameter\n",
            helpCStr,   /*Print out the help message*/
            inputErrCSStr     /*Print out the error*/
        ); /*Let user know about the invalid parameter*/

        exit(1);
    } /*If have an error*/

    if(!(useMajConBl | useRaconConBl | useMedakaConBl))
    { /*If the user said to ingore all consensus building steps*/
        printf("Current settings have turned off all consensus");
        printf(" building methods.\nSelect a consensus step by");
        printf(" removing: -enable-medaka or -enable-racon\n");
    } /*If the user said to ingore all consensus building steps*/

    /*******************************************************************
    # Main Sec-3 Sub-2: Set up minimap2 & racon commands
    *******************************************************************/

    tmpMiniMapCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
    tmpMiniMapCStr = cpParmAndArg(tmpMiniMapCStr, "-t", threadsCStr);

    tmpRaconCStr = cStrCpInvsDelm(raconCmdCStr, raconCMD);
    tmpRaconCStr = cpParmAndArg(tmpRaconCStr, "-t", threadsCStr);

    /*******************************************************************
    # Main Sec-3 Sub-3: Add prefix to file names
    *******************************************************************/

    strcpy(logFileCStr, prefCStr);
    strcpy(binFileCStr, prefCStr);
    strcpy(statFileCStr, prefCStr);
    strcpy(readCntFileCStr, prefCStr);

    tmpCStr = prefCStr;

    while(*tmpCStr != 0)
    { /*While need to find the length of the prefix*/
        ++tmpCStr;
        ++lenPrefUChar;
    } /*While need to find the length of the prefix*/

    strcpy(logFileCStr + lenPrefUChar, "--log.txt"); /*Finsh log name*/
    strcpy(readCntFileCStr + lenPrefUChar, "--read-counts.tsv");

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

    if(useRaconConBl & 1)
    { /*If using racon, get the version used*/
        /*Set up racon check*/
        tmpCStr = cpParmAndArg(tmpCmdCStr, "racon", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);

        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Racon version: ");
        fclose(logFILE); /*Closing to avoid system appending to file*/

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
    } /*Else if not using racon*/

    /******************************************************************\
    * Main Sec-4 Sub-4: Check if medaka exists
    \******************************************************************/

    if(useMedakaConBl & 1)
    { /*If using medaka, check version*/
        /*Set up non-miniconda command*/
        tmpCStr = cStrCpInvsDelm(tmpCmdCStr, medakaCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "medaka", "--version");
        tmpCStr = cpParmAndArg(tmpCStr, ">>", logFileCStr);
        tmpCStr = cpParmAndArg(tmpCStr, medakaCMDEnd, "");

        logFILE = fopen( logFileCStr, "a");
        fprintf(logFILE, "Medaka version: ");
        fclose(logFILE);

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

            fprintf(
                logFILE,
                "    - Using medaka installed by miniconda\n"
            );
            condaBl = 1; /*Using medaka from miniconda*/
        } /*If python virtual enviorment medaka does not exist*/

        else
        { /*Else if found the python virtual enviorment medaka*/
            condaBl = 0; /*Not using medaka from miniconda*/
            fprintf(logFILE, "    - Using medaka installed by python\n");
        } /*Else if found the python virtual enviorment medaka*/
    } /*If using medaka, check version*/

    else
    { /*Else if not using medaka*/
        logFILE = fopen(logFileCStr, "a");
        fprintf(logFILE, "Not using Medaka\n");
        fclose(logFILE);
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

    stdinFILE = fopen(refsPathCStr, "r"); /*Open the reference file*/

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
    fprintf(logFILE, "    -model %s \\\n", medakaModelCStr);

    fprintf(
        logFILE,
        "    -min-reads-per-bin %ju \\\n",
        minReadsPerBinULng
    );

    fprintf(
        logFILE,
        "    -max-reads-per-con %ju \\\n",
        numReadConsULng
    );

    fprintf(logFILE, "    -extra-consensus-steps %u \\\n", numPolishUI);
    fprintf(logFILE, "    -rounds-racon %u \\\n", rndsRaconUC);
    fprintf(logFILE,"    -enable-majority-consensus %u \\\n", useMajConBl);
    fprintf(logFILE,"    -maj-con-min-bases %f \\\n", minBasePerFlt);
    fprintf(logFILE,"    -maj-con-min-base-q %u \\\n", minBaseQUC);
    fprintf(logFILE,"    -maj-con-min-ins %f \\\n", minInsPercFlt);
    fprintf(logFILE,"    -maj-con-min-ins-q %u \\\n", minInsQUC);
    fprintf(logFILE,"    -enable-racon %u \\\n", useRaconConBl);
    fprintf(logFILE,"    -enable-medaka %u \\\n", useMedakaConBl);

    /******************************************************************\
    * Main Sec-5 Sub-2: Print out Comparision settings for binning
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-ref-snps %f \\\n",
        
        minReadRefDiff.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-diff %f \\\n",
        
        minReadRefDiff.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-dels %f \\\n",
        
        minReadRefDiff.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-inss %f \\\n",
        
        minReadRefDiff.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-ref-indels %f \\\n",
        
        minReadRefDiff.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-3: Print out Comparision settings for read mapping
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-read-snps %f \\\n",
        
        minReadReadDiff.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-read-diff %f \\\n",
        
        minReadReadDiff.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-read-dels %f \\\n",
        
        minReadReadDiff.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-read-inss %f \\\n",
        
        minReadReadDiff.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-read-indels %f \\\n",
        
        minReadReadDiff.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-3: Print out Comparision settings for clustering
    \******************************************************************/

    fprintf(
        logFILE,
        "    -read-con-snps %f \\\n",
        
        minReadConDiff.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -read-con-diff %f \\\n",
        
        minReadConDiff.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -read-con-dels %f \\\n",
        
        minReadConDiff.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -read-con-inss %f \\\n",
        
        minReadConDiff.minInssFlt
    );

    fprintf(
        logFILE,
        "    -read-con-indels %f \\\n",
        
        minReadConDiff.minIndelsFlt
    );

    /******************************************************************\
    * Main Sec-5 Sub-4: Print Comparision settings for consensus compare
    \******************************************************************/

    fprintf(
        logFILE,
        "    -con-con-snps %f \\\n",
        
        minConConDiff.minSNPsFlt
    );

    fprintf(
        logFILE,
        "    -con-con-diff %f \\\n",
        
        minConConDiff.minDiffFlt
    );

    fprintf(
        logFILE,
        "    -con-con-dels %f \\\n",
        
        minConConDiff.minDelsFlt
    );

    fprintf(
        logFILE,
        "    -con-con-inss %f \\\n",
        
        minConConDiff.minInssFlt
    );

    fprintf(
        logFILE,
        "    -con-con-indels %f \\\n",
        
        minConConDiff.minIndelsFlt
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
         "    -read-ref-min-mapqq %u \\\n",
         readToRefMinStats.minMapqUInt
     );

     fprintf(
         logFILE,
         "    -read-ref-min-median-q %f \\\n",
         readToRefMinStats.minMedianQFlt
     );

     fprintf(
         logFILE,
         "    -read-ref-min-mean-q %f \\\n",
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
         "    -read-ref-min-read-length %u \\\n",
         readToRefMinStats.minReadLenULng
     );

     fprintf(
         logFILE,
         "    -read-ref-max-read-length %u \\\n",
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
         readToReadMapMinStats.minQChar
     );

     fprintf(
         logFILE,
         "    -read-read-min-mapqq %u \\\n",
         readToReadMapMinStats.minMapqUInt
     );

     fprintf(
         logFILE,
         "    -read-read-max-a-ins-homo %u \\\n",
         readToReadMapMinStats.maxHomoInsAry[0]
     );

     fprintf(
         logFILE,
         "    -read-read-max-t-ins-homo %u \\\n",
         readToReadMapMinStats.maxHomoInsAry[10]
     );

     fprintf(
         logFILE,
         "    -read-read-max-c-ins-homo %u \\\n",
         readToReadMapMinStats.maxHomoInsAry[1]
     );

     fprintf(
         logFILE,
         "    -read-read-max-g-ins-homo %u \\\n",
         readToReadMapMinStats.maxHomoInsAry[3]
     );

     fprintf(
         logFILE,
         "    -read-read-max-a-del-homo %u \\\n",
         readToReadMapMinStats.maxHomoDelAry[0]
     );

     fprintf(
         logFILE,
         "    -read-read-max-t-del-homo %u \\\n",
         readToReadMapMinStats.maxHomoDelAry[10]
     );

     fprintf(
         logFILE,
         "    -read-read-max-c-del-homo %u \\\n",
         readToReadMapMinStats.maxHomoDelAry[1]
     );

     fprintf(
         logFILE,
         "    -read-read-max-g-del-homo %u \\\n",
         readToReadMapMinStats.maxHomoDelAry[3]
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

    fflush(logFILE); /*Flush output to the log file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Main Sec-6: Find initial bins with references
    ^    main sec-6 sub-1: Run minimap2 & read in first line of output
    ^    main sec-6 sub-2: Check if first line is valid
    ^    main sec-6 sub-3: Check if is a header 
    ^    main sec-6 sub-4: Check if query has multiple sequence entries
    ^    main sec-6 sub-5: Trim and score sam file alignments
    ^    main sec-6 sub-6: Set up bin file names
    ^    main sec-6 sub-7: Check if can open the bin & stats file
    ^    main sec-6 sub-8: Append reads to file & add to tree
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-6 Sub-1: Run minimap2 & read in first line of output
    \******************************************************************/

    /*Running minimap2 with on thread so that I can detect duplicate
      entries. Otherwise minimap2 will have no order for ouput
      mappings*/
    tmpCStr = cStrCpInvsDelm(tmpCmdCStr, minimap2CMD);
    tmpCStr = cpParmAndArg(tmpCStr, "-t", "1");
    tmpCStr = cpParmAndArg(tmpCStr, refsPathCStr, fqPathCStr);

    stdinFILE = popen(tmpCmdCStr, "r"); /*run minimap2*/

    blankSamEntry(&samStruct); /*Remove old stats in sam file*/
    errUChar = readSamLine(&samStruct, stdinFILE);
        /*get the first line from minimap2*/

    /******************************************************************\
    * Main Sec-6 Sub-2: Check if first line is valid
    \******************************************************************/

    if(!(errUChar & 1))
    { /*If an error occured*/
        pclose(stdinFILE);    /*No longer need open (due to error*/
        freeStackSamEntry(&samStruct);
        freeStackSamEntry(&refStruct);

        if(errUChar & 2)
        { /*If minimap2 did not run*/
            /*Let user know minimap2 errored out*/
            fprintf(stderr, "Minimap2 errored out\n");

            fprintf(
                stderr,
                "Check input files using minimap2 --eqx -a -x map-ont\n"
            ); /*Let user know how to replicated the error*/

            /*Let user know minimap2 errored out by log*/
            fprintf(logFILE, "Minimap2 errored out\n");

            fprintf(
                logFILE,
                "Check input files using minimap2 --eqx -a -x map-ont\n"
            ); /*Let user know how to replicated the error*/

            fclose(logFILE);
            exit(1);
        } /*If minimap2 did not run*/

        else if(errUChar & 64)
        { /*If was a memory allocation error*/
            fprintf(logFILE, "Ran out of memory when mapping reads to");
            fprintf(logFILE, " references using minimap2 (binning)\n");
            fclose(logFILE);
            exit(1); /*Error printed to user already*/
        } /*If was a memory allocation error*/

        else
        { /*Else, let user know something happened*/
            fprintf(stderr, "minimap2 error during binning step\n");
            fprintf(logFILE, "minimap2 error during binning step\n");
            fclose(logFILE);
            exit(1);
        } /*Else, let user know something happened*/
    } /*If an error occured*/

    if(*(samStruct.samEntryCStr) != '@')
    { /*If their is no header line, minimap2 likely errored out*/
        fprintf(
            stderr,
            "error: minimap2 output missing header (binning)\n"
        );
        fprintf(
            logFILE,
            "error: minimap2 output missing header (binning)\n"
        );

        fclose(logFILE);
        exit(1);
    } /*If minimap2 did not produce a header, it likely error out*/

    /******************************************************************\
    * Main Sec-6 Sub-3: Get past the header
    \******************************************************************/

    newSam = &samStruct;
    oldSam = 0;         /*Set up for first round*/

    while(errUChar & 1)
    { /*While not past the first header*/
        if(*newSam->samEntryCStr == '@')
        { /*If was a header*/
            blankSamEntry(&samStruct); /*Remove old stats in sam file*/
            errUChar = readSamLine(newSam, stdinFILE); /*read new line*/
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        /*Convert & print out sam file entry*/
        errUChar = trimSamEntry(newSam);

        if(errUChar >> 2)
        { /*If entry did not have a sequence, discard*/
            /*make the new alignment (not dupicate) the old alignment*/
            tmpSam = newSam;
            newSam = oldSam;
            oldSam = tmpSam;

            /*Read the next sam entry*/
            errUChar = readSamLine(newSam, stdinFILE);

            continue;
        } /*If entry did not have a sequence, discard*/

        findQScores(newSam); /*Find the Q-scores*/

        scoreAln(
            &readToRefMinStats, /*thesholds for read to reference map*/
            newSam,
            samZeroStruct, /*Not using reference for scoring*/ 
            &zeroUChar,    /*Not using reference, so no Q-score*/
            &zeroUChar     /*Not using reference, so no deletions*/
        );

        oldSam = newSam;    /*Swap pointers around*/
        newSam = &refStruct;

        blankSamEntry(newSam); /*Remove old stats in sam file*/
        /*Read in the next line*/
        errUChar = readSamLine(newSam, stdinFILE);
        break; /*Found the second entry*/
    } /*While not past the first header*/

    /**************************************************************\
    * Main Sec-6 Sub-4: Check if for mulitple entries for same query
    *    - Sometimes minimap2 will split a read into two different
    *      alignments when the read has duplicate regions. This
    *      splitting will mess up some of my later functions.
    \**************************************************************/

    while(errUChar & 1)
    { /*While their is a samfile entry to read in*/

        /*Putting this here so I can detect & remove chimeras*/
        dupBool = 0; /*start off assuming not duplicate*/

        while(isSamAlnDup(newSam, oldSam) & 1)
        { /*While the entries are duplicates*/
            if(*newSam->seqCStr != '*')
                dupBool = 1;

            blankSamEntry(newSam); /*Remove previous reads stats*/
            errUChar = readSamLine(newSam, stdinFILE); /*get next line*/

            if(!(errUChar & 1))
                break;      /*end of file or other error*/
        } /*While the entries are duplicates*/

        if(!(errUChar & 1))
            break;         /*end of file or other error*/

        if(dupBool & 1)
        { /*If was a duplicate*/
            /*make the new alignment (not dupicate) the old alignment*/
            tmpSam = newSam;
            newSam = oldSam;
            oldSam = tmpSam;
    
            errUChar = readSamLine(newSam, stdinFILE); /*get next line*/
            continue; /*Restart with next line with new old seq*/
        } /*If was a duplicate*/

        /**************************************************************\
        * Main Sec-6 Sub-5: Trim and score sam file alignment
        \**************************************************************/

        /*Convert & print out sam file entry*/
        errUChar = trimSamEntry(newSam);

        if(errUChar >> 2)
        { /*If entry did not have a sequence, discard*/
            /*make the new alignment (not dupicate) the old alignment*/
            tmpSam = newSam;
            newSam = oldSam;
            oldSam = tmpSam;

            /*Read the next sam entry*/
            errUChar = readSamLine(newSam, stdinFILE);
            continue;
        } /*If entry did not have a sequence, discard*/

        findQScores(newSam); /*Find the Q-scores*/

        scoreAln(
            &readToRefMinStats, /*thesholds for read to reference map*/
            newSam,
            samZeroStruct, /*Not using reference for scoring*/ 
            &zeroUChar,    /*Not using reference, so no Q-score*/
            &zeroUChar     /*Not using reference, so no deletions*/
        );

        if(
          checkRead(&readToRefMinStats, oldSam) == 0 ||
          !(checkIfKeepRead(&minReadRefDiff, oldSam) & 1)
        ) { /*If the read is under the min quality, discard*/
            /*make the new alignment (not dupicate) the old alignment*/
            tmpSam = newSam;
            newSam = oldSam;
            oldSam = tmpSam;

            blankSamEntry(newSam); /*Remove old stats in sam file*/
            /*Read in the next line*/
            errUChar = readSamLine(newSam, stdinFILE);

            continue;
        } /*If the read is under the min quality, discard*/

        /**************************************************************\
        * Main Sec-6 Sub-5: Set up bin file names
        \**************************************************************/

        /*Build the bin file name*/
        tmpCStr = binFileCStr + lenPrefUChar;
        strcpy(tmpCStr, "--");
        tmpCStr += 2;
        cpTmpCStr = oldSam->refCStr;

        while(*cpTmpCStr > 16)
        { /*Copy over reference id*/
            *tmpCStr = *cpTmpCStr;
            ++tmpCStr;
            ++cpTmpCStr;
        } /*Copy over reference id*/

        strcpy(tmpCStr, ".fastq"); /*Add in fastq ending*/

        tmpCStr = statFileCStr + lenPrefUChar;
        strcpy(tmpCStr, "--");
        tmpCStr += 2;
        cpTmpCStr = oldSam->refCStr;

        while(*cpTmpCStr > 16)
        { /*Copy over reference id*/
            *tmpCStr = *cpTmpCStr;
            ++tmpCStr;
            ++cpTmpCStr;
        } /*Copy over reference id*/

        strcpy(tmpCStr, "--stats.tsv"); /*Add in stats file ending*/

        /**************************************************************\
        * Main Sec-6 Sub-6: Check if can open the bin & stats file
        \**************************************************************/

        statFILE = fopen(statFileCStr, "a");
        fqBinFILE = fopen(binFileCStr, "a");

        if(statFILE == 0)
        { /*If can not open the stats file*/
            fprintf(
                stderr,
                "Can not create the file to hold stats (1st binning)\n"
            ); /*Let user know an error occured*/

            fprintf(
                logFILE,
                "Can not create the file to hold stats (1st binning)\n"
            ); /*Let user know an error occured*/
       
             freeStackSamEntry(&samStruct);
             freeStackSamEntry(&refStruct);
             fclose(logFILE);

             if(fqBinFILE != 0)
                 fclose(fqBinFILE);
 
             exit(1);
        } /*If can not open the stats file*/

        if(fqBinFILE == 0)
        { /*If can not open the fastq binning file*/
            fprintf(
                stderr,
                "Can not create fastq to write reads to (1st binning)\n"
            ); /*Let user know an error occured*/

            fprintf(
                logFILE,
                "Can not create fastq to write reads to (1st binning)\n"
            ); /*Let user know an error occured*/
       
             freeStackSamEntry(&samStruct);
             freeStackSamEntry(&refStruct);
             fclose(logFILE);
             fclose(statFILE);

             exit(1);
        } /*If can not open the stats file*/

        /**************************************************************\
        * Main Sec-6 Sub-7: Append reads to file & add to tree
        \**************************************************************/

        *cpTmpCStr = '\0'; /*Make ref id into a c-string*/

        tmpBin =
            insBinIntoTree(
                oldSam->refCStr,
                binFileCStr,       /*Fastq file for the bin*/
                statFileCStr,      /*Stats file for the bin*/
                &binTree, /*Root of bin tree*/
                binStack  /*Stack to use in rebalencing tree*/
        ); /*Find or add bin to tree*/

        *cpTmpCStr = '\t'; /*Convert ref id back to sam entry*/

        if(tmpBin == 0)
        { /*If a memory error occured*/
            fprintf(
                stderr,
                "Ran out of memory when binning reads first time\n"
            );

            fprintf(
                logFILE,
                "Ran out of memory when binning reads first time\n"
            );

            pclose(stdinFILE);
            fclose(logFILE);
            fclose(statFILE);
            fclose(fqBinFILE);
            freeStackSamEntry(&samStruct);
            freeStackSamEntry(&refStruct);

            exit(1);
        } /*If a memory error occured*/

        if(tmpBin->numReadsULng == 1) /*This is a new bin*/
            printStatsHeadUChar = 1;  /*Add header to stats file*/

        /*Print out the old sam entry (is not a duplicate)*/
        /*Add sequence and stats to their files*/
        samToFq(oldSam, fqBinFILE); /*Print sequence to fastq file*/

        /*Print the stats to its bin file*/
        printSamStats(oldSam, &printStatsHeadUChar, statFILE);

        fclose(fqBinFILE);
        fclose(statFILE);

        fqBinFILE = 0;  /*So program knows that no file is open*/
        statFILE = 0;   /*So program knows that no file is open*/
            
        blankSamEntry(oldSam); /*Remove old stats in sam file*/

        /*Swap pointers around (the new alignent becomes the old)*/
        tmpSam = newSam;
        newSam = oldSam;
        oldSam = tmpSam;

        /*Read in the next line*/
        errUChar = readSamLine(newSam, stdinFILE);
    } /*While their is a samfile entry to read in*/

    if(!(dupBool & 1))
    { /*If need to print out the stats for the last read*/
        if(checkRead(&readToRefMinStats, newSam) == 0)
        { /*If keeping the read*/

            if(fqBinFILE != 0)
                fclose(fqBinFILE);
            if(statFILE != 0)
                fclose(statFILE);

            /*Build the bin file name*/
            tmpCStr = binFileCStr + lenPrefUChar;
            strcpy(tmpCStr, "--");
            tmpCStr += 2;
            tmpCStr = cStrCpInvsDelm(tmpCStr, oldSam->refCStr);
            cpTmpCStr = cStrCpInvsDelm(statFileCStr, binFileCStr);

            strcpy(tmpCStr, ".fastq");        /*Add in fastq ending*/
            strcpy(cpTmpCStr, "--stats.tsv"); /*Add in fastq ending*/

            fqBinFILE = fopen(binFileCStr, "a");
            samToFq(newSam, fqBinFILE);
            fclose(fqBinFILE);
            fqBinFILE = 0;

            statFILE = fopen(statFileCStr, "a");
            printSamStats(newSam, &printStatsHeadUChar, statFILE);
            fclose(statFILE);
            statFILE = 0;
        } /*If keeping the read*/
    } /*If need to print out the stats for the last read*/

    pclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Main Sec-7: Non-reference based binning steps
    ^    main sec-7 sub-1: Check if can do memory allocation
    ^    main sec-7 sub-2: Extract reads to build consensus
    ^    main sec-7 sub-3: Build first consensuses
    ^    main sec-7 sub-4: Copy the best read to it's cluster
    ^    main sec-7 sub-5: Build the 2nd consensus & map reads
    ^    main sec-7 sub-6: Compare new consensus to old consensuses
    ^    main sec-7 sub-7: Update list of clusters in bin
    ^    main sec-7 sub-8: Add clusters to bin & move to next bin
    ^        - Before going in the bin list is an balanced tree
    ^        - after this the bins are in a list, with the left pointer
    ^          pointing towards the bins & the right pointer pointing
    ^          towards clusters (unless user did not want clustering)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Main Sec-7 Sub-1: Check if can do memory allocation
    \******************************************************************/

    /*Convert our bin tree to a list (No longer need AVL tree speed)*/
    cnvtBinTreeToList(&binTree);
    clustOn = binTree;
    lastBin = binTree; /*so I can reset pointers when removing bin*/
    tmpBin = 0;

    /*Blank the rad count file*/
    statFILE = fopen(readCntFileCStr, "w");/*File to recored counts to*/
    fprintf(statFILE, "\nBins\t*\t0\t*\t*\n");

    while(clustOn != 0)
    { /*While have reads to cluster*/
        clustUChar = 0;
        lastClust = clustOn; /*Head of cluster list*/

        fprintf(
            statFILE,
            "%s\t%lu",
            clustOn->refIdCStr,
            clustOn->numReadsULng
        ); /*Recording the number of reads in the bin*/

        if(clustOn->numReadsULng < minReadsPerBinULng)
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
 
        while(clustOn->numReadsULng >= minReadsPerBinULng)
        { /*While have reads to bin*/
            if(tmpBin == 0)
                tmpBin = malloc(sizeof(struct readBin));

            /*Make sure all bin values set to defaults/blanked*/
            tmpBin->numReadsULng = 0;
            tmpBin->refIdCStr[0] = '\0';
            tmpBin->fqPathCStr[0] = '\0';
            tmpBin->statPathCStr[0] = '\0';
            tmpBin->bestReadCStr[0] = '\0';
            tmpBin->topReadsCStr[0] = '\0';
            tmpBin->consensusCStr[0] = '\0';
            tmpBin->rightChild = 0;
            tmpBin->leftChild = 0;
            tmpBin->balUChar = 1; /*To mark keeping*/

            if(tmpBin == 0)
            { /*If errored out*/
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
            } /*If errored out*/

            /**********************************************************\
            * Main Sec-7 Sub-2: Extract reads to build consensus
            \**********************************************************/

            extractBestRead(clustOn); /*Find the best read*/

            useMapqBl = 0; /*Do not use mapping quality in extraction*/

            findBestXReads(
                &numReadConsULng,/*max number reads to build consensus*/
                &numReadsKeptULng,/*Number reads extracted*/
                &minReadReadDiff,/*Min percent simularity to keep read*/
                threadsCStr,     /*Number threads to use with minimap2*/
                &useMapqBl,         /*tells if extacting reads by mapq*/
                &readToReadMapMinStats,/*read to read map thresolds*/
                &samStruct,       /*Struct to use for reading sam file*/
                &refStruct,       /*holds the reference (0 to ignore)*/
                clustOn           /*Bin working on*/
            ); /*Find the reads to bin with the reference*/

            if(numReadsKeptULng < minReadsPerBinULng)
                continue; /*This is not a cluster*/

            strcpy(tmpBin->fqPathCStr, clustOn->topReadsCStr);

            /**********************************************************\
            * Main Sec-7 Sub-3: Build first consensuses
            \**********************************************************/

            errUChar =
                buildCon(
                    threadsCStr,    /*Number threads to use*/
                    &clustUChar,    /*Cluter number*/
                    useMajConBl,  /*Use majority (maj) consensus step?*/
                    &minBasePerFlt, /*For majority to call deleletions*/
                    &minInsPercFlt, /*% of reads to keep ins in major*/
                    minBaseQUC,     /*min q to keep base in majority*/
                    minInsQUC,      /*min q to keep ins in majority*/
                    useRaconConBl,  /*Use Racon in building consnesus?*/ 
                    &rndsRaconUC,   /*Number of rounds to run racon*/
                    useMedakaConBl, /*Use Medaka?*/ 
                    medakaModelCStr, /*Model to use with medaka*/
                    &condaBl,  /*1: Use conda medaka install 0: python*/
                    clustOn,       /*fastq file to build the consensus*/
                    &samStruct /*For reading in sequences or sam file*/
            ); /*Build the consensus (calling my wrapper)*/

            if(errUChar & 16) /*to few reads mapped to the good read*/
                continue; /*This is not a cluster*/

            /**********************************************************\
            * Main Sec-7 Sub-4: Copy the best read to it's cluster
            \**********************************************************/

            ++tmpBin->numReadsULng; /*Account for the best read*/
            tmpFILE = fopen(clustOn->bestReadCStr, "r");
            fqBinFILE = fopen(tmpBin->fqPathCStr, "a");
            tmpULng = 1; /*Reset for while loop*/

            do
            { /*While have the best read to copy over*/
                tmpULng =
                    fread(
                        refStruct.samEntryCStr,
                        sizeof(uint8_t),
                        refStruct.lenBuffULng,
                        tmpFILE
                );

                fwrite(
                    refStruct.samEntryCStr,
                    sizeof(uint8_t),
                    tmpULng,
                    fqBinFILE
                );
            } while(tmpULng != 0); /*While have best read to copy over*/
    
            fclose(tmpFILE);
            fclose(fqBinFILE);

            /**********************************************************\
            * Main Sec-7 Sub-5: Build the 2nd consensus & map reads
            \**********************************************************/

            for(unsigned int rndUI = 0; rndUI < numPolishUI; ++rndUI)
            { /*Loop till have done all the users requested polishing*/
                useMapqBl = 1;/*Use mapping quality for read selection*/

                /*Set up the consnesus as the next best read*/
                remove(clustOn->bestReadCStr);
                tmpCStr = clustOn->bestReadCStr;

                while(*tmpCStr != '\0')
                    ++tmpCStr;
                *(tmpCStr - 1) = 'a'; /*convert .fastq to .fasta*/

                rename(clustOn->consensusCStr, clustOn->bestReadCStr);
                clustOn->consensusCStr[0] = '\0'; /*Remove old name*/

                /*Rebiuld the consensus*/
                findBestXReads(
                    &numReadConsULng,/*max # reads to build consensus*/
                    &numReadsKeptULng,/*Number reads extracted*/
                    &minReadConDiff, /*Min % simularity to keep read*/
                    threadsCStr,     /*# threads to use with minimap2*/
                    &useMapqBl,      /*select best reads by mapq*/
                    &readToConMinStats,/*read to consensus min stats*/
                    &samStruct,       /*for reading sam file entries*/
                    samZeroStruct,   /*Not using consensus for scoring*/
                    clustOn           /*Bin working on*/
                ); /*Find the reads to bin with the reference*/

                if(numReadsKeptULng < minReadsPerBinULng)
                { /*If I could not map enough reads to the consensus*/
                    /*Remove the former best read*/
                    remove(tmpBin->fqPathCStr);
                    tmpBin->fqPathCStr[0] = '\0';
                    errUChar = 16;/*Let next part know problem happend*/
                    break; /*This is not a cluster*/
                } /*If I could not map enough reads to the consensus*/

                /*Rebuild the consensus using the new best reads*/
                errUChar =
                    buildCon(
                        threadsCStr,    /*Number threads to use*/
                        &clustUChar,    /*Cluter number*/
                        useMajConBl,    /*Use majority consensus step?*/
                        &minBasePerFlt,  /*For maj to call deleletions*/
                        &minInsPercFlt, /*% of reads to keep ins*/
                        minBaseQUC,     /*min q to keep base*/
                        minInsQUC,      /*min q to keep ins*/
                        useRaconConBl,  /*Use Racon?*/ 
                        &rndsRaconUC,   /*# of rounds to run racon*/
                        useMedakaConBl, /*Use Medaka?*/ 
                        medakaModelCStr, /*Model to use with medaka*/
                        &condaBl,  /*1: conda medaka install 0: python*/
                        clustOn,    /*fastq file to build consensus*/
                        &samStruct /*reading in sequences or sam file*/
                ); /*Build the consensus (calling my wrapper)*/

                if(errUChar & 16)
                { /*If could not map enough reads to build consensus*/
                    /*Remove the former best read*/
                    remove(tmpBin->fqPathCStr);
                    tmpBin->fqPathCStr[0] = '\0';
                    break; /*This is not a cluster*/
                } /*If could not map enough reads to build consensus*/
            } /*Loop till have done all the users requested polishing*/

            if(errUChar & 16)
                continue;

            remove(clustOn->bestReadCStr);/*Remove best read consensus*/
            /*Copy the consensus name to the clusters bin*/
            strcpy(tmpBin->consensusCStr, clustOn->consensusCStr);

            binReadToCon(
                &minReadConDiff,    /*Min score to keep a read*/
                &clustUChar,        /*Cluster on*/
                clustOn,            /*Bin working on*/
                tmpBin,             /*Bin to hold the cluster*/
                &samStruct,         /*To hold temporary input*/
                &readToConMinStats, /*Settings to keep read to con*/
                minimap2CmdCStr,    /*Minimap2 command*/
                tmpMiniMapCStr
            ); /*Find reads that mapp to the consensus*/

            /**********************************************************\
            * Main Sec-7 Sub-6: Compare new consensus to old consensuses
            \**********************************************************/

            if(errUChar & 16)
                continue; /*Consensus was not built*/

            *clustOn->consensusCStr = '\0';

            bestBin =
                cmpCons(
                    &minConConDiff,/*Min difference to keep consensus*/
                    tmpBin,       /*Bin with consensus to compare*/
                    clustOn,     /*Other clusters to comapre to*/
                    &samStruct,  /*Struct to hold input from minimap2*/
                    &refStruct,  /*Struct to hold input from minimap2*/
                    &conToConMinStats, /*Cons to consensus thresholds*/
                    minimap2CmdCStr,  /*Minimap2 command*/
                    tmpMiniMapCStr
            ); /*Compares a consenses to other consensuses*/

            /**********************************************************\
            * Main Sec-7 Sub-7: Update list of clusters in bin
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
            ++clustUChar; /*Move to the next cluster*/
        } /*While have reads to bin*/

        /**************************************************************\
        * Main Sec-7 Sub-8: Add clusters to bin & move to next bin
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

            if(tmpBin->numReadsULng < minReadsPerBinULng)
            { /*If not enough reads to keep*/
                binDeleteFiles(tmpBin); /*Remove its files*/
                tmpBin->balUChar = -1;
                tmpBin = tmpBin->rightChild;
                continue;
            } /*If not enough reads to keep*/

            bestBin =
                cmpCons(
                    &minConConDiff, /*Holds score of closest consensus*/
                    tmpBin,         /*Consensus to check*/
                    clustOn->leftChild, /*Other bins & their clusters*/
                    &samStruct,   /*Struct to hold input from minimap2*/
                    &refStruct,   /*Struct to hold input from minimap2*/
                    &conToConMinStats,  /*Cons to consensus thresholds*/
                    minimap2CmdCStr,   /*Minimap2 command*/
                    tmpMiniMapCStr
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
                        &minConConDiff, /*score of closest consensus*/
                        tmpBin,         /*Consensus to check*/
                        clustOn->leftChild, /*Other bins*/
                        &samStruct,         /*hold minimap2 output*/
                        &refStruct,         /*hold input from minimap2*/
                        &conToConMinStats,  /*min thresholds*/
                        minimap2CmdCStr,   /*Minimap2 command*/
                        tmpMiniMapCStr
                ); /*Compares a consenses to other consensuses*/
            } /*While have clusters with highly similar consensuses*/

            ++clustUChar;
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

    fclose(logFILE);
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
    unsigned int *numPolishUI, /*Number to times to rebuild consensus*/
    char *medakaModelCStr, /*Model to use with medaka*/
    char *useMajConBl,   /*Use majority consensus to build a consensus*/
    float *minBasePerFlt, /*Majority con min # reads to keep base*/
    float *minInsPerFlt, /*Majority con min # reads to keep insertion*/
    unsigned char *majConBaseQUC, /*Maj con min Q-score to keep base*/
    unsigned char *majConInsQUC, /*Maj con min Q-score to keep insert*/
    char *useRaconConBl, /*Use racon to build a consensus*/
    char *useMedakaConBl,/*Use medaka to build a consensus*/
    struct minDiff *minReadRefDiff, /*min difference for read to ref
                                      alignment to be different*/
    struct minDiff *minReadReadDiff, /*minReadRefDiff, for read read*/
    struct minDiff *minReadConDiff,/*minReadRefDiff, for read consenus*/
    struct minDiff *minConConDiff,/*minReadRef, for consensus consenus*/
    uint64_t *maxReadsPerConULng, /*Max reads for consensus building*/
    uint64_t *minReadsPerBinULng, /*min reads to keep a bin*/
    unsigned char *rndsRaconUC,   /*rounds to run racon*/
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
            strcpy(threadsCStr, inputCStr); /*Number of threads*/

        else if(strcmp(parmCStr, "-model") == 0)
            strcpy(medakaModelCStr, inputCStr);

        else if(strcmp(parmCStr, "-min-reads-per-bin") == 0)
            *minReadsPerBinULng = strtoul(inputCStr, &tmpCStr, 10);

        else if(strcmp(parmCStr, "-max-reads-per-con") == 0)
            *maxReadsPerConULng = strtoul(inputCStr, &tmpCStr, 10);

        else if(strcmp(parmCStr, "-rounds-racon") == 0)
            cStrToUChar(inputCStr, rndsRaconUC);
 
        else if(strcmp(parmCStr, "-extra-consensus-steps") == 0)
            cStrToUInt(inputCStr, numPolishUI);

        else if(strcmp(parmCStr, "-enable-majority-consensus") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            *useMajConBl = !*useMajConBl;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/

        else if(strcmp(parmCStr, "-enable-racon") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            *useRaconConBl = !useRaconConBl;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/

        else if(strcmp(parmCStr, "-enable-medaka") == 0)
        { /*Else if user is ussing the best read instead of consensus*/
            *useMedakaConBl = !*useMedakaConBl;
            --intArg; /*Account for this being a true or false*/
        } /*Else if user is ussing the best read instead of consensus*/
            
        else if(strcmp(parmCStr, "-maj-con-min-bases") == 0)
            sscanf(inputCStr, "%f",  minBasePerFlt);

        else if(strcmp(parmCStr, "-maj-con-min-ins") == 0)
            sscanf(inputCStr, "%f",  minInsPerFlt);

        else if(strcmp(parmCStr, "-maj-con-min-base-q") == 0)
            cStrToUChar(inputCStr, majConBaseQUC);

        else if(strcmp(parmCStr, "-maj-con-min-ins-q") == 0)
            cStrToUChar(inputCStr, majConInsQUC);
            
        /**************************************************************\
        * Fun-1 Sec-2 Sub-2: Percent difference settings
        \**************************************************************/

        else if(strcmp(parmCStr, "-read-ref-snps") == 0)
            sscanf(inputCStr, "%f", &minReadRefDiff->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-ref-diff") == 0)
            sscanf(inputCStr, "%f", &minReadRefDiff->minDiffFlt);
        else if(strcmp(parmCStr, "-read-ref-dels") == 0)
            sscanf(inputCStr, "%f", &minReadRefDiff->minDelsFlt);
        else if(strcmp(parmCStr, "-read-ref-inss") == 0)
            sscanf(inputCStr, "%f", &minReadRefDiff->minInssFlt);
        else if(strcmp(parmCStr, "-read-ref-indels") == 0)
            sscanf(inputCStr, "%f", &minReadRefDiff->minIndelsFlt);

        else if(strcmp(parmCStr, "-read-read-snps") == 0)
            sscanf(inputCStr, "%f", &minReadReadDiff->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-read-diff") == 0)
            sscanf(inputCStr, "%f", &minReadReadDiff->minDiffFlt);
        else if(strcmp(parmCStr, "-read-read-dels") == 0)
            sscanf(inputCStr, "%f", &minReadReadDiff->minDelsFlt);
        else if(strcmp(parmCStr, "-read-read-inss") == 0)
            sscanf(inputCStr, "%f", &minReadReadDiff->minInssFlt);
        else if(strcmp(parmCStr, "-read-read-indels") == 0)
            sscanf(inputCStr, "%f", &minReadReadDiff->minIndelsFlt);

        else if(strcmp(parmCStr, "-read-con-snps") == 0)
            sscanf(inputCStr, "%f", &minReadConDiff->minSNPsFlt);
        else if(strcmp(parmCStr, "-read-con-diff") == 0)
            sscanf(inputCStr, "%f", &minReadConDiff->minDiffFlt);
        else if(strcmp(parmCStr, "-read-con-dels") == 0)
            sscanf(inputCStr, "%f", &minReadConDiff->minDelsFlt);
        else if(strcmp(parmCStr, "-read-con-inss") == 0)
            sscanf(inputCStr, "%f", &minReadConDiff->minInssFlt);
        else if(strcmp(parmCStr, "-read-con-indels") == 0)
            sscanf(inputCStr, "%f", &minReadConDiff->minIndelsFlt);

        else if(strcmp(parmCStr, "-con-con-snps") == 0)
            sscanf(inputCStr, "%f", &minConConDiff->minSNPsFlt);
        else if(strcmp(parmCStr, "-con-con-diff") == 0)
            sscanf(inputCStr, "%f", &minConConDiff->minDiffFlt);
        else if(strcmp(parmCStr, "-con-con-dels") == 0)
            sscanf(inputCStr, "%f", &minConConDiff->minDelsFlt);
        else if(strcmp(parmCStr, "-con-con-inss") == 0)
            sscanf(inputCStr, "%f", &minConConDiff->minInssFlt);
        else if(strcmp(parmCStr, "-con-con-indels") == 0)
            sscanf(inputCStr, "%f", &minConConDiff->minIndelsFlt);

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

        /**************************************************************\
        * Fun-1 Sec-2 Sub-4: scoreReads read to read settings
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

        /**************************************************************\
        * Fun-1 Sec-2 Sub-5: scoreReads read to consensus settings
        \**************************************************************/

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

/*---------------------------------------------------------------------\
| Output:
|    Modifies: minMap2CmdCStr to have the contents of refCStr and fqCStr
\---------------------------------------------------------------------*/
void addMinimap2Files(
    char *minMap2CmdCStr, /*Pointer to file portion of minimap2 cmd*/
    char *refFileCStr,    /*Has name of the file with references*/
    char *readFileCStr   /*Has name of fastq file with reads*/
) /*Adds the reference and read file names to the minimap2 command*/
{ /*addMinimap2Files*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: Fun-2 Sec-1 Sub-1: addMinimap2Files
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    *minMap2CmdCStr = ' '; /*Add space between previous commands*/
    ++minMap2CmdCStr;      /*Move to point to copy*/

    while(*refFileCStr != '\0')
    { /*While have a file name to copy over*/
        *minMap2CmdCStr = *refFileCStr;
        ++minMap2CmdCStr;
        ++refFileCStr;
    } /*While have a file name to copy over*/

    *minMap2CmdCStr = ' '; /*Add space between reference and read*/
    ++minMap2CmdCStr;

    while(*readFileCStr != '\0')
    { /*While their is the read file name to copy over*/
        *minMap2CmdCStr = *readFileCStr;
        ++minMap2CmdCStr;
        ++readFileCStr;
    } /*While their is a read file name to copy over*/

    *minMap2CmdCStr = '\0'; /*Make into c-string*/
    return;
} /*addMinimap2Files*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: raconCmdCStr to hold the file names input
\---------------------------------------------------------------------*/
void addRaconFiles(
    char *raconCmdCStr,  /*Points to end of racon command*/
    char *refFileCStr,   /*Holds file with references*/
    char *samFileCStr,   /*Holds sam file to use with racon*/
    char *readsFileCStr, /*Holds reads to use*/
    char *outFileCStr    /*Holds file name to output consensus to*/
) /*Adds the input files to the racon command*/
{ /*addRaconFiles*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: addRaconFiles
    '    fun-3 sec-1: Copy the read file name
    '    fun-3 sec-2: Copy the sam file name
    '    fun-3 sec-3: Copy the reference file name
    '    fun-3 sec-4: Copy the out file name into the command
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-1: Copy the read file name into the command
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *raconCmdCStr = ' ';
    ++raconCmdCStr;

    while(*readsFileCStr != '\0')
    { /*While have to copy the read file name*/
        *raconCmdCStr = *readsFileCStr;
        ++raconCmdCStr;
        ++readsFileCStr;
    } /*While have to copy the rad file name*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-2: Copy the sam file name
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *raconCmdCStr = ' '; /*Add space between ref and sam file*/
    ++raconCmdCStr;

    while(*samFileCStr != '\0')
    { /*While have to copy the sam file name*/
        *raconCmdCStr = *samFileCStr;
        ++raconCmdCStr;
        ++samFileCStr;
    } /*While have to copy the sam file name*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-3: Copy the reference file name
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *raconCmdCStr = ' '; /*Add space for file commands*/
    ++raconCmdCStr;      /*Move off space*/

    while(*refFileCStr != '\0')
    { /*While have to copy the reference file*/
        *raconCmdCStr = *refFileCStr;
        ++raconCmdCStr;
        ++refFileCStr;
    } /*While have to copy the reference file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-4: Copy the out file name into the command
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Add redirect and space to racon command*/
    *raconCmdCStr = ' ';
    ++raconCmdCStr;
    *raconCmdCStr = '>';
    ++raconCmdCStr;
    *raconCmdCStr = ' ';
    ++raconCmdCStr;

    while(*outFileCStr != '\0')
    { /*While have to copy the read file name*/
        *raconCmdCStr = *outFileCStr;
        ++raconCmdCStr;
        ++outFileCStr;
    } /*While have to copy the rad file name*/

    *raconCmdCStr = '\0';  /*Add null to ensure is a c-string*/

    return;
} /*addRaconFiles*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        -0: If does not meet min stats
|        -1: If meets min stats
\---------------------------------------------------------------------*/
uint8_t checkRead(
    struct minAlnStats *minStats, /*Structer holding min requirments*/
    struct samEntry *samStruct    /*Structer with the alignemtn stats*/
) /*Checks if the sam entry meets the min user requirements*/
{ /*checkRead*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-4 TOC: Fun-4 Sec-1 Sub-1: checRead
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(samStruct->mapqUChar < minStats->minMapqUInt)
       return 0; /*move onto the next entry (has to low mapq)*/
           /*Minimap2 only gives MAPQ for best match*/

   if(
      samStruct->medianQFlt < minStats->minMedianQFlt &&
      *samStruct->qCStr !='*'
   )
       return 0; /*move onto the next entry (median Q-score to low)*/

   if(
      samStruct->meanQFlt < minStats->minMeanQFlt &&
      *samStruct->qCStr !='*'
   )
       return 0; /*move onto the next entry (mean Q-score to low)*/

   if(
       samStruct->readLenUInt > minStats->maxReadLenULng &&
       minStats->maxReadLenULng > 0
   )
            return 0; /*move onto the next entry (read is to long)*/

    if(
        samStruct->medianAligQFlt < minStats->minAlignedMedianQFlt &&
        *samStruct->qCStr !='*'
    )
        return 0; /*aligned median Q-score low*/

    if(
        samStruct->meanAligQFlt < minStats->minAlignedMeanQFlt &&
        *samStruct->qCStr != '*'
    )
        return 0; /*aligned mean Q-score to low*/

    if(samStruct->readAligLenUInt < minStats->minReadLenULng)
        return 0; /*Read to to short*/

    return 1;
} /*checkRead*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - readToBlank to have all stats set to 0 & c-strins to start 
|          with '\0'
\---------------------------------------------------------------------*/
void blankReadList(
    struct readStat *readToBlank
) /*Blanks a read list struct*/
{ /*blankReadList*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-5 TOC: Fun-5 Sec-1 Sub-1: blankReadList
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Blank Ids and mapping quality*/
    readToBlank->mapqUChar = 0;
    readToBlank->queryIdCStr[0] = '\0';
    readToBlank->refIdCStr[0] = '\0';

    /*Blank the Q-score stats*/
    readToBlank->medianQFlt = 0;
    readToBlank->medianAligQFlt = 0;
    readToBlank->meanQFlt = 0;
    readToBlank->meanAligQFlt = 0;

    /*Blank the read length stats*/
    readToBlank->readLenUInt = 0;
    readToBlank->readAligLenUInt = 0;

    /*Blank the total error stats*/
    readToBlank->numMatchUInt = 0;
    readToBlank->numSNPUInt = 0;
    readToBlank->numDelUInt = 0;
    readToBlank->numInsUInt = 0;

    /*Blank the kept error stats*/
    readToBlank->numKeptMatchUInt = 0;
    readToBlank->numKeptSNPUInt = 0;
    readToBlank->numKeptDelUInt = 0;
    readToBlank->numKeptInsUInt = 0;

    return;
} /*blankReadList*/

/*----------------------------------------------------------------------
# Output:
#    Modifies: newReadList to have the same values as oldReadList
----------------------------------------------------------------------*/
void cpReadList(
    struct readStat *newReadList, /*Read to copy stats to*/
    struct readStat *oldReadList /*Read to copy stats from*/
) /*Copies stats from one read list to another*/
{ /*cpReadList*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ' Fun-6 TOC: Fun-6 Sec-1 Sub-1: cpReadList
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Blank Ids and mapping quality*/
    newReadList->mapqUChar = oldReadList->mapqUChar;
    strcpy(newReadList->queryIdCStr, oldReadList->queryIdCStr);
    strcpy(newReadList->refIdCStr, oldReadList->refIdCStr);

    /*Blank the Q-score stats*/
    newReadList->medianQFlt = oldReadList->medianQFlt;
    newReadList->medianAligQFlt = oldReadList->medianAligQFlt;
    newReadList->meanQFlt = oldReadList->meanQFlt;
    newReadList->meanAligQFlt = oldReadList->meanAligQFlt;

    /*Blank the read length stats*/
    newReadList->readLenUInt = oldReadList->readLenUInt;
    newReadList->readAligLenUInt = oldReadList->readAligLenUInt;

    /*Blank the total error stats*/
    newReadList->numMatchUInt = oldReadList->numMatchUInt;
    newReadList->numSNPUInt = oldReadList->numSNPUInt;
    newReadList->numDelUInt = oldReadList->numDelUInt;
    newReadList->numInsUInt = oldReadList->numInsUInt;

    /*Blank the kept error stats*/
    newReadList->numKeptMatchUInt = oldReadList->numKeptMatchUInt;
    newReadList->numKeptSNPUInt = oldReadList->numKeptSNPUInt;
    newReadList->numKeptDelUInt = oldReadList->numKeptDelUInt;
    newReadList->numKeptInsUInt = oldReadList->numKeptInsUInt;

    return;
} /*cpReadList*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Prints:                                                         |
|        - The read with the highest medain Q-score to the fastq file  |
|          name stored in binIn->bestReadCStr                          |
|    - Modifies:                                                       |
|        - Fastq file binIn-fqPathCStr to not have the best read       |
|        - Stats file binIn-statsPathCStr to not have the best read    |
|    - Returns:                                                        |
|        - 1: If sucessfull                                            |
|        - 2: For blank structer                                       |
|        - 4: For no fastq file                                        |
|        - 8: For no stats file                                        |
|        - 16: For error when extracting stats                         |
|        - 32: If could not open a temporary fastq file                |
|        - 64: If could not copy the fastq reads to their bins         |
\---------------------------------------------------------------------*/
uint8_t extractBestRead(
    struct readBin *binIn        /*Bin to extract best read from*/
) /*Splits bin into read with highest Q-score & other all other reads*/
{ /*extractBestRead*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-7 TOC: extractBestRead
    '    fun-7 sec-1: variable declerations                            \
    '    fun-7 sec-2: Check if Bin exists & set up best read file name /
    '    fun-7 sec-3: Check if fastq & stats file can be opened        \
    '    fun-7 sec-4: Open temporary stat file & write header          /
    '    fun-7 sec-5: Copy old stat file, except for best read         \
    '    fun-7 sec-6: Make temporay stat file the new bin stat file    /
    '    fun-7 sec-7: Extract the best read                            \
    '    fun-7 sec-8: Clean up, close files & make tmp file bin fq file/
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-1: variable declerations                               v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpFqCStr ="tmp-202301110827-reads-09876543211234567890.fastq",
        *tmpStatCStr ="tmp-202301110827-stats-09876543211234567890.tsv",
        *tmpCStr = 0;

    int8_t ignoreC = 0;

    uint8_t
        errUChar = 0,     /*Holds error messages*/
        onHeaderBool = 1; /*Tells if need to ignore the first line*/

    struct readStat
        tmpRead,
        bestRead;            /*Holds best read to extract*/

    FILE 
        *inFILE = 0,
        *outFILE = 0,      /*Holds the read to polish with*/
        *otherOutFILE = 0; /*Holds every read except the polish read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-7 Sec-2: Check if Bin exists & set up best read file name    v
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(binIn == 0)
        return 2;    /*The structure has nothing*/

    strcpy(binIn->bestReadCStr, binIn->fqPathCStr); /*Copy the file name*/

    tmpCStr = binIn->bestReadCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    tmpCStr -= 6; /*Get on '.' part of ".fastq"*/

    strcpy(
        tmpCStr,
        "--best-read.fastq"
    ); /*Add in the best read ending*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-7 Sec-3: Check if can open fastq file & stats file           v
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    inFILE = fopen(binIn->statPathCStr, "r"); /*Open the stats file*/

    if(inFILE == 0)
        return 4;    /*The fastq file has nothing*/

    outFILE = fopen(binIn->fqPathCStr, "r"); /*Open the fastq file*/

    if(outFILE == 0)
    { /*If can not open the fastq file*/
        fclose(outFILE);
        return 8;
    } /*If can not open the fastq file*/

    fclose(outFILE); /*Was just a quick test to see if could open*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-4: Open temporary stat file & write header             v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the temporary stat file to write non-best read to*/
    outFILE = fopen(tmpStatCStr, "w");

    fprintf(outFILE, "Read\tRef\tMAPQ\treadLength\talignedLength");
    fprintf(outFILE, "\tmatches\tkeptMatches\tmismatches\tinsertions"); 
    fprintf(outFILE, "\tdeletions\tmedianQ\tmeanQ\talignedMedianQ"); 
    fprintf(outFILE, "\talignedMeanQ\tTotalMismatches"); 
    fprintf(outFILE, "\tTotalInsertions\tTotalDeletions\n");

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-5: Copy old stat file, except for best read            v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /*Read in the frist line (assume is best read*/
    errUChar = readStatsFileLine(inFILE, &onHeaderBool, &bestRead);
        /*onHeaderBool tells if have read in the header*/

    if(!(errUChar & 1))
        return 16;     /*File error, return 16 so user knows*/

    /*Read in next line, so have something to compare*/
    errUChar = readStatsFileLine(inFILE, &onHeaderBool, &tmpRead);

    if(!(errUChar & 1))
        return 16;     /*File error, return 16 so user knows*/
   
    while(errUChar & 1)
    { /*While not at end of file or no problems*/
        ignoreC =
            (int8_t) (
              ((short) tmpRead.mapqUChar - (short) bestRead.mapqUChar)
              >> (sizeof(short) << 3)
            );
            /*tmpRead - bestRead is negative if bestRead > tmpRead
                - converting to hsort, so can have negatives
              (tmp - best) >> bitsInSht: Keeps only the negative bit
            */

        ignoreC |=
            (int8_t) (
                (
                  (int32_t) tmpRead.medianQFlt -
                  (int32_t) bestRead.medianQFlt
                )
              >> ((sizeof(int32_t) << 3) - 1) /*One off max size*/
            ) & ignoreC;
            /* tmp - best: is negative if best > tmp
               (tmp - best) >> 32: 1 if best > tmp; 0 if tmp>best
                   - change to sizeof(.* to deal with complier warnings
               & ingnoreUC: sets to 0 if mapq was best
            */

        ignoreC |=
            (int8_t) (
                (
                  (int32_t) tmpRead.readLenUInt -
                  (int32_t) bestRead.readLenUInt
                )
              >> ((sizeof(int32_t) << 3) - 1) /*One of max size*/
            ) & ignoreC;
            /* (int64_t) tmp - (int64_t) best: is negative if best > tmp
               (tmp - best) >> 64: 1 if negative, 0 if positive
               & ingnoreUC: sets to 0 if mapq was best
            */

        /*At this point ignoreC is 1 (discard) read or 0 (keep)*/
        if(ignoreC & 1)
            printReadStat(&tmpRead, outFILE);
        else
        { /*else have a new best mapq*/
            printReadStat(&bestRead, outFILE);
            cpReadList(&bestRead, &tmpRead); /*copy the read*/
        } /*else have a new best mapq*/

        errUChar =
            readStatsFileLine(
                inFILE,     /*File with line to grab*/
                &onHeaderBool, /*Tells if their is a header line*/
                &tmpRead   /*Will hold the stats from the stats file*/
        ); /*Read in the next line*/
    } /*While not at end of file or no problems*/

    if(errUChar != 2) 
        return 16; /*File error*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-6: Make the temporay stat file the new bin stat file   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

     /*Close the open files (no longer need to read from)*/
     fclose(inFILE);
     fclose(outFILE);

     remove(binIn->statPathCStr); /*Remove the old stats file*/

     /*Make temporary stats file without best read the stats file*/
     rename(tmpStatCStr, binIn->statPathCStr);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-7: Extract the best read                               v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the fastq file (already check if could open)*/
    inFILE = fopen(binIn->fqPathCStr, "r");

    /*Open the temporary file to hold the best read*/
    outFILE = fopen(binIn->bestReadCStr, "w");

    if(outFILE == 0)
        return 32;

    /*File to hold all reads except the best read*/
    otherOutFILE = fopen(tmpFqCStr, "w");

    if(otherOutFILE == 0)
        return 32;

    if(
        fqOneIdExtract(
            bestRead.queryIdCStr,
            inFILE,     /*fastq file to extract read from*/
            outFILE,
            otherOutFILE
        ) != 1
    ) /*If could not extract the best read*/
        return 64;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-8: Clean up, close files & make tmp file bin fq file   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    --binIn->numReadsULng; /*Account for the extracted best read*/

    /*No longer need files open*/
    fclose(inFILE);
    fclose(outFILE);
    fclose(otherOutFILE);

    remove(binIn->fqPathCStr); /*Remove the old file*/

    /*Make the tempoaray fastq with all but teh best read to
      the bin fastq file*/
    rename(tmpFqCStr, binIn->fqPathCStr);

    return 1; /*Success*/
} /*extractBestRead*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|    - Modifies:                                                       |
|        - readStruct to have stats from the next line in statsFILE    |
\---------------------------------------------------------------------*/
uint8_t readStatsFileLine(
    FILE *statsFILE,            /*File with line to grab*/
    uint8_t *onHeaderBool,      /*1: skip one line, 0: grab first line*/
    struct readStat *readStruct /*Holds the stats from the stats file*/
) /*Reads single line from printSamStats function in samEntryStruct.c*/
{ /*readStatsFileLine*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-8 TOC: readStatsFileLine
    '    fun-8 sec-1: variable declerations
    '    fun-8 sec-2: Initalize and read in first line
    '    fun-8 sec-3: If need to, read past header
    '    fun-8 sec-4: Copy read and reference id's
    '    fun-8 sec-5: Get mapq and read lengths
    '    fun-8 sec-6: Get number of matches and mismatches
    '    fun-8 sec-7: Get number of insertions and deletions
    '    fun-8 sec-8: Get Q-scores
    '    fun-8 sec-9: Get aligned Q-scores
    '    fun-8 sec-10: Get total number of SNPs, insertions, & deletions
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint16_t lenBuffUSht = 1024; /*is greater than one full line*/

    char
        buffCStr[lenBuffUSht],
        *tmpCStr = 0,   /*String to copy from*/
        *cpTmpCStr = 0, /*String to copy to*/
        *eofCStr = 0;   /*Tells me if at end of file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-2: Initalize and read in first line
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankReadList(readStruct); /*Make sure no leftover data*/ 

    if(statsFILE == 0)
        return 8;

    eofCStr =
        fgets(
            buffCStr,
            lenBuffUSht,
            statsFILE
    ); /*Read in the line*/

    if(eofCStr == 0)
        return 2;         /*End of file, so no line to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-3: If need to, read past header
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*onHeaderBool == 1)
    { /*If on the header, move past*/
        /*Read in the next line*/
        eofCStr = fgets(buffCStr, lenBuffUSht, statsFILE);

        if(eofCStr == 0)
            return 2;   /*End of file, so only a header*/

        *onHeaderBool = 0; /*No longer on the header*/
    } /*If on the header, move past*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-4: Copy read and reference id's
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = buffCStr;
    cpTmpCStr = readStruct->queryIdCStr;

    while(*tmpCStr > 31)
    { /*While have a read id to copy*/
        *cpTmpCStr = *tmpCStr;
        ++cpTmpCStr;
        ++tmpCStr;

        if(*tmpCStr == '\0')
            return 16; /*Early end of file or line*/
    } /*While have a read id to copy*/

    *cpTmpCStr = '\0'; /*make sure a c-string*/
    ++tmpCStr; /*Get off the tab*/

    cpTmpCStr = readStruct->refIdCStr;

    while(*tmpCStr > 31)
    { /*While have a reference id to copy*/
        *cpTmpCStr = *tmpCStr;
        ++cpTmpCStr;
        ++tmpCStr;

        if(*tmpCStr == '\0')
            return 16; /*Early end of file or line*/
    } /*While have a reference id to copy*/

    *cpTmpCStr = '\0'; /*make sure a c-string*/
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-5: Copy mapq and read lengths
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = cStrToUChar(tmpCStr, &readStruct->mapqUChar); /*Get mapq*/
    ++tmpCStr; /*Get off the tab*/

    /*Get read length*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->readLenUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the algined read length*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->readAligLenUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-6: Get number of matches and mismatches
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of matches*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numMatchUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept matches*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numKeptMatchUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept SNPs*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numKeptSNPUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-7: Get number of insertions and deletions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of kept insertions*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numKeptInsUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of kept deletions*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numKeptDelUInt);
    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-8 Sec-8: Get mean and median aligned Q-scores
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    sscanf(tmpCStr, "%f", &readStruct->meanQFlt); /*Get mean Q-score*/

    /*Move to the next entry*/
    while(*tmpCStr > 32)
        ++tmpCStr;

    ++tmpCStr; /*Get off the tab*/

    sscanf(tmpCStr, "%f", &readStruct->medianQFlt); /*Get median Q*/

    /*Move to the next entry*/
    while(*tmpCStr > 32)
        ++tmpCStr;

    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-9: Get median and mean aligned Q-scores
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the mean aligned Q-score*/
    sscanf(tmpCStr, "%f", &readStruct->meanAligQFlt);

    /*Move to the next entry*/
    while(*tmpCStr > 32)
        ++tmpCStr;

    ++tmpCStr; /*Get off the tab*/

    /*Get the median aligned Q-score*/
    sscanf(tmpCStr, "%f", &readStruct->medianAligQFlt);

    /*Move to the next entry*/
    while(*tmpCStr > 32)
        ++tmpCStr;

    ++tmpCStr; /*Get off the tab*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-10: Get total number of SNPs, insertions, & deletions
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of SNPs*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numSNPUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of insertions*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numInsUInt);
    ++tmpCStr; /*Get off the tab*/

    /*Get the number of deletions*/
    tmpCStr = cStrToUInt(tmpCStr, &readStruct->numDelUInt);
    ++tmpCStr; /*Get off the tab*/

    return 1; /*Sucess*/
} /*readStatsFileLine*/

/*----------------------------------------------------------------------
# Output:
#    Prints: fastq entry to outFILE if finds
#    Returns:
#        - 1: if found and printed the id
#        - 2: If could not find the id (no printing)
#        - 4: If the fqFILE does not exist
#        - 8: If the outFILE does not exist
#        - 16: If an incomplete entry (EOF, but missing lines)
#        - 32: If the idCStr (id looking for) entry is incomplete
----------------------------------------------------------------------*/
uint8_t fqOneIdExtract(
    char *idCStr,  /*'\0' terminated read id to extract*/
    FILE *fqFILE,     /*fastq file to extract read from*/
    FILE *keptFILE,   /*File with the target read*/
    FILE *outFILE     /*fastq file to write read to*/
) /*Extracts one read id from a fastq file*/
{ /*fqOneIdExtract*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-9 TOC: fqOneIdExtract
    '    fun-9 sec-1: variable declerations
    '    fun-9 sec-2: Check if valid file and set up for read in
    '    fun-9 sec-3: Read in the header and compare read ids
    '    fun-9 sec-4: If not a match, move past sequence & spacer line
    '    fun-9 sec-5: If not a match, move past q-score lines
    '    fun-9 sec-6: Is match, print out header
    '    fun-9 sec-7: Is match, print sequence and spacer line
    '    fun-9 sec-8: Is match, print Q-score line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint16_t
        numSeqLineUSht = 0, /*Number of sequence lines in fastq file*/
        lenBuffUSht = 1024;

    char
        *tmpCStr = 0,
        *idIterCStr = 0,
        lineCStr[lenBuffUSht];

    FILE
        *tmpFILE = outFILE; /*File to write read to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Check if valid file and set up for read in
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(fqFILE == 0)
        return 4;
    if(outFILE == 0)
        return 8;
    if(keptFILE == 0)
        return 16;

    /*Make sure have marks set to determine if read in full line*/
    lineCStr[lenBuffUSht - 1] = '\0';
    lineCStr[lenBuffUSht - 2] = '\0';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: Read in the header and compare read ids
    #    fun-9 sec-3 sub-1: Read in header and compare ids
    #    fun-9 sec-3 sub-2: Finsh reading in header if not a match
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*******************************************************************
    # Fun-9 Sec-3 Sub-1: Read in header and compare ids
    *******************************************************************/

    while(fgets(lineCStr, lenBuffUSht, fqFILE))
    { /*While have a line to read in*/
        tmpCStr = lineCStr;
        idIterCStr = idCStr;

        /*Make sure off @ symbols for header*/
        if(*tmpCStr == '@')
            ++tmpCStr;

        if(*idIterCStr == '@')
            ++idIterCStr;

        while(*tmpCStr == *idIterCStr)
        { /*While the two strings are the same*/
            ++tmpCStr;
            ++idIterCStr;

            if(*idIterCStr == '\0') /*idIterCStr always ends with null*/
                break; /*Done with comparision*/
        } /*While the two strings are the same*/

        if(*idIterCStr == '\0')
            tmpFILE = keptFILE;
        else
            tmpFILE = outFILE;

        /*Print out the read in part of the header*/
        fprintf(tmpFILE, "%s", lineCStr);

        /***************************************************************
        # Fun-9 Sec-3 Sub-2: Finsh reading in header if not a match
        ***************************************************************/

        /*The header could have tabs, so can not do < 16*/
        /*I can not just seq till +\n, because Illumina fastq's break 
          into two lines, So I need to find the number of sequence lines
        */
        while(
            !(lineCStr[lenBuffUSht - 2] == '\0' ||
            lineCStr[lenBuffUSht - 2] == '\n')
        ) { /*While not on the next line*/
            lineCStr[lenBuffUSht - 2] = '\0';

            /*Grab next part of header from file*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/

            /*Print out next part of the header*/
            fprintf(tmpFILE, "%s", lineCStr);
        } /*While not on the next line*/

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Fun-9 Sec-4: If not a match, move past sequence & spacer line
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        lineCStr[lenBuffUSht - 2] = '\0'; /*Set up the maker*/

        /*Grab the first part of the sequence*/
        tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);
        numSeqLineUSht = 0; /*reset to 0 (number sequence lines)*/

        do { /*While not on the spacer entry*/

            /*Print out the read in part of the header*/
            fprintf(tmpFILE, "%s", lineCStr);

            /*I can get away with < 16 here because I know their are 
              no tabs in the sequence entry*/
            if(lineCStr[lenBuffUSht - 2] < 16)
                ++numSeqLineUSht; /*Record the number of lines read in*/

            lineCStr[lenBuffUSht - 2] = '\0'; /*Reset EOL marker*/

            /*Get next part of entry before spacer*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/

        } while(lineCStr[0] != '+'); /*While not on spacer*/
        /*fastq sequence line ends with \n+\nq-score line*/

        fwrite("+\n", sizeof(uint8_t), 2, tmpFILE);

        while(lineCStr[lenBuffUSht - 2] > 16)
        { /*While have extra heade entries*/
            lineCStr[lenBuffUSht - 2] = '\0'; /*reset EOL marker*/
            tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

            if(tmpCStr == 0)
                 return 32;    /*Premature end of file*/
        } /*While have extra heade entries*/
        /*Just to remove extra user formating on spacer*/

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Fun-9 Sec-5: If not a match, move past q-score lines
        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        for(uint16_t uShtCnt = 0; uShtCnt < numSeqLineUSht; ++uShtCnt)
        { /*Loop through all lines in the q-socre entry*/

             /*the Q-score entry will not have tabs, so < 16 works*/
        
             do { /*While not on the next line*/
                lineCStr[lenBuffUSht - 2] = '\0'; /*Reset marker*/

                /*Get next part of Q-score entry*/
                tmpCStr = fgets(lineCStr, lenBuffUSht, fqFILE);

                if(tmpCStr == 0)
                     return 32;    /*Premature end of file*/

                /*Print out the next part of the q-score entry*/
                fprintf(tmpFILE, "%s", lineCStr);
            } while(lineCStr[lenBuffUSht - 2] > 16); /*On next line?*/
        } /*Loop through all lines in the q-socre entry*/
    } /*While have a line to read in*/ 

    return 1; /*Sucess*/
} /*fqOneIdExtract*/

/*---------------------------------------------------------------------\
| Output:                                                              |
|   Prints: Prints out variables in readToPrint structer               |
\---------------------------------------------------------------------*/
void printReadStat(
    struct readStat *readToPrint, /*Read to print stats out for*/
    FILE *outFILE                  /*File to print read to*/
) /*Prints out the stats in a readStat structure to a file*/
{ /*printReadStat*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-10 TOC: Sec-1 Sub-1: printReadStat                           /
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(readToPrint == 0)
        return;

    if(outFILE == 0)
        return;

    /*Print out the entry stats*/
    fprintf(
        outFILE,
        "%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%f\t%f\t%f",
        readToPrint->queryIdCStr, 
        readToPrint->refIdCStr,
        readToPrint->mapqUChar,
        readToPrint->readLenUInt,
        readToPrint->readAligLenUInt,
        readToPrint->numMatchUInt,
        readToPrint->numKeptMatchUInt,
        readToPrint->numKeptSNPUInt,
        readToPrint->numKeptInsUInt,
        readToPrint->numKeptDelUInt,
        readToPrint->medianQFlt,
        readToPrint->meanQFlt,
        readToPrint->medianAligQFlt,
        readToPrint->meanAligQFlt
    ); /*1st printf: print out stats*/

    fprintf(
        outFILE,
        "\t%u\t%u\t%u\n",
        readToPrint->numSNPUInt,
        readToPrint->numInsUInt,
        readToPrint->numDelUInt
    ); /*2nd fprintf: print out stats*/

    return;
} /*printReadStat*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|     - 1: if succeeded or no reads mapped (numReadsKept=0 for no reads)
|     - 4: if could not read reference
|     - 8: if could not open the fastq file
|     - 16: if minimap2 errored out or returned nothing
| Note:
|     - Score: percMult * (keptSNPs + keptIns + keptDels) / read length
|     - minSimUSht ranges from 1 (0.01%) to precMult (100%)            
\---------------------------------------------------------------------*/
uint8_t findBestXReads(
    const uint64_t *numReadConsULng, /*# reads for bulding a consensus*/
    uint64_t *numReadsKeptULng,  /*Number of reads binned to con*/
    struct minDiff *maxDifference,
    char *threadsCStr,           /*Number threads to use with minimap2*/
    const char *useMapqBl,       /*1: use mapping quality in selection*/
    struct minAlnStats *minStats,/*Min stats to cluster reads together*/
    struct samEntry *samStruct,  /*Struct to use for reading sam file*/
    struct samEntry *refStruct,  /*holds the reference (0 to ignore)*/
    struct readBin *binTree      /*Bin working on*/
) /*Extract the top reads that mapped to the selected best read*/
{ /*findBestXReads*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   \\ Fun-11 TOC: findBestXReads                                       /
   //    fun-11 sec-1: variable declerations                           \
   \\    fun-11 sec-2: Check if can open files                         /
   //    fun-11 sec-3: Map reads to best read & select top x           \
   \\    fun-11 sec-4: Trim, score, & select best mapped reads         /
   //    fun-11 sec-5: Set up the best x read file name                \
   \\    fun-11 sec-6: Set up hash table to extract reads              /
   //    fun-11 sec-7: Extract reads with fastq greps hash extract     \
   \\    fun-11 sec-8: Clean up                                        /
   //                                                                  \
   \\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-1: variable declerations                              v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    int32_t 
        lenBuffUInt = 1 << 15; /*wil make a 65536 byte array*/

    uint8_t
        errUChar = 0, /*Holds error output*/
        oneUChar = 1,
        digPerKeyUChar = 0,
        zeroUChar = 0;

    char
        minimapCmdCStr[1024],
        *tmpCStr = 0,
        buffCStr[lenBuffUInt];  /*Buffer to extract reads with*/

    int32_t
        lenIdUInt = 100;  /*number of characters allowed for read id*/

    uint16_t
        topScoreUSht = 128 + (2020 * *useMapqBl), /*Max score*/
        lowScoreUSht = topScoreUSht, /*Lowest scoring kept read*/
        scoreUSht = 0;    /*Score of a single read*/

    uint64_t
        hashSizeULng = 0;     /*Size of hash table for read extraction*/

    unsigned long 
        majicNumULng = 0;     /*Number to use with hashing*/

    struct readInfo
        **hashTbl = 0,     /*Hash table of read ids for extraction*/
        *tmpRead = 0,
        readScoreAry[*numReadConsULng], /*Stack for easer work*/
        *readOn = readScoreAry,
        *scoresArray[topScoreUSht]; /*look up table for scores*/

    struct readNodeStack
        searchStack[200];  /*Used for fqGetIds*/

    FILE 
        *testFILE = 0,
        *bestReadsFILE = 0,  /*Holds the read to polish with*/
        *stdinFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-11 Sec-2: Check if can open files and copy reference
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(binTree == 0)
        return 2;    /*The structure has nothing*/

    /*Open the temporary file to hold the best read*/
    testFILE = fopen(binTree->bestReadCStr,"r");

    if(testFILE == 0)
        return 4;

    if(refStruct != 0)
    { /*If using the reference for deletion scoring*/
        blankSamEntry(refStruct);

        /*Read in the reference sequence*/
        errUChar = readRefFqSeq(testFILE, refStruct);

        if(!(errUChar & 1))
            return 4;
    } /*If using the reference for deletion scoring*/

    fclose(testFILE); /*No longer need open*/

    /*See if the fastq file can be opened*/
    testFILE = fopen(binTree->fqPathCStr, "r");

    if(testFILE == 0)
        return 8;     /*Can not open the fastq file*/

    fclose(testFILE); /*Do not need open right know*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-3: Map reads to best read & select top x              v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *numReadsKeptULng = 0; /*Make sure start at 0 reads*/

    /*Build the command to run minimap2*/
    tmpCStr = cStrCpInvsDelm(minimapCmdCStr, minimap2CMD);
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    tmpCStr =
       cpParmAndArg(tmpCStr,binTree->bestReadCStr, binTree->fqPathCStr);
    
    stdinFILE = popen(minimapCmdCStr, "r"); /*Get the minimap2 results*/

    blankSamEntry(samStruct); /*Make sure start with blank*/

    /*Read in a single sam file line to check if valid (header)*/
    errUChar = readSamLine(samStruct, stdinFILE);

    if(!(errUChar & 1))
    { /*If an error occured*/
        pclose(stdinFILE);
        return 16;
    } /*If an error occured*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-4: Trim, score, & select best mapped reads            v
    ^    fun-11 sec-4 sub-1: Trim read & remove if no sequence         v
    ^    fun-11 sec-4 sub-2: Score read                                v
    ^    fun-11 sec-4 sub-3: Check if Score meets requirements         v
    ^    fun-11 sec-4 sub-4: Check if at max number of reads to keep   v
    ^    fun-11 sec-4 sub-5: Check if read beats lowest kept read scorev
    ^    fun-11 sec-4 sub-6: Check if decided to keep the read         v
    ^    fun-11 sec-4 sub-7: Add kept read to the score list           v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    | Fun-11 Sec-4 Sub-1: Trim read & remove read if no sequence        |
    \******************************************************************/

    while(errUChar & 1)
    { /*While their is a samfile entry to read in*/
        if(*samStruct->samEntryCStr == '@')
        { /*If was a header*/
            errUChar =
                readSamLine(
                    samStruct,
                    stdinFILE
            ); /*If header read new line*/

            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        /*Convert & print out sam file entry*/
        errUChar = trimSamEntry(samStruct);

        if(errUChar >> 2)
        { /*If entry did not have a sequence, discard*/
            errUChar =
                readSamLine(
                    samStruct,
                    stdinFILE
            ); /*Read the next sam entry*/

            continue;
        } /*If entry did not have a sequence, discard*/

        /**************************************************************\
        * Fun-11 Sec-4 Sub-2: Score read
        \**************************************************************/

        findQScores(samStruct); /*Find the Q-scores*/

        scoreAln(
            minStats,
            samStruct,
            refStruct,   /*Reference struct to score deletions with*/ 
            &oneUChar,   /*Mapped read has Q-score*/
            &oneUChar    /*Make sure set to one if using reference*/
        ); /*Score the alignment*/

        /**************************************************************\
        * Fun-11 Sec-4 Sub-3: Check if Score meets requirements
        \**************************************************************/

           
        if(
            samStruct->flagUSht & 4 || /*Is an unmapped read*/
            samStruct->mapqUChar < minStats->minMapqUInt ||
            !(checkIfKeepRead(maxDifference, samStruct) & 1)
        ) { /*If the read Is to different from the reference*/
            /*Move to the next entry*/
            blankSamEntry(samStruct); /*Make sure start with blank*/
            errUChar = readSamLine(samStruct, stdinFILE);
            continue;
        } /*If the read Is to different from the reference*/

        /**************************************************************\
        * Fun-11 Sec-4 Sub-4: Check if at max number of reads to keep
        \**************************************************************/

        tmpRead = 0;

        /*Discard decimal part of the Q-score*/
        scoreUSht =
           (((uint16_t) samStruct->mapqUChar << 2)
              & (uint16_t) *useMapqBl
           ) + (uint16_t) samStruct->medianQFlt;
        /*Logic:
            mapping quality goes to 0 (ignored) if useMapqBl is 0.
            mapqUChar is 8 bits (1111 1111), mapqUChar << 2 max is 2027.
            The median Q-score is 7 bits (11 1111), with a max of 127).
            mapq << 2:  11 1111 1100
            median Q:       111 1111
                        || |||| ||||
                              16 4 1
                             32 8 2
                            64
            A median Q-score >= 4 will outweigh mapqs < 2
            A median Q-score >= 8 will outwight a mapq < 4
            A median Q-score >= 16 will outwight a mapq < 8
            A median Q-score >= 32 will outwight a mapq < 16
            A median Q-score >= 64 will outwight a mapq < 32
            A mapq >= 32 will outweigh all median Q-scores*/

        if(*numReadsKeptULng < *numReadConsULng)
        { /*If still accepting new reads*/
            readOn->leftChild = 0;
            readOn->balanceChar = 0;

            if(scoreUSht < lowScoreUSht) /*Check if new lowest score*/
                lowScoreUSht = scoreUSht;

            if(scoresArray[lowScoreUSht] == 0)
            { /*If this is the first time I got this score*/
                scoresArray[lowScoreUSht] = readOn;
                readOn->rightChild = 0; /*Mark end of lifo*/
            } /*If this is the first time I got this score*/

            else
            { /*Else I already have reads with the same score*/
                tmpRead = scoresArray[lowScoreUSht];
                scoresArray[lowScoreUSht] = readOn;
                readOn->rightChild = tmpRead;
            } /*Else I already have reads with the same score*/

            readOn->idBigNum =
                    makeBigNumStruct(samStruct->queryCStr, &lenIdUInt);

            ++readOn; /*Move to next open read*/
            ++(*numReadsKeptULng);
        } /*If still accepting new reads*/

        /**************************************************************\
        | Fun-11 Sec-4 Sub-5: Check if read beets lowest kept read score|
        \**************************************************************/

        else if(scoreUSht > lowScoreUSht)
        { /*Else if only keeping better reads*/

            tmpRead = scoresArray[lowScoreUSht];

            /*Check if still have other reads at the low score*/
            if(tmpRead != 0 && tmpRead->rightChild != 0)
                scoresArray[lowScoreUSht] = tmpRead->rightChild;
                     /*Mark remove the first low scoring read*/
            else
            { /*else need to find the next lowest score*/
                scoresArray[lowScoreUSht] = tmpRead->rightChild;
                ++lowScoreUSht; /*Move to the next score*/

                while(scoresArray[lowScoreUSht] == 0)
                    ++lowScoreUSht;  /*Move to read with higher score*/
            } /*else need to find the next lowest score*/

            /*Check if need to build the stack*/
            if(scoresArray[scoreUSht] != 0)
                tmpRead->rightChild = scoresArray[scoreUSht];
            else
                tmpRead->rightChild = 0;

            strToBackwardsBigNum(
                tmpRead->idBigNum,
                samStruct->queryCStr,
                &lenIdUInt
            ); /*Convert query id to a big number*/

            scoresArray[scoreUSht] = tmpRead;
        } /*Else if only keeping better reads*/

        /*Move to the next read*/
        blankSamEntry(samStruct); /*Make sure start with blank*/
        errUChar = readSamLine(samStruct, stdinFILE);
    } /*While their is a samfile entry to read in*/

    pclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-5: Set up the best x read file name                   v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Free the big numbers stored in the readInfo structs*/
    if(*numReadsKeptULng == 0)
        return 1;

    strcpy(binTree->topReadsCStr, binTree->fqPathCStr);
    tmpCStr = binTree->topReadsCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    tmpCStr -= 6; /*Get to . in .fastq*/
    strcpy(tmpCStr, "--top-reads.fastq"); /*Copy in top read ending*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-6: Set up hash table to extract reads                 v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpRead = readScoreAry;
    tmpRead->rightChild = 0; /*Ensure no circular lists*/
    readOn = tmpRead + 1;

    for(uint32_t intRead = 1; intRead < *numReadsKeptULng; ++intRead)
    { /*For all kept reads, set built the list for the hash table*/
       readOn->rightChild = tmpRead;
       tmpRead = readOn;
       ++readOn; /*move to the next read in the array*/
    } /*For all kept reads, set built the list for the hash table*/

    hashTbl = 
        readListToHash(
            tmpRead,
            numReadsKeptULng,
            searchStack,        /*Used for searching the hash table*/
            &hashSizeULng,      /*Will hold Size of hash table*/
            &digPerKeyUChar,    /*Power of two hash size is at*/
            &majicNumULng       /*Holds majick number for kunths hash*/
    ); /*Build the hash table*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-7: Extract reads with fastq greps hash extract        v
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open fastq file and file to store the best reads in*/
    testFILE = fopen(binTree->fqPathCStr, "r");
    bestReadsFILE = fopen(binTree->topReadsCStr, "w");

    tmpRead = 0; /*So extract reads knows not doing AVL tree search*/
 
    extractReads(
        testFILE,       /*File with reads to extract*/
        bestReadsFILE,  /*File to write reads to*/
        buffCStr,
        lenBuffUInt,
        majicNumULng,   /*Holds majick number for kunths hash*/
        digPerKeyUChar, /*Power of two hash size is at*/
        &zeroUChar,      /*Print the matches*/
        tmpRead,
        hashTbl
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-11 Sec-8: Clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(uint64_t uLngCnt = 0; uLngCnt < *numReadsKeptULng; ++uLngCnt)
        freeBigNumStruct(&readScoreAry[uLngCnt].idBigNum);

    free(hashTbl); /*Free the hash table pointers*/

    fclose(testFILE);
    fclose(bestReadsFILE);

    return 1;
} /*findBestXReads*/

/*---------------------------------------------------------------------\
| Output:
|    Uses: Racon to build a consensus (file name in bin->consensusCStr)
\---------------------------------------------------------------------*/
void buildConWithRacon(
    const uint8_t *clustUChar,         /*Cluster on*/
    const unsigned char *rndsRaconUC,  /*Number of rounds to run racon*/
    const char *useFqConBl,    /*1 consensus is fastq, 0 use best read*/
    char *threadsCStr,             /*Number threads to use with racon*/
    struct readBin *conBin            /*Bin working on*/
) /*Builds a consensus using racon*/
{ /*buildConWithRacon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ' Fun-12 TOC: buildConWithRacon
    '    fun-12 sec-1: Variable declerations
    '    fun-12 sec-2: Set up the consensus name
    '    fun-12 sec-3: Build the consensus using racon
    '    fun-12 sec-4: Rename the consensus if needed
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-12 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *refFileCStr = conBin->bestReadCStr,
        minimapCmdCStr[1024], /*Holds command to run minimap2*/
        raconCmdCStr[1024],   /*Holds command to run racon*/
        tmpConCStr[128], /*Hold the consensus name for fq consensus*/
        *tmpCStr = 0,
        *tmpFileCStr = 0,    /*For swapping consensus file names*/
        *tmpFastaFileCStr = "tmp-2023-12-01-con-1563577017123412.fasta",
        *tmpSamFileCStr = "tmp-2023-12-01-map-15635770171234122342.sam";

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-12 Sec-2: Set up the consensus name
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if need to build the consensus name*/
    if(!(*useFqConBl & 1))
    { /*If using the best read for the first round of racon*/

        tmpCStr =
            cStrCpInvsDelm(conBin->consensusCStr, conBin->fqPathCStr);
        tmpCStr -= 6; /*Get to end of .fastq*/
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--cluster-");

        /*Copy cluster number into file name; tmpCStr will point to end*/
        tmpCStr = uCharToCStr(tmpCStr, *clustUChar);
        tmpCStr=cStrCpInvsDelm(tmpCStr, "--consensus.fasta");

        tmpFileCStr = conBin->consensusCStr;
        refFileCStr = conBin->bestReadCStr;
    } /*If using the best read for the first round of racon*/

    else
    { /*I am using a consensus*/
        tmpFileCStr = tmpFastaFileCStr;
        tmpCStr = cStrCpInvsDelm(tmpConCStr, conBin->consensusCStr);

        if(*(tmpCStr - 1) == 'q')
        { /*if need to change fastq to fasta*/
            tmpCStr = conBin->consensusCStr;

            while(*tmpCStr != '\0')
                ++tmpCStr;
            --tmpCStr;

            *tmpCStr = 'a'; /*change fastq to fasta*/
            refFileCStr = tmpConCStr;
        } /*if need to change fastq to fasta*/

        else
            refFileCStr = conBin->consensusCStr;
    } /*I am using a consensus*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-12 Sec-3: Build the consensus using racon
    ^    fun-12 sec-3 sub-1: Build sam file with minimap2
    ^    fun-12 sec-3 sub-1: Build consensus with racon
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-12 Sec-3 Sub-1: Build sam file with minimap2
    \******************************************************************/

    for(
        uint8_t uCharRound = 0;
        uCharRound < *rndsRaconUC;
        ++uCharRound
    ) { /*Loop until have met user requirment for racon*/

        /*Build and run the minimap2 command*/
        tmpCStr = cStrCpInvsDelm(minimapCmdCStr, minimap2CMD);
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
        tmpCStr =
            cpParmAndArg(tmpCStr, refFileCStr, conBin->topReadsCStr);
        cpParmAndArg(tmpCStr, ">", tmpSamFileCStr);

        system(minimapCmdCStr);

        /**************************************************************\
        * Fun-12 Sec-3 Sub-1: Build consensus with racon
        \**************************************************************/

        tmpCStr = cStrCpInvsDelm(raconCmdCStr, raconCMD);
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
        tmpCStr = 
            cpParmAndArg(tmpCStr, conBin->topReadsCStr, tmpSamFileCStr);

        tmpCStr = cpParmAndArg(tmpCStr, refFileCStr, "> ");
        tmpCStr = cStrCpInvsDelm(tmpCStr, tmpFileCStr);

        system(raconCmdCStr); /*Run racon*/

        if(tmpFileCStr == conBin->consensusCStr)
        { /*If named the last consensus after the final name*/
            tmpFileCStr = tmpFastaFileCStr;
            refFileCStr = conBin->consensusCStr;
        } /*If named the last consensus after the final name*/

        else
        { /*Else named the last consensus after the temporary name*/
            tmpFileCStr = conBin->consensusCStr;
            refFileCStr = tmpFastaFileCStr;
        } /*Else named the last consensus after the temporary name*/
    } /*Loop until have met user requirment for racon*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-12 Sec-4: Rename consensus if needed
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    remove(tmpSamFileCStr); /*No longer need this*/

    /*If next round would have output to the consensus file*/
    if(tmpFileCStr != conBin->consensusCStr)
    { /*If need to rename the final consensus*/
        remove(conBin->consensusCStr);
        rename(tmpFileCStr, conBin->consensusCStr);
    } /*If need to rename the final consensus*/

    else                    /*Last consensus saved to consensus file*/
        remove(tmpFastaFileCStr); /*No longer need*/

    if(*useFqConBl & 1)
        remove(tmpConCStr); /*If used a fastq consensus*/
    return;
} /*buildConWithRacon*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 1: if succeded
|    Modifies:
|        - fastq in binClust->fqPathCStr to be the fastq for the cluster
|        - fastq in binTree->fqPathCStr to not have clustered reads
\---------------------------------------------------------------------*/
uint8_t binReadToCon(
    struct minDiff *maxDifference,  /*Has min thresholds to keep reads*/
    const uint8_t *clustUChar,      /*Cluster on*/
    struct readBin *binTree,        /*Bin working on*/
    struct readBin *binClust,       /*Bin to assign reads & consensus*/
    struct samEntry *samStruct,     /*To hold temporary input*/
    struct minAlnStats *minStats,   /*Min stats needed to keep a read*/
    char *minimapCmdCStr,        /*Minimap2 command to run*/
    char *addFileToMiniCmdCStr   /*Ponts to position to add file
                                      names to in minimap2 command*/
) /*Maps reads to consensus and keeps reads that meet user criteria*/
{ /*binReadToCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-13 TOC: binReadToCon
    '    fun-13 sec-1: Variable declerations
    '    fun-13 sec-2: Set defaults & run minimap2
    '    fun-13 sec-3: Make temporary files & open the old stats file
    '    fun-13 sec-4: Check each minimap2 alignment to see if keeping
    '    fun-13 sec-5: Clean up and rename files
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t
        errUChar = 0,
        oneUChar = 1,
        zeroUChar = 0,
        headBool = 0; /*Tells if frist round in stats file*/

    char
        *tmpCStr = 0,
        *tmpStatsCStr = "2023-01-17-1239-stats-tmp-019283274561234.tsv",
        *tmpFqCStr = "2023-01-17-1239-fastq-tmp-019283274561234.fastq";

    struct samEntry
        *zeroSam = 0; /*Just to tell no reference struct*/

    FILE
        *tmpStatsFILE = 0, /*Stats keeping*/
        *clustFILE = 0,/*Holds reads that mapped to the consensuses*/
        *otherBinFILE = 0,/*Holds reads that did not map*/
        *stdinFILE = 0;   /*Holds minimap2 output*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-2: Set defaults & run minimap2
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    binTree->numReadsULng = 0;  /*Reseting size after binning*/
    binClust->numReadsULng = 0; /*For counting number reads in bin*/

    addMinimap2Files(
        addFileToMiniCmdCStr,
        binClust->consensusCStr,
        binTree->fqPathCStr
    ); /*Add the file names to the minimap2 command*/

    stdinFILE = popen(minimapCmdCStr, "r");

    /*Remove the old stats data in the structures*/
    blankSamEntry(samStruct);

    /*Read First line so can check if errored out*/
    errUChar = readSamLine(samStruct, stdinFILE);

    if(*samStruct->samEntryCStr != '@')
    { /*If their is no header*/
        pclose(stdinFILE);
        return 2; /*Minimap2 failed*/        
    } /*If their is no header*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-3: Make temporary files & open the old stats file
    ^     fun-13 sec-3 sub-1: open the temporary files & bin stat file
    ^     fun-13 sec-3 sub-2: Build the clusters fastq name
    ^     fun-13 sec-3 sub-3: read the first line of the stats file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-13 Sec-3 Sub-1: open the temporary files & bin stat file
    \******************************************************************/

    tmpStatsFILE = fopen(tmpStatsCStr, "w"); /*Open the temp file*/

    otherBinFILE = fopen(tmpFqCStr, "w"); /*file for discarded reads*/

    /******************************************************************\
    * Fun-13 Sec-3 Sub-2: Build the clusters fastq name
    \******************************************************************/

    /*Copy bin fastq name to cluster fastq name*/
    strcpy(binClust->fqPathCStr, binTree->fqPathCStr);

    tmpCStr = binClust->fqPathCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    tmpCStr -= 6; /*get to '.' in ".fastq"*/

    /*Add --*/
    *tmpCStr = '-';
    ++tmpCStr;
    *tmpCStr = '-';
    ++tmpCStr;

    strcpy(tmpCStr, "cluster-");
    tmpCStr += 8;

    /*Add the cluster number to the fastq file name*/
    tmpCStr = uCharToCStr(tmpCStr, *clustUChar);
    strcpy(tmpCStr, ".fastq"); /*Add fastq ending to cluster fastq*/

    clustFILE = fopen(binClust->fqPathCStr, "w"); /*make cluster fastq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-4: Check each minimap2 alignment to see if keeping
    ^    fun-13 sec-4 sub-1: Check if is a header
    ^    fun-13 sec-4 sub-2: Score read
    ^    fun-13 sec-4 sub-3: Check if Score meets requirements
    ^    fun-13 sec-4 sub-4: Read next line from minimap2 & stats file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-13 Sec-4 Sub-1: Check if is a header
    \******************************************************************/

    while(errUChar & 1)
    { /*While their is a samfile entry to read in*/

        if(*samStruct->samEntryCStr == '@')
        { /*If was a header*/
            /*Read in next entry*/
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        if(samStruct->flagUSht & 4)
        { /*Make sure the read mapped to something*/
            /*Read in next entry*/
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*Make sure the read mapped to something*/

        /**************************************************************\
        * Fun-13 Sec-4 Sub-2: Score read
        \**************************************************************/

       findQScores(samStruct); /*Find the Q-scores*/

        scoreAln(
            minStats,
            samStruct,
            zeroSam,    /*Do not use reference for dels (no q-score)*/ 
            &oneUChar,   /*Mapped read has Q-score*/
            &zeroUChar   /*Reference has no q-score entry*/
        ); /*Score the alignment*/

        /**************************************************************\
        * Fun-13 Sec-4 Sub-3: Check if Score meets requirements
        \**************************************************************/

        if(
            samStruct->mapqUChar < minStats->minMapqUInt ||
            !(checkIfKeepRead(maxDifference, samStruct) & 1)
        ) { /*If read does not belong in this cluster*/
            /*Print out fastq & stats to the temporary files
              These will be made into the bins fastq files later*/
            samToFq(samStruct, otherBinFILE);
            printSamStats(samStruct, &headBool, tmpStatsFILE);
                /*Need the Q-scores for future clustering steps 
                  Other stats not big deal
                  I would like to save the orginal stats, but the order
                    is different from the fastq file
                */

            ++binTree->numReadsULng; /*Update total scores in bin*/
        } /*If read does not belong in this cluster*/

        else
        { /*else teh read belongs to the cluster*/
            samToFq(samStruct, clustFILE); /*Print read to cluster fq*/
            ++binClust->numReadsULng; /*Adding another read to the bin*/
        } /*else teh read belongs to the cluster*/

        /**************************************************************\
        * Fun-13 Sec-4 Sub-4: Read next line from minimap2 & stats file
        \**************************************************************/

        blankSamEntry(samStruct);

         /*Read in the next line*/
        errUChar = readSamLine(samStruct, stdinFILE);
    } /*While their is a samfile entry to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-13 Sec-5: Clean up and rename files
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    pclose(stdinFILE);
    fclose(clustFILE);
    fclose(otherBinFILE);
    fclose(tmpStatsFILE);

    /*Remove the bins old fastq & stats files*/
    remove(binTree->fqPathCStr);
    remove(binTree->statPathCStr);

    /*Assign the temporary fastq & stats files to the bin*/
    rename(tmpFqCStr, binTree->fqPathCStr);
    rename(tmpStatsCStr, binTree->statPathCStr);

    return 1; /*No errors*/
} /*binReadToCon*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        - readBin Struct with the closest consensus
|    - Modifies:
|        - closesConUint to hold the score of the most similar consensus
| Note:
|    - clusters with ->balUChar < 0 will be ignored
\---------------------------------------------------------------------*/
struct readBin * cmpCons(
    struct minDiff *maxDifference,  /*Has min thresholds to keep reads*/
    struct readBin *conBin,       /*Bin with consensus to compare*/
    struct readBin *conBinTree,   /*Tree of consensus to compare to*/
    struct samEntry *samStruct,   /*Struct to hold input from minimap2*/
    struct samEntry *refStruct,   /*Struct to hold input from minimap2*/
    struct minAlnStats *minStats, /*Min stats needed to keep a error*/
    char *minimapCmdCStr,      /*Minimap2 command to run*/
    char *addFileToMiniCmdCStr /*Pionts to position to add file names
                                    in minimap2 command*/
) /*Compares a consenses to a another consensus. This will do a
    recursive call if conBinTree has children*/
{ /*cmpCons*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-14 TOC: cmpCons
    '     fun-14 sec-1: Variable declerations
    '     fun-14 sec-2: Check if have valid user input
    '     fun-14 sec-3: Run minimap2 and check first line
    '     fun-14 sec-5: Score the consensus to the reference
    '     fun-14 sec-6: Compare this score to the other consensuses
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-14 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t
        errUChar = 0, /*For holding error returns from functions*/
        zeroUChar = 0;

    char
        *tmpCStr = 0;

    uint32_t
        incBuffUInt = 10000; /*Amount to increase buff size each time*/

    FILE
        *stdinFILE = 0; /*File to see if input files are valid*/

    struct readBin
        *refBin = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-14 Sec-2: Check if have valid user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(conBin == 0)
        return 0;

    if(conBin->leftChild != 0)
        return 0; /*This is a bin, not a cluster*/

    if(conBin->balUChar < 0)
        return 0;

    if(conBinTree == 0)
        return 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-14 Sec-4: Read in target consensus sequence
    ^    fun-14 sec-4 sub-1: Read past the header
    ^    fun-14 sec-4 sub-2: Read in the sequence
    ^    fun-14 sec-4 sub-3: Count the sequence length
    ^ Note:
    ^    - I should change this to the consensus I am checking & then
    ^      have a marker that tells this function to use the input
    ^      reference
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-14 Sec-4 Sub-1: Read past the header
    \******************************************************************/

    blankSamEntry(refStruct);

    /*Set up null endings for lines*/
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Open the reference file for reading*/
    stdinFILE = fopen(conBin->consensusCStr, "r");

    /*Read in the header, I know it worked, because of minimap2*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

    while(
     !(*(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) !='\0' ||
       *(refStruct->samEntryCStr + refStruct->lenBuffULng -2)!='\n')
    ) { /*While have a header to read in*/
        /*Read in the next part of the header*/
        fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

        /*Resetup markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';
    } /*While have a header to read in*/

    /******************************************************************\
    * Fun-14 Sec-4 Sub-2: Read in the sequnence
    \******************************************************************/

    tmpCStr = refStruct->samEntryCStr;
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Read in the first part of the sequence*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, stdinFILE);

    while(*(refStruct->samEntryCStr+refStruct->lenBuffULng-2) > 16)
    { /*While on the sequence line*/
        refStruct->samEntryCStr =
            realloc(
                refStruct->samEntryCStr,
                refStruct->lenBuffULng + incBuffUInt
            );

        if(refStruct->samEntryCStr == 0)
        { /*memory allocation error*/
            fclose(stdinFILE);
            fclose(stdinFILE);
            return 0;
        } /*Memory allocation error*/

        /*Set pointer to new buffer*/
        tmpCStr = refStruct->samEntryCStr + refStruct->lenBuffULng;
        refStruct->lenBuffULng += incBuffUInt; /*Update buff size*/

        /*Reset new line markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng-2) ='\0';

        /*Read in next part of reference sequence*/
        fgets(tmpCStr, refStruct->lenBuffULng, stdinFILE);
    } /*While on the sequence line*/

    fclose(stdinFILE);

    /******************************************************************\
    * Fun-14 Sec-4 Sub-3: Count the sequence length
    \******************************************************************/

    refStruct->seqCStr = refStruct->samEntryCStr;

    /*scoreAln assumes the Q-score entry is not null, so I am setting
      it to something to avoid errors. This is ok since I am telling
      scoreAln that their is not Q-score for the reference
    */
    refStruct->qCStr = refStruct->seqCStr;
    tmpCStr = refStruct->seqCStr;
    refStruct->readLenUInt = 0;

    while(*tmpCStr != '\0')
    { /*While have bases to count*/
        ++refStruct->readLenUInt;
        ++tmpCStr;
    } /*While have bases to count*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ^ Fun-14 Sec-5: Score the consensus to the reference
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(conBinTree != 0)
    { /*While have consensuses to compare*/
        blankSamEntry(samStruct);
        refBin = conBinTree->rightChild;

        while(refBin != 0)
        { /*While have another clusters consensus to compare*/
            if(refBin->balUChar < 0)
            { /*If this cluster has been marked to be skipped*/
                refBin = refBin->rightChild;
                continue;
            } /*If this cluster has been marked to be skipped*/

            blankSamEntry(samStruct);

            addMinimap2Files(
                addFileToMiniCmdCStr,
                conBin->consensusCStr,
                refBin->consensusCStr
            ); /*Add the file names to the minimap2 command*/

            stdinFILE = popen(minimapCmdCStr, "r"); /*run minimap2*/
            errUChar = readSamLine(samStruct, stdinFILE); /*1st line*/

            if(*samStruct->samEntryCStr != '@')
            { /*If no header*/
                fclose(stdinFILE);
                return 0;
            } /*If no header*/

            while(errUChar & 1)
            { /*While on the haeder lines*/
                errUChar = readSamLine(samStruct, stdinFILE);

                if(*samStruct->samEntryCStr != '@')
                    break; /*If not a header*/
            } /*While on the haeder lines*/

            scoreAln(
                minStats,
                samStruct,
                refStruct,   /*Reference struct to score deletions with*/ 
                &zeroUChar,  /*Mapped consensus has no Q-score*/
                &zeroUChar   /*Mapped consensus has no Q-score*/
            ); /*Score the alignment*/

            pclose(stdinFILE); /*Done with this file*/

            if(checkIfKeepRead(maxDifference, samStruct) & 1)
                return refBin; /*If consensus look the same*/

            refBin = refBin->rightChild;
        } /*While have another clusters consensus to compare*/

        conBinTree = conBinTree->leftChild;
    } /*While have consensuses to compare*/

    return 0;
} /*cmpCons*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - binToKeep->fqPathCStr to hold binToMerge reads
|    Deletes:
|        - File named after binToMerge->fqPathCStr
\---------------------------------------------------------------------*/
void mergeBins(
    struct readBin *binToKeep,  /*Bin to merge into*/
    struct readBin *binToMerge /*Bin to merge into binToKeep*/
) /*Merge two bins togetehr*/
{ /*mergeBins*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/
    ' Fun-15 TOC: mergeBins
    '    fun-15 sec-1: Variable declerations
    '    fun-15 sec-2: Open and check files
    '    fun-15 sec-3: Copy mergeFILE contents to end of keepFILE
    '    fun-15 sec-4: Close files and clean up binToMerge
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint32_t
        lenBuffUInt = 1 << 2^15; /*2^16, which is about 64kb*/

    char
        buffCStr[lenBuffUInt];

    FILE
        *readsKeepFILE = 0,
        *readsMergeFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-2: Open and check files
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readsKeepFILE = fopen(binToKeep->fqPathCStr, "a");

    if(readsKeepFILE == 0)
        return;

    readsMergeFILE = fopen(binToMerge->fqPathCStr, "r");

    if(readsMergeFILE == 0)
    { /*If the merge file does not exist*/
        fclose(readsKeepFILE);
        return;
    } /*If the merge file does not exist*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-3: Copy mergeFILE contents to end of keepFILE
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Copy reads from one bin to the other*/
    while(fread(buffCStr, sizeof(uint8_t), lenBuffUInt, readsMergeFILE))
        fwrite(buffCStr, sizeof(uint8_t), lenBuffUInt, readsKeepFILE);
            
    binToKeep->numReadsULng += binToMerge->numReadsULng;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-15 Sec-4: Close files and clean up binToMerge
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fclose(readsKeepFILE);
    fclose(readsMergeFILE);

    binDeleteFiles(binToMerge); /*Remove all files in the bin*/

    return;
} /*mergeBins*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 0: Read does not meet minimum thresholds in maxDifference
|        - 1: Keep the read (meets minimum thresholds)
\---------------------------------------------------------------------*/
uint8_t checkIfKeepRead(
    struct minDiff *maxDifference, /*Has thresholds to keep reads*/
    struct samEntry *samStruct     /*Has stats for read*/
) /*Checks if read does not meet one threshold in maxDifference*/
{ /*checkIfKeepRead*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-16 TOC: Sec-1 Sub-1: checkIfKeepRead
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    float
        percSnpFlt = 0,
        percInsFlt = 0,
        percDelFlt = 0;

    /*This might be better if I found these separately and had if 
          statements between each calcuation, however this is clearer
      Need * 1.0 to ensure it does not truncate*/
    percSnpFlt =
        samStruct->numKeptSNPUInt / (samStruct->readLenUInt * 1.0);

    percInsFlt =
        samStruct->numKeptInsUInt / (samStruct->readLenUInt * 1.0);

    percDelFlt =
        samStruct->numKeptDelUInt / (samStruct->readLenUInt * 1.0);

    if(maxDifference->minSNPsFlt < percSnpFlt) /*percent SNPs to high*/
        return 0;
    /*The sequences in the sam alignment are to different*/
    if(maxDifference->minDiffFlt < percSnpFlt + percDelFlt + percInsFlt)
        return 0;
    /*To many indels is to hight*/
    if(maxDifference->minIndelsFlt < percInsFlt + percDelFlt)
        return 0;
    if(maxDifference->minInssFlt < percInsFlt) /*To many insertions*/
        return 0;
    if(maxDifference->minDelsFlt < percDelFlt) /*To many deletions*/
        return 0;

    return 1;
} /*checkIfKeepRead*/

/*---------------------------------------------------------------------\
| Output:
|    - Modifies:
|        - newBin to have stats in samStruct
\---------------------------------------------------------------------*/
void samEntryToReadStat(
    struct readStat *newBin,   /*Read bin to hold stats from samStruct*/
    struct samEntry *samStruct /*copy stats from this struct*/
) /*Copies stats from a samEntry struct to a readStat struct*/
{ /*samEntryToReadBin*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-17 TOC: samEntryToReadStat
    '    fun-17 sec-2: Copy stats
    '    fun-17 sec-1: Variable declerations
    '    fun-17 sec-3: Copy the query id
    '    fun-17 sec-4: Copy the reference id
    '
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-17 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpCStr = 0,
        *tmpCpCStr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-17 Sec-2: Copy stats
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    newBin->mapqUChar = samStruct->mapqUChar;

    newBin->medianQFlt = samStruct->medianQFlt;
    newBin->medianAligQFlt = samStruct->medianAligQFlt;
    newBin->meanQFlt = samStruct->meanQFlt;
    newBin->meanAligQFlt = samStruct->meanAligQFlt;

    
    newBin->readLenUInt = samStruct->readLenUInt;
    newBin->readAligLenUInt = samStruct->readAligLenUInt;

    newBin->numMatchUInt = samStruct->numMatchUInt;
    newBin->numSNPUInt = samStruct->numSNPUInt;
    newBin->numKeptSNPUInt = samStruct->numKeptSNPUInt;
    newBin->numDelUInt = samStruct->numDelUInt;
    newBin->numKeptDelUInt = samStruct->numKeptSNPUInt;
    newBin->numInsUInt = samStruct->numInsUInt;
    newBin->numKeptInsUInt = samStruct->numKeptSNPUInt;

    tmpCpCStr = newBin->queryIdCStr;
    tmpCStr = samStruct->queryCStr;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-17 Sec-3: Copy the query id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpCStr > 16)
    { /*While have the query id to copy over*/
        *tmpCpCStr = *tmpCStr;
        ++tmpCStr;
        ++tmpCpCStr;
    } /*While have the query id to copy over*/

    tmpCpCStr = '\0'; /*Mark end of c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-17 Sec-4: Copy the reference id
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCpCStr = newBin->refIdCStr;
    tmpCStr = samStruct->refCStr;

    while(*tmpCStr > 16)
    { /*While have the reference id to copy over*/
        *tmpCpCStr = *tmpCStr;
        ++tmpCStr;
        ++tmpCpCStr;
    } /*While have the reference id to copy over*/

    tmpCpCStr = '\0'; /*Mark end of c-string*/
    return;
} /*samStructToReadBin*/

/*---------------------------------------------------------------------\
| Output:
|   Returns:
|     - 1: if both alignments are for the same query
|     - 0: if alignments are for different querys
\---------------------------------------------------------------------*/
uint8_t isSamAlnDup(
    struct samEntry *samAln,    /*sam alignemnt to check if duplicate*/
    struct samEntry *lastSamAln /*last sam alignment to compare to*/
) /*Checks if both sam entrys (alignments) have the same querys*/
{ /*isSamAlnDup*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-18 TOC: Sec-1 Sub-1: isSamAlnDup
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char
        *oldQueryCStr = lastSamAln->queryCStr,
        *newQueryCStr = samAln->queryCStr;

    while(*oldQueryCStr == *newQueryCStr)
    { /*While both querys are the same*/
        ++oldQueryCStr;
        ++newQueryCStr;

        if(*oldQueryCStr < 33)
            return 1; /*Alignments are the same*/
    } /*While both querys are the same*/

    return 0;   
} /*isSamAlnDup*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        binTree->consensusCStr to hold the polished consensus
|    Returns:
|        - 1 If built a consensus
|        - 2 If input file does not exists
|        - 4 IF consensus not built
\---------------------------------------------------------------------*/
unsigned char medakaPolish(
    const char *condaBl,    /*Tells if using miniconda medaka*/
    char *medakaModelCStr,  /*Model to use with medaka*/
    char *threadsCStr,      /*Number threads to use*/
    struct readBin *conBin  /*bin with consensus & top reads*/
) /*Polish a consensus with medaka using the best reads*/
{ /*medakaPolish*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-19 TOC: Sec-1 Sub-1: medakaPolish
    '    fun-19 sec-1: Variable declerations
    '    fun-19 sec-2: See if the best reads and consensus files exist
    '    fun-19 sec-3: Build the command to run medaka
    '    fun-19 sec-4: Run medaka & clean up extra files
    '    fun-19 sec-5: rename consensus & delete directory medaka made
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-19 Sec-1: Variable declerations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char
        medakaConPathCStr[256],   /*Path to consesnsus medaka made*/
        medDirCStr[256],           /*For temporary c-string building*/
        medakaCmdCStr[1024],      /*Holds command to run medaka*/
        *tmpCStr = medakaCmdCStr; /*Temporary, for manipulating cStrs*/

    unsigned long
        fileLenULng = 0;         /*See if files have something in them*/

    FILE 
        *testFILE = 0; /*Test if files exist*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-19 Sec-2: See if the best reads and consensus files exist
    ^    fun-19 sec-2 sub-1: See if the conssus file exists
    ^    fun-19 sec-2 sub-2: See if the best reads file exists
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*******************************************************************
    * Fun-19 Sec-2 Sub-1: See if the conssus file exists
    *******************************************************************/

    /*Open file*/
    testFILE = fopen(conBin->consensusCStr, "r");

    if(testFILE == 0)
        return 2;    /*No file was made*/

    /*Finde the length of the file*/
    fseek(testFILE, 0L, SEEK_END);   /*Find end of file*/
    fileLenULng = ftell(testFILE);   /*Find offset (length) of end*/
    fclose(testFILE);                /*No longer need the file open*/

    if(fileLenULng < 10)
        return 2;    /*Nothing or little in the file*/

    /*******************************************************************
    * Fun-19 Sec-2 Sub-2: See if the best reads file exists
    *******************************************************************/

    testFILE = fopen(conBin->bestReadCStr, "r");

    if(testFILE == 0)
        return 2;    /*No file was made*/

    /*Finde the length of the file*/
    fseek(testFILE, 0L, SEEK_END);   /*Find end of file*/
    fileLenULng = ftell(testFILE);   /*Find offset (length) of end*/
    fclose(testFILE);                /*No longer need the file open*/

    if(fileLenULng < 10)
        return 2;    /*Nothing or little in the file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-19 Sec-3: Build the command to run medaka
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Make temporary director for medaka*/
    tmpCStr = cStrCpInvsDelm(medDirCStr, conBin->fqPathCStr);
    tmpCStr -= 6; /*move to "." in ".fastq"*/
    tmpCStr = cStrCpInvsDelm(tmpCStr, "--medaka");

    /*Copy the path to the consensus medaka will build*/
    tmpCStr = cStrCpInvsDelm(medakaConPathCStr, medDirCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, "/consensus.fasta");

    /*Copy the enviroment activate command*/
    if(*condaBl & 1)
        tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medCondCMD);
    else
        tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medakaCMD);

    *tmpCStr = ' ';
    ++tmpCStr;

    tmpCStr = cStrCpInvsDelm(tmpCStr, "medaka_consensus");
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-i", conBin->topReadsCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-d", conBin->consensusCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-m", medakaModelCStr);
    tmpCStr = cpParmAndArg(tmpCStr, "-o", medDirCStr);

    /*Copy the deactivate command in*/
    if(*condaBl & 1)
        tmpCStr = cStrCpInvsDelm(tmpCStr, medCondCMDEnd);
    else
        tmpCStr = cStrCpInvsDelm(tmpCStr, medakaCMDEnd);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-19 Sec-4: Run medaka & clean up extra files
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    system(medakaCmdCStr); /*Run medaka command*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr, "/calls_to_draft.bam");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr, "/calls_to_draft.bam.bai");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr,"/consensus.fasta.gaps_in_draft_coords.bed");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, medDirCStr);
    cStrCpInvsDelm(tmpCStr,"/consensus_probs.hdf");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*Remove the temporary mapping files made by medaka*/
    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, conBin->consensusCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".fai");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*Remove the temporary mapping files made by medaka*/
    tmpCStr = cStrCpInvsDelm(medakaCmdCStr, conBin->consensusCStr);
    tmpCStr = cStrCpInvsDelm(tmpCStr, ".map-ont.mmi");
    testFILE = fopen(medakaCmdCStr, "r");

    if(testFILE != 0)
    { /*If need to delete the file*/
        fclose(testFILE);
        remove(medakaCmdCStr);
    } /*If need to delete the file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-19 Sec-5: rename consensus & delete directory medaka made
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    testFILE = fopen(medakaConPathCStr, "r");

    if(testFILE == 0)
        return 4;    /*Medaka did not build a consensus*/

    /*See if the test file has anything in it*/
    fseek(testFILE, 0L, SEEK_END);   /*Find end of file*/
    fileLenULng = ftell(testFILE);   /*Find offset (length) of end*/
    fclose(testFILE);                /*No longer need the file open*/

    if(fileLenULng < 10)
        return 4;    /*Nothing or little in the file*/

    remove(conBin->consensusCStr); /*remove the old consensus*/
    rename(medakaConPathCStr, conBin->consensusCStr);
    remove(medDirCStr); /*Delete directory made by medaka*/

    return 1;
} /*medakaPolish*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold the copied C-string
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cStrCpInvsDelm(
    char *cpToCStr,  /*C-string to copy values to*/
    char *cpFromCStr /*C-string to copy*/
) /*Copy one c-string till an tab, newline, or '\0' (keeps spaces)*/
{ /*cStrCpInvsDelm*/


    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-20 TOC: Sec-1 Sub-1: cStrCpInvsDelim
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    while(*cpFromCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpFromCStr;
        ++cpToCStr;
        ++cpFromCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cStrCpInvsDelim*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold space, parameter, space, & argument
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cpParmAndArg(
    char *cpToCStr,   /*Holds copied parameter and argement*/
    char *cpParmCStr, /*Paramater to copy*/
    char *cpArgCStr   /*Argument to copy*/
) /*Copies adds a space and copies a paramater and a agrugment*/
{ /*cpParmAndArg*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-21 TOC: Sec-1 Sub-1: cpSpaceCStr
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    *cpToCStr = ' ';
    ++cpToCStr;
    
    while(*cpParmCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpParmCStr;
        ++cpToCStr;
        ++cpParmCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = ' ';
    ++cpToCStr;

    while(*cpArgCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpArgCStr;
        ++cpToCStr;
        ++cpArgCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cpParmAndArg*/

/*---------------------------------------------------------------------\
| Output:
|    - Returns:
|        o 1 If succeded
|        o 2 if could not open or read in the reference sequence
|        o 4 if could not open the bestReads file
|        o 16 if had to few sequences map to the selected read
|        o 32 if minimap2 crashed
|        o 64 for memory allocation errors
\---------------------------------------------------------------------*/
unsigned char simpleMajCon(
    const unsigned char *clustUChar, /*Number of cluster working on*/
    struct readBin *binStruct, /*Has best read & top reads files*/
    struct samEntry *samStruct,/*For reading in sam file entries*/
    unsigned char minBaseQUC,  /*Has min mapq to keep an base*/
    unsigned char minInsQUC,   /*Has min mapq to keep an insertion*/
    char *threadsCStr,         /*Number threads to use with minimap2*/
    const float *minBasePerFlt, /*Min percentage of bases needed to
                                       not be a deletion (counts all
                                       bases that pass the quality
                                       threshold)*/
    const float *minInsPercFlt  /*Min percentage of supporting reads
                                  needed to keep an insertion*/
) /*Builds a majority consensus from the best reads & top read*/
{ /*simpleMajCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-22 TOC: simpleMajCon
    '    fun-22 sec-1: Variable declerations
    '    fun-22 sec-2: Check if the bestRead and topReads files exist
    '        - Also reads in the reference sequence
    '    fun-22 sec-3: Prepare the minimap2 command
    '    fun-22 sec-4: Make consensus array & initalize with reference
    '    fun-22 sec-5: Map reads to the reference read
    '    fun-22 sec-7: Print out cosensus & do clean up
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char qEntryBl = 0;           /*Marks if reference has Q-core entry*/
    char minimap2CmdCStr[1024];  /*Holds minimap2 command to run*/
    char *tmpCStr = 0;           /*Temp ptr for c-string manipulations*/
    char *cigCStr = 0;           /*Reading the cigar entry*/
    char *seqCStr = 0;           /*Manipulating/reading the sequence*/
    char *qCStr = 0;             /*Manipulating/reading q-score entry*/

    unsigned char qScoreUChar = 0; /*Holds the Q-score for a base*/
    uint8_t errUChar = 0;        /*Holds error messages from functions*/

    uint32_t cigEntryUInt = 0;  /*Holds number of bases in cigar entry*/

    uint64_t minNumBasesUL = 0; /*Min number of bases to keep a base*/
    uint64_t minInsUL = 0;      /*Min number of insertions to keep ins*/
    uint64_t numBasesUL = 0;    /*Number of bases at a position*/
    uint64_t numSeqUL = 0;      /*Number of mapped sequences*/
    uint64_t numMisSeqUL = 0;   /*Number of mapped sequences*/

    struct baseStruct *headBase = 0; /*Head of the list of bases*/
    struct baseStruct *incBase = 0;  /*First base at a position*/
    struct baseStruct *lastBase = 0; /*Base before tmpBase*/
    struct baseStruct *tmpBase = 0;  /*Base position working on*/

    FILE *stdinFILE = 0;        /*For reading and writing files*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-2: Check if the bestRead and topReads files exist
    ^    - Also reads in the reference sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    stdinFILE = fopen(binStruct->bestReadCStr, "r");

    if(stdinFILE == 0)
        return 2;

    blankSamEntry(samStruct); /*Make sure start with blank*/

    tmpCStr = binStruct->bestReadCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    if(*(tmpCStr - 1) == 'q')
    { /*If is a fastq file*/
        errUChar = readRefFqSeq(stdinFILE, samStruct); /*fastq file*/
        qEntryBl = 1;
    } /*If is a fastq file*/

    else if(*(tmpCStr - 1) == 'a')
    { /*Else the best read is a fasta file*/
        errUChar = readInConFa(binStruct->bestReadCStr, samStruct);
        qEntryBl = 0;
    } /*Else the best read is a fasta file*/

    else
    { /*Else reference is the consensus*/
        errUChar = readInConFa(binStruct->consensusCStr, samStruct);
        qEntryBl = 0;
    } /*Else reference is the consensus*/

    if(errUChar & 64)
        return 64; /*Memory allocation error*/
    if(!(errUChar & 1))
        return 2; /*issue with consensus*/

    /*Make sure have doube the sequence + q-score size for insertions*/
    if(samStruct->readLenUInt * 2 > samStruct->lenBuffULng)
    { /*If need to increase the size of the buffer*/
        samStruct->samEntryCStr =
            realloc(
                samStruct->samEntryCStr,
                sizeof(char) * samStruct->readLenUInt * 2
        ); /*Resize the structs buffer*/

        if(samStruct->samEntryCStr == 0)
            return 64;
    } /*If need to increase the size of the buffer*/

    /*Make the consensus array*/
    fclose(stdinFILE);

    if(errUChar & 64)
        return 64;
    if(!(errUChar & 1))
        return 4; /*Invalide reference*/

    /*Check if can open the top reads file*/
    stdinFILE = fopen(binStruct->topReadsCStr, "r");

    if(stdinFILE == 0)
        return 4;
     /*Close file in Fun-22 Sec-4*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-3: Prepare the minimap2 command
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr =
        cStrCpInvsDelm(binStruct->consensusCStr, binStruct->fqPathCStr);
    tmpCStr -= 6; /*Get to end of .fastq*/
    tmpCStr = cStrCpInvsDelm(tmpCStr, "--cluster-");
    tmpCStr = uCharToCStr(tmpCStr, *clustUChar); /*Add cluster number*/
    tmpCStr=cStrCpInvsDelm(tmpCStr, "--consensus.fasta");

    /*Prepare the minimap2 command*/
    tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    cpParmAndArg(
        tmpCStr,
        binStruct->bestReadCStr,
        binStruct->topReadsCStr
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-4: Make the consensus array & initalize with reference
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    seqCStr = samStruct->seqCStr;
    qCStr = samStruct->qCStr;
    tmpBase = 0;

    /*Initalize base struct using the reference read*/

    while(*seqCStr > 16)
    { /*While have bases to read in*/
        if(*seqCStr < 33)
            continue;

        tmpBase = malloc(sizeof(struct baseStruct));

        if(tmpBase == 0)
        { /*If had a memory allocation error*/
            lastBase = 0;       /*So function knows no previous bases*/

            /*Free all the made bases*/
            while(headBase != 0)
                headBase = freeBaseStruct(&headBase, lastBase);

            return 64;
        } /*If had a memory allocation error*/

        if(headBase == 0)
            headBase = tmpBase;           /*If is the first base*/
        else
            lastBase->nextBase = tmpBase; /*Set up list*/

        tmpBase->nextBase = 0; /*Only matches at this point*/
        tmpBase->altBase = 0;  /*No inserts at this point*/

        qScoreUChar = *qCStr - Q_ADJUST; /*Get the Q-score*/

        /*qEntryBl ensures if only files if their was a Q-score entry*/
        if(qEntryBl && qScoreUChar < minBaseQUC)
        { /*If base is to low of quality to keep, make a blank struct*/
            tmpBase->baseChar = 0;     /*Base is beinging discared*/
            tmpBase->errTypeChar = 0;  /*match=0, snp=1, ins=2*/
            tmpBase->numBasesUL = 0;
        } /*If base is to low of quality to keep, make a blank struct*/

        else
        { /*Else keeping the base*/
            tmpBase->baseChar = *seqCStr;
            tmpBase->errTypeChar = 0;       /*match = 0, snp=1, ins=2*/
            tmpBase->numBasesUL = 1;        /*match = 0, snp=1, ins=2*/
        } /*Else keeping the base*/

        lastBase = tmpBase;
        ++seqCStr;
        ++qCStr;
    } /*While have bases to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-5: Map reads to the reference read
    ^    fun-22 sec-5 sub-1: run miniamp2 and read in first line
    ^    fun-22 sec-5 sub-2: Read in sequence and position on first base
    ^    fun-22 sec-5 sub-3: Add matches & SNPs to consensus array
    ^    Fun-22 Sec-5 sub-4: Get base matching cigar (match/snp or ins)
    ^    fun-22 sec-5 sub-5: Check if should keep base or discard
    ^    fun-22 sec-5 sub-6: Find matching base at position
    ^    fun-22 sec-5 sub-7: Make new base (no matching) or update count
    ^    fun-22 sec-5 sub-8: Make sure on an insertion
    ^    fun-22 sec-5 sub-9: Check if should keep ins
    ^    fun-22 sec-5 sub-10: Find matching ins at pos
    ^    fun-22 sec-5 sub-11: check if need to make ins
    ^    fun-22 sec-5 sub-12: Ingnore deletions
    ^    fun-22 sec-5 sub-13: Ingnore soft masking
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-22 Sec-5 Sub-1: run miniamp2 and read in first line
    \******************************************************************/

    stdinFILE = popen(minimap2CmdCStr, "r");
    blankSamEntry(samStruct); /*Make sure start with blank*/

    /*Read in a single sam file line to check if valid (header)*/
    errUChar = readSamLine(samStruct, stdinFILE);

    if(!(errUChar & 1))
    { /*If an error occured*/
        pclose(stdinFILE);
        return 32;
    } /*If an error occured*/

    /******************************************************************\
    * Fun-22 Sec-5 Sub-2: Read in sequence and position on first base
    \******************************************************************/

    while(errUChar & 1)
    { /*While have alignments to read in from the sam file*/
        if(*samStruct->samEntryCStr == '@')
        { /*If on a header entry, read in next entry*/
            errUChar = readSamLine(samStruct, stdinFILE);
            continue;
        } /*If on a header entry, read in next entry*/

        cigCStr = samStruct->cigarCStr;
        seqCStr = samStruct->seqCStr;
        qCStr = samStruct->qCStr;

        if(*seqCStr == '*' || (*qCStr == '*' && *(qCStr + 1) == '\t'))
        { /*If no entry to check*/
            errUChar = readSamLine(samStruct, stdinFILE);
            continue;
        } /*If no entry to check*/

        if(samStruct->flagUSht & 4)
        { /*If was an unampped read*/
            errUChar = readSamLine(samStruct, stdinFILE);
            ++numMisSeqUL;
            ++numSeqUL;
            continue;
        } /*If was an unampped read*/

        /*Get first base in the list*/
        incBase = headBase;
        cigEntryUInt = samStruct->posOnRefUInt - 1;
            /*-1 for 1 index for posOnRef, but 0 index fo incBase*/

        while(cigEntryUInt > 0)
        { /*While not on the starting base*/
            incBase = incBase->nextBase;
            --cigEntryUInt;
        } /*While not on the starting base*/

        while(*cigCStr != '\t')
        { /*While not at the end of the sam alignment sequence*/
            /*Get the cigar entry*/
            readCigEntry(&cigCStr, &cigEntryUInt);

        /**************************************************************\
        * Fun-22 Sec-5 Sub-3: Add matches & SNPs to consensus array
        \**************************************************************/

            switch(*cigCStr)
            { /*switch: check the error type & add bases to consensus*/
                case 'X':              /*snp, similar loop to match*/
                case '=':              /*Match, similar loop to snp*/
                case '\t':             /*Match at end of cigar*/
                /*Switch: Match, end of cigar match, snp, or insertion*/
                    while(cigEntryUInt > 0)
                    { /*While have bases to add to the bases list*/
                        --cigEntryUInt;

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-4: Find first non-insertion
                        \**********************************************/

                        while(incBase->errTypeChar & 2)
                        { /*Get off the insertion*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;
                        } /*Get off the insertion*/

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-5: Check if should keep base
                        \**********************************************/

                        tmpBase = incBase; /*Base position on*/
                        qScoreUChar = *qCStr - Q_ADJUST; /*get Q-score*/

                        if(qScoreUChar < minBaseQUC)
                        { /*If the base is low quality, skip*/
                            /*Move to the next base in the list*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;

                            /*Move to the next base in the sam entry*/
                            ++seqCStr;
                            ++qCStr;
                            continue;
                        } /*If the base is low quality, skip*/

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-6: Find alt matching base
                        \**********************************************/

                        while(tmpBase->baseChar != *seqCStr)
                        { /*While the bases are not equal*/
                            if(tmpBase->baseChar == 0)
                                break; /*Empty base*/

                            lastBase = tmpBase;
                            tmpBase = tmpBase->altBase;

                            if(tmpBase == 0)
                                break; /*Need to create a new base*/
                        } /*While the bases are not equal*/

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-7: Make or update alt base
                        \**********************************************/

                        if(tmpBase == 0)
                        { /*If need to create a new base*/
                            lastBase->altBase =
                                malloc(sizeof(struct baseStruct));

                            tmpBase = lastBase->altBase;
                            tmpBase->nextBase = 0; /*Is an alterantive*/
                            tmpBase->altBase = 0;  /*No alternatives*/

                            tmpBase->numBasesUL = 0;
                            tmpBase->baseChar = *seqCStr;
                            tmpBase->errTypeChar = 0;
                        } /*If need to create a new base*/

                        else if(tmpBase->baseChar == 0)
                            tmpBase->baseChar = *seqCStr;

                        ++tmpBase->numBasesUL;

                        /*Move to the next base*/
                        lastBase = incBase;
                        incBase = incBase->nextBase;
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to add to the bases list*/
            
                    break;
                /*Switch: Match, end of cigar match, snp, or insertion*/


                case 'I':   /*Is an insertion*/
                /*Switch: for insertions*/

                    while(cigEntryUInt > 0)
                    { /*While have bases to add to the bases list*/
                        --cigEntryUInt;

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-8: Make sure on an insertion
                        \**********************************************/

                        /*Insertions at the ends will likely be treated
                          as softmasks by minimap2. This means that
                          incBase will != 0 for insertions
                        */
                        if(!(incBase->errTypeChar & 2))
                        { /*If next base is not an insertion*/
                            tmpBase=malloc(sizeof(struct baseStruct));

                            if(lastBase->errTypeChar & 2)
                            { /*If the last base was an ins*/ 
                                tmpBase->nextBase = incBase->nextBase;
                                incBase->nextBase = tmpBase;
                                incBase = tmpBase;
                            } /*If the last base was an ins*/ 

                            else
                            { /*else the last base was not an ins*/ 
                                tmpBase->nextBase = incBase;
                                lastBase->nextBase = tmpBase;
                                incBase = tmpBase;
                            } /*else the last base was not an ins*/ 

                            /*Set defualts for the new base*/
                            incBase->altBase = 0;
                            incBase->errTypeChar = 2;
                            incBase->baseChar = 0;
                            incBase->numBasesUL = 0;
                        } /*If next base is not an insertion*/

                        /**********************************************\
                        * Fun-22 Sec-5 Sub-9: Check if should keep ins
                        \**********************************************/
    
                        tmpBase = incBase; /*Base position on*/
                        qScoreUChar = *qCStr - Q_ADJUST; /*get Q-score*/
    
                        if(qScoreUChar < minInsQUC)
                        { /*If the base is low quality, skip*/
                            /*Move to the next base in the list*/
    
                            /*This avoids scattered insertions*/
                            if(
                                cigEntryUInt > 0 &&  /*More insertions*/
                                !(incBase->nextBase->errTypeChar & 2)
                                    /*^^Next base is not an insertion*/
                            ) { /*If have other insertions to process*/
                                /*Make a new base for next insertion*/
                                tmpBase =
                                    malloc(sizeof(struct baseStruct));
    
                                tmpBase->nextBase = incBase->nextBase;
                                incBase->nextBase = tmpBase;
                                lastBase = incBase;
                                incBase = tmpBase;
    
                                /*Set the defaulst for the new base*/
                                incBase->altBase = 0;
                                incBase->errTypeChar = 2;
                                incBase->baseChar = 0;
                                incBase->numBasesUL = 0;
                            } /*If have other insertions to process*/
    
                            else
                            { /*Else is safe to move to the next base*/
                                lastBase = incBase;
                                incBase = incBase->nextBase;
                            } /*Else is safe to move to the next base*/
    
                            /*Move to the next base in the sam entry*/
                            ++seqCStr;
                            ++qCStr;
                            continue;
                        } /*If the base is low quality, skip*/
    
                        /**********************************************\
                        * Fun-22 Sec-5 Sub-10: Find matching ins at pos
                        \**********************************************/
    
                        while(tmpBase->baseChar != *seqCStr)
                        { /*While the bases are not equal*/
                            if(tmpBase->baseChar == 0)
                                break; /*Empty base*/
    
                            lastBase = tmpBase;
                            tmpBase = tmpBase->altBase;
    
                            if(tmpBase == 0)
                                 break; /*Need to create a new base*/
                        } /*While the bases are not equal*/
    
                        /**********************************************\
                        * Fun-22 Sec-5 Sub-11: check if need to make ins
                        \**********************************************/
    
                        if(tmpBase == 0)
                        { /*If need to create a new base*/
                            lastBase->altBase =
                                    malloc(sizeof(struct baseStruct));
    
                            tmpBase = lastBase->altBase;
                            tmpBase->nextBase = 0; /*Is an alterantive*/
                            tmpBase->altBase = 0;  /*No alternatives*/
    
                            tmpBase->errTypeChar = 2; /*Insertion*/
                            tmpBase->baseChar = *seqCStr;
                            tmpBase->numBasesUL = 0;
                        } /*If need to create a new base*/
    
                        else if(tmpBase->baseChar == 0)
                            tmpBase->baseChar = *seqCStr;
    
                        ++tmpBase->numBasesUL;
    
                        /*This avoids scattered deletions*/
                        if(
                            cigEntryUInt > 0 &&  /*More insertions*/
                            !(incBase->nextBase->errTypeChar & 2)
                                    /*^^Next base is not an insertion*/
                        ) { /*If have other insertions to process*/
                            /*Make a new base for next insertion*/
                            tmpBase =
                                malloc(sizeof(struct baseStruct));
    
                            tmpBase->nextBase = incBase->nextBase;
                            incBase->nextBase = tmpBase;
                            lastBase = incBase;
                            incBase = tmpBase;
    
                            /*Set the defaulst for the new base*/
                            incBase->altBase = 0;
                            incBase->errTypeChar = 2;
                            incBase->baseChar = 0;
                            incBase->numBasesUL = 0;
                        } /*If have other insertions to process*/
    
                        else /*Can move to the next base*/
                        { /*Else is safe to move to the next base*/
                            lastBase = incBase;
                            incBase = incBase->nextBase;
                        } /*Else is safe to move to the next base*/
    
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to add to the bases list*/
                
                    break;
                /*Switch: for insertions*/

                /******************************************************\
                * Fun-22 Sec-5 Sub-12: Ignore deletions
                \******************************************************/

                case 'D':
                    /*Deletions will be marked as empty and will be 
                      caught if their are to few bases, so can ignore.*/

                    while(cigEntryUInt > 0)
                    { /*While have bases to ignore*/
                        --cigEntryUInt;

                        /*Move past other sequences insertions*/
                        while(incBase->errTypeChar & 2)
                            incBase = incBase->nextBase; 

                        lastBase = incBase;
                        incBase = incBase->nextBase;
                    } /*While have bases to ignore*/

                    break;

                /******************************************************\
                * Fun-22 Sec-5 Sub-13: Ignore soft masking
                \******************************************************/

                case 'S':
                /*Switch: ingnore soft maskes 'S'*/
                    while(cigEntryUInt > 0)
                    { /*While have bases to ignore*/
                        --cigEntryUInt;
                        ++seqCStr;
                        ++qCStr;
                    } /*While have bases to ignore*/

                    break;
                /*Switch: ingnore soft maskes 'S'*/
            } /*switch: check the error type & add bases to consensus*/
        } /*While not at the end of the sam alignment sequence*/

        ++numSeqUL;
        errUChar = readSamLine(samStruct, stdinFILE);
    } /*While have alignments to read in from the sam file*/

    pclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-6: Merge bases into a single majority consensus
    ^    fun-22 sec-6 sub-1: Set up for deciding bases to keep
    ^    fun-22 sec-6 sub-2: Find best alternate base for each position
    ^    fun-22 sec-6 sub-3: Check if should keep base
    ^    fun-22 sec-6 sub-4: Add the base & Q-score to the sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-22 Sec-6 Sub-1: Set up for deciding bases to keep
    \******************************************************************/

    /*Minimum number of bases needed keep an SNP, match, or insterion*/
    minNumBasesUL = numSeqUL * *minBasePerFlt;
    minInsUL = numSeqUL * *minInsPercFlt;

    seqCStr = samStruct->samEntryCStr;
    samStruct->readLenUInt = 0;
    lastBase = 0;                     /*For freeing the list*/

    /******************************************************************\
    * Fun-22 Sec-6 Sub-2: Find the best alternate base for each position
    \******************************************************************/

    while(headBase != 0)
    { /*For all positions in the consensus, remove non-majority bases*/
        incBase = headBase;         /*Move to 1st alternative base*/
        tmpBase = headBase->altBase;
        numBasesUL = incBase->numBasesUL;
       
        while(tmpBase != 0)
        { /*While have bases to remove*/
            numBasesUL += tmpBase->numBasesUL;

            if(incBase->numBasesUL < tmpBase->numBasesUL)
                headBase = freeBaseStruct(&incBase, lastBase);
                /*Freeing puts tmpBase were incBase is*/

            else
            { /*If removing the current base in the list*/
                /*Need to make sure freeBaseStruct does nothing funny*/
                incBase->altBase = tmpBase->altBase;
                tmpBase->altBase = 0;
                freeBaseStruct(&tmpBase, lastBase);
            } /*If removing the current base in the list*/

            tmpBase = incBase->altBase; /*Move to next alternate base*/
        } /*While have bases to remove*/

        /**************************************************************\
        * Fun-22 Sec-6 Sub-3: Check if should keep base
        \**************************************************************/

        if(headBase->errTypeChar & 2)
        { /*If is an insertion*/
            if(numBasesUL < minInsUL)
            { /*If have to few insertions to have support*/
                headBase = freeBaseStruct(&headBase, lastBase);
                continue;                /*Move on to the next base*/       
            } /*If have to few insertions to have support*/
        } /*If is an insertion*/

        else
        { /*If is a match or SNP*/
            if(numBasesUL < minNumBasesUL)
            { /*If have to few SNPs or matches to have support*/
                headBase = freeBaseStruct(&headBase, lastBase);
                continue;                /*Move on to the next base*/       
            } /*If have to few SNPs or matches to have support*/
        } /*If is a match or SNP*/

        /**************************************************************\
        * Fun-22 Sec-6 Sub-4: Add the base & Q-score to the sequence
        \**************************************************************/

        /*Put the base into the sequence*/
        *seqCStr = headBase->baseChar; /*Set the base*/
        ++seqCStr; /*Move to next sequence entry*/

        /*Free base and move to the next base*/
        headBase = freeBaseStruct(&headBase, lastBase);
    } /*For all positions in the consensus, remove non-majority bases*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-22 Sec-7: Print out cosensus & do clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(numMisSeqUL > minNumBasesUL)
        return 16;                        /*If had to few mapped reads*/

    *seqCStr = '\0';

    /*Reset the sequence & q-score line pointers (for clairity)*/
    seqCStr = samStruct->samEntryCStr;

    /*Write consensus as fasta file*/
    stdinFILE = fopen(binStruct->consensusCStr, "w");
    fprintf(stdinFILE, ">%s\n%s\n", binStruct->consensusCStr, seqCStr);
    fclose(stdinFILE);

    return 1;
} /*simpleMajCon*/

/*---------------------------------------------------------------------\
| Output:
|   o Creates:
|      - Fasta file with reads from the fastq file
|   o Returns:
|      - 1: Success
|      - 2: No fastq file (or invalid fastq file)
|      - 4: No fasta file
|      - 64: memory allocation error
\---------------------------------------------------------------------*/
unsigned char fqToFa(
    char *fqToCnvtCStr,         /*File name of fastq to convert*/
    char *outFaCStr,            /*File name of new fasta file*/
    struct samEntry *samStruct  /*Holds fastq file entries*/
) /*Converts a fastq file to a fasta file*/
{ /*fqToFa*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-23 TOC: Sec-1 Sub-1: fqToFq
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char *tmpCStr = 0;
    unsigned long lenFqUL = 0;
    FILE *fqFILE = 0;
    FILE *faFILE = 0;

    fqFILE = fopen(fqToCnvtCStr, "r");

    if(fqFILE == 0)
        return 2;

    fseek(fqFILE, 0, SEEK_END); /*Go to end of file*/
    lenFqUL = ftell(fqFILE);  /*Get offset (gives file size)*/
    fseek(fqFILE, 0, SEEK_SET); /*Go back to start of file*/

    faFILE = fopen(outFaCStr, "w");

    if(faFILE == 0)
    { /*If could not open the fasta file*/
        fclose(fqFILE);
        return 4;
    } /*If could not open the fasta file*/
    
    while(ftell(fqFILE) < lenFqUL)
    { /*While have more file to read in*/
        if(!(readRefFqSeq(fqFILE, samStruct) & 1))
            return 64; /*Report the memory error (may not be)*/

        *samStruct->samEntryCStr = '>'; /*replace @ with >*/
        tmpCStr = samStruct->qCStr;

        while(*tmpCStr != '+')
            --tmpCStr;

        *tmpCStr = '\0'; /*Make the sequence into a c-string*/
        fprintf(faFILE, "%s", samStruct->samEntryCStr);
    } /*While have more file to read in*/

    fclose(fqFILE);
    fclose(faFILE);
    return 1;
} /*fqToFa*/

/*---------------------------------------------------------------------\
| Output:
|    o Creates:
|      - Fasta file with the consensus
|    o Modifies:
|      - clustOn->consensusCStr to have file name of created consensus 
|    o Returns:
|      - 1 for success
|      - 2 or 4 for file errors
|      - 16 for minimap2 error
|      - 64 for memory allocation error
\---------------------------------------------------------------------*/
unsigned char buildCon(
    char *threadsCStr,     /*Number threads to use with system calls*/
    unsigned char *clustUC,/*Cluter number*/
    char useMajConBl,      /*Use the majority (maj) consensus step*/
    float *minBasePerFlt,   /*For majority call del for bases < %*/
    float *minInsPercFlt,  /*for majority call min bases to keep ins*/
    unsigned char minBaseQUC, /*min Q-score to keep base (majority)*/
    unsigned char minInsQUC, /*Has min mapq to keep an insertion*/
    char useRaconConBl,    /*Use Racon in building a consnesus*/ 
    const unsigned char *rndsRaconUC, /*Number of rounds to run racon*/
    char useMedakaConBl,     /*Use the Medaka in building a consensus*/ 
    char *modelCStr,         /*Model to use with medaka*/
    const char *condaBl,     /*1: Use conda medaka install 0: python*/
    struct readBin *clustOn, /*Has fastq file to build the consensus*/
    struct samEntry *samStruct /*For reading in sequences or sam file*/
) /*Builds a consensus using the best read & top reads in clustOn*/
{ /*buildCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-24 TOC: buildCon
    '     fun-24 sec-1: variable declerations
    '     fun-24 sec-2: Build majority consensus if asked for
    '     fun-24 sec-3: Build consensus with racon if asked for
    '     fun-24 sec-4: Build consensus with medaka if asked for
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-24 Sec-1: variable declerations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char nameCStr[1024]; 
    char *tmpCStr = 0;
    uint8_t errUC = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-24 Sec-2: Build majority consensus if asked for
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(useMajConBl & 1)
    { /*If need to build a simple majority consensus first*/

        errUC = 
            simpleMajCon(
                clustUC,     /*Cluster on*/
                clustOn,     /*Holds fastq file and top read*/
                samStruct,   /*For reading in alignments*/
                minBaseQUC,  /*Min q-score needed to keep an base*/
                minInsQUC,   /*Has min mapq to keep an insertion*/
                threadsCStr, /*Number of threads to use*/
                minBasePerFlt, /*% missing bases needed for del*/
                minInsPercFlt /*% of supporting reads to keep ins*/
        );

        if(!(errUC & 1))
            return errUC;
    } /*If need to build a simple majority consensus first*/
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-24 Sec-3: Build consensus with racon if asked for
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(useRaconConBl & 1)
        buildConWithRacon(
            clustUC, /*Cluster on*/
            rndsRaconUC,/*Number rounds to run racon*/
            &useMajConBl,        /*use best read or consensus*/
            threadsCStr,
            clustOn
        ); /*Builds a consensus using racon*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-24 Sec-4: Build consensus with medaka if asked for
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(useMedakaConBl & 1)
    { /*If using medaka to polish*/

        if(useRaconConBl == 0 && useMajConBl == 0)
        { /*If using the best read*/
            /*Set up the consensus name*/
            tmpCStr =
                cStrCpInvsDelm(
                    clustOn->consensusCStr,
                    clustOn->fqPathCStr
            ); /*Copy over consensus name to convert to fastq*/

            /*Details to add for the final part of the name*/
            tmpCStr -= 6; /*Get to end of .fastq*/
            tmpCStr = cStrCpInvsDelm(tmpCStr, "--cluster-");
            tmpCStr = uCharToCStr(tmpCStr, *clustUC);
            tmpCStr=cStrCpInvsDelm(tmpCStr,"--consensus.fasta");

            /*Make a fasta of the fastq file*/
            fqToFa(clustOn->consensusCStr, nameCStr, samStruct);
        } /*If using the best read*/

        /*Build the consensus with medaka*/
        errUC = medakaPolish(condaBl, modelCStr, threadsCStr, clustOn);

        if(!(errUC & 1))
            return errUC;
    } /*If using medaka to polish*/

    return 1;
} /*buildCon*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies
|    - refStruct: To hold the sequence (no header)
|  - Returns
|    - 1 if succeeded
|    - 2 if file does not exist
|    - 4 invalid file
|    - 64 memory allocation error
| Note:
|  - Fasta file should only have one sequence
\---------------------------------------------------------------------*/
unsigned char readInConFa(
    char *conFaToReadCStr, /*Name of fasta file with the consensus*/
    struct samEntry *refStruct /*Sam struct to hold consensus*/
) /*Reads in reference sequence in fasta file*/
{ /*readInConFa*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-25 TOC: readInConFa
    '   fun-25 sec-1: variable declerations
    '   fun-25 sec-2: Check if file exists and read in header
    '   fun-25 sec-2: Read in the sequence
    '   fun-25 sec-3: Set up q-score null entry, get length, and return
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-25 Sec-1: variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    uint32_t incBuffUInt = 1000;
    FILE *faFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-25 Sec-2: Check if file exists and move past header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the reference file for reading*/
    faFILE = fopen(conFaToReadCStr, "r");

    if(faFILE == 0)
        return 2;

    blankSamEntry(refStruct);

    /*Set up null endings for lines*/
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Read in the header, I know it worked, because of minimap2*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, faFILE);

    while(
     !(*(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) !='\0' ||
       *(refStruct->samEntryCStr + refStruct->lenBuffULng -2)!='\n')
    ) { /*While have a header to read in*/
        /*Read in the next part of the header*/
        fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, faFILE);

        if(*refStruct->samEntryCStr == '\0')
            return 4; /*Failed to read in anything*/

        /*Resetup markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 1) = '\0';
        *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';
    } /*While have a header to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-25 Sec-2: Read in the sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    tmpCStr = refStruct->samEntryCStr;
    *(refStruct->samEntryCStr + refStruct->lenBuffULng - 2) = '\0';

    /*Read in the first part of the sequence*/
    fgets(refStruct->samEntryCStr, refStruct->lenBuffULng, faFILE);

    while(*(refStruct->samEntryCStr+refStruct->lenBuffULng-2) > 16)
    { /*While on the sequence line*/
        refStruct->samEntryCStr =
            realloc(
                refStruct->samEntryCStr,
                refStruct->lenBuffULng + incBuffUInt
            );

        if(refStruct->samEntryCStr == 0)
        { /*memory allocation error*/
            fclose(faFILE);
            fclose(faFILE);
            return 64;
        } /*Memory allocation error*/

        /*Set pointer to new buffer*/
        tmpCStr = refStruct->samEntryCStr + refStruct->lenBuffULng;
        refStruct->lenBuffULng += incBuffUInt; /*Update buff size*/

        /*Reset new line markers*/
        *(refStruct->samEntryCStr + refStruct->lenBuffULng-2) ='\0';

        /*Read in next part of reference sequence*/
        fgets(tmpCStr, refStruct->lenBuffULng, faFILE);


        if(*refStruct->samEntryCStr == '\0')
            return 4; /*Failed to read in anything*/
    } /*While on the sequence line*/

    fclose(faFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-25 Sec-3: Set up q-score null entry, get length, and return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    refStruct->seqCStr = refStruct->samEntryCStr;

    /*Find the end of the sequence*/
    tmpCStr = refStruct->seqCStr;
    refStruct->readLenUInt = 0;    /*So can find the length as well*/

    while(*tmpCStr > 16)
    { /*Find hte end of the sequence*/
        ++tmpCStr;
        ++refStruct->readLenUInt;
    } /*Find hte end of the sequence*/

    /*Set up Q-score entry, so score reads knows is empty*/
    refStruct->qCStr = tmpCStr;
    strcpy(refStruct->qCStr, "\t*\t");     /*marking no Q-score entry*/

    return 1;
} /*readInConFa*/

/*---------------------------------------------------------------------\
| Output:
|   o Returns:
|     - The next base in the list (if altBase != 0, returns altBase)
|   o Frees:
|     - The baseToFree from memory
|   o Modifies:
|     - If their is an alternative base (baseToFree->altBase != 0)
|         o lastBase->nextBase is set to baseToFree->altBase
|         o bastToFree is set to bastToFree->altBase
|     - If their is not alternative base (bastToFree->altBase == 0)
|         o lastBase->nextBase is set to baseToFree->nextBase
|         o bastToFree is set to 0
| Note:
|    o This function assumes that the next base pointer (nextBase) is
|      is alwasy 0 (not set) for teh alternate base pointer (altBase==0)
\---------------------------------------------------------------------*/
struct baseStruct * freeBaseStruct(
    struct baseStruct **baseToFree, /*Insertion list to free*/
    struct baseStruct *lastBase     /*Base to assign pointers to*/
) /*Frees an base from a linked list of bases*/
{ /*freeBaseStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-26 TOC: Sec-1 Sub-1: freeBaseStruct
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(*baseToFree == 0)
        return 0;
    
    if((*baseToFree)->altBase != 0)
    { /*If have an alternative base to keep as the next base*/
        (*baseToFree)->altBase->nextBase = (*baseToFree)->nextBase;

        if(lastBase != 0)
        { /*If have a previous base*/
            lastBase->nextBase = (*baseToFree)->altBase;
            free(*baseToFree);
            *baseToFree = lastBase->nextBase;
        } /*If have a previous base*/

        else
        { /*Else this was the first base*/
            lastBase = (*baseToFree)->altBase;
            free(*baseToFree);
            *baseToFree = lastBase;
        } /*Else this was the first base*/

        return *baseToFree;
    } /*If have an alternative base to keep as the next base*/

    else
    { /*Else if their is only a next base*/
        if(lastBase != 0)
        { /*If have a previous base in the list*/
            lastBase->nextBase = (*baseToFree)->nextBase;
            free(*baseToFree);
            *baseToFree = 0; /*Tell that their are no alternate bases*/
            return lastBase->nextBase;
        } /*If have a previous base in the list*/

        else
        { /*Else baseToFree was the first base*/
            lastBase = (*baseToFree)->nextBase;
            free(*baseToFree);
            *baseToFree = 0; /*Tell that their are no alternate bases*/
            return lastBase;
        } /*Else baseToFree was the first base*/
    } /*Else if their is only a next base*/
} /*freeBaseStruct*/
