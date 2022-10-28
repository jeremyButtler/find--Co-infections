/*##############################################################################
# Name: clustGraph.c
# Use: puts reads into a graph by mapq and then prints out identified clusters
# Input:
#    -f: File to build clusters for (See command to get file beneath)
#        - minimap2 --dual=yes -ax ava-ont reads.fastq reads.fastq |
#          awk 'BEGIN{OFS="\t"}; {if($1 !~ /^@/){print $1, $3, $5}}' > file.tsv;
#        - See note for details on what each option does
#        - Requried if stdout not specified
#    -min-shared-edges: Min number of edges to assume two reads are in a cluster
#        - Default
#    -min-num-mapped-reads: Min number of mapped reads needed to keep a read
#    -min-reads-per-cluster: Min number of reads to keep a cluster
#        - Default
#    -min-mapq: Min mapq to consider two reads form an edge (in a cluster)
#        - Default
#    -stdin: Take alignment file from stdout (NEED TO SET UP)
#        - minimap2 --dual=yes -ax ava-ont reads.fastq reads.fastq |
#          awk 'BEGIN{OFS="\t"}; {if($1 !~ /^@/){print $1, $3, $5}}' |
#          clustGraph.c -stdout [other options..]
# Output: stdout: prints out each read by cluster
# Note: All versus all alignment command:
#       minimap2 --dual=yes -ax ava-ont reads.fastq reads.fastq
#           --dual=yes keeps duplicates alignemnts
#               This is needed for my code to maintain a low memory footprint
#               but comes at the cost of minimap2 speed 
#           -a tells minimap2 to output a sam file
#           -x ava-ont tells minimap2 to set settings for an an all versus all
#              nanopore sequencer read alignemnt
# Note: Use awk 'BEGIN{OFS="\t"}; {if($1 !~ /^@/){print $1, $3, $5}}' > file.tsv
#       to format the sam file into the correct format for this program
#       - BEGIN{OFS="\t"}: Make sure the ouput is tab deliminated
#       - if($1 !~ /^@/): gets rid of the sam header sections
#       - print $1: prints out the subject read name
#       - print $3: prints out the reference read name
#       - print $5: prints out the mapping quality for the alignemnt
# Note: Input will look like: read\treference\tmapq\n
# Note: Line size of single input should be < 10kb
##############################################################################*/

/*
Sam file table for first 11 columns (all sam files have)
| Col | Field |  Type  |        Brief description              |
|:---:|:-----:|:------:|:-------------------------------------:|
|  1  | QNAME | String |       Query template NAME             |
|  2  | FLAG  |  Int   |          bitwise FLAG                 |
|  3  | RNAME | String |     Reference sequence NAME           |
|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
|  5  | MAPQ  |  Int   |          MAPping Quality              |
|  6  | CIGAR | String |            CIGAR string               |
|  7  | RNEXT | String | Reference name of the mate/next read  |
|  8  | PNEXT |  Int   |   Position of the mate/next read      |
|  9  | TLEN  |  Int   |      observed Template LENgth         |
| 10  |  SEQ  | String |          segment SEQuence             |
| 11  | QUAL  | String | ASCII of Phred-scaled base QUALity+33 |
*/


/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    fun-1: main: function that runs everything
#    fun-2: checkInput: check and process the user input (TO BE WRITTEN)
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <stdio.h> /*For input and output operations*/
#include "clusterGraph.h" /*Structs and funs for read (cluster) graph*/
    /*
      Includes: <stdlib.h>
                clustGraphReadTree.c              build tree with read names
                clustGraphStructsReadInfo.c       nodes in read tree
                clustGraphStructReadInfoStack.c   to build stack of readInfo
                clustGraphStructGraphNode.c       nodes in graph of reads
                clustGraphStructGraphNodeStack.c  to build stack of graphNodes
                clustGraphStructBridgeNode.c      marks bridges in graph
                clustGraphStructBridgeNodeStack.c make bridgeNode stacks
                clustGraphStructMinValues.c       holds user input
    */

char * checkInput(int *lenArgsInt,        /*Number of arugments user input*/
                  char *argsCStr[],       /*Argumenas & parameters input*/
                  char **fileCStr,        /*Holds path/name of file working on*/
                  struct minValues *minStats,/*Holds thresholds user provided*/
                  char *stdinChar         /*1: use stdin input, 0: use file*/
); /*Checks user input & puts input into variables for later use*/


int main(int lenArgsInt, char *argsCStr[])
{ /*main function*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: main function
    #    fun-1 sec-1: variable declerations
    #    fun-1 sec-2: Check user input and index file (NEED TO BUILD)
    #    fun-1 sec-3: Build my cluster graph
    #    fun-1 sec-4: Print out clusters, free graph, free readInfo tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *lineInCStr = 0,              /*Holds the input line*/
         *fileCStr = 0,                /*file to open*/
         stdinChar = 0,                /*If 1 taking input from stdin*/
         *inputChar = 0;               /*Holds arguemnt that had input error*/

    int buffSizeInt = 10000;           /*Size of buffer for reading a file*/

    unsigned long mappedReadsULng = 0; /*number reads mapped to a single read*/

    FILE *readFile = 0;                /*Input file*/

    struct readInfo *readTree = 0,     /*readInfo node tree with read names*/
                    *lastRead = 0;     /*Used to detect when change reads*/ 

    struct graphNode *graphRoot = 0;   /*Root of my read graph*/

    struct graphNodeStack *clustStack = 0;/*Keep track of clusters sharing edges
                                            with a single read*/
    struct minValues minStats;        /*min stats to keep alignment or cluster*/    

    char *helpMesgCStr = "clustGraph -f file.tsv [options...] \
             \n  -f: File to build clusters for \
             \n      - Format needs to be read\trefernce\tstat \
             \n      - Stat: any number, I used mapping quality \
             \n      - Required, unless -stdin is specified \
             \n  -stdin: Take alignment file from stdout of other program \
             \n  -min-shared-edges: Min number of shared edges needed to put read in cluster \
             \n      - Default: 50 \
             \n  -min-num-mapped-reads: Min number of edges needed to keep a read \
             \n      - Default: 100 \
             \n  -min-reads-per-cluster: Min number of reads needed to keep a cluster \
             \n      - Default: 100 \
             \n  -min-mapq: Min maping quality to keep an alignment \
             \n      - Note: This is not limited to mapping quality \
             \n      - Default: 13 \
             \n Output: tsv format to stdout: Read\\tclusterNumber \
         "; /*Help message*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check user input and index file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initMinValuesStruct(&minStats);   /*set my default values*/
    /*CHECK USER INPUT FUNCTION HERE*/

    inputChar = checkInput(&lenArgsInt,
                           argsCStr,
                           &fileCStr,
                           &minStats,
                           &stdinChar 
    ); /*Get the user input*/

    if(inputChar == 0)
    { /*Input ok, echoud user input*/
        if(stdinChar == 0)
            fprintf(stderr, "Input settings:\n    file: %s\n", fileCStr);
        else
            fprintf(stderr, "Input settings:\n    Input from stdin\n");

        fprintf(stderr,
                "    Min number of edges to assign read to cluster: %lu\n",
               minStats.minSharedEdgesUInt
        );
        fprintf(stderr,
                "    Min number of alignments to keep a read: %lu\n",
                minStats.minReadsULng
        );
        fprintf(stderr,
                "    Min number of reads to keep a cluster: %lu\n",
                minStats.minReadsPerClustULng
        );
        fprintf(stderr, 
                "    Min mapping quality to keep an aligment: %f\n",
                minStats.minMapqDbl
        );
    } /*Input ok, echoud user input*/

    else if(strcmp(inputChar, "-h") == 0)
    { /*If user wanted the help message*/
         printf("%s\n", helpMesgCStr);
         return 1;
    } /*If user wanted the help message*/

    else if(inputChar != 0)
    { /*If user wanted the help message*/
         printf("%s\n%s is invalid\n", helpMesgCStr, inputChar);
         return 1;
    } /*If user wanted the help message*/


    if(stdinChar == 0)
    { /*If taking input from file*/
        readFile = fopen(fileCStr, "r"); /*input ok, open file*/

        if(readFile == 0)
        { /*If no file was oppened*/
            printf("Input file (%s) could not be opened\n", fileCStr);
            return 1;
        } /*If no file was oppened*/
    } /*If taking input from file*/

    lineInCStr = malloc(sizeof(char) * (buffSizeInt + 1));
    *(lineInCStr + buffSizeInt) = '\0'; /*Make sure always a c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Build my cluster graph
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(stdinChar == 0) /*Input from file*/
    { /*if taking output from the file*/
        while(fgets(lineInCStr, buffSizeInt, readFile))
        { /*While there are lines in the file*/
            lastRead = insertGraphNode(lineInCStr,
                                       &readTree,
                                       &graphRoot,
                                       lastRead,
                                       &minStats,
                                       &mappedReadsULng,
                                       &clustStack
            ); /*Insert reads in graph with input from file*/

            if(lastRead == 0)
            { /*If malloc errored out*/
                fprintf(stderr, "main: cluster.c:214\n");
                printAndFreeGraph(&graphRoot, &minStats, 0);
                freeReadTree(&readTree);
                free(lineInCStr);
                fclose(readFile);
                exit(-1);
            } /*If malloc errored out*/
        } /*While there are lines in the file*/
    } /*if taking output from the file*/

    else
    { /*if taking output from stdin*/
        while(fgets(lineInCStr, buffSizeInt, stdin))
        { /*While stdin still has lines*/
            lastRead = insertGraphNode(lineInCStr,
                                       &readTree,
                                       &graphRoot,
                                       lastRead,
                                       &minStats,
                                       &mappedReadsULng,
                                       &clustStack
            ); /*insert nodes in graph using stdin input*/

            if(lastRead == 0)
            { /*If malloc errored out*/
                fprintf(stderr, "main: cluster.c:239\n");
                printAndFreeGraph(&graphRoot, &minStats, 0);
                freeReadTree(&readTree);
                free(lineInCStr);
                exit(-1);
            } /*If malloc errored out*/
        } /*While stdin still has lines*/
    } /*if taking output from stdin*/

     /*Check if should keep the last cluster*/
     if(mappedReadsULng < minStats.minReadsULng ||
        mappedReadsULng < minStats.minSharedEdgesUInt)
        { /*If read is not worth keeping*/
            if(breakCluster(&(lastRead)->nodeInGraph, &graphRoot) == 0)
            { /*If malloc failed to allocate memory*/
                fprintf(stderr, "main: cluster.c:263\n");
                return 0;
            } /*If malloc failed to allocate memory*/
         
            lastRead->doneChar = 4;  /*Ignoring read for futre*/

            while(clustStack != 0)
            { /*Resest and free the cluster stacks*/
                clustStack->graphStruct->readCntULng = 0;
                popGraphStackStruct(&clustStack);
            } /*Resest and free the cluster stacks*/
        } /*If read is not worth keeping*/

     else 
	     assignReadToCluster(lastRead->nodeInGraph,
                             &clustStack,
                             &minStats,
                             &graphRoot
         ); /*Assign the last read in read to a bridge*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Print out clusters, free graph, free readInfo tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    printAndFreeGraph(&graphRoot, &minStats, 1);
    freeReadTree(&readTree);
    free(lineInCStr);

    if(readFile != 0)
        fclose(readFile);

    exit(0);
} /*main function*/

/*##############################################################################
# Output:
#    Returns: 0 if no errors, pointer to argumet errored on for errors
#    Modifies: fileCStr, variables in minStats, & stdinChar to hold input
##############################################################################*/
char * checkInput(int *lenArgsInt,           /*Number of arugments user input*/
                 char *argsCStr[],          /*Argumenas & parameters input*/
                 char **fileCStr,         /*Holds path/name of file working on*/
                                            /*Needs to be array or allocated*/
                 struct minValues *minStats,/*Holds thresholds user provided*/
                 char *stdinChar            /*1: use stdin input, 0: use file*/
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
            *fileCStr = singleArgCStr;
        else if(strcmp(tmpCStr, "-stdin") == 0)
        { /*If if taking input from stdin*/
            *stdinChar = 1;
            intArg--;                  /*Account for incurment at end of loop*/
        } /*If taking input from stdin*/
        else if(strcmp(tmpCStr, "-h") == 0)
            return tmpCStr;            /*Help message asked for*/

        else if(strcmp(tmpCStr, "-min-shared-edges") == 0)
            (*minStats).minSharedEdgesUInt = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-min-num-mapped-reads") == 0)
            (*minStats).minReadsULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-min-reads-per-cluster") == 0)
            (*minStats).minReadsPerClustULng = strtoul(singleArgCStr, NULL, 10);
        else if(strcmp(tmpCStr, "-min-mapq") == 0)
            (*minStats).minMapqDbl = strtoul(singleArgCStr, NULL, 10);
        else
            return tmpCStr;
        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
