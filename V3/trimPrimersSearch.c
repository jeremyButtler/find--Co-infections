/*######################################################################
# Name: trimPrimersSearch
# Use:
#   o Sets up hash table of reads with primer mappings, searchs for
#     the reads in the hash table, & trims off the primers from the 
#     target reads.
# Output
#   o Trimmed reads are output to a fastq file
# Includes:
#   - "trimPrimersHash.h"
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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimPrimersSerach SOF: Start Of Functions
'  - fun-1 trimPrimers:
'    o Wrapper functions for a series of functions that map the reads
'      to primers, builds a hash table for mappings with primers,
'      & extracts & trims reads to their primer mappings.
'   - fun-2 primReadsExtReads:
'     o Extract target reads from fastq file with hash table or tree
'   - fun-3 trimAndPrintRead
'     o Trims off primer regions and prints out the untrimmed region.
'       Multiple fastq entrries are printed out for reads with multiple
'       primer targets. 
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "trimPrimersSearch.h"

/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 2 if could not open the fasta file
|     o 4 if could not open the fastq file
|     o 8 if could not open the output file
|     o 16 if both fasta and fastq file coming from stdin
|     o 32 for invalid fastq file
|     o 64 for memory allocation error
\---------------------------------------------------------------------*/
unsigned char trimPrimers(
    char *faPathCStr,  /*Path to fasta file with primers to map*/
    char *pafFileCStr,     /*Paf with read mappings; skips minimap2*/
    char stdinPafBl,       /*Paf with read mappings from stdin*/
    char *fqPathCStr,      /*Path to fastq file with reads to trim*/
    char *outPathCStr,     /*Path to fastq file to write trimmed reads*/
    char *threadsCStr,     /*Number of threads to use with minimap2*/
    char hashSearchBl      /*1: do hash search, 0: do Tree search*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-1 TOC: trimPrimers
   '   - Wrapper functions for a series of functions that map the reads
   '     to primers, builds a hash table for mappings with primers,
   '     & extracts & trims reads to their primer mappings.
   '   o fun-1 sec-1: Variable declerations
   '   o fun-1 sec-2: Check if files exist
   '   o fun-1 sec-3: Build the tree or hash table for reads
   '   o fun-1 sec-4: extract/trim reads, handle errors, clean up, &exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Majic number for kunth multiplicative hashing*/
    unsigned char errUC = 0;       /*Tells if error in fastq entry*/
    struct readPrimHash hashST;
    FILE *pafFILE = 0;           /*For skipping minimap2*/
    FILE *fqFILE = 0;            /*For getting reads from fastq file*/
    FILE *outFILE = 0;           /*File to write extracted reads to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Check if files exist
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(pafFileCStr == 0 && !(stdinPafBl & 1))
    { /*If running minimap2*/
        /*Check if the primer fasta file exists*/
        outFILE = fopen(faPathCStr, "r");

        if(outFILE == 0)
            return 2;      /*No input file with read ids*/

        fclose(outFILE); /*Minimap2 will handel this*/
        outFILE = 0;
    } /*If running minimap2*/

    else if(!(stdinPafBl & 1))
    { /*Else if user provided a paf file*/
        pafFILE = fopen(pafFileCStr, "r");

        if(pafFILE == 0)
            return 2; /*File error*/
    } /*Else if user provided a paf file*/

    else pafFILE = stdin; /*taking input from stdin*/

    /*Check if the fastq file exits*/
    outFILE = fopen(fqPathCStr, "r");

    if(outFILE == 0)
        return 2;      /*No input file with read ids*/

    fclose(outFILE); /*Minimap2 will handel this*/
    outFILE = 0;

    /*Check if I can open the outupt file*/
    if(outPathCStr != 0) 
    { /*If not using stdout*/
        outFILE = fopen(outPathCStr, "w");

        if(outFILE == 0) return 2;

        fclose(outFILE); /*Will open when needed*/
        outFILE = 0;
    } /*If not using stdout*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Build the tree or hash table for reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make a list of read ids. This is not very efficent for the AVL
      tree, but is needed for the hash table. I am doing this in a more
      simple way to reduce the number of functions I need to make*/
    errUC =
        makeReadPrimList(
            faPathCStr,
            pafFILE,
            fqPathCStr,
            threadsCStr,
            &hashST
    ); /*Make the linked list of read ids & primer coordinates used to
         make a hash table or tree*/

    if(pafFILE != 0)
    { /*If need to close the paf file*/
        fclose(pafFILE);
        pafFILE = 0;
    } /*If need to close the paf file*/

    if(!(errUC & 1))
        return errUC;
    
    /*Make the AVL tree*/
    if(hashSearchBl == 0) readPrimListToTree(1, &hashST.readTree); 
    else errUC = readPrimListToHash(1, &hashST);

    if(!(errUC & 1))
        return errUC; /*Something errored out, likely memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-4: extract/trim reads, handle errors, clean up, & exit
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Open the outfile for the trimmed reads*/
    if(outPathCStr == 0) 
        outFILE = stdout;
    else
        outFILE = fopen(outPathCStr, "w");

    fqFILE = fopen(fqPathCStr, "r");
    errUC = extracAndTrimReads(fqFILE, outFILE, &hashST);

    fflush(outFILE); /*Make sure nothing in buffer*/
    fclose(fqFILE); /*No longer need open*/
    fclose(outFILE); /*No longer need open*/
    freeReadPrimHashST(0, 1, &hashST);

    if(errUC & 64) return 64; /*memory allocation error*/
    if(errUC & 32) return 32; /*Not a valide fastq file*/

    return 1; /*Success*/
} /*fastqExtract*/

/*---------------------------------------------------------------------\
| Output:
|    stdout: Prints out reads in hash table
|    Returns:
|      o 1 for success
|      o 32 if was not a valid fastq file
|      o 64 for memory allocation errors
\---------------------------------------------------------------------*/
unsigned char extracAndTrimReads(
    FILE *fqFILE,               /*fastq file to search through*/
    FILE *outFILE,              /*File to write extracted reads to*/
    struct readPrimHash *hashST /*Holds tree/hash table variables*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-2 TOC: primReadsExtReads
   '  - Extract target reads from fastq file with hash table or tree
   '  o fun-2 sec-1: Variable declerations
   '  o fun-2 sec-2: Extract target reads from fastq file
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char dummyConvertUChar = '0'; /*blank number to initalize bignum*/

    unsigned char EOFUC = 0;    /*Marks if at the end of the file*/
    int32_t lenIdInt = 0;  /*Holds length of read id*/

    struct samEntry samST;  /*For reading in fastq entries*/
    struct bigNum *idBigNum = 0;

    /*Holds node of read id found in tree search*/
    struct readPrim *lastRead = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Extract target reads from fastq file
    ^    fun-2 sec-3 sub-1: Get the read name of a single fastq entry
    ^    fun-2 sec-3 sub-2: Determine if read is in tree
    ^    fun-2 sec-3 sub-3: Decide if should keep read
    ^    fun-2 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Intalize the bigNum struct*/
    lenIdInt = 256;
    idBigNum = makeBigNumStruct(&dummyConvertUChar, &lenIdInt);

    if(idBigNum == 0) return 64; /*memory allocation error*/

    initSamEntry(&samST);
    EOFUC = readRefFqSeq(fqFILE, &samST, 0);
            /*Inputing 0, so that new lines are removed*/

    while(EOFUC & 1) /*0 is the readSamLin EOF flag*/
    { /*While there are lines in the file*/
        /*Convert the query id to a big number*/
        strToBackwardsBigNum(idBigNum, samST.queryCStr, &lenIdInt);

        /*Checking if doing an avl tree or hash table search*/
        if(hashST->hashTbl == 0)
            lastRead = searchReadPrimTree(idBigNum, hashST->readTree);
        else lastRead = findReadPrimInHash(idBigNum, hashST);

        /*If printing the read*/
        if(lastRead != 0) trimAndPrintRead(&samST, lastRead, outFILE);

        /*Get the next line*/
        EOFUC = readRefFqSeq(fqFILE, &samST, 0);
            /*Inputing 0, so that new lines are removed*/
        /*Normally I would blank samST, but in this case I am not
          worried, since I am only reading in a fastq entry and not
          finding any stats*/
    } /*While there are lines in the file*/

    freeBigNumStruct(&idBigNum);
    freeStackSamEntry(&samST);

    if(EOFUC == 0) return 1; /*End of file*/
    if(EOFUC & 64) return 64;

    return 32; /*No errors*/
} /*primReadsExtReads*/

/*---------------------------------------------------------------------\
| Output: Prints trimmed reads to outFILE
\---------------------------------------------------------------------*/
void trimAndPrintRead(
    struct samEntry *samST,  /*Has buffer and sequence to process*/
    struct readPrim *readIn, /*Has sorted primer coordinates to cut at*/
    FILE *outFILE            /*File to output everything to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-3 TOC: trimAndPrintRead
   '  - Trims off primer regions and prints out the untrimmed region.
   '    Multiple fastq entrries are printed out for reads with multiple
   '    primer targets. 
   '  o fun-3 sec-1: Variable declerations
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char readIdCStr[256];

   char *seqIterCStr = 0;        /*ending point to print read out*/
   char *qIterCStr = 0;        /*ending point to print read out*/
   char *seqCStr = 0;  /*Starting point to print read out*/
   char *qCStr = 0;  /*Starting point to print read out*/

   unsigned int numBasesToPrintUI = 0; /*Number of bases to print out*/
   unsigned short dupUS = 0; /*Duplicate read on*/
   unsigned int posOnUI = 0; /*Holds Position at in read*/
   struct primCord *primST = readIn->primCordST;

   seqIterCStr = readIdCStr;
   *seqIterCStr = '@';           /*header for sam file start*/
   ++seqIterCStr;

   qIterCStr = samST->queryCStr;

   while(*qIterCStr > 32)
   { /*While have a read id to copy over*/
       *seqIterCStr = *qIterCStr;
       ++seqIterCStr;
       ++qIterCStr;
   } /*While have a read id to copy over*/

   /*Copy the read id for printing*/
   cStrCpInvsDelm(seqIterCStr, "--D1_19");

   seqIterCStr = samST->seqCStr; /*Sequence entry*/
   qIterCStr = samST->qCStr;  /*Q-score entry*/

   while(primST != 0)
   { /*While have primer coordinates to trim*/
      numBasesToPrintUI = 0;

      seqCStr = seqIterCStr;
      qCStr = qIterCStr;

      while(posOnUI < primST->startUI)
      { /*While not at the start of the primer*/
          ++seqIterCStr; /*Move to the next base*/
          ++qIterCStr; /*Move to the next base*/
          ++posOnUI;
          ++numBasesToPrintUI;
      } /*While not at the start of the primer*/

      /*Get off the first primer base*/
      if(numBasesToPrintUI > 0) --numBasesToPrintUI;

      if(numBasesToPrintUI > 0)
      { /*If I need to print out an entry*/
          /*Deal with the header*/
          fprintf(outFILE, "%s-%u\n", readIdCStr, dupUS);
          ++dupUS; /*Count that I am moving to the next duplicate*/

          /*Print out the rest of the fastq entry*/
          fwrite(seqCStr, sizeof(char), numBasesToPrintUI, outFILE);
          fwrite("\n+\n", sizeof(char), 3, outFILE);
          fwrite(qCStr, sizeof(char), numBasesToPrintUI, outFILE);
          fwrite("\n", sizeof(char), 1, outFILE);
      } /*If I need to print out an entry*/

      /*Get off the primer part of the alignment*/
      while(posOnUI <= primST->endUI)
      { /*While not at the end of the primer*/
          ++seqIterCStr; /*Move to the next base*/
          ++qIterCStr; /*Move to the next base*/
          ++posOnUI;
      } /*While not at the end of the primer*/

      primST = primST->nextCord;
   } /*While have primer coordinates to trim*/

   /*Write out the remaing part of the sequence*/
   if(posOnUI < samST->readLenUInt)
   { /*If have more sequence to write out*/
       numBasesToPrintUI = samST->readLenUInt - posOnUI;

       /*Deal with the header*/
       fprintf(outFILE, "%s-%u\n", readIdCStr, dupUS);
       ++dupUS; /*Count that I am moving to the next duplicate*/

       /*Print out the rest of the fastq entry*/
       fwrite(seqIterCStr, sizeof(char), numBasesToPrintUI, outFILE);
       fwrite("\n+\n", sizeof(char), 3, outFILE);
       fwrite(qIterCStr, sizeof(char), numBasesToPrintUI, outFILE);
       fwrite("\n", sizeof(char), 1, outFILE);
   } /*If have more sequence to write out*/

   return;
} /*trimAndPrintRead*/
