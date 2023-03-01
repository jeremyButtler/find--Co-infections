/*##############################################################################
# Name: fastqGrepSearchFastq
# Use:
#    Extracts target reads from fastq file
# Requires:
#    fastqGrepAVLTree
#    fastqGrepStructs (called by fastqGrepAVLTree)
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepSearchFastq TOC:
#    fun-1: fastqExtractTree: 
#    fun-2: buildAvlTree: Builds an balanced tree of read id's from filter file
#    fun-3: extractReads: Extract reads from fastq file with tree/hash
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fqGetIdsSearchFq.h"

/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 2 if invalid filter file
|     o 4 if invalid input fastq file
|     o 8 if could not open the output file
|     o 16 if both filter and fastq file coming from stdin
\---------------------------------------------------------------------*/
uint8_t fastqExtract(
    char *filtPathCStr,        /*Path to file with read ids to extract*/
    char *fqPathCStr,          /*Path to fastq file to extract reads*/
    char *outPathCStr,         /*Path to fastq file to to write reads*/
    uint8_t sizeReadStackUChar,/*Number of elements to use in stack*/
    uint32_t lenBuffUI,     /*Size of buffer to read input with*/
    uint8_t hashSearchChar,    /*1: do hash search, 0: do Tree search*/
    uint8_t printReverseChar   /*1: Keep reads in filter file
                                 0: ingore reads in filter file*/
) /*Searches and extracts reads from a fastq file using read id's*/
{ /*fastqExtract*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Pull out reads of interest
    '    fun-1 sec-1: Variable declerations
    '    fun-1 sec-2: Check if files exist
    '    fun-1 sec-3: Build the tree for reads
    '    fun-1 sec-4: Call tree or hash table function to search file
    '    fun-1 sec-5: Handle errors, clean up, & exit
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char buffCStr[lenBuffUI + 1];  /*Holds the input line*/

    /*Majic number for kunth multiplicative hashing*/
    unsigned long majicNumULng = 0;
    uint8_t digPerKeyUChar = 0;  /*Number digitis used per key in hash*/

    /*Assume hash tabel failed, makes AVL error out when fail*/
    uint8_t hashFailedBool = 1;
    uint64_t hashSizeULng = 0;       /*Will hold Size of hash table*/

    uint64_t fastqErrULng = 0;       /*Tells if error in fastq entry*/

    /*Stack to use for searching tree. (depth = 73 = 10^18 nodes)*/
    struct readNodeStack readStack[sizeReadStackUChar + 2];
    struct readInfo *readTree = 0;
    struct readInfo **hashTbl = 0;

    FILE *filtFILE = 0;          /*to file with ID's to search for*/
    FILE *fqFILE = 0; /*fastq file to search*/
    FILE *outFILE = 0;           /*File to write extracted reads to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Check if files exist
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(filtPathCStr == 0 && fqPathCStr == 0)
        return 16;

    if(filtPathCStr == 0)
         filtFILE = stdin;
    else
        filtFILE = fopen(filtPathCStr, "r");

    if(filtFILE == 0)
        return 2;      /*No input file with read ids*/

    if(fqPathCStr == 0)
        fqFILE = stdin;
    else
        fqFILE = fopen(fqPathCStr, "r");

    if(fqFILE == 0)
    { /*If could not open the fastq file*/
        fclose(filtFILE);
        return 4;      /*No input file with reads to extract*/
    } /*If could not open the fastq file*/

    if(outPathCStr == 0) 
        outFILE = stdout;
    else
        outFILE = fopen(outPathCStr, "w");

    if(outFILE == 0)
    { /*if I could not open the output file*/
        fclose(filtFILE);
        fclose(fqFILE);
    } /*if I could not open the output file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Build the tree or hash table for reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make sure start & end of my stacks are marked*/
    readStack[0].readNode = 0;
    readStack[sizeReadStackUChar + 1].readNode = 0;

    if(hashSearchChar == 0)
    { /*If just using the avl tree for searching*/
        readTree =
            buildAvlTree(
                filtFILE,    /*File with target read ids*/
                readStack,   /*Stack for searching trees*/
                buffCStr,  /*Buffer to hold one line from file*/
                lenBuffUI /*Size of buffer*/
        ); /*Build the tree of reads to search*/
    } /*If just using the avl tree for searching*/

    else
    { /*Else I am searching using a hash function*/
        hashTbl =
            makeReadHash(
                filtFILE,     /*Name of file with read ids to hash*/
                buffCStr,     /*Buffer to hold one line from file*/
                lenBuffUI,   /*Size of buffer to hold each line of the file*/
                readStack,      /*Stack for searching trees in hash table*/
                &hashSizeULng,  /*Will hold Size of hash table*/
                &digPerKeyUChar, /*Number digitis used per key in hash*/
                &majicNumULng,   /*Will hold the majic number*/
                &hashFailedBool /*Holds if manged to make hash table*/
        ); /*Build the hash table*/
    } /*Else I am searching using a hash function*/

    fclose(filtFILE); /*No longer need open*/

    if(readTree == 0 && hashFailedBool == 1)
    { /*If calloc errored out in making the tree*/
        fprintf(
            stderr,
            "calloc failed: fastqGrepSearchFastq.c: Fun-1: 99\n"
        ); /*Warn user calloc failed*/

        return 0;
    } /*If calloc errored out in making the tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Call the tree or hash table function to search the file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fastqErrULng =
        extractReads(
            fqFILE,         /*to fastq file to search through*/
            outFILE,           /*File to write extracted reads to*/
            buffCStr,        /*Buffer to hold input from fastq file*/
            lenBuffUI,      /*Size of buffCStr*/
            majicNumULng,      /*Holds majick number for kunths hash*/
            digPerKeyUChar,    /*Digits needed to get a key*/
            &printReverseChar,
            readTree,          /*AVL tree to search if hashTbl == 0*/
            hashTbl            /*hash table to search*/
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Handle errors, clean up, & exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fclose(fqFILE); /*No longer need open*/
    fclose(outFILE); /*No longer need open*/

    /*Check if freeing tree or hash table with tree*/
    if(hashSearchChar == 0)
        freeReadTree(&readTree, readStack);
    else
    { /*If used hashing, free the hashing variables*/
        freeHashTbl(&hashTbl, &hashSizeULng, readStack);
    } /*If used hashing, free the hashing variables*/

    if(fastqErrULng == 0)
        return 0;         /*Not a valid fastq file*/

    return 1; /*Success*/
} /*fastqExtract*/

/*##############################################################################
# Output:
#    Returns: balanced readInfo tree with all read id's in filtFILE
#    Returns: 0 if calloc errored out
##############################################################################*/
struct readInfo * buildAvlTree(
    FILE *filtFILE,                 /*File with read ids to keep or ignore*/
    struct readNodeStack *readStack,  /*Stack to use in building AVL tree*/
    char *buffCStr,        /*Buffer to hold one line from file*/
    uint32_t lenBuffUI         /*Size of buffer to read each line*/
) /*Builds a readInfo tree with read id's in filtFILE*/
{ /*buildAvlTree function*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC:
    '    fun-2 Sec-1: Variable declerations
    '    fun-2 Sec-2: Read id's in filter file and build tree
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpIdCStr = 0;
    unsigned char maxHexChar = 1;  /*Size of bigNum array*/
    uint64_t lenInputULng = lenBuffUI; /*length of read id*/

    struct readInfo *readTree = 0;
    struct readInfo *lastRead = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Read id's in filter file and build tree
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Intalize read in loop*/
   *(buffCStr + lenBuffUI) = '\0';
   tmpIdCStr = buffCStr + lenBuffUI; /*Force initial read in*/
   lenInputULng = lenBuffUI; /*Make sure read in first read*/


   do { /*While ids to read in*/
       lastRead = 
           cnvtIdToBigNum(
               buffCStr,
               lenBuffUI,
               &tmpIdCStr,
               &lenInputULng,
               &maxHexChar,
               filtFILE
       ); /*Read in id and convert to big number*/

       if(lastRead == 0)
       { /*If was a falied read*/
           if(lenInputULng == 0)
           { /*If was a memory allocation error (message already printed)*/
               freeReadTree(&readTree, readStack);
               return 0;
           } /*If was a memory allocation error (message already printed)*/

           break; /*end of file*/
       } /*If was a falied read*/

       if(insertNodeIntoReadTree(lastRead, &readTree, readStack) == 0)
       { /*If id is in tree, need to free*/
           freeReadInfoStruct(&lastRead); /*If id already in tree*/
           lastRead = readTree;           /*Prevent loop ending early*/
       } /*If id is in tree, need to free*/

        while(*tmpIdCStr != '\n')
        { /*While not on next entry*/
            tmpIdCStr++;

            if(*tmpIdCStr == '\0')
            { /*If at the end of the buffer, but not at start of read*/
                if(lenInputULng < lenBuffUI)
                    break; /*At end of file*/

                lenInputULng = fread(
                                   buffCStr,
                                   sizeof(char),
                                   lenBuffUI,
                                   filtFILE
                ); /*Read in more of the file*/

                *(buffCStr + lenInputULng) = '\0';/*make sure a c-string*/
                tmpIdCStr = buffCStr;
            } /*If at the end of the buffer, but not at start of read*/
        } /*While not on next entry*/
   } while(lastRead != 0);/*ids to read in*/

   return readTree;
} /*buildAvlTree function*/

/*##############################################################################
# Output:
#    stdout: Prints out reads in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
uint8_t extractReads(
    FILE *fqFILE,            /*fastq file to search through*/
    FILE *outFILE,              /*File to write extracted reads to*/
    char *buffCStr,        /*Buffer to hold input from fastq file*/
    uint32_t lenBuffUI,      /*Size of buffCStr*/
    unsigned long majicNumULng, /*Holds majick number for kunths hash*/
    uint8_t digPerKeyUChar,     /*Digits needed to get a key*/
    uint8_t *printNonMatchBool, /*1: print non-match, 0: print match*/
    struct readInfo *readTree,  /*For AVL search (hashTbl == 0)*/
    struct readInfo **hashTbl   /*Hash table to search for ids in*/
) /*Extract target reads from fastq file with hash table or tree*/
{ /*extractReadsInHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: Variable declerations
    #    fun-3 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t
       endOfFileChar = 0;       /*Marks if at the end of the file*/

    char
       *readPosCStr = 0,        /*Position in buffCStr*/
       dummyConvertUChar = '0',  /*blanck number to initalize bignum*/
       *startReadCStr = 0;      /*Start of read id for a fastq entry*/

    int32_t
        lenIdInt = 0;    /*Holds length of read id*/

    uint64_t
        lenInputULng = 0; /*Holds number char fread grabbed from file*/

    struct bigNum
        *idBigNum = 0;

    struct readInfo
        *lastRead = 0;          /*Holds node of read id found in tree search*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Extract target reads from fastq file
    #    fun-3 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-3 sec-3 sub-2: Determine if read is in tree
    #    fun-3 sec-3 sub-3: Decide if should keep read
    #    fun-3 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    buffCStr[lenBuffUI] = '\0';          /*Make c-string for fread*/
    readPosCStr = buffCStr + lenBuffUI; /*So will read in buffer*/
    lenInputULng = lenBuffUI;               /*So does not exit early*/

    /*Intalize the bigNum struct*/
    lenIdInt = 256;
    idBigNum = makeBigNumStruct(&dummyConvertUChar, &lenIdInt);

    if(idBigNum == 0)
    { /*If could not allocatem memory*/
        fprintf(stderr, "fastqGrepSearchFastq.c Fun-3 extractReadsInHash\n");
        return 0;
    } /*If could not allocatem memory*/

    while(endOfFileChar != 4)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        endOfFileChar =
            parseFastqHeader(
                buffCStr,     /*Buffer to hold file input*/
                &startReadCStr, /*Will hold start of read name*/
                &readPosCStr,   /*Start of read name, wil hold end*/
                &lenInputULng,
                lenBuffUI,
                &lenIdInt,     /*Will hold the lenght of read id*/
                idBigNum,
                fqFILE
        ); /*Get read name from file*/

         if(endOfFileChar == 0)
         { /*If ended to early*/
             freeBigNumStruct(&idBigNum);
             return 0; /*If fastq file ends on a header*/
         } /*If ended to early*/

         else if(endOfFileChar == 4)
             continue;            /*At end of file, return*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        if(hashTbl == 0)
            lastRead = searchTree(idBigNum, readTree); /*Use avl tree*/
        else 
        { /*Else doing a hash table*/
            lastRead = 
                findReadInHashTbl(
                    idBigNum,
                    &majicNumULng,
                    &digPerKeyUChar,
                    hashTbl
            );  /*See if read id is in the hash table*/
        } /*Else doing a hash table*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read (!! converts address to 1)
          read is match, priting matches only = 1 ^ 0 = 1
          read is match, priting non-matches only = 1 ^ 1 = 0
          read is not match, priting matches only = 0 ^ 0 = 0
          read is not match, priting non-matches only = 0 ^ 1 = 1
        */
        if(((!!lastRead) ^ *printNonMatchBool) == 0)
        { /*If is a read I am not printing out*/

             endOfFileChar =
                 moveToNextFastqEntry(
                     buffCStr,
                     &readPosCStr,   /*Start of read name, wil hold end*/
                     lenBuffUI,
                     &lenInputULng,   /*Holds how many characters fread got*/
                     fqFILE
             ); /*Move to next fastq entry*/

             if(endOfFileChar == 0)
             { /*If ended to early*/
                 freeBigNumStruct(&idBigNum);
                 return 0; /*If fastq file ends on a header*/
             } /*If ended to early*/

            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         endOfFileChar =
             printFastqEntry(
                 buffCStr,
                 &readPosCStr,   /*Start of read name, wil hold end*/
                 &startReadCStr, /*points to start of read name*/
                 lenBuffUI,
                 &lenInputULng,
                 outFILE,
                 fqFILE
         ); /*Print read & move to next read*/

         if(endOfFileChar == 0)
         { /*If ended to early*/
             freeBigNumStruct(&idBigNum);
             return 0; /*If fastq file ends on a header*/
         } /*If ended to early*/
    } /*While there are lines in the file*/

    freeBigNumStruct(&idBigNum);
    return 1; /*No errors*/
} /*extractReadsInHash*/
