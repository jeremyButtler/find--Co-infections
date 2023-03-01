/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start of Functions
'   - fun-1 fastqThreadExtractThread:
'     o Calls function to set up hash table/AVL tree and extract reads
'   - fun-2 extractReadsThread:
'     o Extracts reads by id with multi-thread support
'   - fun-3 getFileLen:
'     o Get length of file (set up to be called by a separate thread)
'   - fun-4 findStartPos:
'     o Find starting position of a thread in a read id extraction
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "fqGetIdsSearchThread.h"

pthread_mutex_t lockerMutex;

/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 2 if invalid filter file
|     o 4 if invalid input fastq file
|     o 8 if could not open the output file
|     o 16 if both filter and fastq file coming from stdin
\---------------------------------------------------------------------*/
uint8_t fastqThreadExtract(
    char *filtPathCStr,        /*Path to file with read ids to extract*/
    char *fqPathCStr,          /*Path to fastq file to extract reads*/
    char *outPathCStr,         /*Path to fastq file to to write reads*/
    unsigned char threadsUC,   /*Number of threads to use*/
    uint8_t sizeReadStackUC,   /*Number of elements to use in stack*/
    uint32_t lenBuffUI,        /*Size of buffer to read input with*/
    uint8_t hashSearchC,    /*1: do hash search, 0: do Tree search*/
    uint8_t printReverseC   /*1: Keep reads in filter file
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
    unsigned long majicNumUL = 0;
    uint8_t digPerKeyUC = 0;  /*Number digitis used per key in hash*/

    /*Assume hash tabel failed, makes AVL error out when fail*/
    uint8_t hashFailedBl = 1;
    uint64_t hashSizeUL = 0;       /*Will hold Size of hash table*/
    uint64_t fastqErrUL = 0;       /*Tells if error in fastq entry*/

    /*For multi-threading*/
    pthread_t threadsAry[threadsUC];
    struct getFileLenStruct fileLenCall;
    unsigned long startUL[threadsUC + 1];/*Starting point each thread*/
    struct extReadsST extAryST[threadsUC];

    /*Stack to use for searching tree. (depth = 73 = 10^18 nodes)*/
    struct readNodeStack readStack[sizeReadStackUC + 2];
    struct readInfo *readTree = 0;
    struct readInfo **hashTbl = 0;

    FILE *filtFILE = 0;          /*to file with ID's to search for*/
    FILE *fqAryFILE[threadsUC]; /*fastq file to search*/
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
        fqAryFILE[0] = stdin;
    else
        fqAryFILE[0] = fopen(fqPathCStr, "r");

    if(fqAryFILE[0] == 0)
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
        fclose(fqAryFILE[0]);
    } /*if I could not open the output file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Build the tree or hash table for reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make sure start & end of my stacks are marked*/
    readStack[0].readNode = 0;
    readStack[sizeReadStackUC + 1].readNode = 0;

    if(threadsUC > 1)
    { /*If I am using multiple threads*/
        fileLenCall.lenFileUL = 0;
        fileLenCall.inFILE = fqAryFILE[0];
                       /*Thread to use, ?, function, parameters*/
        pthread_create(&threadsAry[0], 0, getFileLen, &fileLenCall);
    } /*If I am using multiple threads*/

    if(hashSearchC == 0)
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
                filtFILE,    /*Name of file with read ids to hash*/
                buffCStr,    /*Buffer to hold one line from file*/
                lenBuffUI,   /*Size of buffer*/
                readStack,   /*Stack for searching trees in hash table*/
                &hashSizeUL,  /*Will hold Size of hash table*/
                &digPerKeyUC, /*Number digitis used per key in hash*/
                &majicNumUL,   /*Will hold the majic number*/
                &hashFailedBl /*Holds if manged to make hash table*/
        ); /*Build the hash table*/
    } /*Else I am searching using a hash function*/

    fclose(filtFILE); /*No longer need open*/

    if(threadsUC > 1)
        pthread_join(threadsAry[0], 0); /*Wait till 2nd thread finshes*/

    if(readTree == 0 && hashFailedBl == 1)
    { /*If calloc errored out in making the tree*/
        fprintf(
            stderr,
            "calloc failed: fastqGrepSearchFastq.c: Fun-1: 99\n"
        ); /*Warn user calloc failed*/

        return 0;
    } /*If calloc errored out in making the tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Assing file sections to each thread
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(threadsUC > 1)
    { /*If have more than one thread*/
       /*Get number of bytes for each thread to handel*/
       startUL[threadsUC] = fileLenCall.lenFileUL; /*Ending piont*/

       for(unsigned char ucThread = 0; ucThread < threadsUC; ++ucThread)
       { /*Initalize the thread values*/
           if(ucThread == 0)
               startUL[ucThread] = 0;
           else 
           { /*Else if not on the first thread*/
               if(ucThread == 1)
                   startUL[1]= fileLenCall.lenFileUL / threadsUC;
               else
                   startUL[ucThread] +=startUL[ucThread - 1]+startUL[1];

               fqAryFILE[ucThread] = fopen(fqPathCStr, "r");
               fseek(fqAryFILE[ucThread], startUL[ucThread], SEEK_SET);
           } /*Else if not on the first thread*/

           extAryST[ucThread].startPosUL = &startUL[ucThread];
           extAryST[ucThread].endPosUL = &startUL[ucThread + 1];
           extAryST[ucThread].lenBuffUI = lenBuffUI;
           extAryST[ucThread].majicNumUL = majicNumUL;
           extAryST[ucThread].digPerKeyUC = digPerKeyUC;
           extAryST[ucThread].printNonMatchBl = printReverseC;

           extAryST[ucThread].readTree = readTree;
           extAryST[ucThread].hashTbl = hashTbl;

           extAryST[ucThread].fqFILE = fqAryFILE[ucThread];
           extAryST[ucThread].outFILE = outFILE;
       } /*Initalize the thread values*/
    } /*If have more than one thread*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-4: Call tree or hash table function to search fastq
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(threadsUC < 2)
    { /*If working with only one thread*/
        fastqErrUL =
            extractReads(
                fqAryFILE[0], /*to fastq file to search through*/
                outFILE,      /*File to write extracted reads to*/
                buffCStr,     /*Buffer to hold input from fastq file*/
                lenBuffUI,    /*Size of buffCStr*/
                majicNumUL, /*Holds majick number for kunths hash*/
                digPerKeyUC, /*Digits needed to get a key*/
                &printReverseC,
                readTree,       /*AVL tree to search if hashTbl == 0*/
                hashTbl         /*hash table to search*/
        );
    } /*If working with only one thread*/

    else
    { /*Else I am working with multiple threads*/
        pthread_mutex_init(&lockerMutex,0); /*intiate my mutext*/

        for(unsigned char ucThread = 0; ucThread <threadsUC; ++ucThread)
           pthread_create(
               &threadsAry[ucThread],
               0,
               extractReadsThread,
               &extAryST[ucThread]
           ); /*Launch threads*/

        /*Join all threads so have not loose threads*/
        for(unsigned char ucThread = 0; ucThread <threadsUC; ++ucThread)
            pthread_join(threadsAry[ucThread], 0);

        pthread_mutex_destroy(&lockerMutex);
        fastqErrUL = 1;

        for(unsigned char ucThread=0; ucThread < threadsUC; ++ucThread)
        { /*loop to close all my opened files*/
            fclose(fqAryFILE[ucThread]);
            fqAryFILE[ucThread] = 0;

            if(extAryST[ucThread].retValUC == 0)
                fastqErrUL = 0;
        } /*loop to close all my opened files*/
    } /*Else I am working with multiple threads*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-5: Handle errors, clean up, & exit
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fclose(outFILE); /*No longer need open*/

    /*Check if freeing tree or hash table with tree*/
    if(hashSearchC == 0)
        freeReadTree(&readTree, readStack);
    else
    { /*If used hashing, free the hashing variables*/
        freeHashTbl(&hashTbl, &hashSizeUL, readStack);
    } /*If used hashing, free the hashing variables*/

    if(fastqErrUL == 0)
        return 0;         /*Not a valid fastq file*/

    return 1; /*Success*/
} /*fastqExtract*/

/*---------------------------------------------------------------------\
| Output:
|    stdout: Prints out reads in hash table
\---------------------------------------------------------------------*/
void * extractReadsThread(
    void *parmST /*extReadsST Structer with parameters*/
) /*Extract target reads from fastq file with hash table or tree*/
{ /*extractReadsInHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC:
    #    fun-2 sec-1: Variable declerations
    #    fun-2 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Cast input as the required struct*/
    struct extReadsST *extParmST = (struct extReadsST *) parmST;
    uint8_t endOfFileC = 0;       /*Marks if at the end of the file*/

    char *buffCStr = malloc(sizeof(char) * extParmST->lenBuffUI + 1);
    char *readPosCStr = 0;      /*Position in buffCStr*/
    char dummyConvertUC = '0';  /*blank number to initalize bignum*/
    char *startReadCStr = 0;    /*Start of read id for a fastq entry*/

    int32_t lenIdInt = 0;    /*Holds length of read id*/

    uint64_t lenInputUL = 0;/*Holds number char fread grabbed from file*/
    uint64_t posInFileUL = 0; /*Position at in the file*/

    struct bigNum *idBigNum = 0;

    struct readInfo *lastRead = 0; /*Holds node of read id found*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Extract target reads from fastq file
    ^    fun-2 sec-3 sub-1: Get the read name of a single fastq entry
    ^    fun-2 sec-3 sub-2: Determine if read is in tree
    ^    fun-2 sec-3 sub-3: Decide if should keep read
    ^    fun-2 sec-3 sub-4: Keeping read, print out read
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(buffCStr == 0)
    { /*If had a memory allocation error*/
        free(buffCStr);
        extParmST->retValUC = 64;
        pthread_exit(0);
    } /*If had a memory allocation error*/

    /*Make sure the buffer is a c-string for fread*/
    *(buffCStr + extParmST->lenBuffUI - 1) = '\0';
    lenInputUL =
            fread(
                buffCStr,
                sizeof(char),
                extParmST->lenBuffUI,
                extParmST->fqFILE
    ); /*Read in the first few bytes*/

    readPosCStr = buffCStr;

    if(*extParmST->startPosUL != 0)
    { /*If this is not the first thread*/
        findStartPos(
            &buffCStr,
            &readPosCStr,
            &lenInputUL,
            extParmST
        ); /*Find starting position of first read*/
    } /*If this is not the first thread*/

    else
        extParmST->retValUC = 1; /*Else at the start of the file*/

    if(!(extParmST->retValUC & 1))
        return 0; /*Something when wrong*/

    /*Intalize the bigNum struct*/
    lenIdInt = 256;
    idBigNum = makeBigNumStruct(&dummyConvertUC, &lenIdInt);

    if(idBigNum == 0)
    { /*If could not allocatem memory*/
        free(buffCStr);
        extParmST->retValUC = 64;
        pthread_exit(0);
    } /*If could not allocatem memory*/

    while(endOfFileC != 4)
    { /*While there are lines in the file*/

        /**************************************************************\
        * Fun-2 Sec-3 Sub-1: Get the read name of a single fastq entry
        \**************************************************************/

        endOfFileC =
            parseFastqHeader(
                buffCStr,       /*Buffer to hold file input*/
                &startReadCStr, /*Will hold start of read name*/
                &readPosCStr,   /*Start of read name, wil hold end*/
                &lenInputUL,
                extParmST->lenBuffUI,
                &lenIdInt,     /*Will hold the lenght of read id*/
                idBigNum,
                extParmST->fqFILE
        ); /*Get read name from file*/

         if(endOfFileC == 0)
         { /*If ended to early*/
             free(buffCStr);
             freeBigNumStruct(&idBigNum);
             extParmST->retValUC = 0;
             return 0; /*If fastq file ends on a header*/
         } /*If ended to early*/

         else if(endOfFileC == 4)
             continue;            /*At end of file, return*/

        /**************************************************************\
        * Fun-2 Sec-3 Sub-2: Determine if read is in tree
        \**************************************************************/

        /*Check if using hash table or AVL tree*/
        if(extParmST->hashTbl == 0)
            lastRead = searchTree(idBigNum, extParmST->readTree);
        else 
        { /*Else doing a hash table*/
            lastRead = 
                findReadInHashTbl(
                    idBigNum,
                    &(extParmST->majicNumUL),
                    &(extParmST->digPerKeyUC),
                    extParmST->hashTbl
            );  /*See if read id is in the hash table*/
        } /*Else doing a hash table*/

        /**************************************************************\
        * Fun-2 Sec-3 Sub-3: Decide if should keep read
        \**************************************************************/

        /*Check if should print read (!! converts address to 1)
          read is match, priting matches only = 1 ^ 0 = 1
          read is match, priting non-matches only = 1 ^ 1 = 0
          read is not match, priting matches only = 0 ^ 0 = 0
          read is not match, priting non-matches only = 0 ^ 1 = 1
        */
        if(((!!lastRead) ^ extParmST->printNonMatchBl) == 0)
        { /*If is a read I am not printing out*/
             endOfFileC =
                 moveToNextFastqEntry(
                     buffCStr,
                     &readPosCStr, /*Start of read name, wil hold end*/
                     extParmST->lenBuffUI,
                     &lenInputUL,/*Holds how many characters fread got*/
                     extParmST->fqFILE
             ); /*Move to next fastq entry*/

             if(endOfFileC == 0)
             { /*If ended to early*/
                 free(buffCStr);
                 freeBigNumStruct(&idBigNum);
                 extParmST->retValUC = 0;
                 pthread_exit(0);
             } /*If ended to early*/

             if(ftell(extParmST->fqFILE) >= *extParmST->endPosUL)
             { /*If on the last buffer read in*/
                 /*Find the positoin in the file*/
                 posInFileUL =
                     1 +                        /*Account for \n or @*/
                     ftell(extParmST->fqFILE) - /*piont in file*/
                     lenInputUL -          /*Number of bytes in buffer*/
                     (readPosCStr - buffCStr);  /*byte of position on*/


                 if(posInFileUL >= *extParmST->endPosUL)
                 { /*If finshed all entries assigned to this thread*/
                     free(buffCStr);
                     freeBigNumStruct(&idBigNum);
                     extParmST->retValUC = 1;
                     pthread_exit(0);
                 } /*If finshed all entries assigned to this thread*/
             } /*If on the last buffer read in*/

            continue;
         } /*If is a read I am not printing out*/

        /***************************************************************
        * Fun-2 Sec-3 Sub-4: Keeping read, print out read
        \**************************************************************/

         /*allow only one thread print fastq entries at a time*/
         pthread_mutex_lock(&lockerMutex);

         endOfFileC =
             printFastqEntry(
                 buffCStr,
                 &readPosCStr,   /*Start of read name, wil hold end*/
                 &startReadCStr, /*points to start of read name*/
                 extParmST->lenBuffUI,
                 &lenInputUL,
                 extParmST->outFILE,
                 extParmST->fqFILE
         ); /*Print read & move to next read*/

         /*Allow another thread to print out a fastq entry*/
         pthread_mutex_unlock(&lockerMutex);

         if(endOfFileC == 0)
         { /*If ended to early*/
             free(buffCStr);
             freeBigNumStruct(&idBigNum);
             extParmST->retValUC = 0;
             pthread_exit(0);
         } /*If ended to early*/

         if(ftell(extParmST->fqFILE) >= *extParmST->endPosUL)
         { /*If on the last buffer read in*/
             /*Find the positoin in the file*/
             posInFileUL =
                 1 +                        /*Account for \n or @*/
                 ftell(extParmST->fqFILE) - /*piont in file*/
                 lenInputUL -              /*Number of bytes in buffer*/
                 (readPosCStr - buffCStr);  /*byte of position on*/


             if(posInFileUL >= *extParmST->endPosUL)
             { /*If finshed all entries assigned to this thread*/
                 free(buffCStr);
                 freeBigNumStruct(&idBigNum);
                 extParmST->retValUC = 1;
                 pthread_exit(0);
             } /*If finshed all entries assigned to this thread*/
         } /*If on the last buffer read in*/
    } /*While there are lines in the file*/

    free(buffCStr);
    freeBigNumStruct(&idBigNum);
    extParmST->retValUC = 1;
    pthread_exit(0);
} /*extractReadsInHash*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o lenFileUL in parmStruct to hold the file length
\---------------------------------------------------------------------*/
void * getFileLen(
   void *parmST /*getFileLenStruct having file & return variable*/ 
) /*Get the length of the file in parmStruct. Here for multi threading*/
{ /*getFileLen*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: Sec-1 Sub-1: getFileLen
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    struct getFileLenStruct *parmStruct =
        (struct getFileLenStruct *) parmST;

    fseek(parmStruct->inFILE, 0, SEEK_END);
    parmStruct->lenFileUL = ftell(parmStruct->inFILE);
    fseek(parmStruct->inFILE, 0, SEEK_SET);

    pthread_exit(0);
} /*getFileLen*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 if have more file to read into the buffer
|     o 0 if read remaning part of file into buffer
|   - Modifies:
|     o startOfBuffCStr to piont to start of the next read
|     o If neeed to read in more bytes, will modify buffCStr
|     o lenInUL to hold the number of bytes read into buffCStr
|     o extParmST->startPosUL to point to header of the next read
\---------------------------------------------------------------------*/
unsigned char findStartPos(
    char **buffCStr,             /*Pionts to the end of the buffer*/
    char **startOfBuffCStr,      /*Pionts to the start of the buffer*/
    uint64_t *lenInUL,           /*Number bytes read from fread*/
    struct extReadsST *extParmST /*Holds parameters to use*/
) /*Finds the next read after the starting position*/
{ /*findStartPos*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: findStarPos
    '   o fun-4 sec-1: Variable declerations
    '   o fun-4 sec-2: Find a point were I know I read past a header
    '   o fun-4 sec-3: Update file location variables
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t lastHeaderUL = 0; /*Postion of header to start at*/
    unsigned long cntUL = *extParmST->startPosUL;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-2: Find a point were I know I have read past a header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(lastHeaderUL == 0)
    { /*While I have to find the starting point*/
        while(**startOfBuffCStr != ' ' || **startOfBuffCStr != '\t')
        { /*While I have not found a new entry*/
             if(**startOfBuffCStr == '+')
             { /*If I may be on the spacer entry*/
                 if(lastHeaderUL != 0 &&
                    *(*startOfBuffCStr - 1) == '\n' &&
                    *(*startOfBuffCStr + 1) == '\n'
                 ) break; /*If on the spacer entry*/
             } /*If I may be on the spacer entry*/

             if(**startOfBuffCStr == '@')
             { /*If this may be the header*/
                 /*Check if I may have found a header*/
                 if(*startOfBuffCStr != *buffCStr && 
                    *(*startOfBuffCStr - 1) == '\n'
                 ) { /*If I may have found the header line*/
                     lastHeaderUL = cntUL;
                     *extParmST->startPosUL = lastHeaderUL;
                 } /*If I may have found the header line*/
             } /*If this may be the header*/

             if(**startOfBuffCStr == '\0')
             { /*If need to read in more buffer*/
                 *lenInUL =
                     fread(
                         *buffCStr,
                         sizeof(char),
                         extParmST->lenBuffUI,
                         extParmST->fqFILE
                 ); /*read in next part of file*/

                 if(*lenInUL == 0)
                 { /*If at the end of the file*/
                     extParmST->retValUC = 2;
                     return 0;
                 } /*If at the end of the file*/
 
                 *startOfBuffCStr = *buffCStr;
                 continue;
             } /*If need to read in more buffer*/

             /*Adjust start position so next read knows this reads end*/
             ++cntUL;
             ++(*startOfBuffCStr);
        } /*While I have not found a new entry*/
    } /*While I have to find the starting point*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-3: Update file location variables
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    fseek(extParmST->fqFILE, lastHeaderUL, SEEK_SET);
    *lenInUL = 
        fread(
            *buffCStr,
            sizeof(char),
            extParmST->lenBuffUI,
            extParmST->fqFILE
    ); /*re-read in the header*/

    *startOfBuffCStr = *buffCStr;
    *(*buffCStr + *lenInUL) = '\0';
    /*Make sure the previous thread knows were to end*/
    extParmST->retValUC = 1; /*So user knows succeded*/

    return 1;
} /*findStartPos*/
