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
#    fun-3: extractReadsInTree: Only extract reads in the AVL tree
#    fun-4: extractReadsNotInTree: Only extract reads not in the AVL tree
#    fun-5: extractReadsInHash: Only extract reads in hash table        ADDD
#    fun-6: extractReadsNotInHash: Only extract reads not in hash table ADDD
# Note: I could have written functions 3-6 as a single function, with if 
#       statements. I wrote it as four functions to reduce the number of if
#       statements called in my code. This should give a minor, but not very
#       noticiable performance boost. The logic is I only have to check the
#       user settings once and then launch the correct function.
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fastqGrepSearchFastq.h"

/*##############################################################################
# Output:
#    Stdout: Prints out kept reads
#    Returns: 0 if not a valid fastq file
##############################################################################*/
char fastqExtract(
    FILE *filterFile,   /*FILE object pointing to file with ID's to search for*/
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    unsigned char sizeReadStackUChar, /*Number of elements to use in stack*/
    unsigned long buffSizeULng,       /*Size of buffer to read input with*/
    char hashSearchChar,        /*1: do a hash search, 0: only do AVL Tree*/
    char printReverseChar       /*1: Keep reads in filter file
                                  0: ingore reads in filter file*/
) /*Searches and extracts reads from a fastq file using read id's*/
{ /*fastqExtract*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: Pull out reads of interest
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: Build the tree for reads
    #    fun-1 sec-3: Call the tree or hash table function to search the file
    #    fun-1 sec-4: Handle errors, clean up, & exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        lineInCStr[buffSizeULng + 1],  /*Holds the input line*/
        hashFailedBool = 1;  /*Assume hash tabel failed, makes AVL errors out
                               when fail*/

    unsigned long
        hashSizeULng = 0,       /*Will hold Size of hash table*/
        fastqErrULng = 0;       /*Tells if error in fastq entry*/

    struct readNodeStack
        readStack[sizeReadStackUChar + 2];
           /*Stack to use for searching tree. (tree depth of 73 = 10^18 items)*/

    struct readInfo
        *readTree = 0,
        **hashTbl = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Build the tree or hash table for reads
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make sure start & end of my stacks are marked*/
    readStack[0].readNode = 0;
    readStack[sizeReadStackUChar + 1].readNode = 0;

    if(hashSearchChar == 0)
    { /*If just using the avl tree for searching*/
        readTree =
            buildAvlTree(
                filterFile,          /*File with read ids to keep or ignore*/
                readStack,
                lineInCStr,         /*Buffer to hold one line from file*/
                buffSizeULng   /*Size of buffer to hold each line of the file*/
        ); /*Build the tree of reads to search*/
    } /*If just using the avl tree for searching*/

    else
    { /*Else I am searching using a hash function*/
        hashTbl =
            makeReadHash(
                filterFile,              /*Name of file with read ids to hash*/
                lineInCStr,         /*Buffer to hold one line from file*/
                buffSizeULng,   /*Size of buffer to hold each line of the file*/
                readStack,
                &hashSizeULng,      /*Will hold Size of hash table*/
                &hashFailedBool
        ); /*Build the hash table*/
    } /*Else I am searching using a hash function*/

    fclose(filterFile); /*No longer need open*/

    if(readTree == 0 && hashFailedBool == 1)
    { /*If malloc errored out in making the tree*/
        fprintf(
            stderr,
            "Malloc failed: fastqGrepSearchFastq.c: Fun-1: 108\n"
        ); /*Warn user malloc failed*/

        return 0;
    } /*If malloc errored out in making the tree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Call the tree or hash table function to search the file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(hashSearchChar == 0 && printReverseChar == 0)
    { /*If extracting reads with AVL tree*/
        fastqErrULng =
            extractReadsInTree(
               fastqFile,       /*fastq file to extract reads from*/
               lineInCStr,     /*Buffer to hold input from fastq file*/
               buffSizeULng,   /*Size of lineInCStr*/
               readTree /*root of readInfo tree of read id's to search*/
        ); /*Extract reads using AVL tree*/
    } /*If extracting reads with AVL tree*/

    else if(hashSearchChar == 0 && printReverseChar == 1)
    { /*Else if extracting with AVL tree and keeping reads not in tree*/
        fastqErrULng =
            extractReadsNotInTree(
               fastqFile,       /*fastq file to extract reads from*/
               lineInCStr,     /*Buffer to hold input from fastq file*/
               buffSizeULng,   /*Size of lineInCStr*/
               readTree /*root of readInfo tree of read id's to search*/
        ); /*Extract reads not in AVL tree*/
    } /*Else if extracting with AVL tree and keeping reads not in tree*/

    else if (printReverseChar == 0)
    { /*Else extracting reads in a hash table*/
        fastqErrULng =
            extractReadsInHash(
                fastqFile,
                lineInCStr,   /*Buffer to hold input from fastq file*/
                hashSizeULng, /*holds Size of hash table*/
                buffSizeULng, /*Size of lineInCStr*/
                hashTbl /*Hash table to search for read ids*/
        );
    } /*Else extracting reads in a hash table*/

    else
    { /*Else extracting reads that are missing from a hash table*/
        fastqErrULng =
            extractReadsNotInHash(
                fastqFile,
                lineInCStr,   /*Buffer to hold input from fastq file*/
                hashSizeULng, /*holds Size of hash table*/
                buffSizeULng, /*Size of lineInCStr*/
                hashTbl /*Hash table to search for read ids*/
        );
    } /*Else extracting reads that are missing from a hash table*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Handle errors, clean up, & exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if freeing tree or hash table with tree*/
    if(hashSearchChar == 0)
        freeReadTree(&readTree, readStack);
    else
        freeHashTbl(&hashTbl, hashSizeULng, readStack);

    if(fastqErrULng == 0) /*CHANGE TO IF EQUALS PROBLEM LINE*/
    { /*If this was not a valid fastq file*/
        fclose(fastqFile);
        return 0;  /*main function reports warning*/
    } /*If this was not a valid fastq file*/

    fclose(fastqFile);

    return 1; /*Success*/
} /*fastqExtract*/

/*##############################################################################
# Output:
#    Returns: balanced readInfo tree with all read id's in filterFile
#    Returns: 0 if malloc errored out
##############################################################################*/
struct readInfo * buildAvlTree(
    FILE *filterFile,                 /*File with read ids to keep or ignore*/
    struct readNodeStack *readStack,  /*Stack to use in building AVL tree*/
    char *lineInCStr,                 /*Buffer to hold one line from file*/
    unsigned int buffSizeULng         /*Size of buffer to read each line*/
) /*Builds a readInfo tree with read id's in filterFile*/
{ /*buildAvlTree function*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC:
    #    fun-2 Sec-1: Variable declerations
    #    fun-2 Sec-2: Read id's in filter file and build tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char 
        *readPosCStr = 0;

    unsigned char
        lenNameUInt = 0; /*Length of read id*/

    struct readInfo
        *readTree = 0,
        *lastRead = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Read id's in filter file and build tree
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(
        fgets(
          lineInCStr,
          buffSizeULng,
          filterFile
        ) /*Read a single read name from the file*/
    ) { /*While there are read names to read in*/

        readPosCStr = lineInCStr;
        lenNameUInt = 0;            /*Reset for new read*/

        /*Find the length of the read name*/
        while(
            *readPosCStr != '\0' &&         /*End of lineInCStr*/
            *readPosCStr != ' ' &&          /*End of read name*/
            *readPosCStr != '\n'            /*End of line*/
        ) { /*While still on the read name part of header*/
            readPosCStr++;                  /*Move to next character in id*/
            lenNameUInt++;                  /*Count number of characters in id*/
        } /*While still on the read name part of header*/

        *readPosCStr = '\0';       /*Only want first part, so can ignore rest*/

        lastRead =
            findAddNodeToReadTree(
                lineInCStr,        /*Line is the read name*/
                lenNameUInt,      /*length of readNameCStr*/
                &readTree,         /*read tree to search for read name*/
                readStack          /*Stack to use in searching tree*/
        ); /*Add a read name to my tree of read names*/

        if(lastRead == 0)
        { /*If malloc failed to find memory*/
            fprintf(
                stderr,
                "Malloc failed: fastqGrepSearchFastq.c: Fun-2: 55\n"
            ); /*Warn user malloc failed*/

            return 0;
        } /*If malloc failed to find memory*/
    } /*While there are read names to read in*/

    return readTree;
} /*buildAvlTree function*/

/*##############################################################################
# Output:
#    stdout: Prints extracted reads
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsInTree(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo *readTree  /*root of readInfo tree of read id's to search*/
) /*Extract target reads from fastq file*/
{ /*extractReadsInTree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: Variable declerations
    #    fun-3 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
       *readPosCStr = 0,        /*Position in lineInCStr*/
       *startReadCStr = 0,      /*Start of read id for a fastq entry*/
       endOfFileChar = 0,       /*Marks if at the end of the file*/
       tmpChar = 0;             /*Holds a single character temporarly*/

    int
        lenInputInt = 0;        /*Holds number char fread grabbed from file*/

    struct readInfo
        *lastRead = 0;          /*Holds node of read id found in tree search*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Extract target reads from fastq file
    #    fun-3 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-3 sec-3 sub-2: Determine if read is in tree
    #    fun-3 sec-3 sub-3: Decide if should keep read
    #    fun-3 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    lineInCStr[buffSizeULng] = '\0';          /*Make c-string for fread*/
    readPosCStr = lineInCStr + buffSizeULng; /*So will read in buffer*/
    lenInputInt = buffSizeULng;               /*So does not exit early*/

    while(endOfFileChar != 4)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        endOfFileChar = parseFastqHeader(
                           lineInCStr,     /*Buffer to hold file input*/
                           &startReadCStr, /*Will hold start of read name*/
                           &readPosCStr,   /*Start of read name, wil hold end*/
                           &lenInputInt,
                           buffSizeULng,
                           fastqFile
        ); /*Get read name from file*/

         if(endOfFileChar == 0)
             return 0; /*If fastq file ends on a header*/
         else if(endOfFileChar == 4)
             continue;            /*At end of file, return*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        tmpChar = *readPosCStr;            /*Is a '\n', ' ', or '\t'*/
        *readPosCStr = '\0';               /*So is a c-string (for strcmp)*/

        lastRead =
          searchTree(
            startReadCStr,       /*Line is the read name*/
            readTree             /*read tree to search*/
        ); /*Check if read is one I want to keep*/

        *readPosCStr = tmpChar;            /*Remove c-string*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read*/
        if(lastRead == 0)
        { /*If is a read I am not printing out*/
             endOfFileChar = moveToNextFastqEntry(
                               lineInCStr,
                               &readPosCStr, /*wil hold end of read name*/
                               buffSizeULng,
                               &lenInputInt, /*Holds how many char fread got*/
                               fastqFile
             ); /*Move to next fastq entry*/

             if(endOfFileChar == 0)
                 return 0; /*If there ws */

            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         endOfFileChar = printFastqEntry(
                             lineInCStr,
                             &readPosCStr,   /*Start of read name, wil hold end*/
                             &startReadCStr, /*points to start of read name*/
                             buffSizeULng,
                             &lenInputInt,
                             fastqFile
         ); /*Print read & move to next read*/

         if(endOfFileChar == 0)
             return 0; /*If ended to early*/
    } /*While there are lines in the file*/

    return 1; /*No errors*/
} /*extractReadsInTree*/

/*##############################################################################
# Output:
#    stdout: Prints out reads not in tree
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsNotInTree(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo *readTree  /*root of readInfo tree of read id's to search*/
) /*Extract reads not given as targets from a fastq file*/
{ /*extractReadsNotInTree*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: Variable declerations
    #    fun-3 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
       *readPosCStr = 0,        /*Position in lineInCStr*/
       *startReadCStr = 0,      /*Start of read id for a fastq entry*/
       endOfFileChar = 0,       /*Marks if at the end of the file*/
       tmpChar = 0;             /*Holds a single character temporarly*/

    int
        lenInputInt = 0;        /*Holds number char fread grabbed from file*/

    struct readInfo
        *lastRead = 0;          /*Holds node of read id found in tree search*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Extract target reads from fastq file
    #    fun-4 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-4 sec-3 sub-2: Determine if read is in tree
    #    fun-4 sec-3 sub-3: Decide if should keep read
    #    fun-4 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    lineInCStr[buffSizeULng] = '\0';          /*Make c-string for fread*/
    readPosCStr = lineInCStr + buffSizeULng; /*So will read in buffer*/
    lenInputInt = buffSizeULng;               /*So does not exit early*/

    while(endOfFileChar != 4)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-4 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        endOfFileChar = parseFastqHeader(
                           lineInCStr,     /*Buffer to hold file input*/
                           &startReadCStr, /*Will hold start of read name*/
                           &readPosCStr,   /*Start of read name, wil hold end*/
                           &lenInputInt,
                           buffSizeULng,
                           fastqFile
        ); /*Get read name from file*/

         if(endOfFileChar == 0)
             return 0; /*If fastq file ends on a header*/
         else if(endOfFileChar == 4)
             continue;            /*At end of file, return*/

        /***********************************************************************
        # Fun-4 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        tmpChar = *readPosCStr;            /*Is a '\n', ' ', or '\t'*/
        *readPosCStr = '\0';               /*So is a c-string (for strcmp)*/

        lastRead =
          searchTree(
            startReadCStr,       /*Line is the read name*/
            readTree             /*read tree to search*/
        ); /*Check if read is one I want to keep*/

        *readPosCStr = tmpChar;            /*Remove c-string*/

        /***********************************************************************
        # Fun-4 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read*/
        if(lastRead != 0)
        { /*If is a read I am not printing out*/

             endOfFileChar = moveToNextFastqEntry(
                               lineInCStr,
                               &readPosCStr,   /*Start of read name, wil hold end*/
                               buffSizeULng,
                               &lenInputInt,   /*Holds how many characters fread got*/
                               fastqFile
             ); /*Move to next fastq entry*/

             if(endOfFileChar == 0)
                 return 0; /*If there ws */

            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-4 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         endOfFileChar = printFastqEntry(
                             lineInCStr,
                             &readPosCStr,   /*Start of read name, wil hold end*/
                             &startReadCStr, /*points to start of read name*/
                             buffSizeULng,
                             &lenInputInt,
                             fastqFile
         ); /*Print read & move to next read*/

         if(endOfFileChar == 0)
             return 0; /*If ended to early*/
    } /*While there are lines in the file*/

    return 1; /*No errors*/
} /*extractReadsNotInTree*/

/*##############################################################################
# Output:
#    stdout: Prints out reads in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsInHash(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long hashSizeULng, /*Will hold Size of hash table*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo **hashTbl  /*root of readInfo tree of read id's to search*/
) /*Extract target reads from fastq file with hash table*/
{ /*extractReadsInHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 TOC:
    #    fun-5 sec-1: Variable declerations
    #    fun-5 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
       *readPosCStr = 0,        /*Position in lineInCStr*/
       *startReadCStr = 0,      /*Start of read id for a fastq entry*/
       endOfFileChar = 0,       /*Marks if at the end of the file*/
       tmpChar = 0;             /*Holds a single character temporarly*/

    int
        lenInputInt = 0;        /*Holds number char fread grabbed from file*/

    struct readInfo
        *lastRead = 0;          /*Holds node of read id found in tree search*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Extract target reads from fastq file
    #    fun-5 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-5 sec-3 sub-2: Determine if read is in tree
    #    fun-5 sec-3 sub-3: Decide if should keep read
    #    fun-5 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    lineInCStr[buffSizeULng] = '\0';          /*Make c-string for fread*/
    readPosCStr = lineInCStr + buffSizeULng; /*So will read in buffer*/
    lenInputInt = buffSizeULng;               /*So does not exit early*/

    while(endOfFileChar != 4)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-5 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        endOfFileChar = parseFastqHeader(
                           lineInCStr,     /*Buffer to hold file input*/
                           &startReadCStr, /*Will hold start of read name*/
                           &readPosCStr,   /*Start of read name, wil hold end*/
                           &lenInputInt,
                           buffSizeULng,
                           fastqFile
        ); /*Get read name from file*/

         if(endOfFileChar == 0)
             return 0; /*If fastq file ends on a header*/
         else if(endOfFileChar == 4)
             continue;            /*At end of file, return*/

        /***********************************************************************
        # Fun-5 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        tmpChar = *readPosCStr;            /*Is a '\n', ' ', or '\t'*/
        *readPosCStr = '\0';               /*So is a c-string (for strcmp)*/

        lastRead = 
            findReadInHashTbl(
                startReadCStr, /*Read to search for*/
                hashTbl,     /*Hash table to search for read in*/
                hashSizeULng  /*Size of array holding hash*/
        );  /*See if read id is in the hash table*/

        *readPosCStr = tmpChar;            /*Remove c-string*/

        /***********************************************************************
        # Fun-5 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read*/
        if(lastRead == 0)
        { /*If is a read I am not printing out*/

             endOfFileChar = moveToNextFastqEntry(
                               lineInCStr,
                               &readPosCStr,   /*Start of read name, wil hold end*/
                               buffSizeULng,
                               &lenInputInt,   /*Holds how many characters fread got*/
                               fastqFile
             ); /*Move to next fastq entry*/

             if(endOfFileChar == 0)
                 return 0; /*If there ws */

            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-5 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         endOfFileChar = printFastqEntry(
                             lineInCStr,
                             &readPosCStr,   /*Start of read name, wil hold end*/
                             &startReadCStr, /*points to start of read name*/
                             buffSizeULng,
                             &lenInputInt,
                             fastqFile
         ); /*Print read & move to next read*/

         if(endOfFileChar == 0)
             return 0; /*If ended to early*/
    } /*While there are lines in the file*/

    return 1; /*No errors*/
} /*extractReadsInHash*/

/*##############################################################################
# Output:
#    stdout: Prints out reads not in hash table
#    Returns: 0 if was not a valid fastq file
##############################################################################*/
char extractReadsNotInHash(
    FILE *fastqFile,    /*FILE object pointing to fastq file to search through*/
    char *lineInCStr,   /*Buffer to hold input from fastq file*/
    unsigned long hashSizeULng, /*Will hold Size of hash table*/
    unsigned long buffSizeULng, /*Size of lineInCStr*/
    struct readInfo **hashTbl  /*root of readInfo tree of read id's to search*/
) /*Extract non-target reads from fastq file using a hash table*/
{ /*extractReadsNotInHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 TOC:
    #    fun-6 sec-1: Variable declerations
    #    fun-6 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
       *readPosCStr = 0,        /*Position in lineInCStr*/
       *startReadCStr = 0,      /*Start of read id for a fastq entry*/
       endOfFileChar = 0,       /*Marks if at the end of the file*/
       tmpChar = 0;             /*Holds a single character temporarly*/

    int
        lenInputInt = 0;        /*Holds number char fread grabbed from file*/

    struct readInfo
        *lastRead = 0;          /*Holds node of read id found in tree search*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Extract target reads from fastq file
    #    fun-6 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-6 sec-3 sub-2: Determine if read is in tree
    #    fun-6 sec-3 sub-3: Decide if should keep read
    #    fun-6 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    lineInCStr[buffSizeULng] = '\0';          /*Make c-string for fread*/
    readPosCStr = lineInCStr + buffSizeULng; /*So will read in buffer*/
    lenInputInt = buffSizeULng;               /*So does not exit early*/

    while(endOfFileChar != 4)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-6 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        endOfFileChar = parseFastqHeader(
                           lineInCStr,     /*Buffer to hold file input*/
                           &startReadCStr, /*Will hold start of read name*/
                           &readPosCStr,   /*Start of read name, wil hold end*/
                           &lenInputInt,
                           buffSizeULng,
                           fastqFile
        ); /*Get read name from file*/

         if(endOfFileChar == 0)
             return 0; /*If fastq file ends on a header*/
         else if(endOfFileChar == 4)
             continue;            /*At end of file, return*/

        /***********************************************************************
        # Fun-6 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        tmpChar = *readPosCStr;            /*Is a '\n', ' ', or '\t'*/
        *readPosCStr = '\0';               /*So is a c-string (for strcmp)*/

        lastRead = 
            findReadInHashTbl(
                startReadCStr, /*Read to search for*/
                hashTbl,     /*Hash table to search for read in*/
                hashSizeULng  /*Size of array holding hash*/
        );  /*See if read id is in the hash table*/

        *readPosCStr = tmpChar;            /*Remove c-string*/

        /***********************************************************************
        # Fun-6 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read*/
        if(lastRead != 0)
        { /*If is a read I am not printing out*/

             endOfFileChar = moveToNextFastqEntry(
                               lineInCStr,
                               &readPosCStr,   /*Start of read name, wil hold end*/
                               buffSizeULng,
                               &lenInputInt,   /*Holds how many characters fread got*/
                               fastqFile
             ); /*Move to next fastq entry*/

             if(endOfFileChar == 0)
                 return 0; /*If there ws */

            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-6 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         endOfFileChar = printFastqEntry(
                             lineInCStr,
                             &readPosCStr,   /*Start of read name, wil hold end*/
                             &startReadCStr, /*points to start of read name*/
                             buffSizeULng,
                             &lenInputInt,
                             fastqFile
         ); /*Print read & move to next read*/

         if(endOfFileChar == 0)
             return 0; /*If ended to early*/
    } /*While there are lines in the file*/

    return 1; /*No errors*/
} /*extractReadsNotInHash*/
