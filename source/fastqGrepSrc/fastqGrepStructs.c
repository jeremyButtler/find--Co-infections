/*##############################################################################
# Name:
# Use: Holds structers needed to make a tree of read names
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphGraphAndTreeStructs
#    fun-1: makeReadInfoStruct: makes a readInfo struction with default values
#    fun-2: freeReadInfoStruct: frees a readInfo structer
#    fun-3: pushReadNodeStack: pushes a readNodeStack structer onto stack
#    fun-4: popReadNodeStack: frees readNodeStack & returns next node in stack
#    fun-5: makeBigNumStruct: Converts hex elements in c-string to big number
#    fun-6: strToBackwardsBigNum: Filps c-string & converts to big number
#    fun-7: freeBigNumStruct: Frees bigNum structer
#    Fun-8: cmpBigNums: compare two big numbers
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fastqGrepStructs.h"

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    char *readIdCStr,          /*c-string with read name to copy*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
) /*Allocates memomory and makes a readInfo structer (variables set to 0)*/
{ /*makeReadInfoStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: Make a readInfo structer
    #     fun-1 sec-1: Variable declerations
    #     fun-1 sec-2: Copy the input string
    #     fun-1 sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo
        *readInfoStruct = malloc(sizeof(readInfo));

    if(readInfoStruct == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(
            stderr,
            "Failed to allocate memory: Fun-1 makeReadInfoStruct %s",
            "fastqGrepStructs.c Line 40\n"
        ); /*Let user know that memory allocation failed*/
        return 0;
    } /*If mallac failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Get numeric id of intput string
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Convert hex elements in read id to big number*/
    readInfoStruct->idBigNum = makeBigNumStruct(readIdCStr, numElmUChar);

    if(readInfoStruct->idBigNum == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(
            stderr,
            "Failed to allocate memory for big number: Fun-1 %s",
            "makeReadInfoStruct fastqGrepStructs.c Line 57\n"
        ); /*Let user know that memory allocation failed*/
        return 0;
    } /*If mallac failed to allocate memory*/


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Set other variables to 0
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->balanceChar = 0;
    readInfoStruct->leftChild = 0;
    readInfoStruct->rightChild = 0;

    return readInfoStruct;
} /*makeReadInfoStruct*/

/*##############################################################################
# Output: frees readInforStruct and sets pointer to 0
# Note: Does not free nodeInGraph.
##############################################################################*/
void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
) /*frees a readInfo structer*/
{ /*freeReadInfoStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: free readInfoStruct
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    freeBigNumStruct(&((*readInfoStruct)->idBigNum));
    free(*readInfoStruct);               /*User handles nodeInGraph*/
    *readInfoStruct = 0;

    return;
} /*freeReadInfoStruct*/

/*##############################################################################
# Output:
#    Modifes readStack to piont to last element
##############################################################################*/
void pushReadNodeStack(
    struct readNodeStack **readStackAry, /*Array of read info nodes*/
    struct readInfo *readNode            /*readInfo to assing to next node*/
) /*pushes a readNodeStack structer onto a readNodeStack stack*/
{ /*makeReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: makes a readNodeStack structer
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Move to next node*/
    *readStackAry = (*readStackAry) + 1; /*+ sizeof(readNodeStack);*/
    (*readStackAry)->readNode = readNode;

    return;
} /*makeReadNodeStack*/

/*##############################################################################
# Output: Modifies: readNodeStack to point to next readInfo node in stack
##############################################################################*/
void popReadNodeStack(
    struct readNodeStack **readStackAry /*readInfo Array (stack) to pop*/
) /*frees readNodeStack & sets readNodeStack to next readInfo node in stack*/
{ /*popReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 TOC: frees readNodeStack & returns next node in stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    *readStackAry = (*readStackAry) - 1;/* - sizeof(readNodeStack);*/
    return;
} /*popReadNodeStack*/

/*##############################################################################
# Output:
#    returns: bigNum structer with the converted big number
#    Modifies: numElmUChar if a larger number of unsigned longs are needed
#    returns: 0 if memory alloaction failed
##############################################################################*/
struct bigNum * makeBigNumStruct(
    char *cStrToCnvt,          /*C-string to convert hex elements to big number*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
) /*Converts hex characters in c-string to a bitNum struct with a big number*/
{ /*makeBigNumStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 Sub-1 TOC: makeBigNumStruct
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct bigNum *idBigNum = malloc(sizeof(struct bigNum));

    if(idBigNum == 0)
    { /*If memory allocation failed*/
        fprintf(
            stderr,
            "Memory allocation failed: Fun-5 makeBigNumStruct %s",
            "fastqGrepStructs.c Line 153\n"
        ); /*Print error message to user*/

        return 0; 
    } /*If memory allocation failed*/

    idBigNum->lenAllElmChar = 0; /*So function knows to reallocate memory*/
    idBigNum->bigNumAryULng = 0;  /*Just in case realloc needs*/

    /*Convert the c-string to a big number*/
    strToBackwardsBigNum(idBigNum, cStrToCnvt, numElmUChar); 

    if(idBigNum->lenUsedElmChar == 0)
    { /*If memory allocation failed*/
        freeBigNumStruct(&idBigNum); /*Errored out, so no longer need*/

        fprintf(
            stderr,
            "Memory allocation failed in c-strin to number conversion: %s",
            "Fun-5 makeBigNumStruct fastqGrepStructs.c line 170\n"
        ); /*Print error message to user*/

        return 0; 
    } /*If memory allocation failed*/

    return idBigNum;
} /*makeBigNumStruct*/

/*##############################################################################
# Output:
#    Modifies: idBigNum to hold the converted big number.
#        - Sets idBigNum->lenUsedULng to 0 if memory reallocation failed
#    Modifies: numElmUChar if a larger number of unsigned longs are needed
##############################################################################*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,   /*Holds the output big number*/
    char *cStrToCnvt,          /*C-string to convert to large number*/
    unsigned char *numElmUChar /*Number of unsinged longs needed*/
) /*Flips c-string & converts to big number*/
{ /*strToBackwardsBigNum*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 TOC: strToGackwardsBigNum
    #    fun-6 sec-1: Variable declerations
    #    fun-6 sec-2: Resize big number array if needed
    #    fun-6 sec-3: Convert string to big number (work backwards)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: Variable declerations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    unsigned long *elmOnPtrULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Resize big number array if needed
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    idBigNum->lenUsedElmChar = 0; /*Make sure starts at 0 (re-finding)*/

    /*Make sure have an array large enough to store number*/
    if(idBigNum->lenAllElmChar < *numElmUChar)
    { /*If I need to make the array biger*/
        if(idBigNum->bigNumAryULng == 0)
            idBigNum->bigNumAryULng =
                malloc(sizeof(unsigned long) * (*numElmUChar));
        else
            idBigNum->bigNumAryULng =
                realloc(
                    idBigNum->bigNumAryULng,
                    sizeof(unsigned long) * (*numElmUChar)
            ); /*Need to reallocate memory*/
 
        if(idBigNum->bigNumAryULng == 0)
        { /*If memory reallocation failed*/
            idBigNum->lenUsedElmChar = 0;
            fprintf(
                stderr,
               "Memory allocation failed: Fun-6 fastqGrepStructs.c line 225\n"
            ); /*Print error message to user*/

            return;
        } /*If memory reallocation failed*/

        idBigNum->lenAllElmChar = (*numElmUChar);
    } /*If I need to make the array biger*/

    else
        *numElmUChar = idBigNum->lenAllElmChar;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Convert string to big number (work backwards)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(*cStrToCnvt == '@')
        cStrToCnvt++; /*Move off header for fastq entry*/

    while(*cStrToCnvt > 32)
    { /*While there are array elements to fill*/

        /*This allows me to work on strings that have no lengths provided.
          The idea is that fastq id's are generally the same length, so I can
          just keep providing the same numElmUChar again. If they are not I
          will eventaully hit the max length, and have it set
        */
        if(idBigNum->lenUsedElmChar >= idBigNum->lenAllElmChar)
        { /*If need to resize the unsigned long array*/
            (*numElmUChar)++; /*Number of elements needed are bigger*/

            idBigNum->bigNumAryULng =
                realloc(
                    idBigNum->bigNumAryULng,
                    sizeof(unsigned long) * (*numElmUChar)
            ); /*Need to make the unsigned long array biger*/

            idBigNum->lenAllElmChar++; /*More elements added to structer*/
            (*numElmUChar)++;
        } /*If need to resize the unsigned long array*/

        /*Graph unsigned long element working on*/
        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        *elmOnPtrULng = 0;

        for(
            unsigned char charBit = 0;       /*Using ULng to make easiy*/
            charBit < (sizeof(unsigned long) << 3); /*While bits to fill in*/
            charBit += 4                     /*Bits used per hex character*/
        ) { /*For empty bits in the current big number unsigned long element*/
            if(*cStrToCnvt < 33)
                break; /*If have finshed converting the hex string*/

            if(*cStrToCnvt > 47 && *cStrToCnvt < 71) /*0-9 or A-F, (covers 0-15)*/
                *elmOnPtrULng = *elmOnPtrULng + (((*cStrToCnvt)-48) << charBit);

            else if(*cStrToCnvt > 96 && *cStrToCnvt < 103) /*a-f, (covers 10-15)*/
                *elmOnPtrULng = *elmOnPtrULng + (((*cStrToCnvt)-87) << charBit);
            else
                charBit -= 4;   /*make sureo only recored conversion*/

            cStrToCnvt++; /*move to next character in id*/
        } /*For empty bits in the current big number unsigned long element*/

        (idBigNum->lenUsedElmChar)++; /*Track number ULngs acctualy used*/
    } /*While there are array elements to fill*/

    return;
} /*strToBackwardsBigNum*/

/*##############################################################################
# Output:
#    frees: idBigNum & sets pointer to 0
##############################################################################*/
void freeBigNumStruct(
    struct bigNum **idBigNum  /*Address to bigNum structer to free*/
) /*Frees a bigNum struct from memory*/
{ /*freeBigNumStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1 Sub-1 TOC: freeBigNumStruct
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    free((*idBigNum)->bigNumAryULng); /*Should be 0 or on heap*/
    free(*idBigNum);

    *idBigNum = 0;  /*Make sure user can never use again*/

    return;
} /*freeBigNumStruct*/

/*##############################################################################
# Output:
#    Returns: 0 if both equal, < 0 if first is smaller, > 0 if first is bigger
##############################################################################*/
unsigned long cmpBigNums(
    struct bigNum *bigNumOne,
    struct bigNum *bigNumTwo
) /*Compares bigNumOne to bigNumTwo to see if equal, >, or <*/
{ /*cmpBigNums*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1 Sub-1 TOC: cmpBigNums
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*If one stuct has more elements*/
    if(bigNumOne->lenUsedElmChar != bigNumTwo->lenUsedElmChar)
        return bigNumOne->lenUsedElmChar - bigNumTwo->lenUsedElmChar;

    for(char charElm = bigNumOne->lenUsedElmChar - 1; charElm > -1; charElm--)
    { /*For all unsinged longs in the big number*/
        if(
             *(bigNumOne->bigNumAryULng + charElm) !=
             *(bigNumTwo->bigNumAryULng + charElm)
        ) { /*If one number is bigger*/
             return
                 *(bigNumOne->bigNumAryULng + charElm) -
                 *(bigNumTwo->bigNumAryULng + charElm);
        } /*If one number is bigger*/
    } /*For all unsinged longs in the big number*/

    return 0;
} /*cmpBigNums*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fun-9 Sec-1 Sub-1 TOC:
# Output:
#    Returns: unsigned char with the number of unsigned longs needed to hold
#             the hex c-string
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
unsigned char cnvtStrLenToNumHexULng(
    const unsigned long *lenHexCStrULng /*Number of characters in hex c-string*/
) /*Converts the number of unsinged longs needed to store the hex characters*/
{ /*cnvtStrLenToNumHexULng*/
    return 1 + (*lenHexCStrULng >> (sizeof(unsigned long) >> 1));
} /*cnvtStrLenToNumHexULng*/

