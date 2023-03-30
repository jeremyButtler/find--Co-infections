/*##############################################################################
# Name:
# Use: Holds structers needed to make a tree of read names
##############################################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: fqGetsIdsStructs
'    fun-1 makeReadInfoStruct:
'      o Makes a readInfo struction with default values
'    fun-2 freeReadInfoStruct:
'      o Frees a readInfo structer
'    fun-3 pushReadNodeStack:
'      o Pushes a readNodeStack structer onto stack
'    fun-4 popReadNodeStack:
'      o Frees readNodeStack & returns next node in stack
'    fun-5 makeBigNumStruct:
'      o Converts hex elements in c-string to big number
'    fun-6 strToBackwardsBigNum:
'      o Filps c-string & converts to big number
'    fun-7 freeBigNumStruct:
'      o Frees bigNum structer
'    fun-8 cmpBigNums:
'      o Compare two big numbers
'    fun-9 cnvtIdToBigNum:
'      o Read in read id line & convert to big num
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "fqGetIdsStructs.h"

/*Make look up table to look if character is valid hex character
    64 is invisivle character
    32 is non-hex character (printable)
*/
char hexTblCharAry[] =
    {
     /*0-32 (invisible)*/
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
     64,

     // Special character : (need so I can distinguish Illumina reads)
     17,

     // Special characters : ; < = > ? @ (58 to 64)
     53, 54, 55, 56, 57, 58,

     /*Numbers: 48 (0) to 57 (9)*/
     16, 1, 2, 3, 4, 5, 6, 7, 8, 9, /*48-57 Numbers*/

     /*Special characters : ; < = > ? @ (58 to 64)*/
     52, 53, 54, 55, 56, 57, 58,

     /*Hext A-F (65-70)*/
     10, 11, 12, 13, 14, 15,

     /*5 bit G-V (71 to 85) (my liimit for a 5bit with ':' support*/
     18,19,20,21,22,23,24,25,26,27,28,29,30,31,

     /*W,X,Y,Z (86 to 90)*/
     32,32,33,34,35,36,

     /*specialcharacters [ \ ] ^ _ ` (91 to 96)*/
     59,60,61,62,63,49,

     /*a-f (97 to 102)*/
     10, 11, 12, 13, 14, 15,

     /*5 bit g-v (71 to 85) (my liimit for a 5bit with ':' support*/
     18,19,20,21,22,23,24,25,26,27,28,29,30,31,

     /*u, w,x,y,z (86 to 90)*/
     32,32,33,34,35,36,

     /*special characters { | } ~ (123 to 126)*/
     32, 32, 32, 32,

     /*127 del*/
     32,

     /*unsigned extended range (just in case)*/
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32
    };

/*---------------------------------------------------------------------\
| Output: Modifies: readInfoStruct to have default values (all 0's)
\---------------------------------------------------------------------*/
struct readInfo * makeReadInfoStruct(
    char *readIdCStr,         /*c-string with read name to copy*/
    const int32_t *lenCStrUInt /*Length of cString to convert*/
) /*Allocates memomory & makes a readInfo struct (variables set to 0)*/
{ /*makeReadInfoStruct*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Make a readInfo structer
    '     fun-1 sec-1: Variable declerations
    '     fun-1 sec-2: Copy the input string
    '     fun-1 sec-3: Set other variables to 0
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Get numeric id of intput string
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Convert hex elements in read id to big number*/
    readInfoStruct->idBigNum = makeBigNumStruct(readIdCStr,lenCStrUInt);

    if(readInfoStruct->idBigNum == 0)
    { /*If mallac failed to allocate memory*/
        fprintf(
            stderr,
            "Failed to allocate memory for big number: Fun-1 %s",
            "makeReadInfoStruct fastqGrepStructs.c Line 57\n"
        ); /*Let user know that memory allocation failed*/
        return 0;
    } /*If mallac failed to allocate memory*/


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Set other variables to 0
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readInfoStruct->balanceChar = 0;
    readInfoStruct->leftChild = 0;
    readInfoStruct->rightChild = 0;

    return readInfoStruct;
} /*makeReadInfoStruct*/

/*---------------------------------------------------------------------\
| Output: frees readInforStruct and sets pointer to 0
| Note: Does not free nodeInGraph.
\---------------------------------------------------------------------*/
void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
) /*frees a readInfo structer*/
{ /*freeReadInfoStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1 TOC: free readInfoStruct
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    freeBigNumStruct(&((*readInfoStruct)->idBigNum));
    free(*readInfoStruct);               /*User handles nodeInGraph*/
    *readInfoStruct = 0;

    return;
} /*freeReadInfoStruct*/

/*---------------------------------------------------------------------\
| Output:
|    Modifes readStack to piont to last element
\---------------------------------------------------------------------*/
void pushReadNodeStack(
    struct readNodeStack **readStackAry, /*Array of read info nodes*/
    struct readInfo *readNode            /*readInfo to assing to next node*/
) /*pushes a readNodeStack structer onto a readNodeStack stack*/
{ /*makeReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-3 Sec-1 TOC: makes a readNodeStack structer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Move to next node*/
    *readStackAry = (*readStackAry) + 1; /*+ sizeof(readNodeStack);*/
    (*readStackAry)->readNode = readNode;

    return;
} /*makeReadNodeStack*/

/*---------------------------------------------------------------------\
| Output: Modifies: readNodeStack to point to next node in stack
\---------------------------------------------------------------------*/
void popReadNodeStack(
    struct readNodeStack **readStackAry /*readInfo Array (stack) to pop*/
) /*Sets readNodeStack to next readInfo node in stack*/
{ /*popReadNodeStack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-4 Sec-1 TOC: frees readNodeStack & returns next node in stack
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    *readStackAry = (*readStackAry) - 1;/* - sizeof(readNodeStack);*/
    return;
} /*popReadNodeStack*/

/*---------------------------------------------------------------------\
| Output:
|    returns: bigNum structer with the converted big number
|    returns: 0 if memory alloaction failed
\---------------------------------------------------------------------*/
struct bigNum * makeBigNumStruct(
    char *cnvtCStr,/*C-string to convert hex elements to big number*/
    const int32_t *lenCStrUInt /*Length of cString to convert*/
) /*Converts hex characters in c-string to a bigNum struct*/
{ /*makeBigNumStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-5 Sec-1 Sub-1 TOC: makeBigNumStruct
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct bigNum *idBigNum = malloc(sizeof(struct bigNum));

    if(idBigNum == 0) /*Memory allocation error*/
        return 0; 

    idBigNum->lenAllElmChar = 0;  /*tell to reallocate memory*/
    idBigNum->bigNumAryIOrL = 0;  /*Just in case realloc needs*/

    /*Convert the c-string to a big number*/
    strToBackwardsBigNum(idBigNum, cnvtCStr, lenCStrUInt); 

    if(idBigNum->lenUsedElmChar == 0)
    { /*If memory allocation failed*/
        freeBigNumStruct(&idBigNum); /*Errored out, so no longer need*/
        return 0; 
    } /*If memory allocation failed*/

    return idBigNum;
} /*makeBigNumStruct*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: idBigNum to hold the converted big number.
|        - Sets idBigNum->lenUsedULng to 0 if memory reallocation failed
\---------------------------------------------------------------------*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,         /*Holds the output big number*/
    char *cnvtCStr,       /*C-string to convert to large number*/
    const int32_t *lenCStrUInt
) /*Flips c-string & converts to big number*/
{ /*strToBackwardsBigNum*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-6 TOC: strToGackwardsBigNum
    '    fun-6 sec-1: Variable declerations
    '    fun-6 sec-2: Resize big number array if needed
    '    fun-6 sec-3: Convert string to big number (work backwards)
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char charBit = 0; /*Number bits to shift*/

    /*compiler settings for memory (MEM) or speed compile*/
    #ifndef MEM
        #if defOSBit == 64
            int *elmILPtr = 0;
        #else
            short *elmILPtr = 0;
        #endif
    #else
        long *elmILPtr = 0;
    #endif

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-2: Resize big number array if needed
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Get the number of limbs needed to store the number*/ 
    idBigNum->lenUsedElmChar = 1 + (*lenCStrUInt / defBitsPerLimb);
       /*+ 1 to handle truncation by division*/

    /*Make sure have an array large enough to store number*/
    if(idBigNum->lenAllElmChar < idBigNum->lenUsedElmChar)
    { /*If I need to make the array biger*/
        if(idBigNum->bigNumAryIOrL == 0)
            idBigNum->bigNumAryIOrL =
                #ifndef MEM
                     #if defOSBit == 64
                        malloc(sizeof(int) *idBigNum->lenUsedElmChar);
                    #else
                        malloc(sizeof(short) *idBigNum->lenUsedElmChar);
                    #endif
                #else
                    malloc(sizeof(long) *idBigNum->lenUsedElmChar);
                #endif
        else
            idBigNum->bigNumAryIOrL =
                realloc(
                    idBigNum->bigNumAryIOrL,
                    #ifndef MEM
                         #if defOSBit == 64
                            sizeof(int) * idBigNum->lenUsedElmChar
                        #else
                            sizeof(short) * idBigNum->lenUsedElmChar
                        #endif
                    #else
                        sizeof(long) * idBigNum->lenUsedElmChar
                    #endif
            ); /*Need to reallocate memory*/

 
        if(idBigNum->bigNumAryIOrL == 0)
        { /*If memory reallocation failed*/
            idBigNum->lenUsedElmChar = 0;
            return;
        } /*If memory reallocation failed*/

        idBigNum->lenAllElmChar = idBigNum->lenUsedElmChar;
    } /*If I need to make the array biger*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-6 Sec-3: Convert string to big number (work backwards)
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    idBigNum->lenUsedElmChar = 0;

    #ifndef MEM
        idBigNum->totalL = 0;
    #endif

    do { /*While there are array elements to fill*/
        /*Graph unsigned long element working on*/
        elmILPtr =idBigNum->bigNumAryIOrL+ idBigNum->lenUsedElmChar;
        (idBigNum->lenUsedElmChar)++;/*Track number ULngs used*/
        *elmILPtr = 0;
        charBit = 0;

        while(charBit < defMaxDigPerLimb)
        { /*while their are empty bits in the current big number limb*/
            if(hexTblCharAry[(unsigned char) *cnvtCStr] & 64)
                break; /*If have finshed converting the hex string*/

            if(!(hexTblCharAry[(unsigned char) *cnvtCStr] & 32))
            { /*If is a hex character*/
                *elmILPtr = *elmILPtr << defBitsPerChar;
                *elmILPtr+=hexTblCharAry[(unsigned char) *cnvtCStr];
                ++charBit;
            } /*If is a hex character*/

            ++cnvtCStr; /*move to next character in id*/
        } /*while their are empty bits in the current big number limb*/

        #ifndef MEM
            /*add up all limbs for faster comparisons with Illumina*/
            idBigNum->totalL += *elmILPtr;
        #endif
    } while((unsigned char) *cnvtCStr > 32);
    /*While still on the read name part of header*/

    return;
} /*strToBackwardsBigNum*/

/*---------------------------------------------------------------------\
| Output:
|    frees: idBigNum & sets pointer to 0
\---------------------------------------------------------------------*/
void freeBigNumStruct(
    struct bigNum **idBigNum  /*Address to bigNum structer to free*/
) /*Frees a bigNum struct from memory*/
{ /*freeBigNumStruct*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-7 Sec-1 Sub-1 TOC: freeBigNumStruct
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    free((*idBigNum)->bigNumAryIOrL); /*Should be 0 or on heap*/
    free(*idBigNum);
    *idBigNum = 0;  /*Make sure user can never use again*/

    return;
} /*freeBigNumStruct*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|      - 0 if both equal,
|      - < 0 if query is smaller than subject
|      - > 0 if qeury is bigger than subject
\---------------------------------------------------------------------*/
long cmpBigNums(
    struct bigNum *bigNumQ, /*Query to compare against*/
    struct bigNum *bigNumS  /*Subject to compare to*/
) /*Compares bigNumQ to bigNumS to see if equal, >, or <*/
{ /*cmpBigNums*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-8 Sec-1 Sub-1 TOC: cmpBigNums
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*A speed setting*/ 
    #ifndef MEM
        if(bigNumQ->totalL != bigNumS->totalL)
            return bigNumQ->totalL - bigNumS->totalL;
    #endif

    /*If one stuct has more elements*/
    if(bigNumQ->lenUsedElmChar != bigNumS->lenUsedElmChar)
        return bigNumQ->lenUsedElmChar - bigNumS->lenUsedElmChar;

    for(short sElm = bigNumQ->lenUsedElmChar - 1; sElm > -1; --sElm)
    { /*For all unsinged longs in the big number*/
        if(*(bigNumQ->bigNumAryIOrL +sElm)
           != *(bigNumS->bigNumAryIOrL +sElm)
        ) return
             *(bigNumQ->bigNumAryIOrL + sElm) -
             *(bigNumS->bigNumAryIOrL + sElm);
    } /*For all unsinged longs in the big number*/

    return 0;
} /*cmpBigNums*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|      - EndNameCStr to point to '\n', ' ', & '\t' at end of read name
|      - LenIdULng to hold the length of the read id
|      - LenInputInt to hold length of input buffer
|    Returns:
|      - 0 if fails or end of file (lenIdULng < buffSizeInt)
|      - pointer to struct with bigNum struct having converted read id
\---------------------------------------------------------------------*/
struct readInfo * cnvtIdToBigNum(
    char *bufferCStr,  /*buffer to hold fread input (can have data)*/
    uint32_t buffSizeInt,   /*Size of buffer to work on*/
    char **endCStr, /*Points to start of id, will point to end*/
    uint64_t *lenInputULng, /*Length of input from fread*/
    unsigned char *lenBigNumChar, /*Holds size to make bigNumber*/
    FILE *idFILE         /*Fastq file to get data from*/
) /*Converts read id to bigNum read id, will grab new file input*/
{ /*cnvtIdToBigNum*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-9 cnvtIdToBigNum TOC:
    '    fun-9 sec-1: Variable declerations
    '    fun-1 sec-2: Initalize readInfo & bigNum structs
    '    fun-1 sec-3: Move to first character on header
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned char charBit = 0;
    struct bigNum *idBigNum = malloc(sizeof(struct bigNum));
    struct readInfo *readNode = malloc(sizeof(struct readInfo));

    #ifndef MEM
        #if defOSBit == 64
            int *elmILPtr = 0;
        #else
            short *elmILPtr = 0;
        #endif
    #else
        long *elmILPtr = 0;
    #endif

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-2: Initalize readInfo & bigNum structs
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(idBigNum == 0 || readNode == 0)
    { /*If memory allocation failed*/
        if(idBigNum != 0) free(idBigNum);
        if(readNode != 0) free(readNode);

        *lenInputULng = 0; /*Make sure user detects failure*/
        return 0; 
    } /*If memory allocation failed*/

    readNode->balanceChar = 0;
    readNode->leftChild = 0;
    readNode->rightChild = 0;
    readNode->idBigNum = idBigNum;

    idBigNum->lenUsedElmChar = 0;

    #ifndef MEM
        idBigNum->totalL = 0;

         #if defOSBit == 64
           idBigNum->bigNumAryIOrL=malloc(sizeof(int)*(*lenBigNumChar));
        #else
         idBigNum->bigNumAryIOrL=malloc(sizeof(short)*(*lenBigNumChar));
        #endif
    #else
        idBigNum->bigNumAryIOrL = malloc(sizeof(long)*(*lenBigNumChar));
    #endif

    idBigNum->lenAllElmChar = *lenBigNumChar;

    if(idBigNum->bigNumAryIOrL == 0)
    { /*If memory reallocation failed*/
        *lenInputULng = 0; /*Make sure user detects failure*/
        free(idBigNum);
        free(readNode);
        return 0;
    } /*If memory reallocation failed*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-3: Move to first character on header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(**endCStr == '\n') (*endCStr)++;

    if(**endCStr == '\0')
    { /*If at the end of the buffer, but not at start of read*/

        if(*lenInputULng < buffSizeInt)
        { /*If at end of file*/
          free(idBigNum->bigNumAryIOrL);
          free(idBigNum);
          free(readNode);
          return 0;                    /*Done with file*/
        } /*If at end of file*/

        *lenInputULng = fread(
                           bufferCStr,
                           sizeof(char),
                           buffSizeInt,
                           idFILE
        ); /*Read in more of the file*/

      *endCStr = bufferCStr;
    } /*If at the end of the buffer, but not at start of read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-9 Sec-4: Convert string to big number
    ^    fun-9 sec-4 sub-1: Check if need to resize the long array
    ^    fun-9 sec-4 sub-2: Convert characters until long is full
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-9 Sec-4 Sub-1: Check if need to resize the long array
    \******************************************************************/
      
    do { /*While still on the read name part of header*/

        if(idBigNum->lenUsedElmChar >= idBigNum->lenAllElmChar)
        { /*If need to reallocate memory*/
            (idBigNum->lenAllElmChar)++;
            (*lenBigNumChar)++;
            idBigNum->bigNumAryIOrL =
                realloc(
                    idBigNum->bigNumAryIOrL,
                    #ifndef MEM
                         #if defOSBit == 64
                            sizeof(int) * idBigNum->lenAllElmChar
                        #else
                            sizeof(short) * idBigNum->lenAllElmChar
                        #endif /*Check if 32 bit or 64 bit*/
                    #else
                        sizeof(long) * idBigNum->lenAllElmChar
                    #endif
            ); /*Rellocate memory for the array*/

            if(idBigNum == 0 || readNode == 0)
            { /*If memory allocation failed*/
                free(idBigNum->bigNumAryIOrL);
                free(idBigNum);
                free(readNode);
                *lenInputULng = 0; /*Make sure user detects failure*/
                return 0; 
            } /*If memory allocation failed*/
        } /*If need to reallocate memory*/

        /*Graph unsigned long element working on*/
        elmILPtr =idBigNum->bigNumAryIOrL +idBigNum->lenUsedElmChar;
        (idBigNum->lenUsedElmChar)++; /*Track number of used longs*/
        *elmILPtr = 0;
        charBit = 0;

       /***************************************************************\
       * Fun-9 Sec-4 Sub-1: Convert characters until long is full
       \***************************************************************/

        while(charBit < defMaxDigPerLimb)
        { /*while their are empty bits in the current big number limb*/
            if(hexTblCharAry[(unsigned char) **endCStr] & 64)
                break; /*If finshed converting the c-string*/

            if(!(hexTblCharAry[(unsigned char) **endCStr] & 32))
            { /*If was a hex character*/
                *elmILPtr = *elmILPtr << defBitsPerChar;
                *elmILPtr += hexTblCharAry[(unsigned char) **endCStr];
                ++charBit;
            } /*If was a hex character*/

            ++(*endCStr);

            if(**endCStr == '\0')
            { /*If ran out of buffer & need to read in more of file*/
                if(*lenInputULng < buffSizeInt)
                    break; /*at end of file*/

                *lenInputULng = fread(bufferCStr,
                                     sizeof(char),
                                     buffSizeInt,
                                     idFILE
                ); /*Read in more of the file*/

                *(bufferCStr + *lenInputULng) = '\0';/*make a c-string*/
                *endCStr = bufferCStr;
                continue;
            } /*If ran out of buffer & need to read in more of file*/
        } /*while their are empty bits in the current big number limb*/

        #ifndef MEM
            /*add up all limbs for faster comparisons with Illumina*/
            idBigNum->totalL += *elmILPtr;
        #endif
    } while((unsigned char) **endCStr > 32);
    /*While still on the read name part of header*/

    return readNode; /*Copied name sucessfully*/
} /*cnvtIdToBigNum*/

/* Ascii table
Dec  Char                           Dec  Char     Dec  Char     Dec Char
---------                           ---------     ---------     --------
  0  NUL (null)                      32  SPACE     64  @         96  `
  1  SOH (start of heading)          33  !         65  A         97  a
  2  STX (start of text)             34  "         66  B         98  b
  3  ETX (end of text)               35  #         67  C         99  c
  4  EOT (end of transmission)       36  $         68  D        100  d
  5  ENQ (enquiry)                   37  %         69  E        101  e
  6  ACK (acknowledge)               38  &         70  F        102  f
  7  BEL (bell)                      39  '         71  G        103  g
  8  BS  (backspace)                 40  (         72  H        104  h
  9  TAB (horizontal tab)            41  )         73  I        105  i
 10  LF  (NL line feed, new line)    42  *         74  J        106  j
 11  VT  (vertical tab)              43  +         75  K        107  k
 12  FF  (NP form feed, new page)    44  ,         76  L        108  l
 13  CR  (carriage return)           45  -         77  M        109  m
 14  SO  (shift out)                 46  .         78  N        110  n
 15  SI  (shift in)                  47  /         79  O        111  o
 16  DLE (data link escape)          48  0         80  P        112  p
 17  DC1 (device control 1)          49  1         81  Q        113  q
 18  DC2 (device control 2)          50  2         82  R        114  r
 19  DC3 (device control 3)          51  3         83  S        115  s
 20  DC4 (device control 4)          52  4         84  T        116  t
 21  NAK (negative acknowledge)      53  5         85  U        117  u
 22  SYN (synchronous idle)          54  6         86  V        118  v
 23  ETB (end of trans. block)       55  7         87  W        119  w
 24  CAN (cancel)                    56  8         88  X        120  x
 25  EM  (end of medium)             57  9         89  Y        121  y
 26  SUB (substitute)                58  :         90  Z        122  z
 27  ESC (escape)                    59  ;         91  [        123  {
 28  FS  (file separator)            60  <         92  \        124  |
 29  GS  (group separator)           61  =         93  ]        125  }
 30  RS  (record separator)          62  >         94  ^        126  ~
 31  US  (unit separator)            63  ?         95  _        127  DEL
*/

