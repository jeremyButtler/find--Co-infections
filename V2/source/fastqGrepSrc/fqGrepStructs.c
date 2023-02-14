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

#include "fqGrepStructs.h"

/*Make look up table to look if character is valid hex character
    64 is invisivle character
    32 is non-hex character (printable)
*/
char hexTblCharAry[] =
    {
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,/*0-32 (invisible)*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, /*non-#*/

     0, 1, 2, 3, 4, 5, 6, 7, 8, 9, /*48-57 Numbers*/

     32, 32, 32, 32, 32, 32, 32, /*Between numbers & uppcase letters*/

     10, 11, 12, 13, 14, 15,     /*A-F (65-70)*/

     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,/*G-Z(71-96)*/
     32,32,32,32,32,32,/*specialcharacters91to96*/

     10, 11, 12, 13, 14, 15,     /*a-f (97 to 102)*/

     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32, /*g-z 98-122*/
     32, 32, 32, 32, 32, /*Final special characters (123 to 127)*/
     /*unsigned extended range (just in case)*/
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,
     32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32
    };

/*##############################################################################
# Output: Modifies: readInfoStruct to have default values (all 0's)
##############################################################################*/
struct readInfo * makeReadInfoStruct(
    unsigned char *readIdCStr,         /*c-string with read name to copy*/
    const unsigned long *lenCStrULng /*Length of cString to convert*/
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
    readInfoStruct->idBigNum = makeBigNumStruct(readIdCStr, lenCStrULng);

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
#    returns: 0 if memory alloaction failed
##############################################################################*/
struct bigNum * makeBigNumStruct(
    unsigned char *cStrToCnvt,/*C-string to convert hex elements to big number*/
    const unsigned long *lenCStrULng /*Length of cString to convert*/
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
    strToBackwardsBigNum(idBigNum, cStrToCnvt, lenCStrULng); 

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
##############################################################################*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,         /*Holds the output big number*/
    unsigned char *cStrToCnvt,       /*C-string to convert to large number*/
    const unsigned long *lenCStrULng
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

    char charBit = 0; /*Number bits to shift*/

    unsigned long *elmOnPtrULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Resize big number array if needed
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    idBigNum->lenUsedElmChar =
        1 + (*lenCStrULng >> (sizeof(unsigned long) >> 1));
        /* Each hex number takes up 4 bits
           sizeof returns number of bytes so multipy by 8 (<< 3)
           (sizeof(unsigned long) << 3) >> 1 = sizeof(unsigned long) >> 1)
        */

    /*Make sure have an array large enough to store number*/
    if(idBigNum->lenAllElmChar < idBigNum->lenUsedElmChar)
    { /*If I need to make the array biger*/
        if(idBigNum->bigNumAryULng == 0)
            idBigNum->bigNumAryULng =
                malloc(sizeof(unsigned long) * idBigNum->lenUsedElmChar);
        else
            idBigNum->bigNumAryULng =
                realloc(
                    idBigNum->bigNumAryULng,
                    sizeof(unsigned long) * idBigNum->lenUsedElmChar
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

        idBigNum->lenAllElmChar = idBigNum->lenUsedElmChar;
    } /*If I need to make the array biger*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Convert string to big number (work backwards)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    idBigNum->lenUsedElmChar = 0;

    while(*cStrToCnvt > 32)
    { /*While there are array elements to fill*/

        /*Graph unsigned long element working on*/
        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        *elmOnPtrULng = 0;
        charBit = 0;

        while(charBit < (sizeof(unsigned long) << 3))
        { /*while empty bits in the current big number unsigned long element*/
            if(hexTblCharAry[*cStrToCnvt] & 64)
                break; /*If have finshed converting the hex string*/

            if(!(hexTblCharAry[*cStrToCnvt] & 32))
            { /*If is a hex character*/
                (*elmOnPtrULng) += (hexTblCharAry[*cStrToCnvt] << charBit);
                charBit += 4;
            } /*If is a hex character*/

            cStrToCnvt++; /*move to next character in id*/
        } /*while empty bits in the current big number unsigned long element*/

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

    for(short shtElm = bigNumOne->lenUsedElmChar - 1; shtElm > -1; shtElm--)
    { /*For all unsinged longs in the big number*/
        if(
             *(bigNumOne->bigNumAryULng + shtElm) !=
             *(bigNumTwo->bigNumAryULng + shtElm)
        ) { /*If one number is bigger*/

            if(
                *(bigNumOne->bigNumAryULng + shtElm) >
                *(bigNumTwo->bigNumAryULng + shtElm)
            )
                return 1;
            else
                return -1;
             /*ISSUE WITH ULNG NOT GOING NEGATIVE*/
             /*return
                 *(bigNumOne->bigNumAryULng + shtElm) -
                 *(bigNumTwo->bigNumAryULng + shtElm);*/
        } /*If one number is bigger*/
    } /*For all unsinged longs in the big number*/

    return 0;
} /*cmpBigNums*/

/* Ascii table
Dec  Char                           Dec  Char     Dec  Char     Dec  Char
---------                           ---------     ---------     ----------
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

