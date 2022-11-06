/*##############################################################################
# Name: fastqGrepHash.c
#   Use:
#     Uses an hash + AVL tree to find if reads is in user suplied list
##############################################################################*/

#include "fastqGrepHash.h"

char /*First 392 digitits (defing as array so can modify)*/
   GOLDEN_RATIO_CSTR[] = "618033988749894848204586834365638117720309179805762862135448622705260462818902449707207204189391137484754088075386891752126633862223536931793180060766726354433389086595939582905638322661319928290267880675208766892501711696207032221043216269548626296313614438149758701220340805887954454749246185695364864449241044320771344947049565846788509874339442212544877066478091588460749988712400765";

double
   powTwoPerTen = 3.3333333333333333333333333333333333333333333;
   /*there are 2^(3.33) per 10^1. This is so I can avoid gmp right shift,
     which seems to error out on realloc for O2*/

/*Make look up table to look if character is valid hex character
    64 is invisivle character
    32 is non-hex character (printable)
*/
char hexCharAry[] =
    {
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
     64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,/*0-32 (invisible)*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, /*non-#*/

     0, 1, 2, 3, 4, 5, 6, 7, 8, 9, /*48-57 Numbers*/

     32, 32, 32, 32, 32, 32, 32, /*Between numbers & uppcase letters*/

     10, 11, 12, 13, 14, 15,     /*A-F*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32,  /*G-` (71 to 96)*/

     10, 11, 12, 13, 14, 15,     /*a-f (97 to 102)*/

     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32, 32, 32, 32, 32, 32, 32, 32 /*g to ascii limit (103 to 256)*/
    };

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepHash TOC:
#    fun-3: Insert read into hash table
#    fun-4: findReadInHashTbl: Search hash table for a particler read id
#    fun-5: freeHashTbl: free a hash table
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Output:
#    Returns: pointer to array of readInfo nodes that is the hash
#             0: if failed
#    Modifies: hashSizeUlng to hold hte length of the returned array
#    Modifies: numCharChar to hold the number of characters used in hash
#    Modifies: digPerKeyUChar to hold the hash size (as multiple of two)
#    Modifies: magickNumULng to hold the magick number for this hash table
##############################################################################*/
struct readInfo ** makeReadHash(
    FILE * filtFILE,                  /*file with read id's to filter by*/
    char * buffCStr,                  /*Buffer to hold input from file*/
    unsigned long lenBuffULng,        /*Length of buffer (buffCStr)*/
    struct readNodeStack *readStackAry,  /*Stack, (as array) for searching*/
    unsigned long *hashSizeULng,      /*Will hold Size of hash table*/
    unsigned char *digPerKeyUChar,    /*Power of two hash size is at*/
    unsigned long *majicNumULng,      /*Holds majick number for kunths hash*/
    char *failedChar                  /*Tells if did not make hash table*/
) /*Makes a read hash array using input read ids*/
{ /*makeReadHash function*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 TOC:
   #     fun-1 Sec-1: variable declerations
   #     fun-1 Sec-2: Read in the file
   #     fun-1 sec-3: Find the hash table size & make hash table
   #     fun-1 sec-4: Find the magick number
   #     fun-1 sec-5: Build hash
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-1: variable declerations
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   char *tmpBuffCStr = 0;

   unsigned char maxHexChar = 2; /*Max hex digits to read in*/

   int lenInputInt = 0; /*Number characters in buffer*/

   unsigned long
       lenIdULng = 0; /*Stores read id length*/

   struct readInfo
       **hashTbl = 0,                 /*hash table*/
       *readTree = 0,                 /*Singe node in hash table or tree*/
       *tmpRead = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-2: Read in the file
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *hashSizeULng = 0; /*Intalize*/
   *(buffCStr + lenBuffULng) = '\0';
   tmpBuffCStr = buffCStr + lenBuffULng; /*Force initial read in*/
   lenInputInt = lenBuffULng; /*Make sure read in first read*/

   do { /*While ids to read in*/
       tmpRead = 
           cnvtIdToBigNum(
               buffCStr,
               lenBuffULng,
               &tmpBuffCStr,
               &lenInputInt,
               &lenIdULng,
               &maxHexChar,
               filtFILE
       ); /*Read in id and convert to big number*/

       if(readTree == 0)
            readTree = tmpRead;
        else
        { /*else, adding new node to list*/
            tmpRead->rightChild = readTree;
            readTree = tmpRead;
        } /*else, adding new node to list*/

        while(*tmpBuffCStr != '\n')
        { /*While not on next entry*/
            tmpBuffCStr++;

            if(*tmpBuffCStr == '\0')
            { /*If at the end of the buffer, but not at start of read*/
                if(lenInputInt < lenBuffULng)
                    break; /*At end of file*/

                lenInputInt = fread(
                                   buffCStr,
                                   sizeof(char),
                                   lenBuffULng,
                                   filtFILE
                ); /*Read in more of the file*/

                *(buffCStr + lenInputInt) = '\0';/*make sure a c-string*/
                tmpBuffCStr = buffCStr;
            } /*If at the end of the buffer, but not at start of read*/
        } /*While not on next entry*/

        (*hashSizeULng)++; /*Count number of reads*/
   } while(lenInputInt == lenBuffULng && tmpRead != 0); /*While ids to read in*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-3: Find the hash table size & make hash tabel
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *digPerKeyUChar = 0;

   /*Find power of two to use (will go one bit over to make sure lands at 0)*/
   while(*hashSizeULng > 0)
   { /*While need to find the next power of two*/
       *hashSizeULng = (*hashSizeULng) >> 1;
       (*digPerKeyUChar)++;
   } /*While need to find the next power of two*/

   (*digPerKeyUChar)++; /*Make double size, so is a bit sparse*/

   (*hashSizeULng) = 1 << (*digPerKeyUChar); /*Nearest power of two*/

   hashTbl = calloc((*hashSizeULng + 1), sizeof(struct readInfo *));
      /*Calloc intalizes with 0's*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-4: Find the magick number
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Find how many digits of my golden ratio to use*/
      /*
        The idea here is 2^(x * 3.333...) ~ 10^x, which gives me the number of
          base 10 digits per limb. This only covers the positions were having
          10^x = 2^(x * 3.333...), so will not always cover the last postion.
          For example an 64 bit number has at most 20 digits, but 64 / 3.333...
          is 19. This is due to 10^20 > 2^20 = 18,446,744,073,709,551,615.
      */

   for(
       unsigned char uCharGoldDig = 0;
       uCharGoldDig < ((sizeof(unsigned long) << 3) / powTwoPerTen);
       uCharGoldDig++
   ) /*For each digit in the godlen number I need to convert to a long*/
       *majicNumULng = 10*(*majicNumULng) + GOLDEN_RATIO_CSTR[uCharGoldDig]-48;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-5: Build hash
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(readTree != 0)
   { /*While their are reads to put in hash table*/
       tmpRead = readTree->rightChild;
       readTree->rightChild = 0;

       insertHashEntry(
           hashTbl,           /*Hash table to insert read into*/
           readTree,          /*readNode to insert into hash table*/
           *majicNumULng,  /*Will hold final majick number*/
           *digPerKeyUChar,
           readStackAry       /*Stack, (as array) for searching*/
        ); /*Insert the read into the hash table*/

        readTree = tmpRead;  /*move to the next read*/
   } /*While their are reads to put in hash table*/

   *failedChar = 0;  /*Succeed in making hash tabel*/

   return hashTbl;  /*Return head of hash table*/
} /*makeReadHash function*/

/*##############################################################################
# Output:
#   Returns: integer holding the hashed value
##############################################################################*/
unsigned long calcHash(
    struct bigNum *idBigNum,     /*read id converted to number to hash*/
    unsigned long majicNumULng, /*Majick number to mulitply (size of one limb)*/
    unsigned char digPerKeyUChar  /*Hash table size is 2^digPerKeyUChar power*/
) /*Calculate the hast for a read id*/
{ /*calcHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: Find the hash for a single read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(idBigNum->lenUsedElmChar < 2) /*The number fit in a single limb*/
       return
         *(idBigNum->bigNumAryULng) >>
         ((sizeof(unsigned long) << 3) - digPerKeyUChar); /*<< 3 (byte to bit)*/

   return
       (
        *(idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar - 1) +
        *(idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar - 2)
       ) *
       majicNumULng >>
       ((sizeof(unsigned long) << 3) - digPerKeyUChar);

   /*
     Idea is that that Kunths multiplicative hash is only only keeping the most
     significant digits. So, I can cheat & only multiply the most siginificant
     digits. This will loose some percisions since remainders are not passed,
     but removes big number arithmetic. I add the two most significant limbs,
     beacuse one limb may not be completely filled.
   */
} /*calcHash*/

/*##############################################################################
# Output:
#    Modifies: Inserts readNode into the hashTbl.
##############################################################################*/
void insertHashEntry(
    struct readInfo **hashTbl,     /*Hash table to insert read into*/
    struct readInfo *readNode,     /*readNode to insert into hash table*/
    unsigned long majicNumULng,            /*Majick number to mulitply by*/
    unsigned char digPerKeyUChar,  /*Hash table size is 2^digPerKeyUChar power*/
    struct readNodeStack *readStack /*Stack, (as array) for searching*/
) /*Iinserts a read into a hash table*/
{ /*insertHashEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 Sub-1 TOC: Insert read into hash table
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long hashULng = 0;                  /*Hold the hash value*/

    /*Find the has value*/
    hashULng = calcHash(readNode->idBigNum, majicNumULng, digPerKeyUChar);

   if(*(hashTbl + hashULng) == 0)
       *(hashTbl + hashULng) = readNode;  /*If is first node at hash value*/
   else
   { /*Else already have nodes at hash, insert node into tree*/
       if(
           insertNodeIntoReadTree(
               readNode,             /*node with read to insert in tree*/
               (hashTbl + hashULng), /*readInfo tree to insert read into*/
               readStack             /*Stack, (as array) for searching*/
           ) == 0
       ) { /*If read was a duplicate*/
           freeReadInfoStruct(&readNode); /*Free the duplicate node*/
       } /*If read was a duplicate*/
   } /*Else already have nodes at hash, insert node into tree*/

   return;
} /*insertHashEntry*/

/*##############################################################################
# Output:
#   Returns: readInfo node if found read in hash table, otherwise 0
##############################################################################*/
struct readInfo * findReadInHashTbl(
    struct bigNum *idBigNum,       /*read id converted to number to hash*/
    unsigned long majicNumULng,   /*Majick number to mulitply by*/
    unsigned char digPerKeyUChar, /*Hash table size is 2^digPerKeyUChar power*/
    struct readInfo **hashTbl     /*Hash table to search for read in*/
)  /*find a readInfo node in a hash table using the read id*/
{ /*findReadInHashTbl*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 Sub-1 TOC: Search hash table for a particler read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo *tmpReadTree = 0;

    unsigned long hashULng = 0;

    /*Find hash value*/
    hashULng = calcHash(idBigNum, majicNumULng, digPerKeyUChar);
    tmpReadTree = *(hashTbl + hashULng); /*Get tree at hash table entry*/

    return searchTree(idBigNum, tmpReadTree); /*Search hash tree to find id*/
} /*findReadInHashTbl*/

/*##############################################################################
# Output:
#    Frees: hash table & sets hashTblToFree  to 0
##############################################################################*/
void freeHashTbl(
    struct readInfo ***hashTblToFree,
    unsigned long hashSizeULng,       /*Size of hash table*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
) /*Frees a hash table of read trees*/
{ /*freeeHashTbl*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 Sub-1 TOC: free a hash table
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo
        *readTree = **hashTblToFree;

    unsigned long
        uLngTblElm = 0;

    while(uLngTblElm <= hashSizeULng)
    { /*While there are entries to free in the hash table*/
        if(readTree != 0)
        { /*If have a node to free*/
            freeReadTree(
                &readTree,         /*Tree of sequence ids to remove*/
                readStack              /*Stack to use in freeing*/
            ); /*Free the read tree at each stack*/
        } /*If have a node to free*/

        uLngTblElm++;  /*Set up for the next hashed tree*/
        readTree = *(*hashTblToFree + uLngTblElm);
    } /*While there are entries to free in the hash table*/

    free(*hashTblToFree);
    *hashTblToFree = 0;

    return;
} /*freeeHashTbl*/

/*##############################################################################
# Output:
#    Modifies: endNameCStr to pont to the '\n', ' ', & '\t' at end of read name
#    Modifies: lenIdULng to hold the length of the read id
#    Modifies: lenInputInt to hold length of input buffer
#    Returns: readInfo struct with bigNum struct having read id converted to hex
#        - 0 if fails or end of file (lenIdULng < buffSizeInt)
##############################################################################*/
struct readInfo * cnvtIdToBigNum(
    char *bufferCStr,        /*buffer to hold fread input (can have data)*/
    int buffSizeInt,         /*Size of buffer to work on*/
    char **endNameCStr,      /*Points to start of read name, will point to end*/
    int *lenInputInt,        /*Length of input from fread*/
    unsigned long *lenIdULng,/*Lengtho of the read id*/
    unsigned char *lenBigNumChar, /*Holds size to make bigNumber*/
    FILE *fastqFile          /*Fastq file to get data from*/
) /*Reads input from file & sets pointer to start of read name*/
{ /*parseFastqHeader*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC:
    #    fun-1 sec-1: Variable declerations
    #    fun-1 sec-2: loop through buffer & check if need to get input from file
    #    fun-1 sec-3: Copy header over, recored length, & incurment pointers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char charBit = 0;
    struct bigNum *idBigNum = malloc(sizeof(struct bigNum));
    struct readInfo *readNode = malloc(sizeof(struct readInfo));
    unsigned long *elmOnPtrULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: loop through buffer & check if need to get input from file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(idBigNum == 0 || readNode == 0)
    { /*If memory allocation failed*/
        if(idBigNum != 0)
            free(idBigNum);
        if(readNode != 0)
            free(readNode);

        fprintf(
            stderr,
            "Memory allocation failed: Fun-5 makeBigNumStruct %s",
            "fastqGrepStructs.c Line 153\n"
        ); /*Print error message to user*/

        return 0; 
    } /*If memory allocation failed*/

    readNode->balanceChar = 0;
    readNode->leftChild = 0;
    readNode->rightChild = 0;
    readNode->idBigNum = idBigNum;

    *lenIdULng = 0; /*Make sure starts at 0*/
    idBigNum->lenUsedElmChar = 0;
    idBigNum->bigNumAryULng = malloc(sizeof(unsigned long) * (*lenBigNumChar));
    idBigNum->lenAllElmChar = *lenBigNumChar;

    if(idBigNum->bigNumAryULng == 0)
    { /*If memory reallocation failed*/
        fprintf(
            stderr,
           "Memory allocation failed: Fun-6 fastqGrepStructs.c line 225\n"
        ); /*Print error message to user*/

        free(idBigNum);
        free(readNode);
        return 0;
    } /*If memory reallocation failed*/

    if(**endNameCStr == '\n')
        (*endNameCStr)++;

    if(**endNameCStr == '\0')
    { /*If at the end of the buffer, but not at start of read*/

        if(*lenInputInt < buffSizeInt)
        { /*If at end of file*/
          free(idBigNum->bigNumAryULng);
          free(idBigNum);
          free(readNode);
          return 0;                    /*Done with file*/
        } /*If at end of file*/

        *lenInputInt = fread(
                           bufferCStr,
                           sizeof(char),
                           buffSizeInt,
                           fastqFile
        ); /*Read in more of the file*/

      *endNameCStr = bufferCStr;
    } /*If at the end of the buffer, but not at start of read*/

      
    do { /*While still on the read name part of header*/

        if(idBigNum->lenUsedElmChar >= idBigNum->lenAllElmChar)
        { /*If need to reallocate memory*/
            (idBigNum->lenAllElmChar)++;
            (*lenBigNumChar)++;
            idBigNum->bigNumAryULng =
                realloc(idBigNum->bigNumAryULng, idBigNum->lenAllElmChar);

            if(idBigNum == 0 || readNode == 0)
            { /*If memory allocation failed*/
                free(idBigNum->bigNumAryULng);
                free(idBigNum);
                free(readNode);

                fprintf(
                    stderr,
                    "Memory allocation failed: Fun-5 makeBigNumStruct %s",
                    "fastqGrepStructs.c Line 153\n"
                ); /*Print error message to user*/

                return 0; 
            } /*If memory allocation failed*/
        } /*If need to reallocate memory*/

        /*Graph unsigned long element working on*/
        elmOnPtrULng = idBigNum->bigNumAryULng + idBigNum->lenUsedElmChar;
        *elmOnPtrULng = 0;
        charBit = 0;

        while(charBit < (sizeof(unsigned long) << 3))
        { /*while empty bits in the current big number unsigned long element*/
            if(hexCharAry[**endNameCStr] & 64)
                break; /*If have finshed converting the hex string*/

            else if(!(hexCharAry[**endNameCStr] & 32)) /*Was hex char*/
            { /*Else if is a hex character*/
                *elmOnPtrULng =
                    *elmOnPtrULng +
                    (hexCharAry[**endNameCStr] << charBit);
                charBit += 4;   /*make sureo only recored conversion*/
            } /*Else if is a hex character*/

            (*lenIdULng)++;
            (*endNameCStr)++;

            if(**endNameCStr == '\0')
            { /*If ran out of buffer & need to read in more of the file*/
                if(*lenInputInt < buffSizeInt)
                    return readNode; /*at end of file*/

                *lenInputInt = fread(bufferCStr,
                                     sizeof(char),
                                     buffSizeInt,
                                     fastqFile
                ); /*Read in more of the file*/

                *(bufferCStr + *lenInputInt) = '\0';/*make sure a c-string*/
                *endNameCStr = bufferCStr;
                continue;
            } /*If ran out of buffer & need to read more of the file*/
        } /*while empty bits in the current big number unsigned long element*/

        (idBigNum->lenUsedElmChar)++; /*Track number ULngs acctualy used*/
    } while(**endNameCStr > 32); /*While still on the read name part of header*/

    return readNode; /*Copied name sucessfully*/
} /*parseFastqHeader*/

