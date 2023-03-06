/*##############################################################################
# Name: fastqGrepHash.c
#   Use:
#     Uses an hash + AVL tree to find if reads is in user suplied list
##############################################################################*/

#include "fqGetIdsHash.h"

uint8_t /*First 392 digitits (defing as array so can modify)*/
   GOLDEN_RATIO[] = "618033988749894848204586834365638117720309179805762862135448622705260462818902449707207204189391137484754088075386891752126633862223536931793180060766726354433389086595939582905638322661319928290267880675208766892501711696207032221043216269548626296313614438149758701220340805887954454749246185695364864449241044320771344947049565846788509874339442212544877066478091588460749988712400765";

double
   powTwoPerTen = 3.3333333333333333333333333333333333333333333;
   /*there are 2^(3.33) per 10^1. This is so I can avoid gmp right shift,
     which seems to error out on realloc for O2*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
\ fastqGrepHash TOC:                                                   \
/    fun-1 makeReadHash:                                               /
\        - Make a hash table of read id's from a file of read ids      \
/    fun-2 calcHash:                                                   /
\        - Find the hash from the read id bignumber                    \
/    fun-3 Insert read into hash table                                 /
\        - Insert a read id (is big number) into a hash table          \
/    fun-4 findReadInHashTbl:                                          /
\        -  Search hash table for a particler read id                  \
/    fun-5 freeHashTbl:                                                /
\        - Free a hash table                                           \
/    fun-6 readListToHash:                                             /
\        - Convert a lined readList to a hash table                    \
/        - Only use readList->rightChild pointer in the linked list    /
\    fun-7 findMajicNumber:                                            \
/        - Get the first 64 bits of the golden ratio                   /
\                                                                      \
/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
    FILE * filtFILE,               /*file with read id's to filter by*/
    char * buffCStr,            /*Buffer to hold input from file*/
    uint64_t lenBuffULng,          /*Length of buffer (buffCStr)*/
    struct readNodeStack *readStackAry, /*Stack, (as array) for search*/
    uint64_t *hashSizeULng,       /*Will hold Size of hash table*/
    uint8_t *digPerKeyUChar,      /*Power of two hash size is at*/
    unsigned long *majicNumULng,  /*Holds majick number for kunths hash*/
    uint8_t *failedChar           /*Tells if did not make hash table*/
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

   char
       *tmpBuffCStr = 0;

   uint8_t
       maxHexChar = 1; /*Max hex digits to read in*/

   uint64_t
       lenInputULng = 0; /*Number characters in buffer*/

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
   lenInputULng = lenBuffULng; /*Make sure read in first read*/

   do { /*While ids to read in*/
       tmpRead = 
           cnvtIdToBigNum(
               buffCStr,
               lenBuffULng,
               &tmpBuffCStr,
               &lenInputULng,
               &maxHexChar,
               filtFILE
       ); /*Read in id and convert to big number*/

       if(tmpRead == 0)
       { /*If was a falied read*/
           if(lenInputULng == 0)
           { /*If was a memory allocation error (message already printed)*/
               freeReadTree(&readTree, readStackAry);
               return 0;
           } /*If was a memory allocation error (message already printed)*/

           break; /*end of file*/
       } /*If was a falied read*/

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
                if(lenInputULng < lenBuffULng)
                    break; /*At end of file*/

                lenInputULng = fread(
                                   buffCStr,
                                   sizeof(char),
                                   lenBuffULng,
                                   filtFILE
                ); /*Read in more of the file*/

                *(buffCStr + lenInputULng) = '\0';/*make sure a c-string*/
                tmpBuffCStr = buffCStr;
            } /*If at the end of the buffer, but not at start of read*/
        } /*While not on next entry*/

        (*hashSizeULng)++; /*Count number of reads*/
   } while(tmpRead != 0);/*ids to read*/

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
   *majicNumULng = findMajicNumber();

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
           majicNumULng,  /*Will hold final majick number*/
           digPerKeyUChar,
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
uint64_t calcHash(
    struct bigNum *idBigNum, /*read id converted to number to hash*/
    const unsigned long *majicNumULng, /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar /*Hash table size 2^digPerKeyUChar*/
) /*Calculate the hast for a read id*/
{ /*calcHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: Find the hash for a single read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Hash calucation for speed or memory efficent (MEM) compiles*/ 
    #ifndef MEM
        return
             (idBigNum->totalL * *majicNumULng)
          >> (defBitsInUL- *digPerKeyUChar);
    #else
        if(idBigNum->lenUsedElmChar == 1)
        { /*If only had one element*/
            return
                  (*(idBigNum->bigNumAryIOrL) * *majicNumULng)
               >> (defBitsInUL - *digPerKeyUChar);
                 /*<< 3 (byte to bit)*/
        } /*If only had one element*/
    
        /*Make sure get at least on full size array in hash*/
        return
           ( 
              *(idBigNum->bigNumAryIOrL + idBigNum->lenUsedElmChar - 1)
            + *(idBigNum->bigNumAryIOrL + idBigNum->lenUsedElmChar - 2)
            * *majicNumULng
           )
            >> (defBitsInUL - *digPerKeyUChar);
                 /*<< 3 (byte to bit)*/
    #endif

   /*Using unsinged long, because want to vary with arcitecture*/
   /*
     Idea is that that Kunths multiplicative hash is only keeping the
     most significant digits. So, I can cheat & only multiply the most
     siginificant digits. This will loose some percisions since
     remainders are not passed,
   */
} /*calcHash*/

/*##############################################################################
# Output:
#    Modifies: Inserts readNode into the hashTbl.
##############################################################################*/
void insertHashEntry(
    struct readInfo **hashTbl,   /*Hash table to insert read into*/
    struct readInfo *readNode,   /*readNode to insert into hash table*/
    const unsigned long *majicNumULng,  /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar, /*Hash table size 2^digPerKeyUChar*/
    struct readNodeStack *readStack /*Stack, (as array) for searching*/
) /*Iinserts a read into a hash table*/
{ /*insertHashEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 Sub-1 TOC: Insert read into hash table
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint64_t hashULng = 0; /*Hold the hash value*/

    hashULng =
        calcHash(
            readNode->idBigNum,
            majicNumULng,
            digPerKeyUChar
   ); /*Find the has value*/

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
    const unsigned long *majicNumULng, /*Majick number to mulitply by*/
    const uint8_t *digPerKeyUChar,/*Hash table size is 2^digPerKeyUChar*/
    struct readInfo **hashTbl     /*Hash table to search for read in*/
)  /*find a readInfo node in a hash table using the read id*/
{ /*findReadInHashTbl*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 Sub-1 TOC: Search hash table for a particler read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo *tmpReadTree = 0;

    uint64_t hashULng = 0;

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
    const uint64_t *hashSizeULng,     /*Size of hash table*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
) /*Frees a hash table of read trees*/
{ /*freeeHashTbl*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 Sub-1 TOC: free a hash table
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo
        *readTree = **hashTblToFree;

    uint64_t
        uLngTblElm = 0;

    while(uLngTblElm <= *hashSizeULng)
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

/*---------------------------------------------------------------------\
| Output:                                                              |
|    Returns:                                                          |
|        - readInfo: hash table (as heap array)                        |
|        - 0: if failed                                                |
|    Modifies:                                                         |
|        - hashSizeUlng to hold length of the returned array           |
|        - numCharChar: to hold the number of characters used in hash  |
|        - digPerKeyUChar: to hold the hash size (as multiple of two)  |
|        - magickNumULng: to hold magick number for this hash table    |
\---------------------------------------------------------------------*/
struct readInfo ** readListToHash(
    struct readInfo * readList,   /*List of id's to build hash from.
                                  only use the rightChild ptr in list*/
    const uint64_t *lenListULng,   /*Number of id's in readList*/
    struct readNodeStack *readStackAry, /*Stack (array) for searching*/
    uint64_t *hashSizeULng,  /*Will hold Size of hash table*/
    uint8_t *digPerKeyUChar,       /*Power of two hash size is at*/
    unsigned long *majicNumULng  /*Holds majick number for kunths hash*/
) /*Makes a read hash array using input read ids*/
{ /*makeReadHash function*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   \ Fun-6 TOC:                                                        /
   /     fun-6 sec-1: Find the hash table size & make hash table       \
   \     fun-6 sec-2: Find the magick number                           /
   /     fun-6 sec-3: Build hash                                       \
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-6 Sec-1: variable declerations                                v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct readInfo
       *tmpRead = 0,
       **hashTbl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-6 Sec-2: Find hash table size & make hash tabel               v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *hashSizeULng = *lenListULng;
   *digPerKeyUChar = 0;

   /*Find power of two to use (will go 1 bit over to make sure 0)*/
   while(*hashSizeULng > 0)
   { /*While need to find the next power of two*/
       *hashSizeULng = (*hashSizeULng) >> 1;
       (*digPerKeyUChar)++;
   } /*While need to find the next power of two*/

   (*digPerKeyUChar)++; /*Make double size, so is a bit sparse*/

   (*hashSizeULng) = 1 << (*digPerKeyUChar); /*Nearest power of two*/

   hashTbl = calloc((*hashSizeULng + 1), sizeof(struct readInfo *));
      /*Calloc intalizes with 0's*/

   if(hashTbl == 0)
       return 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-6 Sec-3: Find the magick number                               v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Find how many digits of my golden ratio to use*/
   *majicNumULng = findMajicNumber();

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-6 Sec-2: Build hash                                           v
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(readList != 0)
   { /*While their are reads to put in hash table*/
       tmpRead = readList->rightChild;
       readList->rightChild = 0;

       insertHashEntry(
           hashTbl,           /*Hash table to insert read into*/
           readList,          /*readNode to insert into hash table*/
           majicNumULng,  /*Will hold final majick number*/
           digPerKeyUChar,
           readStackAry       /*Stack, (as array) for searching*/
        ); /*Insert the read into the hash table*/

        readList = tmpRead;  /*move to the next read*/
   } /*While their are reads to put in hash table*/

   return hashTbl;  /*Return head of hash table*/
} /*makeReadHash function*/

/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o Unsigned long with the first 64 bits of the golden number
\---------------------------------------------------------------------*/
unsigned long findMajicNumber()
/*Gets first 64 bits of the goldent ratio*/
{ /*findMajicNumber*/

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-7 TOC: Sec-1 Sub-1: findMajicNumber
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   unsigned long majicNumUL = 0;
   unsigned char numDigUC =((sizeof(unsigned long) << 3)/powTwoPerTen);

      /*
        The idea here is 2^(x * 3.333...) ~ 10^x, which gives me the
        number of base 10 digits per limb. This only covers the
        positions were having 10^x = 2^(x * 3.333...), so will not
        always cover the last postion.  For example an 64 bit number has
        at most 20 digits, but 64 / 3.333...  is 19. This is due to
        10^20 > 2^20 = 18,446,744,073,709,551,615.
      */

   /*Grab the first 64 bits of numbers in golden number into a long*/
   for(unsigned char uCGoldDig = 0; uCGoldDig < numDigUC; ++uCGoldDig)
       majicNumUL = 10*majicNumUL + (GOLDEN_RATIO[uCGoldDig] & ~48);

   return majicNumUL;
} /*findMajicNumber*/
