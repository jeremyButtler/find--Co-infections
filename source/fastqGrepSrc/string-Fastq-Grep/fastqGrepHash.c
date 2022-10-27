/*##############################################################################
# Name: fastqGrepHash.c
#   Use:
#     Uses an hash + AVL tree to find if reads is in user suplied list
##############################################################################*/

#include "fastqGrepHash.h"

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
##############################################################################*/
struct readInfo ** makeReadHash(
    FILE * filtFILE,                  /*file with read id's to filter by*/
    char * buffCStr,                  /*Buffer to hold input from file*/
    unsigned long lenBuffULng,        /*Length of buffer (buffCStr)*/
    struct readNodeStack *readStackAry,  /*Stack, (as array) for searching*/
    unsigned long *hashSizeULng,      /*Will hold Size of hash table*/
    char *failedChar                  /*Tells if did not make hash table*/
) /*Makes a read hash array using input read ids*/
{ /*makeReadHash function*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 TOC:
   #     fun-1 Sec-1: variable declerations
   #     fun-1 Sec-2: Read in the file
   #     fun-1 Sec-3: Build hash
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-1: variable declerations
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char
       *readPosCStr = 0;              /*Incurment through read id*/

   unsigned char
       lenNameUChar = 0;              /*Length of readPosCStr*/

   struct readInfo
       **hashTbl = 0,                 /*hash table*/
       *readTree = 0,                 /*Singe node in hash table or tree*/
       *tmpRead = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-2: Read in the file
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *hashSizeULng = 0; /*Intalize*/

    while(
        fgets(
            buffCStr,                  /*buffer to put file line into*/
            lenBuffULng,               /*Length of the buffer*/
            filtFILE                  /*File to extract read id's from*/
        ) /*Get the next read id (line) in the file*/
    ) { /*While there are read id's to read in*/

        readPosCStr = buffCStr;
        lenNameUChar = 0;            /*Reset for new read*/

        /*Find the length of the read name*/
        while(
            *readPosCStr != '\0' &&          /*End of lineInCStr*/
            *readPosCStr != ' ' &&          /*End of read name*/
            *readPosCStr != '\t' &&         /*End of read name (?)*/
            *readPosCStr != '\n'            /*End of line*/
        ) { /*While still on the read name part of header*/
            readPosCStr++;
            lenNameUChar++;
        } /*While still on the read name part of header*/

        *readPosCStr = '\0';       /*Only want first part, so can ignore rest*/

        tmpRead =
            makeReadInfoStruct(
                buffCStr,           /*buffer with name*/
                lenNameUChar,       /*Length of input name*/
                '@'                 /*Ignore @ symbol at start*/
        ); /*Make a read info struct with read name*/

        if(tmpRead == 0)
        { /*If could not allocate memeory*/
            fprintf(
                stderr,
                "makeReadHash Fun-1 fastqGrepHash:67, Failed memory allocation)"
            ); /*Let user know error happened*/

            *failedChar = 1;
            return 0;
        } /*If could not allocate memeory*/

        if(readTree == 0)
            readTree = tmpRead;
        else
        { /*else, adding new node to list*/
            tmpRead->rightChild = readTree;
            readTree = tmpRead;
        } /*else, adding new node to list*/

        (*hashSizeULng)++; /*Count number of reads*/
    } /*While there are read id's to read in*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # Fun-1 Sec-3: Build hash
   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   hashTbl = malloc(sizeof(struct readInfo *) * (*hashSizeULng + 1));

   /*Intalize hash table & make sure ends with 0*/
   for(unsigned long uLngCnt = 0; uLngCnt <= *hashSizeULng; uLngCnt++)
       *(hashTbl + uLngCnt) = 0; /*Set end of array to null*/

   while(readTree != 0)
   { /*While their are reads to put in hash table*/
       tmpRead = readTree->rightChild;
       readTree->rightChild = 0;

       insertHashEntry(
           hashTbl,           /*Hash table to insert read into*/
           readTree,          /*readNode to insert into hash table*/
           *hashSizeULng,     /*Size of array holding hash*/
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
    char *readIdCStr,             /*read id to hash*/
    unsigned long hashSizeULng    /*Size of array holding hash*/
) /*Calculate the hast for a read id*/
{ /*calcHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1 TOC: Find the hash for a single read id
    #   Fun-2 Sec-1: Find hash for single read id
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *tmpCStr = readIdCStr;

    unsigned long
        hashULng = 1;

    if(*tmpCStr == '@')
        tmpCStr++;        /*Move off the @ symbol, so it is not hashed*/ 

    while(
        *tmpCStr != '\n' &&              /*If at new line (end of id)*/
        *tmpCStr != ' ' &&               /*Space marks end of the read id*/
        *tmpCStr != '\0'                 /*Null is the end of a c-string*/
    ) { /*While there are characters to hash*/

        if(
            *tmpCStr > 64 &&
            *tmpCStr < 91
        ) { /*If is uppercase*/
            hashULng = hashULng * (*tmpCStr + 32); /*convert to lowercase*/
        } /*If is uppercase*/

        else
            hashULng = hashULng * (*tmpCStr);

        tmpCStr++; /*Move to next character*/
    } /*While there are characters to hash*/

    return hashULng % hashSizeULng;     /*Return the final hash*/
} /*calcHash*/

/*##############################################################################
# Output:
#    Modifies: Inserts readNode into the hashTbl.
##############################################################################*/
void insertHashEntry(
    struct readInfo **hashTbl,          /*Hash table to insert read into*/
    struct readInfo *readNode,          /*readNode to insert into hash table*/
    unsigned long hashSizeULng,         /*Size of array holding hash*/
    struct readNodeStack *readStack   /*Stack, (as array) for searching*/
) /*Iinserts a read into a hash table*/
{ /*insertHashEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: Insert read into hash table
    #   Fun-3 Sec-1: Insert read into hash table
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long
        hashULng = 0;                  /*Hold the hash value*/

    hashULng =
        calcHash(
            readNode->idCStr,        /*read id to hash*/
            hashSizeULng             /*Size of array holding hash*/
    ); /*Calculate the hast for a read id*/

   if(*(hashTbl + hashULng) == 0)
       *(hashTbl + hashULng) = readNode;  /*If is first node at hash value*/
   else
   { /*Else already have nodes at hash, insert node into tree*/
       insertNodeIntoReadTree(
           readNode,              /*node with read to insert in tree*/
           (hashTbl + hashULng), /*readInfo tree to insert read into*/
           readStack             /*Stack, (as array) for searching*/
       );
   } /*Else already have nodes at hash, insert node into tree*/

   return;
} /*insertHashEntry*/

/*##############################################################################
# Output:
#   Returns: readInfo node if found read in hash table, otherwise 0
##############################################################################*/
struct readInfo * findReadInHashTbl(
    char *readIdCStr,                   /*Read to search for*/
    struct readInfo **hashTbl,          /*Hash table to search for read in*/
    unsigned long hashSizeULng          /*Size of array holding hash*/
)  /*find a readInfo node in a hash table using the read id*/
{ /*findReadInHashTbl*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1 TOC: Search hash table for a particler read id
    #   Fun-4 Sec-1: Search hash table to see if holds a read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct readInfo
        *tmpReadTree = 0;
    unsigned long 
        hashULng = 0;

    hashULng = 
        calcHash(
            readIdCStr,              /*read id to hash*/
            hashSizeULng             /*Size of array holding hash*/
    ); /*Calculate the hast for a read id*/

   tmpReadTree = *(hashTbl + hashULng); /*Get tree at hash table entry*/

   /*if(tmpReadTree == 0)
       return 0;*/ /*I think this is just EXTRA, but want to make sure*/

   return 
       searchTree(
           readIdCStr,
           tmpReadTree 
   ); /*Search the hash position to see if read exists in tree*/
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
