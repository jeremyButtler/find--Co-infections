#ifndef ALIGNMENTSFUN_H
#define ALIGNMENTSFUN_H

#include "defaultSettings.h"
#include "cStrToNumberFun.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define defClearNonAlph (1 | 2 | 4 | 8 | 16) // clear 64 bit and case
#define defToUper (1 | 2 | 4 | 8 | 16 | 64)
#define defMoveLeft 0
#define defMoveUp 1
#define defMoveDiagnol 2
#define defMoveMatch 3

#define defMatchFlag 12 // match (base | 8)
#define defBaseFlag 4   // match or snp
#define defInsFlag 2    // insertion
#define defDelFlag 1    // deletion
    // Clear 32nd bit (marks lower case)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - st-01 alnSet:
'     o Holds settings for my alignment program
'  - fun-04 cnvtAlnErrToSeq:
'     o Uses an array of error types and a c-string with the a sequence
'       to make one part of an alignment
'  - fun-05 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-06 getElmFromToBitUCAry:
'     o Get an element from a two bit array
'  - fun-07 twoBitAryShiftBytsForNewElm:
'     o Make room in a unit8_t for two more two-bit elements
'  - fun-08 twoBitAryMoveToNextElm:
'     o Moves to the next element in a two-bit array
'  - fun-09 towBitAryMoveBackOneElm:
'     o Moves back one element in a 2-bit array
'  - fun-10 towBitAryMoveBackXElm:
'     o Moves back X elements in a 2-bit array
'  - fun-11 makeAlnSet:
'     o Makes & initalizes an alnSet structer on the heap
'  - fun-12 getBasePairScore:
'     o Get the score for a pair of bases from an alignment structure
'  - fun-13 initAlnSet:
'     o Set values in altSet (alingment settings) structure to defaults
'  - fun-14 setAlnSetSnpPenalty:
'     o Changes SNP/Match penalty for one query/reference combination
'  - fun-15 freeAlnSet:
'     o Frees and alnSet (alignment settings) structure
'  - fun-16 checkIfBasesMatch:
'     o Check if two bases are the same (includes anonymous bases)
'  - fun-17 cnvtAlnErrAryToLetter:
'     o Converts an alignment error array from my Needleman-Wunsch
'       alignment into an array of letters (I = insertion, D = deletion,
'       = = match, X = snp) [These codes are from the eqx cigar entry]
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| ST-01: alnSet
| Use: Holds settings for my alignment program
\---------------------------------------------------------------------*/
typedef struct alnSet
{ /*alnSet*/
   //kmer mapping variables
   double percOverlapD;   // Percentage of overlap between chunks
   int32_t maxChunksI;    // Maximum number of chunks to process at once
   int32_t minChunkLenI;  // Point to Start Needleman-Wunsch alignment
   double minPercKmerD;   // Min % of shared kmers to keep a chunk map

   // Needleman-Wunsch variables
   int16_t snpPenaltyC[26][26];   // Penalty for mismatches (in matrix)
     // Size is due to wanting a look up table that can handle
     // anonymous bases. Most cells will be set to 0.
     // value = snpPenaltyC[(uint8_t) (base1 & defClearNonAlph) - 1 ]
     //                    [(uint8_t) (base2 & defClearNonAlph) - 1 ]
   int32_t gapStartPenaltyI;     // Penalty for starting an indel
   int32_t gapExtendPenaltyI;    // Penalty for extending an indel
   uint32_t minScoreUI;        // Minimum score needed to keep alignment
}alnSet;

/*---------------------------------------------------------------------\
| Output: Heap alloacted C-string with alignment for the input sequence
\---------------------------------------------------------------------*/
char * cnvtAlnErrToSeq(
    char *seqCStr,        // c-string with query sequence
    int32_t seqStartI,    // Were alignment starts on query (index 1)
    char queryBl,         // 1: input sequence is query, 0: is reference
    uint8_t *alnErrUCAry, // Holds error types (ends with 0)
    uint32_t lenErrAryUI  // length of alnErrUCAry (for sequence array)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: cnvtAlnQueryErrToSeq
   '  - Uses an array of error types and a c-string with the a sequence
   '    to make one part of an alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Returns:
|    o array with flags for snp/match, insertion, and deletions at each
|      position. (1 = snp/match, 2 = insertion, 4 = deletion)
|    o 0 for memory allocation errors
|  - Modifies:
|    o lenErrAryUI to hold the length of the returned array
\---------------------------------------------------------------------*/
uint8_t * NeedleManWunschAln(
    char *queryCStr,        // Full query sequence as c-string
    int32_t queryStartI,    // Starting query coordinate for alignment
    int32_t queryEndI,      // Ending query coordinate for alignment
    char *refCStr,          // Full reference sequence as c-string
    int32_t refStartI,      // Starting reference coordinate for aln
    int32_t refEndI,        // Ending reference coordinate for alignment
    struct alnSet *settings,// Settings for the alignment
    uint32_t *lenErrAryUI,  // Will hold the return arrays length
    long *scoreL            // Score for the alignment (bottom right)
    // *startI and *endI paramaters should be index 1
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: NeedleManWunschAln
   '  - Perform a Needleman-Wunsch alignment on input sequences
   '  o fun-06 sec-1: Variable declerations
   '  o fun-06 sec-2: Allocate memory for alignment
   '  o fun-06 sec-3: Fill in the initial negatives for the reference
   '  o fun-06 sec-4: Fill the matrix with scores
   '  o fun-06 sec-5: Find the best path
   '  o fun-06 sec-6: Clean up and invert error array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Two bits of interest from the two bit array
\---------------------------------------------------------------------*/
uint8_t getTwoBitAryElm(
    uint8_t *twoBitAryUC,// Array of bytes that contain 4 2 bit elements
    uint8_t *charElmOnUC // The two bits on in the unit8_t
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
   '  - Get an element from a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAryUC to be shifted to (move last two bits further up) if
|      the unit8_t has room for two more bits
|    o twoBitAryUC to point to next uint8_t element if need to shift to
|      the next set of four two bits
|    o Incurements bitUC to hold number of two bits in the current
|      unit8_t
|    o Sets bitUC to 0 if moving to the next uint8_t element
| Note:
|  - This function does not change the memory size
|  - You have to change the two bit values your self
\---------------------------------------------------------------------*/
void twoBitAryShiftBitsForNewElm(
    uint8_t **twoBitAryUC, // To bit array to add space for a new element
    uint8_t *bitUC         // Which bit am I on in the uint8_t
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: twoBitAryShiftBytsForNewElm
   '  - Make room in a unit8_t for two more two-bit elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAryUC to point to next uint8_t element if the next two bit
|      elements are in the next unit8_t
|      (one uint8_t holds four 2-bit elements)
|    o Incurements bitUC to be on the next 2-bit element in the current
|      uint8_t
|    o Sets bitUC to 0 if moving to the next uint8_t in the array
\---------------------------------------------------------------------*/
void twoBitAryMoveToNextElm(
    uint8_t **twoBitAryUC, // To bit array to move to next element in
    uint8_t *bitUC         // Which bit am I on in the current uint8_t
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAryUC to point to the previous uint8_t element if the r
|      previous two bit elements are in the previous unit8_t
|      (one uint8_t holds four 2-bit elements)
|    o Decurments bitUC to be on the previouls 2-bit element in the
|      current uint8_t
|    o Sets bitUC to 0 if moving to the previous uint8_t in the array
\---------------------------------------------------------------------*/
void twoBitAryMoveBackOneElm(
    uint8_t **twoBitAryUC, // 2-bit array to move to previous element in
    uint8_t *bitUC         // Which bit am I on in the current uint8_t
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: towBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAryUC to point to the target uint8_t element
|    o Sets bitUC to be on the target bit count of the shift
\---------------------------------------------------------------------*/
void twoBitAryMoveBackXElm(
    uint8_t **twoBitAryUC, // 2-bit array to move to previous element in
    uint8_t *bitUC,        // Which bit am I on in the current uint8_t
    int32_t shiftByI       // How many elements to shift back by
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: towBitAryMoveBackXElm
   '  - Moves back X elements in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Returns: Initalized alnSet structer or 0 for memory error
| Note: I prefer using stack, so this is more here for someone else
\---------------------------------------------------------------------*/
struct alnSet * makeAlnSetST(
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: makeAlnSet
   '  - Makes & initalizes an alnSet structer on the heap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Returns score of a single pair of bases
\---------------------------------------------------------------------*/
int16_t getBasePairScore(
    const char *queryBaseC, // Query base of pair to get score for
    const char *refBaseC,   // Reference base of pair to get score for
    struct alnSet *alnSetST // structure with scoring matrix to change
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: getBasePairScore
   '  - Get the score for a pair of bases from an alignment structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings) structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: one score in an snp/match scoring matrix
\---------------------------------------------------------------------*/
void setBasePairScore(
    const char *queryBaseC, // Query base to change score for
    const char *refBaseC,   // Reference base to change score for
    int16_t newScoreC,      // New value for [query][ref] combination
    struct alnSet *alnSetST // structure with scoring matrix to change
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: 1: if bases were a match; 0 if not
\---------------------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC, // Query base to see if same as reference
    char *refBaseC    // Reference base to see if same as query
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Check if two bases are the same (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnErrUCAry to have letters instead of codes
\---------------------------------------------------------------------*/
void cnvtAlnErrAryToLetter(
    uint8_t *alnErrUCAry   // Array array from Needleman (fun-05) to 
                          // convert to legible characters
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
   '  - Converts an alignment error array from my Needleman-Wunsch
   '    alignment into an array of letters (I = insertion, D = deletion,
   '    = = match, X = snp) [These codes are from the eqx cigar entry]
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

unsigned long readInScoreFile(
    struct alnSet *alnSetST,  // structure with scoring matrix to change
    FILE *scoreFILE           // File of scores for a scoring matrix
);  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-18 TOC: readInScoreFile
    '  - Reads in a file of scores for a scoring matrix
    '  o fun-18 sec-1: Variable declerations and buffer set up
    '  o fun-18 sec-2: Read in line and check if comment
    '  o fun-18 sec-3: Get score, update matrix, & move to next line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
