/*######################################################################
# Name alignmentsFun
# Use:
#  o Holds functions for doing pairwise alignments (Needleman/waterman)
# Includes:
#   - "defaultSettings.h"
#   - "cStrToNumberFun.h"
#   - "twoBitArrays.h"
# C Standard libraries:
#   - <stdlib.h>
#   - <stdio.h>
#   o <stdint.h>
######################################################################*/

#ifndef ALIGNMENTSFUN_H
#define ALIGNMENTSFUN_H

#include "defaultSettings.h"
#include "cStrToNumberFun.h"
#include "twoBitArrays.h"

#include <stdlib.h>
#include <stdio.h>

#define defClearNonAlph (1 | 2 | 4 | 8 | 16) // clear 64 bit and case
#define defToUper (1 | 2 | 4 | 8 | 16 | 64)
    // Clear 32nd bit (marks lower case)

#define defMoveStop 0    // Do not move
#define defMoveLeft 1    // Move left (deletion) in alignment matrix
#define defMoveUp 2      // Move up (insertion) in alignment matrix
#define defMoveDiagnol 3 // Move on a diagnol (snp/match) in alignment

#define defBaseFlag 4   // match or snp
#define defInsFlag 2    // insertion
#define defDelFlag 1    // deletion
#define defSoftQueryFlag 8 // Softmask a query base
#define defSoftRefFlag 16  // Softmask a reference base

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  - st-01 alnSet:
'     o Holds settings for my alignment program
'  - fun-03 cnvtAlnErrToSeq:
'     o Uses an array of error types and a c-string with the a sequence
'       to make one part of an alignment
'  - fun-04 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input sequences
'  - fun-05 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-07 checkIfBasesMatch
'     o Check if two bases are the same (includes anonymous bases)
'  - fun-12 makeAlnSet:
'     o Makes & initalizes an alnSet structer on the heap
'  - fun-13 getBasePairScore:
'     o Get the score for a pair of bases from an alignment structure
'  - fun-14 initAlnSet:
'     o Set values in altSet (alingment settings) structure to defaults
'  - fun-15 setAlnSetSnpPenalty:
'     o Changes SNP/Match penalty for one query/reference combination
'  - fun-16 freeAlnSet:
'     o Frees and alnSet (alignment settings) structure
'  - fun-17 cnvtAlnErrAryToLetter:
'     o Converts an alignment error array from my Needleman-Wunsch
'       alignment into an array of letters (I = insertion, D = deletion,
'       = = match, X = snp) [These codes are from the eqx cigar entry]
'  - fun-18 readInScoreFile
'     o Reads in a file of scores for a scoring matrix
'  - fun-19 getIndelScore:
'     o Gets an indel score for the current cell
'  - fun-20 updateDirAndScore:
'     o Picks the best score and direction for the current base pairs
'       being compared in a Needleman Wunsch alignment
'  - fun-21 updateDirAndScoreWater:
'     o Picks the best score and direction for the current base pairs
'       being compared in a Waterman Smith alignment
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
   char diagnolPriorityC; // 0 favor snps; 1 kinda; 2 do not
   char topPriorityC;     // 0 favor insertions; 1 kinda; 2 do not
   char leftPriorityC;    // 0 favor deletions; 1 kinda; 2 do not
   char matchPriorityBl;  // Match scores always beat all other scores

   // Needleman-Wunsch / Waterman Smith variables
   char useNeedleBl;
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
   ' Fun-03 TOC: Sec-1 Sub-1: cnvtAlnQueryErrToSeq
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
uint8_t * WatermanSmithAln(
    char *queryCStr,        // Full query sequence as c-string
    int32_t queryStartI,    // Starting query coordinate for alignment
    int32_t queryEndI,      // Ending query coordinate for alignment
    char *refCStr,          // Full reference sequence as c-string
    int32_t refStartI,      // Starting reference coordinate for aln
    int32_t refEndI,        // Ending reference coordinate for alignment
    struct alnSet *settings,// Settings for the alignment
    uint32_t *lenErrAryUI,  // Will hold the return arrays length
    long *scoreL        // Score for the alignment (bottom right)
    // *startI and *endI paramaters should be index 1
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: WatermanSmithAln
   '  - Perform a Waterman Smith alignment on input sequences
   '  o fun-04 sec-1: Variable declerations
   '  o fun-04 sec-2: Allocate memory for alignment
   '  o fun-04 sec-3: Fill in the initial negatives for the reference
   '  o fun-04 sec-4: Fill the matrix with scores
   '  o fun-04 sec-5: Find the best path
   '  o fun-04 sec-6: Clean up and add softmasking to start
   '  o fun-04 sec-7: Invert the error type array
   '  o fun-04 sec-8: Add softmasking to the end and return the array
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
   '  o fun-05 sec-1: Variable declerations
   '  o fun-05 sec-2: Allocate memory for alignment
   '  o fun-05 sec-3: Fill in the initial negatives for the reference
   '  o fun-05 sec-4: Fill the matrix with scores
   '  o fun-05 sec-5: Find the best path
   '  o fun-05 sec-6: Clean up and invert error array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: 1: if bases were a match; 0 if not
\---------------------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC, // Query base to see if same as reference
    char *refBaseC    // Reference base to see if same as query
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Check if two bases are the same (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



/*---------------------------------------------------------------------\
| Output: Returns: Initalized alnSet structer or 0 for memory error
| Note: I prefer using stack, so this is more here for someone else
\---------------------------------------------------------------------*/
struct alnSet * makeAlnSetST(
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeAlnSet
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
   ' Fun-13 TOC: Sec-1 Sub-1: getBasePairScore
   '  - Get the score for a pair of bases from an alignment structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: initAlnSet
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
   ' Fun-15 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: alnErrUCAry to have letters instead of codes
\---------------------------------------------------------------------*/
void cnvtAlnErrAryToLetter(
    char *refSeqCStr,      // Reference sequence for detecting matches
                           // Input 0 to ignore matches
    char *querySeqCStr,    // query sequence for detecting matches
                           // Input 0 to ignore matches
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

/*---------------------------------------------------------------------\
| Output: Returns the score for an Needleman Wunsch indel
\---------------------------------------------------------------------*/
long getIndelScore(
    uint8_t *lastDirUC,      // Cell gettign last score from
    uint8_t *lastBitUC,      // two bit element on in lastDirUC
    struct alnSet *alnSetST, // Holds the gap open & extension penalties
    long *lastBaseL          // Has score of the last base
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-19 TOC: Sec-1 Sub-1: getIndelScore
   '  - Gets an indel score for the current cell
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: scoreOnL and dirOnUC to hold best score & direction
\---------------------------------------------------------------------*/
void updateDirAndScoreNeedle(
    uint8_t *dirOnUCPtr,    // Direction on with first two bits cleared
    struct alnSet *alnSetST, // Has preference for score selection
    long *scoreTopL,     // Score for an insertion
    long *scoreDiagnolL, // Score for an match/snp
    long *scoreLeftL,    // The score for an deletion
    long *scoreOnL       // Score to update
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-20 TOC: updateDirAndScore
   '  - Picks the best score and direction for the current base pairs
   '    being compared in a Needleman Wunsch alignment
   '  o fun-20 sec-1: Matches->insertions->deletions
   '  o fun-20 sec-2: Matches->deletions->insertions
   '  o fun-20 sec-3: Insertions->matches->deletions
   '  o fun-20 sec-4: Deletions->matches->insertions
   '  o fun-20 sec-5: Insertions->deletions->matches
   '  o fun-20 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: scoreOnL and dirOnUC to hold best score & direction
\---------------------------------------------------------------------*/
void updateDirAndScoreWater(
    uint8_t *dirOnUCPtr,     // Direction on with first two bits cleared
    struct alnSet *alnSetST, // Has preference for score selection
    long *scoreTopL,     // Score for an insertion
    long *scoreDiagnolL, // Score for an match/snp
    long *scoreLeftL,    // The score for an deletion
    long *scoreOnL       // Score to update
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-21 TOC: updateDirAndScoreWater
   '  - Picks the best score and direction for the current base pairs
   '    being compared in a Waterman Smith alignment
   '  o fun-21 sec-1: Matches->insertions->deletions
   '  o fun-21 sec-2: Matches->deletions->insertions
   '  o fun-21 sec-3: Insertions->matches->deletions
   '  o fun-21 sec-4: Deletions->matches->insertions
   '  o fun-21 sec-5: Insertions->deletions->matches
   '  o fun-21 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#endif
