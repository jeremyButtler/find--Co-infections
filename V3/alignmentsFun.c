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

#include "alignmentsFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-03 cnvtAlnErrToSeq:
'     o Uses an array of error types and a c-string with the a sequence
'       to make one part of an alignment
'  - fun-04 WatermanSmithAln:
'     o Perform a Waterman Smith alignment on input sequences
'  - fun-05 NeedleManWunschAln:
'     o Perform a Needleman-Wunsch alignment on input sequences
'  - fun-07 checkIfBasesMatch
'     o Check if two bases are the same (includes anonymous bases)
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
| Output: Heap alloacted C-string with alignment for the input sequence
\---------------------------------------------------------------------*/
char * cnvtAlnErrToSeq(
    char *seqCStr,        // c-string with query sequence
    int32_t seqStartI,    // Were alignment starts on query (index 1)
    char queryBl,         // 1: input sequence is query, 0: is reference
    uint8_t *alnErrUCAry, // Holds error types (ends with 0)
    uint32_t lenErrAryUI  // length of alnErrUCAry (for sequence array)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: Sec-1 Sub-1: cnvtAlnQueryErrToSeq
   '  - Uses an array of error types and a c-string with the a sequence
   '    to make one part of an alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *baseCStr = seqCStr + seqStartI - 1;
   char *tmpBaseCStr = 0;
   char *seqAlnCStr = 0;
   uint8_t *errUCPtr = alnErrUCAry;

   seqAlnCStr = malloc(sizeof(char) * lenErrAryUI + 2);

   if(seqAlnCStr == 0) return 0;  // memory allocation error

   tmpBaseCStr = seqAlnCStr;

   while(*errUCPtr != 0)
   { // While I have bases or an alignment to copy over
       switch(*errUCPtr)
       { // Switch; check the error type
           case defDelFlag:                    // deletion
               if(queryBl & 1) *tmpBaseCStr = '-';
               else
               { // Else dealing with a reference sequence
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // Else dealing with a reference sequence

               break;
           case defInsFlag:                    // insertion
               if(!(queryBl & 1)) *tmpBaseCStr = '-';

               else
               { // Else dealing with a query sequence
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // Else dealing with a query sequence

               break;

           case defSoftQueryFlag:
               if(queryBl != 0)
               { // If softmasking a query region
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // If softmasking a query region

               // Else working on the reference sequence
               else *tmpBaseCStr = '-';
               break;


           case defSoftRefFlag:   // Soft masked a reference base
               if(queryBl == 0)
               { // Else if It is a reference base I am soft masking
                   *tmpBaseCStr = *baseCStr;
                   ++baseCStr;
               } // Else if It is a reference base I am soft masking

               // Else working on a query sequence
               else *tmpBaseCStr = '-';
               break;

           // Is bot a query and reference soft mask
           case defSoftQueryFlag + defSoftRefFlag:
               *tmpBaseCStr = *baseCStr;
               ++baseCStr;
               break;

           case defBaseFlag:                    // (match/snp)
               *tmpBaseCStr = *baseCStr;
               ++baseCStr;
               break;
       } // Switch; check the error type

       // Move to the next error type
       ++tmpBaseCStr;
       ++errUCPtr;
   } // While I have bases or an alignment to copy over

  *tmpBaseCStr = '\0'; // Make into a c-string
  return seqAlnCStr;
} // cnvtAlnQueryErrToSeq

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: WatermanSmithAln
   '  - Perform a Waterman Smith alignment on input sequences
   '  o fun-04 sec-1: Variable declerations
   '  o fun-04 sec-2: Allocate memory for alignment
   '  o fun-04 sec-3: Fill in the initial negatives for the reference
   '  o fun-04 sec-4: Fill the matrix with scores
   '  o fun-04 sec-5: Find the best path
   '  o fun-04 sec-6: Clean up and add softmasking to start
   '  o fun-04 sec-7: Mark end of alignment array & clean up end indels
   '  o fun-04 sec-8: Invert the error type array
   '  o fun-04 sec-9: Add softmasking to the end and return the array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-1: Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char swapBuffBl = 1;   // Swap buffers when finshed every 2nd row

   char *refStartCStr = refCStr + refStartI - 1;
   char *queryStartCStr = queryCStr + refStartI - 1;

   char *tmpQueryCStr = queryStartCStr;
   char *tmpRefCStr = refStartCStr;
   char *bestQueryCStr = 0;  // Best end start for the query sequence
   char *bestRefCStr = 0;    // Best end for the reference sequence

   // Find the length of the reference and query
   unsigned long lenQueryUL = queryEndI - queryStartI + 1;
       // +1 to convert back to 1 index (subtraction makes 0)
   unsigned long lenRefUL = refEndI - refStartI + 1;
       // +1 to convert back to 1 index (subtraction makes 0)
   unsigned long lenMatrixUL = 0;

   long snpScoreL = 0;     // Score for single base pair
   long scoreTopL = 0;     // Score when using the top cell
   long scoreDiagnolL = 0; // Score when using the diagnol cell
   long scoreLeftL = 0;    // Score when using the left cell
   long bestScoreL = 0;    // Best score

   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLPtr = 0;  // Score I am currently working on
   long *lastBaseLPtr = 0; // Pointer to cell with last base

   // Two bit array variables
   uint8_t *dirMatrixUC = 0;  // Directions for each score cell
   uint8_t bitElmUC = 0;      // The value of a single bit element
   uint8_t lastBitElmUC = 0;

   uint8_t *dirOnUCPtr = 0;   // Score working on
   uint8_t bitUC = 0;  // Keep track of if need to change direction elm

   uint8_t *topDirUCPtr = 0; // Direction of last score in the matrix
   uint8_t topBitUC = 0;     // Element on for last direction

   uint8_t *leftDirUCPtr = 0;
   uint8_t leftBitUC = 0;
   uint8_t tmpLeftBitUC = 0;

   // For a more through system I would recored every best score
   uint8_t *bestScoreUCPtr = 0; // Start of path to trace back
   uint8_t bestBitUC = 0;       // Start of path to trace back

   uint8_t *alnErrAryUC = 0;   // Error type array (match/snp, ins, del)
   uint32_t numErrUI = 0;      // Number of error detected
   uint8_t *startUCPtr = 0;    // For inverting the alignment array
   uint8_t *endUCPtr = 0;     // For inverting the alignment array
   uint8_t swapUC = 0;         // For swaping elements in the array

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-2: Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   scoreMatrixL = calloc(2 * (lenRefUL + 1), sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion zero query column
   if(scoreMatrixL == 0) return 0;

   // This gives me an array of 0's (so no extra clearing)
   // This array is used for finding the directions
   lenMatrixUL = 1 + (((lenRefUL + 1) * (lenQueryUL + 1)) >> 2);
   dirMatrixUC = calloc(lenMatrixUL, sizeof(uint8_t));
       // Make the direction array for the scoring array
       // 1 + is to make sure have enough chars
       // lenRefUL + 1 is to account for the insertion  reference row
       // lenQeurI + 1 is to account for insertion zero query column
       // x >> 2 = x / 4 & accounts for taking 2 bits (char = 8 bits)
          // per cell

   if(dirMatrixUC == 0)
   { // If I do not have a direction matrix for each cell
       free(scoreMatrixL);
       return 0;
   } // If I do not have a direction matrix for each cell

   *lenErrAryUI = lenQueryUL + lenRefUL;
   alnErrAryUC = malloc(sizeof(uint8_t) * *lenErrAryUI);

   if(alnErrAryUC == 0) 
   { /*If I had a memory allocation error*/
       free(scoreMatrixL);
       free(dirMatrixUC);
       return 0;
   } /*If I had a memory allocation error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-3: Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Calloc already set everything to done (defMoveStop = 0)
   scoreOnLPtr = scoreMatrixL + lenRefUL + 1;
       // Move to the first column in next row
   dirOnUCPtr = dirMatrixUC;
   twoBitAryMoveForXElm(&dirOnUCPtr, &bitUC, lenRefUL + 1);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-4: Fill the matrix with scores
   ^  o fun-04 sec-4 sub-1: Fill in the indel column
   ^  o fun-04 sec-4 sub-2: Get scores for insertion, deletion, match
   ^  o fun-04 sec-4 sub-3: Move to the next refernce/query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-04 Sec-4 Sub-1: Fill in the indel column
   \*******************************************************************/

   // Marks the last cell in the score matrix (just to make sure)
   lastBaseLPtr = scoreMatrixL; // Set up the last base
   swapBuffBl = 1;   // Swap buffers when finsh every 2nd row

   // Direction of the upper cell
   topDirUCPtr = dirMatrixUC;
   topBitUC = 0;

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; compare one query base against all reference bases
       *scoreOnLPtr = 0;
       *dirOnUCPtr |= defMoveStop;  // Set to move to top
       leftDirUCPtr = dirOnUCPtr; // Previous direction
       leftBitUC = bitUC;         // bit of previous direction

       // Move to the first base comparison
       ++scoreOnLPtr;  // Get of negative column for the new query base
       ++lastBaseLPtr; // Get of negative column for last query base

       twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
       twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);
       tmpRefCStr = refStartCStr;

       /***************************************************************\
       * Fun-04 Sec-4 Sub-2: Get scores for insertion, deletion, match
       \***************************************************************/

       // First reference bases column (indel column already handled)
       while(*tmpRefCStr != '\0')
       { // loop; compare one query to one reference base

           snpScoreL =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for the diagnol cell (snp/match)
           scoreDiagnolL = *(lastBaseLPtr - 1) + snpScoreL;

           switch(settings->matchPriorityBl) 
           { // Switch: check if matches have priority

               case 1:
                   if(checkIfBasesMatch(tmpQueryCStr, tmpRefCStr) != 0)
                   { // If had matching bases
                       *dirOnUCPtr |= defMoveDiagnol;
                       *scoreOnLPtr = scoreDiagnolL;
                       break;
                   } // If had matching bases

               case 0:   // Either not using match priority or not match
                   scoreTopL =
                       getIndelScore(
                           topDirUCPtr,
                           &topBitUC,
                           settings,
                           lastBaseLPtr
                   ); // Get the score for an insertion

                   // If the limb is not complete the last direction
                   // will always be one shift back 
                   if(leftBitUC < 3) tmpLeftBitUC = 2;
                   else tmpLeftBitUC = 3;

                   scoreLeftL =
                       getIndelScore(
                           leftDirUCPtr,   // part of two bit index
                           &tmpLeftBitUC,  // part of two bit index
                           settings,       // Has gap penalties
                           scoreOnLPtr - 1 // Score of the previous base
                   ); // Get the score for an insertion

                   updateDirAndScoreWater(
                       dirOnUCPtr,
                       settings,   // Has preference for score selection
                       &scoreTopL,     // Score for an insertion
                       &scoreDiagnolL, // Score for an match/snp
                       &scoreLeftL,    // The score for an deletion
                       scoreOnLPtr
                   ); // Update the scores

                   if(*scoreOnLPtr >= bestScoreL)
                   { // Else if have a new best score
                       bestScoreL = *scoreOnLPtr;
                       bestScoreUCPtr = dirOnUCPtr;
                       bestBitUC = bitUC;
                       bestQueryCStr = tmpQueryCStr;
                       bestRefCStr = tmpRefCStr;
                   } // Else if have a new best score
           } // Switch: check if matches have priority

           /************************************************************\
           * Fun-04 Sec-4 Sub-3: Move to the next refernce/query base
           \************************************************************/
       
           // Move to the next cell to score
           ++scoreOnLPtr; // Move to next comparison for this query base
           ++lastBaseLPtr; // Move to next element
           ++tmpRefCStr;   // Move to the next reference base

           leftDirUCPtr = dirOnUCPtr; // Previous direction
           leftBitUC = bitUC;         // bit of previous direction

           twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
           twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);
       } // loop; compare one query to one reference base

       // Will end on the corrnor cell
       *scoreL = *(scoreOnLPtr - 1);

       if(swapBuffBl & 1)
       { // If need to swap the buffers
           scoreOnLPtr = scoreMatrixL; // Restarting scoring
           swapBuffBl = 0;
       } // If need to swap the buffers

       else
       { // Else need to reset the score part
           lastBaseLPtr = scoreMatrixL; // Last base is on first row
           swapBuffBl = 1;
       } // Else need to reset the score part

       ++tmpQueryCStr; // Move to the next query base
   } // loop; compare one query base against all reference bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-5: Find the best path
   ^   o fun-04 sec-5 sub-1: Get to the very last score (bottom right)
   ^   o fun-04 sec-5 sub-2: Find the best path
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-04 Sec-5 Sub-1: Get to the very last score (bottom right)
   \*******************************************************************/

   // Make sure the last index is correct
   finishTwoBitAry(dirOnUCPtr, bitUC);

   // Move to the last direction (currently on direction after last)
   dirOnUCPtr = bestScoreUCPtr;
   bitUC = bestBitUC;
   bitElmUC = getTwoBitAryElm(dirOnUCPtr, &bitUC);
   // Move to bottom right base to trace back the path

   /*******************************************************************\
   * Fun-04 Sec-5 Sub-2: Find the best path
   \*******************************************************************/

   while(bitElmUC != defMoveStop) // While not at the ned of the matrix
   { // While I have more bases to add to the path
       lastBitElmUC = bitElmUC; // Keep track of the last bit

       switch(bitElmUC)
       { // Switch: check what the next base is in the sequence
           case defMoveUp:                    // Move to top (insertion)
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 1);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 1 to get to cell above current
               *(alnErrAryUC + numErrUI) = defInsFlag;  // flag for ins
               --bestQueryCStr;
               break;

           case defMoveDiagnol:           // Move to diagnol (match/snp)
               *(alnErrAryUC + numErrUI) = defBaseFlag; // snp
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 2);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 2 to get to diagnol cell
               --bestQueryCStr;
               --bestRefCStr;
               break;

           case defMoveLeft:              // Move to left (deletion)
               twoBitAryMoveBackOneElm(&dirOnUCPtr, &bitUC);
               *(alnErrAryUC + numErrUI) = defDelFlag;  // flag for del
               --bestRefCStr;
               break;
       } // Switch: check what the next base is in the sequence

       // Move back to the selected base

       // Get the next direction to move
       bitElmUC = getTwoBitAryElm(dirOnUCPtr, &bitUC);

       ++numErrUI; // Account for the added error
   } // While I have more bases to add to the path

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-6: Clean up and add softmasking to start
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(scoreMatrixL);               // No longer need
   free(dirMatrixUC);

   switch(lastBitElmUC)
   { // Switch; check which sequence I am on the off base
       case defMoveDiagnol:   // Both sequences are one bit off
           ++bestQueryCStr;
           ++bestRefCStr;
           break;

       case defMoveLeft:   // Last move was a deletion, ref one base off
           ++bestRefCStr;
           break;

       case defMoveUp:   // Last base was insertion, query one base off
           ++bestQueryCStr;
   } // Switch; check which sequence I am on the off base

   endUCPtr = alnErrAryUC;

   tmpQueryCStr = queryCStr;
   tmpRefCStr = refCStr;
   endUCPtr = (alnErrAryUC + numErrUI);
   lenRefUL = 0; // So can softmask the end bases
   lenQueryUL = 0; // So can softmask the end bases

   // Apply softmasking to the start region
   while(1)
   { // While I have softmasking to do
       if(tmpQueryCStr != bestQueryCStr)
       { // IF have a query base to mask
           *endUCPtr |= defSoftQueryFlag;
           ++tmpQueryCStr;
           ++lenQueryUL;

           if(tmpRefCStr != bestRefCStr)
           { // IF have a query base to mask
               *endUCPtr |= defSoftRefFlag;
               ++tmpRefCStr;
               ++lenRefUL;
           } // IF have a query base to mask

           ++endUCPtr;
           ++numErrUI;
       } // IF have a query base to mask

       else if(tmpRefCStr != bestRefCStr)
       { // IF have a query base to mask
           *endUCPtr |= defSoftRefFlag;
           ++tmpRefCStr;
           ++numErrUI;
           ++endUCPtr;
           ++lenRefUL;
       } // IF have a query base to mask

       else break; // Done
   } // While I have softmasking to do

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-7: Mark end of alignment array & clean up indels at end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // This only happens when the gapopen penalty is 0
   // Mark the end of the error array
   *(alnErrAryUC + numErrUI) = 0;
   startUCPtr = alnErrAryUC;

   // Remove any hangin indels at the end (no matches after indel)
   while(!(*startUCPtr & 4))
   { // While I have no bases
       switch(*startUCPtr)
       { // Switch: Check wich soft mask I apply
           case defDelFlag:
               *startUCPtr = defSoftRefFlag;
               break;

           case defInsFlag:
               *startUCPtr = defSoftQueryFlag;
               break;
       } // Switch: Check wich soft mask I apply

       if(*startUCPtr == 0) break; // At end of the array

       ++startUCPtr;
   } // While I have no bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-8: Invert the error type array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // flag the end of the error array
   startUCPtr = alnErrAryUC;
   endUCPtr = (alnErrAryUC + numErrUI - 1);

   while(startUCPtr < endUCPtr)
   { // While I have elements to inver
       switch(*endUCPtr)
       { // Switch; Check if I need to count a base
           case defBaseFlag:
               ++lenQueryUL;
               ++lenRefUL;
               break;

           case defInsFlag:
               ++lenQueryUL;
               break;

           case defDelFlag:
               ++lenRefUL;
               break;
       } // Switch; Check if I need to count a base

       switch(*startUCPtr)
       { // Switch; Check if I need to count a base
           case defBaseFlag:
               ++lenQueryUL;
               ++lenRefUL;
               break;

           case defInsFlag:
               ++lenQueryUL;
               break;

           case defDelFlag:
               ++lenRefUL;
               break;
       } // Switch; Check if I need to count a base

       swapUC = *endUCPtr;
       *endUCPtr = *startUCPtr;
       *startUCPtr = swapUC;
       ++startUCPtr;
       --endUCPtr;
   } // While I have elements to inver

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-9: Add softmasking to the end and return the array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Add softmasking to the end
   tmpQueryCStr = queryCStr + lenQueryUL + 1;
   tmpRefCStr = refCStr + lenRefUL + 1;
   endUCPtr = alnErrAryUC + numErrUI;

   // Apply softmasking to the start region
   while(1)
   { // While I have softmasking to do
       if(*tmpQueryCStr != '\0')
       { // IF have a query base to mask
           *endUCPtr |= defSoftQueryFlag;
           ++tmpQueryCStr;

           if(*tmpRefCStr != '\0')
           { // IF have a query base to mask
               *endUCPtr |= defSoftRefFlag;
               ++tmpRefCStr;
           } // IF have a query base to mask

           ++endUCPtr;
       } // IF have a query base to mask

       else if(*tmpRefCStr != '\0')
       { // IF have a query base to mask
           *endUCPtr |= defSoftRefFlag;
           ++tmpRefCStr;
           ++endUCPtr;
       } // IF have a query base to mask

       else break; // Done
   } // While I have softmasking to do

   endUCPtr = 0; // Add end of sequence

   return alnErrAryUC;
} // WatermanSmithAln

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
    long *scoreL        // Score for the alignment (bottom right)
    // *startI and *endI paramaters should be index 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: NeedleManWunschAln
   '  - Perform a Needleman-Wunsch alignment on input sequences
   '  o fun-05 sec-1: Variable declerations
   '  o fun-05 sec-2: Allocate memory for alignment
   '  o fun-05 sec-3: Fill in the initial negatives for the reference
   '  o fun-05 sec-4: Fill the matrix with scores
   '  o fun-05 sec-5: Find the best path
   '  o fun-05 sec-6: Clean up and invert error array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-1: Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char firstRoundBl = 1;
   char swapBuffBl = 1;   // Swap buffers when finshed every 2nd row

   char *refStartCStr = refCStr + refStartI - 1;
   char *queryStartCStr = queryCStr + refStartI - 1;

   char *tmpQueryCStr = queryStartCStr;
   char *tmpRefCStr = refStartCStr;

   // Find the length of the reference and query
   unsigned long lenQueryUL = queryEndI - queryStartI + 1;
       // +1 to convert back to 1 index (subtraction makes 0)
   unsigned long lenRefUL = refEndI - refStartI + 1;
   unsigned long lenMatrixUL = 0;

       // +1 to convert back to 1 index (subtraction makes 0)

   long snpScoreL = 0;     // Score for single base pair
   long scoreTopL = 0;     // Score when using the top cell
   long scoreDiagnolL = 0; // Score when using the diagnol cell
   long scoreLeftL = 0;    // Score when using the left cell

   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLPtr = 0;  // Score I am currently working on
   long *lastBaseLPtr = 0; // Pointer to cell with last base

   uint8_t *dirMatrixUC = 0;  // Directions for each score cell
   uint8_t bitElmUC = 0;      // The value of a single bit element
   uint8_t *dirOnUCPtr = 0;   // Score working on
   uint8_t bitUC = 0;  // Keep track of if need to change direction elm

   uint8_t *topDirUCPtr = 0; // Direction of last score in the matrix
   uint8_t topBitUC = 0;     // Element on for last direction

   uint8_t *leftDirUCPtr = 0;
   uint8_t leftBitUC = 0;
   uint8_t tmpLeftBitUC = 0;
     // Every 2 bits tells were to move next 11=top,10=diagnol,01=left
     // Find next cell: top: scoreOnLPtr - lenRefUL;
     // Find next cell: diagnol: scoreOnLPtr - lenRefUL - 1;
     // Find next cell: left: --scoreOnLPtr;

   uint8_t *alnErrAryUC = 0;   // Error type array (match/snp, ins, del)
   uint32_t numErrUI = 0;      // Number of error detected
   uint8_t *startUCPtr = 0;    // For inverting the alignment array
   uint8_t *endUCPtr = 0;     // For inverting the alignment array
   uint8_t swapUC = 0;         // For swaping elements in the array

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-2: Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   scoreMatrixL = calloc(2 * (lenRefUL + 1), sizeof(unsigned long));
       // I need two rows to keep track of the scores (2x)
       // lenRefUL + 1 is to account for insertion zero query column
   if(scoreMatrixL == 0) return 0;

   // This gives me an array of 0's (so no extra clearing)
   // This array is used for finding the directions
   lenMatrixUL = 1 + (((lenRefUL + 1) * (lenQueryUL + 1)) >> 2);
   dirMatrixUC = calloc(lenMatrixUL, sizeof(uint8_t));
       // Make the direction array for the scoring array
       // 1 + is to make sure have enough chars
       // lenRefUL + 1 is to account for the insertion  reference row
       // lenQeurI + 1 is to account for insertion zero query column
       // x >> 2 = x / 4 & accounts for taking 2 bits (char = 8 bits)
          // per cell

   if(dirMatrixUC == 0)
   { // If I do not have a direction matrix for each cell
       free(scoreMatrixL);
       return 0;
   } // If I do not have a direction matrix for each cell

   *lenErrAryUI = lenQueryUL + lenRefUL;
   alnErrAryUC = malloc(sizeof(uint8_t) * *lenErrAryUI);

   if(alnErrAryUC == 0) 
   { /*If I had a memory allocation error*/
       free(scoreMatrixL);
       free(dirMatrixUC);
       return 0;
   } /*If I had a memory allocation error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-3: Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Build up the indels for the reference row
   scoreOnLPtr = scoreMatrixL;
   *scoreOnLPtr = 0;         // Top left cell
   ++scoreOnLPtr;            // Move to first reference base cell
   *scoreOnLPtr = settings->gapStartPenaltyI;
       // Second column of the first row holds the first indel
   ++scoreOnLPtr;            // Move to first reference base cell

   // Get direction matrix in sync with scoring matrix
   dirOnUCPtr = dirMatrixUC; // Move left
   *dirOnUCPtr |= defMoveStop;
   twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);

   *dirOnUCPtr |= defMoveLeft;
   twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);

   tmpRefCStr = refStartCStr + 1; // Move off the first base (is done)

   // Already filled in two cells in this row so, it is lenRef - 1
   while(*tmpRefCStr != '\0')
   { // loop; till have initalized the first row
       *scoreOnLPtr = *(scoreOnLPtr - 1) + settings->gapExtendPenaltyI;
       *dirOnUCPtr |= defMoveLeft;

       ++scoreOnLPtr; // Move to the next element
       // Move past elements in the bit array
       // Not worried about specifiying left shift (0), since calloc
       // already set everything to 0.

       twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
       ++tmpRefCStr;
   } // loop; till have initalized the first row

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-4: Fill the matrix with scores
   ^  o fun-05 sec-4 sub-1: Fill in the indel column
   ^  o fun-05 sec-4 sub-2: Get scores for insertion, deletion, match
   ^  o fun-05 sec-4 sub-3: Move to the next refernce/query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-05 Sec-4 Sub-1: Fill in the indel column
   \*******************************************************************/

   // Marks the last cell in the score matrix (just to make sure)
   lastBaseLPtr = scoreMatrixL; // Set up the last base
   firstRoundBl = 1; // Mark need to add a gap start penalty at start
   swapBuffBl = 1;   // Swap buffers when finsh every 2nd row

   // Direction of the upper cell
   topDirUCPtr = dirMatrixUC;
   topBitUC = 0;

   tmpQueryCStr = queryStartCStr;
   tmpRefCStr = refStartCStr;

   // Starting on the first sequence row
   while(*tmpQueryCStr != '\0')
   { // loop; compare one query base against all reference bases

       switch(firstRoundBl)
       { // Switch, check if on the 1st cell of the 2nd row
           case 1:
             *scoreOnLPtr = *lastBaseLPtr + settings->gapStartPenaltyI;
             firstRoundBl = 0; // Never will be reset
             break;

           case 0:
             *scoreOnLPtr = *lastBaseLPtr + settings->gapExtendPenaltyI;
             break;
           // Else is the first indel for the indel column
       } // Switch, check if on the 1st cell of the 2nd row
 
       *dirOnUCPtr |= defMoveUp;  // Set to move to top
       leftDirUCPtr = dirOnUCPtr; // Previous direction
       leftBitUC = bitUC;         // bit of previous direction

       // Move to the first base comparison
       ++scoreOnLPtr;  // Get of negative column for the new query base
       ++lastBaseLPtr; // Get of negative column for last query base

       twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
       twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);

       tmpRefCStr = refStartCStr;

       /***************************************************************\
       * Fun-05 Sec-4 Sub-2: Get scores for insertion, deletion, match
       \***************************************************************/

       // First reference bases column (indel column already handled)
       while(*tmpRefCStr != '\0')
       { // loop; compare one query to one reference base

           snpScoreL =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for the diagnol cell (snp/match)
           scoreDiagnolL = *(lastBaseLPtr - 1) + snpScoreL;

           switch(settings->matchPriorityBl) 
           { // Switch: check if matches have priority

               case 1:
                   if(checkIfBasesMatch(tmpQueryCStr, tmpRefCStr) != 0)
                   { // If had matching bases
                       *dirOnUCPtr |= defMoveDiagnol;
                       *scoreOnLPtr = scoreDiagnolL;
                       break;
                   } // If had matching bases

               case 0:   // Either not using match priority or not match
                   scoreTopL =
                       getIndelScore(
                           topDirUCPtr,
                           &topBitUC,
                           settings,
                           lastBaseLPtr
                   ); // Get the score for an insertion

                   // If the limb is not complete the last direction
                   // will always be one shift back 
                   if(leftBitUC < 3) tmpLeftBitUC = 2;
                   else tmpLeftBitUC = 3;

                   scoreLeftL =
                       getIndelScore(
                           leftDirUCPtr,   // part of two bit index
                           &tmpLeftBitUC,  // part of two bit index
                           settings,       // Has gap penalties
                           scoreOnLPtr - 1 // Score of the previous base
                   ); // Get the score for an insertion

                   updateDirAndScoreNeedle(
                       dirOnUCPtr,
                       settings,   // Has preference for score selection
                       &scoreTopL,     // Score for an insertion
                       &scoreDiagnolL, // Score for an match/snp
                       &scoreLeftL,    // The score for an deletion
                       scoreOnLPtr
                   ); // Update the scores
           } // Switch: check if matches have priority

           /************************************************************\
           * Fun-05 Sec-4 Sub-3: Move to the next refernce/query base
           \************************************************************/
       
           // Move to the next cell to score
           ++scoreOnLPtr; // Move to next comparison for this query base
           ++lastBaseLPtr; // Move to next element
           ++tmpRefCStr;   // Move to the next reference base

           leftDirUCPtr = dirOnUCPtr; // Previous direction
           leftBitUC = bitUC;         // bit of previous direction

           twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
           twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);
       } // loop; compare one query to one reference base

       // Will end on the corrnor cell
       *scoreL = *(scoreOnLPtr - 1);

       if(swapBuffBl & 1)
       { // If need to swap the buffers
           scoreOnLPtr = scoreMatrixL; // Restarting scoring
           swapBuffBl = 0;
       } // If need to swap the buffers

       else
       { // Else need to reset the score part
           lastBaseLPtr = scoreMatrixL; // Last base is on first row
           swapBuffBl = 1;
       } // Else need to reset the score part

       ++tmpQueryCStr; // Move to the next query base
   } // loop; compare one query base against all reference bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-5: Find the best path
   ^   o fun-05 sec-5 sub-1: Get to the very last score (bottom right)
   ^   o fun-05 sec-5 sub-2: Find the best path
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-05 Sec-5 Sub-1: Get to the very last score (bottom right)
   \*******************************************************************/

   // Make sure the last index is correct
   finishTwoBitAry(dirOnUCPtr, bitUC);

   // Move to the last direction (currently on direction after last)
   twoBitAryMoveBackOneElm(&dirOnUCPtr, &bitUC);
   bitElmUC = getTwoBitAryElm(dirOnUCPtr, &bitUC);
   // Move to bottom right base to trace back the path

   /*******************************************************************\
   * Fun-05 Sec-5 Sub-2: Find the best path
   \*******************************************************************/

   while(bitElmUC != defMoveStop) // While not at the ned of the matrix
   { // While I have more bases to add to the path
       switch(bitElmUC)
       { // Switch: check what the next base is in the sequence
           case defMoveUp:                    // Move to top (insertion)
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 1);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 1 to get to cell above current
               *(alnErrAryUC + numErrUI) = defInsFlag;  // flag for ins
               break;

           case defMoveDiagnol:           // Move to diagnol (match/snp)
               *(alnErrAryUC + numErrUI) = defBaseFlag; // snp
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 2);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 2 to get to diagnol cell
               break;

           case defMoveLeft:              // Move to left (deletion)
               twoBitAryMoveBackOneElm(&dirOnUCPtr, &bitUC);
               *(alnErrAryUC + numErrUI) = defDelFlag;  // flag for del
               break;
       } // Switch: check what the next base is in the sequence

       // Move back to the selected base

       // Get the next direction to move
       bitElmUC = getTwoBitAryElm(dirOnUCPtr, &bitUC);

       ++numErrUI; // Account for the added error
   } // While I have more bases to add to the path

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-6: Clean up and invert error array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(scoreMatrixL);               // No longer need
   free(dirMatrixUC);

   // flag the end of the error array
   *(alnErrAryUC + numErrUI) = 0;
       // Not the last direction (last score) is the starting 0 cell

   startUCPtr = alnErrAryUC;
   endUCPtr = (alnErrAryUC + numErrUI - 1);

   while(startUCPtr < endUCPtr)
   { // While I have elements to inver
       swapUC = *endUCPtr;
       *endUCPtr = *startUCPtr;
       *startUCPtr = swapUC;
       ++startUCPtr;
       --endUCPtr;
   } // While I have elements to inver

   return alnErrAryUC;
} // NeeldeManWunschAln

/*---------------------------------------------------------------------\
| Output: 1: if bases were a match; 0 if not
\---------------------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC, // Query base to see if same as reference
    char *refBaseC    // Reference base to see if same as query
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Check if two bases are the same (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
   // The switch should default to a look up table & will be more clear
   switch(*queryBaseC & defToUper)
   { // Switch: Check if bases are the same
       case 'A':
       // Case: Query is an A
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 1;
               case 'W': return 1;
               case 'M': return 1;
               case 'R': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an A

       case 'T':
       // Case: Query is an T
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an T

       case 'U':
       // Case: Query is an U
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 1;
               case 'U': return 1;
               case 'W': return 1;
               case 'K': return 1;
               case 'B': return 1;
               case 'Y': return 1;
               case 'D': return 1;
               case 'H': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an U

       case 'C':
       // Case: Query is an C
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 1;
               case 'S': return 1;
               case 'M': return 1;
               case 'Y': return 1;
               case 'B': return 1;
               case 'H': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case 'G':
       // Case: Query is an G
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 1;
               case 'S': return 1;
               case 'K': return 1;
               case 'R': return 1;
               case 'B': return 1;
               case 'D': return 1;
               case 'V': return 1;
               case 'N': return 1;
               case 'X': return 1;
               default: return 0;
           } // Switch: Check what the reference base was
       // Case: Query is an G

       case 'W':
       // Case: Query is an W
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               case 'G': return 0;
               case 'S': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an W

       case 'S':
       // Case: Query is an S
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'W': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an S

       case 'M':
       // Case: Query is an M
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'K': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an M

       case 'K':
       // Case: Query is an K
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'C': return 0;
               case 'M': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an K

       case 'R':
       // Case: Query is an R
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               case 'T': return 0;
               case 'U': return 0;
               case 'Y': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an R

       case 'Y':
       // Case: Query is an Y
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               case 'G': return 0;
               case 'R': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an Y

       case 'B':
       // Case: Query is an B
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'A': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an B

       case 'D':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'C': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an C

       case 'H':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'G': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an H

       case 'V':
       // Case: Query is an D
           switch(*refBaseC & defToUper)
           { // Switch: Check what the reference base was
               case 'T': return 0;
               case 'U': return 0;
               default: return 1;
           } // Switch: Check what the reference base was
       // Case: Query is an V

       case 'N': return 1;
       case 'X': return 1;
   } // Switch: Check if bases are the same

   return 0; // not a base
} // checkIfBasesMatch

/*---------------------------------------------------------------------\
| Output: Returns: Initalized alnSet structer or 0 for memory error
| Note: I prefer using stack, so this is more here for someone else
\---------------------------------------------------------------------*/
struct alnSet * makeAlnSetST(
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: makeAlnSet
   '  - Makes & initalizes an alnSet structer on the heap
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   struct alnSet *alnSetST = malloc(sizeof(struct alnSet));

   if(alnSetST == 0) return 0;

   initAlnSet(alnSetST);
   return alnSetST;
} // makeAlnSetST
  
/*---------------------------------------------------------------------\
| Output: Returns score of a single pair of bases
\---------------------------------------------------------------------*/
int16_t getBasePairScore(
    const char *queryBaseC, // Query base of pair to get score for
    const char *refBaseC,   // Reference base of pair to get score for
    struct alnSet *alnSetST // structure with scoring matrix to change
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-13 TOC: Sec-1 Sub-1: getBasePairScore
   '  - Get the score for a pair of bases from an alignment structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   return
       alnSetST->snpPenaltyC
           [(uint8_t) (*queryBaseC & defClearNonAlph) - 1]
           [(uint8_t) (*refBaseC & defClearNonAlph) - 1];
} // getBasePairScore

/*---------------------------------------------------------------------\
| Output: Modifies: alnSetST to have default alignment settings values
\---------------------------------------------------------------------*/
void initAlnSet(
    struct alnSet *alnSetST // Alinment settings structure to initialize
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings) structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Kmer mapping settings
   alnSetST->percOverlapD = defPercOverlap;
   alnSetST->maxChunksI = defMaxChunks;
   alnSetST->minChunkLenI = defMinChunkLen;
   alnSetST->minPercKmerD = defMinPercKmer;
   alnSetST->useNeedleBl = defUseNeedle;

   // Aligment Needleman-Wunsch settings
   alnSetST->gapStartPenaltyI = defGapStartPenalty;
   alnSetST->gapExtendPenaltyI = defGapExtendPenalty;
   alnSetST->minScoreUI = defMinScore;

   alnSetST->matchPriorityBl = defMatchPriority;
   alnSetST->diagnolPriorityC = defDiagnolPriority;
   alnSetST->topPriorityC = defTopPriority;
   alnSetST->leftPriorityC = defLeftPriority;

   // Aligment Needleman-Wunsch scoring matrix
   for(uint8_t colUC = 0; colUC < 26; ++colUC)
   { // loop for all columns in the comparison matrix
       for(uint8_t rowUC = 0; rowUC < 26; ++rowUC)
           alnSetST->snpPenaltyC[colUC][rowUC] = 0;
           // Most of these cells do not have possible combinations
   } // loop for all columns in the comparison matrix

   // Set up scores for non-anonmyous base pairs
   setBasePairScore("a", "a", defAToA, alnSetST);
   setBasePairScore("a", "t", defAToT, alnSetST);
   setBasePairScore("a", "u", defAToT, alnSetST);
   setBasePairScore("a", "g", defAToG, alnSetST);
   setBasePairScore("a", "c", defAToC, alnSetST);

   setBasePairScore("t", "a", defTToA, alnSetST);
   setBasePairScore("t", "t", defTToT, alnSetST);
   setBasePairScore("t", "g", defTToG, alnSetST);
   setBasePairScore("t", "c", defTToC, alnSetST);

   setBasePairScore("u", "a", defTToA, alnSetST);
   setBasePairScore("u", "u", defTToT, alnSetST);
   setBasePairScore("u", "g", defTToG, alnSetST);
   setBasePairScore("u", "c", defTToC, alnSetST);

   setBasePairScore("g", "a", defGToA, alnSetST);
   setBasePairScore("g", "t", defGToT, alnSetST);
   setBasePairScore("g", "u", defGToT, alnSetST);
   setBasePairScore("g", "g", defGToG, alnSetST);
   setBasePairScore("g", "c", defGToC, alnSetST);

   setBasePairScore("c", "a", defCToA, alnSetST);
   setBasePairScore("c", "t", defCToT, alnSetST);
   setBasePairScore("c", "u", defCToT, alnSetST);
   setBasePairScore("c", "g", defCToG, alnSetST);
   setBasePairScore("c", "c", defCToC, alnSetST);

   // non-anonymous base and anonymous base pairs
   setBasePairScore("a", "w", defAToW, alnSetST);
   setBasePairScore("a", "s", defAToS, alnSetST);
   setBasePairScore("a", "m", defAToM, alnSetST);
   setBasePairScore("a", "k", defAToK, alnSetST);
   setBasePairScore("a", "r", defAToR, alnSetST);
   setBasePairScore("a", "y", defAToY, alnSetST);
   setBasePairScore("a", "b", defAToB, alnSetST);
   setBasePairScore("a", "d", defAToD, alnSetST);
   setBasePairScore("a", "h", defAToH, alnSetST);
   setBasePairScore("a", "v", defAToV, alnSetST);
   setBasePairScore("a", "n", defAToN, alnSetST);
   setBasePairScore("a", "x", defAToX, alnSetST);

   setBasePairScore("w", "a", defWToA, alnSetST);
   setBasePairScore("s", "a", defSToA, alnSetST);
   setBasePairScore("m", "a", defMToA, alnSetST);
   setBasePairScore("k", "a", defKToA, alnSetST);
   setBasePairScore("r", "a", defRToA, alnSetST);
   setBasePairScore("y", "a", defYToA, alnSetST);
   setBasePairScore("b", "a", defBToA, alnSetST);
   setBasePairScore("d", "a", defDToA, alnSetST);
   setBasePairScore("h", "a", defHToA, alnSetST);
   setBasePairScore("v", "a", defVToA, alnSetST);
   setBasePairScore("n", "a", defNToA, alnSetST);
   setBasePairScore("x", "a", defXToA, alnSetST);

   setBasePairScore("c", "w", defCToW, alnSetST);
   setBasePairScore("c", "s", defCToS, alnSetST);
   setBasePairScore("c", "m", defCToM, alnSetST);
   setBasePairScore("c", "k", defCToK, alnSetST);
   setBasePairScore("c", "r", defCToR, alnSetST);
   setBasePairScore("c", "y", defCToY, alnSetST);
   setBasePairScore("c", "b", defCToB, alnSetST);
   setBasePairScore("c", "d", defCToD, alnSetST);
   setBasePairScore("c", "h", defCToH, alnSetST);
   setBasePairScore("c", "v", defCToV, alnSetST);
   setBasePairScore("c", "n", defCToN, alnSetST);
   setBasePairScore("c", "x", defCToX, alnSetST);

   setBasePairScore("w", "C", defWToC, alnSetST);
   setBasePairScore("s", "C", defSToC, alnSetST);
   setBasePairScore("m", "C", defMToC, alnSetST);
   setBasePairScore("k", "C", defKToC, alnSetST);
   setBasePairScore("r", "C", defRToC, alnSetST);
   setBasePairScore("y", "C", defYToC, alnSetST);
   setBasePairScore("b", "C", defBToC, alnSetST);
   setBasePairScore("d", "C", defDToC, alnSetST);
   setBasePairScore("h", "C", defHToC, alnSetST);
   setBasePairScore("v", "C", defVToC, alnSetST);
   setBasePairScore("n", "C", defNToC, alnSetST);
   setBasePairScore("x", "C", defXToC, alnSetST);

   setBasePairScore("g", "w", defGToW, alnSetST);
   setBasePairScore("g", "s", defGToS, alnSetST);
   setBasePairScore("g", "m", defGToM, alnSetST);
   setBasePairScore("g", "k", defGToK, alnSetST);
   setBasePairScore("g", "r", defGToR, alnSetST);
   setBasePairScore("g", "y", defGToY, alnSetST);
   setBasePairScore("g", "b", defGToB, alnSetST);
   setBasePairScore("g", "d", defGToD, alnSetST);
   setBasePairScore("g", "h", defGToH, alnSetST);
   setBasePairScore("g", "v", defGToV, alnSetST);
   setBasePairScore("g", "n", defGToN, alnSetST);
   setBasePairScore("g", "x", defGToX, alnSetST);

   setBasePairScore("w", "g", defWToG, alnSetST);
   setBasePairScore("s", "g", defSToG, alnSetST);
   setBasePairScore("m", "g", defMToG, alnSetST);
   setBasePairScore("k", "g", defKToG, alnSetST);
   setBasePairScore("r", "g", defRToG, alnSetST);
   setBasePairScore("y", "g", defYToG, alnSetST);
   setBasePairScore("b", "g", defBToG, alnSetST);
   setBasePairScore("d", "g", defDToG, alnSetST);
   setBasePairScore("h", "g", defHToG, alnSetST);
   setBasePairScore("v", "g", defVToG, alnSetST);
   setBasePairScore("n", "g", defNToG, alnSetST);
   setBasePairScore("x", "g", defXToG, alnSetST);

   setBasePairScore("t", "w", defTToW, alnSetST);
   setBasePairScore("t", "s", defTToS, alnSetST);
   setBasePairScore("t", "m", defTToM, alnSetST);
   setBasePairScore("t", "k", defTToK, alnSetST);
   setBasePairScore("t", "r", defTToR, alnSetST);
   setBasePairScore("t", "y", defTToY, alnSetST);
   setBasePairScore("t", "b", defTToB, alnSetST);
   setBasePairScore("t", "d", defTToD, alnSetST);
   setBasePairScore("t", "h", defTToH, alnSetST);
   setBasePairScore("t", "v", defTToV, alnSetST);
   setBasePairScore("t", "n", defTToN, alnSetST);
   setBasePairScore("t", "x", defTToX, alnSetST);

   setBasePairScore("w", "t", defWToT, alnSetST);
   setBasePairScore("s", "t", defSToT, alnSetST);
   setBasePairScore("m", "t", defMToT, alnSetST);
   setBasePairScore("k", "t", defKToT, alnSetST);
   setBasePairScore("r", "t", defRToT, alnSetST);
   setBasePairScore("y", "t", defYToT, alnSetST);
   setBasePairScore("b", "t", defBToT, alnSetST);
   setBasePairScore("d", "t", defDToT, alnSetST);
   setBasePairScore("h", "t", defHToT, alnSetST);
   setBasePairScore("v", "t", defVToT, alnSetST);
   setBasePairScore("n", "t", defNToT, alnSetST);
   setBasePairScore("x", "t", defXToT, alnSetST);

   // Set u and t to same scores (U is RNA version of T)
   setBasePairScore("u", "w", defTToW, alnSetST);
   setBasePairScore("u", "s", defTToS, alnSetST);
   setBasePairScore("u", "m", defTToM, alnSetST);
   setBasePairScore("u", "k", defTToK, alnSetST);
   setBasePairScore("u", "r", defTToR, alnSetST);
   setBasePairScore("u", "y", defTToY, alnSetST);
   setBasePairScore("u", "b", defTToB, alnSetST);
   setBasePairScore("u", "d", defTToD, alnSetST);
   setBasePairScore("u", "h", defTToH, alnSetST);
   setBasePairScore("u", "v", defTToV, alnSetST);
   setBasePairScore("u", "n", defTToN, alnSetST);
   setBasePairScore("u", "x", defTToX, alnSetST);

   setBasePairScore("w", "u", defWToT, alnSetST);
   setBasePairScore("s", "u", defSToT, alnSetST);
   setBasePairScore("m", "u", defMToT, alnSetST);
   setBasePairScore("k", "u", defKToT, alnSetST);
   setBasePairScore("r", "u", defRToT, alnSetST);
   setBasePairScore("y", "u", defYToT, alnSetST);
   setBasePairScore("b", "u", defBToT, alnSetST);
   setBasePairScore("d", "u", defDToT, alnSetST);
   setBasePairScore("h", "u", defHToT, alnSetST);
   setBasePairScore("v", "u", defVToT, alnSetST);
   setBasePairScore("n", "u", defNToT, alnSetST);
   setBasePairScore("x", "u", defXToT, alnSetST);

   // anonymous base and anonymous base pairs
   setBasePairScore("w", "w", defWToW, alnSetST);
   setBasePairScore("w", "s", defWToS, alnSetST);
   setBasePairScore("w", "m", defWToM, alnSetST);
   setBasePairScore("w", "k", defWToK, alnSetST);
   setBasePairScore("w", "r", defWToR, alnSetST);
   setBasePairScore("w", "y", defWToY, alnSetST);
   setBasePairScore("w", "b", defWToB, alnSetST);
   setBasePairScore("w", "d", defWToD, alnSetST);
   setBasePairScore("w", "h", defWToH, alnSetST);
   setBasePairScore("w", "v", defWToV, alnSetST);
   setBasePairScore("w", "n", defWToN, alnSetST);
   setBasePairScore("w", "x", defWToX, alnSetST);

   setBasePairScore("s", "w", defSToW, alnSetST);
   setBasePairScore("s", "s", defSToS, alnSetST);
   setBasePairScore("s", "m", defSToM, alnSetST);
   setBasePairScore("s", "k", defSToK, alnSetST);
   setBasePairScore("s", "r", defSToR, alnSetST);
   setBasePairScore("s", "y", defSToY, alnSetST);
   setBasePairScore("s", "b", defSToB, alnSetST);
   setBasePairScore("s", "d", defSToD, alnSetST);
   setBasePairScore("s", "h", defSToH, alnSetST);
   setBasePairScore("s", "v", defSToV, alnSetST);
   setBasePairScore("s", "n", defSToN, alnSetST);
   setBasePairScore("s", "x", defSToX, alnSetST);

   setBasePairScore("m", "w", defMToW, alnSetST);
   setBasePairScore("m", "s", defMToS, alnSetST);
   setBasePairScore("m", "m", defMToM, alnSetST);
   setBasePairScore("m", "k", defMToK, alnSetST);
   setBasePairScore("m", "r", defMToR, alnSetST);
   setBasePairScore("m", "y", defMToY, alnSetST);
   setBasePairScore("m", "b", defMToB, alnSetST);
   setBasePairScore("m", "d", defMToD, alnSetST);
   setBasePairScore("m", "h", defMToH, alnSetST);
   setBasePairScore("m", "v", defMToV, alnSetST);
   setBasePairScore("m", "n", defMToN, alnSetST);
   setBasePairScore("m", "x", defMToX, alnSetST);

   setBasePairScore("k", "w", defKToW, alnSetST);
   setBasePairScore("k", "s", defKToS, alnSetST);
   setBasePairScore("k", "m", defKToM, alnSetST);
   setBasePairScore("k", "k", defKToK, alnSetST);
   setBasePairScore("k", "r", defKToR, alnSetST);
   setBasePairScore("k", "y", defKToY, alnSetST);
   setBasePairScore("k", "b", defKToB, alnSetST);
   setBasePairScore("k", "d", defKToD, alnSetST);
   setBasePairScore("k", "h", defKToH, alnSetST);
   setBasePairScore("k", "v", defKToV, alnSetST);
   setBasePairScore("k", "n", defKToN, alnSetST);
   setBasePairScore("k", "x", defKToX, alnSetST);

   setBasePairScore("r", "w", defRToW, alnSetST);
   setBasePairScore("r", "s", defRToS, alnSetST);
   setBasePairScore("r", "m", defRToM, alnSetST);
   setBasePairScore("r", "k", defRToK, alnSetST);
   setBasePairScore("r", "r", defRToR, alnSetST);
   setBasePairScore("r", "y", defRToY, alnSetST);
   setBasePairScore("r", "b", defRToB, alnSetST);
   setBasePairScore("r", "d", defRToD, alnSetST);
   setBasePairScore("r", "h", defRToH, alnSetST);
   setBasePairScore("r", "v", defRToV, alnSetST);
   setBasePairScore("r", "n", defRToN, alnSetST);
   setBasePairScore("r", "x", defRToX, alnSetST);

   setBasePairScore("y", "w", defYToW, alnSetST);
   setBasePairScore("y", "s", defYToS, alnSetST);
   setBasePairScore("y", "m", defYToM, alnSetST);
   setBasePairScore("y", "k", defYToK, alnSetST);
   setBasePairScore("y", "r", defYToR, alnSetST);
   setBasePairScore("y", "y", defYToY, alnSetST);
   setBasePairScore("y", "b", defYToB, alnSetST);
   setBasePairScore("y", "d", defYToD, alnSetST);
   setBasePairScore("y", "h", defYToH, alnSetST);
   setBasePairScore("y", "v", defYToV, alnSetST);
   setBasePairScore("y", "n", defYToN, alnSetST);
   setBasePairScore("y", "x", defYToX, alnSetST);

   setBasePairScore("b", "w", defBToW, alnSetST);
   setBasePairScore("b", "s", defBToS, alnSetST);
   setBasePairScore("b", "m", defBToM, alnSetST);
   setBasePairScore("b", "k", defBToK, alnSetST);
   setBasePairScore("b", "r", defBToR, alnSetST);
   setBasePairScore("b", "y", defBToY, alnSetST);
   setBasePairScore("b", "b", defBToB, alnSetST);
   setBasePairScore("b", "d", defBToD, alnSetST);
   setBasePairScore("b", "h", defBToH, alnSetST);
   setBasePairScore("b", "v", defBToV, alnSetST);
   setBasePairScore("b", "n", defBToN, alnSetST);
   setBasePairScore("b", "x", defBToX, alnSetST);

   setBasePairScore("d", "w", defDToW, alnSetST);
   setBasePairScore("d", "s", defDToS, alnSetST);
   setBasePairScore("d", "m", defDToM, alnSetST);
   setBasePairScore("d", "k", defDToK, alnSetST);
   setBasePairScore("d", "r", defDToR, alnSetST);
   setBasePairScore("d", "y", defDToY, alnSetST);
   setBasePairScore("d", "b", defDToB, alnSetST);
   setBasePairScore("d", "d", defDToD, alnSetST);
   setBasePairScore("d", "h", defDToH, alnSetST);
   setBasePairScore("d", "v", defDToV, alnSetST);
   setBasePairScore("d", "n", defDToN, alnSetST);
   setBasePairScore("d", "x", defDToX, alnSetST);

   setBasePairScore("h", "w", defHToW, alnSetST);
   setBasePairScore("h", "s", defHToS, alnSetST);
   setBasePairScore("h", "m", defHToM, alnSetST);
   setBasePairScore("h", "k", defHToK, alnSetST);
   setBasePairScore("h", "r", defHToR, alnSetST);
   setBasePairScore("h", "y", defHToY, alnSetST);
   setBasePairScore("h", "b", defHToB, alnSetST);
   setBasePairScore("h", "d", defHToD, alnSetST);
   setBasePairScore("h", "h", defHToH, alnSetST);
   setBasePairScore("h", "v", defHToV, alnSetST);
   setBasePairScore("h", "n", defHToN, alnSetST);
   setBasePairScore("h", "x", defHToX, alnSetST);

   setBasePairScore("v", "w", defVToW, alnSetST);
   setBasePairScore("v", "s", defVToS, alnSetST);
   setBasePairScore("v", "m", defVToM, alnSetST);
   setBasePairScore("v", "k", defVToK, alnSetST);
   setBasePairScore("v", "r", defVToR, alnSetST);
   setBasePairScore("v", "y", defVToY, alnSetST);
   setBasePairScore("v", "b", defVToB, alnSetST);
   setBasePairScore("v", "d", defVToD, alnSetST);
   setBasePairScore("v", "h", defVToH, alnSetST);
   setBasePairScore("v", "v", defVToV, alnSetST);
   setBasePairScore("v", "n", defVToN, alnSetST);
   setBasePairScore("v", "x", defVToX, alnSetST);

   setBasePairScore("n", "w", defNToW, alnSetST);
   setBasePairScore("n", "s", defNToS, alnSetST);
   setBasePairScore("n", "m", defNToM, alnSetST);
   setBasePairScore("n", "k", defNToK, alnSetST);
   setBasePairScore("n", "r", defNToR, alnSetST);
   setBasePairScore("n", "y", defNToY, alnSetST);
   setBasePairScore("n", "b", defNToB, alnSetST);
   setBasePairScore("n", "d", defNToD, alnSetST);
   setBasePairScore("n", "h", defNToH, alnSetST);
   setBasePairScore("n", "v", defNToV, alnSetST);
   setBasePairScore("n", "n", defNToN, alnSetST);
   setBasePairScore("n", "x", defNToX, alnSetST);

   setBasePairScore("x", "w", defXToW, alnSetST);
   setBasePairScore("x", "s", defXToS, alnSetST);
   setBasePairScore("x", "m", defXToM, alnSetST);
   setBasePairScore("x", "k", defXToK, alnSetST);
   setBasePairScore("x", "r", defXToR, alnSetST);
   setBasePairScore("x", "y", defXToY, alnSetST);
   setBasePairScore("x", "b", defXToB, alnSetST);
   setBasePairScore("x", "d", defXToD, alnSetST);
   setBasePairScore("x", "h", defXToH, alnSetST);
   setBasePairScore("x", "v", defXToV, alnSetST);
   setBasePairScore("x", "n", defXToN, alnSetST);
   setBasePairScore("x", "x", defXToX, alnSetST);

   return;
} // initAlnSet

/*---------------------------------------------------------------------\
| Output: Modifies: one score in an snp/match scoring matrix
\---------------------------------------------------------------------*/
void setBasePairScore(
    const char *queryBaseC,   // Query base to change score for
    const char *refBaseC,     // Reference base to change score for
    int16_t newScoreC,        // New value for [query][ref] combination
    struct alnSet *alnSetST   // structure with scoring matrix to change
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-15 TOC: Sec-1 Sub-1: setBasePairScore
   '  - Changes SNP/Match penalty for one query/reference combination
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   alnSetST->snpPenaltyC
       [(uint8_t) (*queryBaseC & defClearNonAlph) - 1]
       [(uint8_t) (*refBaseC & defClearNonAlph) - 1]
       =
       newScoreC;

   // This is here to account for U is the RNA version of T
   if((*queryBaseC & defToUper) =='T')
       alnSetST->snpPenaltyC
           [(uint8_t) ('U' & defClearNonAlph) - 1]
           [(uint8_t) (*refBaseC & defClearNonAlph) - 1]
           =
           newScoreC;
   else if((*queryBaseC & defToUper) == 'U')
       alnSetST->snpPenaltyC
           [(uint8_t) ('T' & defClearNonAlph) - 1]
           [(uint8_t) (*refBaseC & defClearNonAlph) - 1]
           =
           newScoreC;

   if((*refBaseC & defToUper) =='T')
       alnSetST->snpPenaltyC
           [(uint8_t)(*queryBaseC & defClearNonAlph)-1]
           [(uint8_t) ('U' & defClearNonAlph) - 1]
           =
           newScoreC;
   else if((*refBaseC & defToUper) =='U')
       alnSetST->snpPenaltyC
           [(uint8_t)(*queryBaseC & defClearNonAlph)-1]
           [(uint8_t) ('T' & defClearNonAlph) - 1]
           =
           newScoreC;

   if((*queryBaseC & defToUper) =='T' && (*refBaseC & defToUper) == 'T')
       alnSetST->snpPenaltyC
           [(uint8_t)('U' & defClearNonAlph) - 1]
           [(uint8_t) ('U' & defClearNonAlph) - 1]
           =
           newScoreC;
   else if((*queryBaseC & defToUper)=='U'&&(*refBaseC & defToUper)=='U')
       alnSetST->snpPenaltyC
           [(uint8_t)('T' & defClearNonAlph) - 1]
           [(uint8_t) ('T' & defClearNonAlph) - 1]
           =
           newScoreC;

   return;
} // setBasePairScore

/*---------------------------------------------------------------------\
| Output: Frees the alnSet structure (does not set pointer to 0)
\---------------------------------------------------------------------*/
void freeAlnSet(
    struct alnSet *alnSetST,  // Alignment settings structure to free
    char stackBl              // 1: alnSetSt on stack; 0: on heap
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(stackBl & 1) return; // Nothing to do
   if(alnSetST != 0) free(alnSetST); // No heap variables
   return; 
} // freeAlnSet
    

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-17 TOC: Sec-1 Sub-1: cnvtAlnErrAryToLetter
   '  - Converts an alignment error array from my Needleman-Wunsch
   '    alignment into an array of letters (I = insertion, D = deletion,
   '    = = match, X = snp) [These codes are from the eqx cigar entry]
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    while(*alnErrUCAry != 0)
    { // While have codes to convert
        switch(*alnErrUCAry)
        { // Switch; check the error type
            case defDelFlag:                    // deletion
                *alnErrUCAry = 'D';
                ++refSeqCStr;           // Move to next reference base
                break;
            case defInsFlag:                    // insertion
                *alnErrUCAry = 'I';
                ++querySeqCStr;         // Move to next query base
                break;

            case defSoftQueryFlag:      // Query base was sofmasked
                *alnErrUCAry = 's';     // mark as query masked
                ++querySeqCStr;         // Move to next query base
                break;

            case defSoftRefFlag:        // Reference base was softmasked
                *alnErrUCAry = 'P';     // Mark as reference mask
                ++refSeqCStr;         // Move to next query base
                break;

            case defSoftRefFlag + defSoftQueryFlag: 
                *alnErrUCAry = 'S';    // Mark as both masked
                ++querySeqCStr;         // Move to next query base
                ++refSeqCStr;         // Move to next query base
                break;

            case defBaseFlag:                   // Match or snp
                if(refSeqCStr == 0 || querySeqCStr == 0)
                    *alnErrUCAry = 'X';    // No idea if match or snp
                else if(checkIfBasesMatch(querySeqCStr,refSeqCStr) == 0)
                    *alnErrUCAry = 'X';    // snp
                else *alnErrUCAry = '=';   // match

                ++refSeqCStr;           // Move to next reference base
                ++querySeqCStr;         // Move to next query base
                break;
        } // Switch; check the error type

        ++alnErrUCAry;
    } // While have codes to convert

    return;
} // cnvtAlnErryAryToLetter


unsigned long readInScoreFile(
    struct alnSet *alnSetST,  // structure with scoring matrix to change
    FILE *scoreFILE           // File of scores for a scoring matrix
) { /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-18 TOC: readInScoreFile
    '  - Reads in a file of scores for a scoring matrix
    '  o fun-18 sec-1: Variable declerations and buffer set up
    '  o fun-18 sec-2: Read in line and check if comment
    '  o fun-18 sec-3: Get score, update matrix, & move to next line
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-1: Variable declerations and buffer set up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint16_t lenBuffUS = 1024;
    char buffCStr[lenBuffUS];
    char *tmpCStr = 0;
    int16_t scoreS = 0;

    buffCStr[lenBuffUS - 1] = '\0';
    buffCStr[lenBuffUS - 2] = '\0';

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-18 Sec-2: Read in line and check if comment
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(fgets(buffCStr, 1024, scoreFILE))
    { // While I have scores to read in
        
        if(buffCStr[0] == '/' && buffCStr[1] == '/')
        { // On a comment, move onto the next line
            while(
                buffCStr[lenBuffUS - 2] != '\0' &&
                buffCStr[lenBuffUS - 2] != '\n'
            ) { // While have more buffer to read in
                buffCStr[lenBuffUS - 2] = '\0';
                fgets(buffCStr, 1024, scoreFILE);
            } // While have more buffer to read in

            // Reset the buffer
            buffCStr[lenBuffUS - 2] = '\0';

            continue;
        } // On a comment, move onto the next line

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
        ^ Fun-18 Sec-3: Get score, update matrix, & move to next line
        \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

        if(buffCStr[0] == '\n')
            continue;                        // Blank line

        if(buffCStr[0] < 64 && buffCStr[2] < 64)
            return ftell(scoreFILE);         // Invalid character
        
        tmpCStr = cStrToInt16(&buffCStr[4], &scoreS);
        setBasePairScore(&buffCStr[0], &buffCStr[2], scoreS, alnSetST);

        if(tmpCStr == &buffCStr[3])
            return ftell(scoreFILE);         // No score

        while(
            buffCStr[lenBuffUS - 2] != '\0' &&
            buffCStr[lenBuffUS - 2] != '\n'
        ) { // While have more buffer to read in
            buffCStr[lenBuffUS - 2] = '\0';
            fgets(buffCStr, 1024, scoreFILE);
        } // While have more buffer to read in

        // Reset the buffer
        buffCStr[lenBuffUS - 2] = '\0';
    } // While I have scores to read in

    return 0;
} // readInScoreFile

/*---------------------------------------------------------------------\
| Output: Returns the score for an Needleman Wunsch indel
\---------------------------------------------------------------------*/
long getIndelScore(
    uint8_t *lastDirUC,      // Cell gettign last score from
    uint8_t *lastBitUC,      // two bit element on in lastDirUC
    struct alnSet *alnSetST, // Holds the gap open & extension penalties
    long *lastBaseL          // Has score of the last base
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-19 TOC: Sec-1 Sub-1: getIndelScore
   '  - Gets an indel score for the current cell
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Check if this is the first indel or not
   switch(getTwoBitAryElm(lastDirUC, lastBitUC))
   { // Switch; check if this is the first indel
     case defMoveStop:     // At the end of the matrix
     case defMoveDiagnol:  // Top base was a SNP
         return *lastBaseL + alnSetST->gapStartPenaltyI;

     case defMoveUp:   // Top base was an insertion
     case defMoveLeft: // Top base was an deletion
         return *lastBaseL + alnSetST->gapExtendPenaltyI;
   } // Switch; check if this is the first indel

   return 0; // Somthing went wrong
} // getIndelScore

/*---------------------------------------------------------------------\
| Output: Modifies: scoreOnL and dirOnUC to hold best score & direction
\---------------------------------------------------------------------*/
void updateDirAndScoreNeedle(
    uint8_t *dirOnUCPtr,     // Direction on with first two bits cleared
    struct alnSet *alnSetST, // Has preference for score selection
    long *scoreTopL,     // Score for an insertion
    long *scoreDiagnolL, // Score for an match/snp
    long *scoreLeftL,    // The score for an deletion
    long *scoreOnL       // Score to update
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-20 TOC: updateDirAndScoreNeedle
   '  - Picks the best score and direction for the current base pairs
   '    being compared in a Needleman Wunsch alignment
   '  o fun-20 sec-1: Matches->insertions->deletions
   '  o fun-20 sec-2: Matches->deletions->insertions
   '  o fun-20 sec-3: Insertions->matches->deletions
   '  o fun-20 sec-4: Deletions->matches->insertions
   '  o fun-20 sec-5: Insertions->deletions->matches
   '  o fun-20 sec-6: Deletions->insertions->matches
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // I decide the best path as I score since all equal scores create
   // an equally valid alignment.

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-20 Sec-1: Matches and then insertions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   switch(alnSetST->diagnolPriorityC)
   { // Switch; get an snp/match priority
       case 0:  // Top priority
       // Case; bases or matches are top priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 1:
               // Case: priority is matches/snps and then insertions
                   if(*scoreDiagnolL >= *scoreTopL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-20 Sec-2: Matches->deletions->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is matches/snps and then deletions
                   if(*scoreDiagnolL >= *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions
       // Case; bases or matches are top priority

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-20 Sec-3: Insertions->matches->deletions
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 1:
       // Case; bases or matches are second priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then matches/snps
                   if(*scoreDiagnolL > *scoreTopL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-20 Sec-4: Deletions->matches->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then matches/snps
                   if(*scoreDiagnolL > *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-20 Sec-5: Insertions->deletions->matches
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 2:
       // Case; bases or matches are last priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then deletions
                   if(*scoreTopL >= *scoreLeftL)
                   { // If diagnol beats insertion

                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol is the highest score
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol is the highest score

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions then deletions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-20 Sec-6: Deletions->insertions->matches
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then insertions
                   if(*scoreTopL > *scoreLeftL)
                   { // If insertion beats deletion

                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol beats insertions
                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                       } // Else the insertion is the best score
                   } // If insertion beats deletion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions and then deletions
           } // Priority for insertions
       // Case; bases or matches are last priority
   } // Switch; get an snp/match priority

   return;
} // updateDirAndScoreNeedle

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
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

   // I decide the best path as I score since all equal scores create
   // an equally valid alignment.

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-21 Sec-1: Matches and then insertions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   switch(alnSetST->diagnolPriorityC)
   { // Switch; get an snp/match priority
       case 0:  // Top priority
       // Case; bases or matches are top priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 1:
               // Case: priority is matches/snps and then insertions
                   if(*scoreDiagnolL >= *scoreTopL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           if(*scoreLeftL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                           return;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-21 Sec-2: Matches->deletions->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is matches/snps and then deletions
                   if(*scoreDiagnolL >= *scoreLeftL)
                   { // If diagnol beats deletion
                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                          if(*scoreTopL <= 0)
                          { // If need to make a stopping point
                              *dirOnUCPtr |= defMoveStop;
                              *scoreOnL = 0;
                              return;
                          } // If need to make a stopping point

                          *dirOnUCPtr |= defMoveUp;
                          *scoreOnL = *scoreTopL;
                          return;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions
       // Case; bases or matches are top priority

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-21 Sec-3: Insertions->matches->deletions
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 1:
       // Case; bases or matches are second priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then matches/snps
                   if(*scoreDiagnolL > *scoreTopL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL >= *scoreLeftL)
                       { // If diagnol beats deletions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats deletions

                       else
                       { // Else the deletion is the best score
                           if(*scoreLeftL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveLeft;
                           *scoreOnL = *scoreLeftL;
                           return;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreTopL >= *scoreLeftL)
                   { // Else the insertion is the best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else the insertion is the best score

                   else
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is the best score
                       
                   return;
               // Case: priority is matches/snps and then insertions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-21 Sec-4: Deletions->matches->insertions
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then matches/snps
                   if(*scoreDiagnolL > *scoreLeftL)
                   { // If diagnol beats deletion

                       if(*scoreDiagnolL >= *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else the insertion is the best score
                   } // If diagnol beats deletion

                   else if(*scoreLeftL >= *scoreTopL)
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is the best score

                   else
                   { // Else the insertion is the best score
                      if(*scoreTopL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveUp;
                      *scoreOnL = *scoreTopL;
                      return;
                   } // Else the insertion is the best score
                       
                   return;
               // Case: priority is matches/snps and then deletions
           } // Priority for insertions

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun-21 Sec-5: Insertions->deletions->matches
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       case 2:
       // Case; bases or matches are last priority
           switch(alnSetST->topPriorityC)
           { // Priority for insertions
               case 0:
               // Case: priority is insertions and then deletions
                   if(*scoreTopL >= *scoreLeftL)
                   { // If diagnol beats insertion
                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol is the highest score
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol is the highest score

                       else
                       { // Else the insertion is the best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else the deletion is the best score
                   } // If diagnol beats insertion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      if(*scoreDiagnolL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                      return;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions then deletions

               /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
               ^ Fun-21 Sec-6: Deletions->insertions->matches
               \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

               case 2:
               // Case: priority is deletions and then insertions
                   if(*scoreTopL > *scoreLeftL)
                   { // If insertion beats deletion
                       if(*scoreDiagnolL > *scoreTopL)
                       { // If diagnol beats insertions
                           if(*scoreDiagnolL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveDiagnol;
                           *scoreOnL = *scoreDiagnolL;
                           return;
                       } // If diagnol beats insertions

                       else
                       { // Else the insertion is the best score
                           if(*scoreTopL <= 0)
                           { // If need to make a stopping point
                               *dirOnUCPtr |= defMoveStop;
                               *scoreOnL = 0;
                               return;
                           } // If need to make a stopping point

                           *dirOnUCPtr |= defMoveUp;
                           *scoreOnL = *scoreTopL;
                           return;
                       } // Else the insertion is the best score
                   } // If insertion beats deletion

                   else if(*scoreLeftL >= *scoreDiagnolL)
                   { // Else the deletion is the best score
                      if(*scoreLeftL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveLeft;
                      *scoreOnL = *scoreLeftL;
                      return;
                   } // Else the deletion is the best score

                   else
                   { // Else the match/snp is the best score
                      if(*scoreDiagnolL <= 0)
                      { // If need to make a stopping point
                          *dirOnUCPtr |= defMoveStop;
                          *scoreOnL = 0;
                          return;
                      } // If need to make a stopping point

                      *dirOnUCPtr |= defMoveDiagnol;
                      *scoreOnL = *scoreDiagnolL;
                      return;
                   } // Else the match/snp is the best score
                       
                   return;
               // Case: priority is insertions and then deletions
           } // Priority for insertions
       // Case; bases or matches are last priority
   } // Switch; get an snp/match priority

   return;
} // updateDirAndScoreWater
