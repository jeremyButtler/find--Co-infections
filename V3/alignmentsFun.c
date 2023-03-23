/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
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
'  - fun-09 twoBitAryMoveBackOneElm:
'     o Moves back one element in a 2-bit array
'  - fun-10 twoBitAryMoveBackXElm:
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
'  - fun-0? NeedleManWunschAlnCostly:
'     o Perform a Needleman-Wunsch alignment on input sequences this
'       variation uses a scoring array and direction array and so is
        more costly (this was before I realized I did not need to keep
'       scores)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "alignmentsFun.h"

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
   ' Fun-04 TOC: Sec-1 Sub-1: cnvtAlnQueryErrToSeq
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
               if(!(queryBl & 1))
               { // If I am dealing with a reference sequence
                   *tmpBaseCStr = '-';
                   break;
               } // If I am dealing with a reference sequence
           case defBaseFlag:                    // (match/snp)
               *tmpBaseCStr = *baseCStr;
               ++baseCStr;
               break;
           case defMatchFlag:
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
uint8_t * NeedleManWunschAln(
    char *queryCStr,        // Full query sequence as c-string
    int32_t queryStartI,    // Starting query coordinate for alignment
    int32_t queryEndI,      // Ending query coordinate for alignment
    char *refCStr,          // Full reference sequence as c-string
    int32_t refStartI,      // Starting reference coordinate for aln
    int32_t refEndI,        // Ending reference coordinate for alignment
    struct alnSet *settings,// Settings for the alignment
    uint32_t *lenErrAryUI,  // Will hold the return arrays length
    int32_t *scoreI        // Score for the alignment (bottom right)
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
   unsigned long lenMatrixUL = 1+(((lenRefUL+1) * (lenQueryUL+1)) >> 2);

       // +1 to convert back to 1 index (subtraction makes 0)

   int32_t snpScoreI = 0;     // Score for single base pair
   int32_t scoreTopI = 0;     // Score when using the top cell
   int32_t scoreDiagnolI = 0; // Score when using the diagnol cell
   int32_t scoreLeftI = 0;    // Score when using the left cell

   long *scoreMatrixL = 0; // matrix to use in alignment
   long *scoreOnLPtr = 0;  // Score I am currently working on
   long *lastBaseLPtr = 0; // Pointer to cell with last base

   uint32_t shiftByI = 0;      // How much to shift back when aliging

   uint8_t *dirMatrixUC = 0;  // Directions for each score cell
   uint8_t bitElmUC = 0;      // The value of a single bit element
   uint8_t *dirOnUCPtr = 0;   // Score working on
   uint8_t bitUC = 0;  // Keep track of if need to change direction elm

   uint8_t *topDirUCPtr = 0; // Direction of last score in the matrix
   uint8_t topBitUC = 0;     // Element on for last direction

   uint8_t *leftDirUCPtr = 0;
   uint8_t leftBitUC = 0;
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
   *dirOnUCPtr |= defMoveMatch;
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
   ^  o fun-05 sec-4 sub-3: Find direction from scores (best score)
   ^  o fun-05 sec-4 sub-4: Move to the next refernce/query base
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

           case 0: *scoreOnLPtr = settings->gapExtendPenaltyI;
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

           snpScoreI =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for the diagnol cell (snp/match)
           scoreDiagnolI = *(lastBaseLPtr - 1) + snpScoreI;

           // Check if I need to find the other scores
           if(checkIfBasesMatch(tmpQueryCStr, tmpRefCStr) & 1)
           { //If I have a match
               *dirOnUCPtr |= defMoveMatch;
               *scoreOnLPtr = scoreDiagnolI;
           } //If I have a match

           else
           { // else the bases do not match

               // Find the score for the top cell (insertion)
               //bitElmUC = getTwoBitAryElm(topDirUCPtr, &topBitUC);

               if(bitElmUC & 1)     // If is the first indel
                 scoreTopI = *lastBaseLPtr + settings->gapStartPenaltyI;
               else                 // Else is part of a larger indel
                 scoreTopI = *lastBaseLPtr +settings->gapExtendPenaltyI;

               // find the score for the left cell (deletion)
               if(bitUC > 0) bitElmUC = (*dirOnUCPtr & (4 | 8)) >> 2;
               else bitElmUC = getTwoBitAryElm(leftDirUCPtr,&leftBitUC);
                  // If handles when have an incomplete element

               if(bitElmUC & 1)     // If is the first indel
                 scoreLeftI =
                     *(scoreOnLPtr - 1) + settings->gapStartPenaltyI;
               else                 // Else is part of a larger indel
                 scoreLeftI =
                     *(scoreOnLPtr - 1) + settings->gapExtendPenaltyI;

               /********************************************************\
               * Fun-05 Sec-4 Sub-3: Find direction (best score)
               \********************************************************/
               
               // The logic here is that I decide the best path as I
               // score, since the best of equal alternatives will
               // always follow an arbitary decision or report all
               // possible paths (I only care about 1). This allows me
               // to use a smaller direction matrix, but does add in
               // some more time.
               if(scoreTopI > scoreDiagnolI)
               { // If top score is better than the diagnol (or as good)
                   if(scoreTopI >= scoreLeftI)
                   { // If have a better top score
                       *dirOnUCPtr |= defMoveUp;
                       *scoreOnLPtr = scoreTopI;
                   } // If have a better top score

                   else
                   { // Else if have a better left score
                       *dirOnUCPtr |= defMoveLeft;
                       *scoreOnLPtr = scoreLeftI;
                   } // Else if have a better left score
               } // If top score is better than the diagnol (or as good)

               else if(scoreLeftI <= scoreDiagnolI)
               { //else if have a better or as good diagnol score
                   *dirOnUCPtr |= defMoveDiagnol;
                   *scoreOnLPtr = scoreDiagnolI;
               } //else if have a better or as good diagnol score

               else
               { //Else have a better left score
                   *dirOnUCPtr |= defMoveLeft;
                   *scoreOnLPtr = scoreLeftI;
               } //Else have a better left score
           } // If the bases do not match

           /************************************************************\
           * Fun-05 Sec-4 Sub-4: Move to the next refernce/query base
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
       *scoreI = *(scoreOnLPtr - 1);

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

   // Move to bottom right base to trace back the path
   switch(bitUC)
   { // Switch; check if need to move to new direction element
       case 0: // moved to new char, so last score is 1st & 2nd bits
          --dirOnUCPtr;  // Move back to the char for the last score
          bitElmUC = *dirOnUCPtr & (1 | 2);
          break;
       case 1:           // Is the last bit 
       case 2:
       case 3:                                // The first two bitsk
          bitElmUC = *dirOnUCPtr >> 2; //Get off score after final score
          bitElmUC = *dirOnUCPtr & (1 | 2);
          break;
   } // Switch; check if need to move to new direction element

   /*******************************************************************\
   * Fun-05 Sec-5 Sub-2: Find the best path
   \*******************************************************************/

   while(bitUC != 0 || dirOnUCPtr != dirMatrixUC)
   { // While I have more bases to add to the path
       switch(bitElmUC)
       { // Switch: check what the next base is in the sequence
           case defMoveUp:                    // Move to top (insertion)
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 1);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 1 to get to cell above current
               *(alnErrAryUC + numErrUI) = defInsFlag;  // flag for ins
               break;

           case defMoveMatch:
               *(alnErrAryUC + numErrUI) = defMatchFlag; // match
               twoBitAryMoveBackXElm(&dirOnUCPtr, &bitUC, lenRefUL + 2);
                   // have lenRefUL (index 1) cells per row, so need to
                   // do lenRefUL + 2 to get to diagnol cell
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
   *(alnErrAryUC + numErrUI - 1) = 0;
       // Not the last direction (last score) is the starting 0 cell

   startUCPtr = alnErrAryUC;
   endUCPtr = (alnErrAryUC + numErrUI - 2);

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
| Output: Two bits of interest from the two bit array
\---------------------------------------------------------------------*/
uint8_t getTwoBitAryElm(
    uint8_t *twoBitAryUC,// Array of bytes that contain 4 2 bit elements
    uint8_t *charElmOnUC // The two bits on in the unit8_t
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
   '  - Get an element from a two bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    switch(*charElmOnUC)
    { // Switch: find the element on
       case 0: return *twoBitAryUC >> 6;          // Last two bits (7,8)
       case 1: return (*twoBitAryUC & (16 | 32)) >> 4; // 5th & 6th bits
       case 2: return (*twoBitAryUC & (4 | 8)) >> 2;   // 3rd & 4th bits
       case 3: return *twoBitAryUC & (1 | 2);      // The first two bits
    } // Switch: find the element on

    return 0;
} // getTwoBitAryElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: twoBitAryShiftBytsForNewElm
   '  - Make room in a unit8_t for two more two-bit elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(*bitUC)
   { // Switch; check if need to move to new direction element
       case 0:
       case 1:
       case 2:
           ++(*bitUC);          // Still have more bits in char
           **twoBitAryUC = (**twoBitAryUC) << 2;
                // Move to next cell in 2 bit array
           return;
       case 3:
           *bitUC = 0;
           ++(*twoBitAryUC);    // Move to next uint8_t in two bit array
           return;
   } // Switch; check if need to move to new direction element
} // twoBitAryShiftBytsForNewElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move the bit above 
   switch(*bitUC)
   { // Switch; check if need to move to a new element
       case 0:
       case 1:
       case 2:
          ++(*bitUC);      // Move to next two bits in uint8_t
          return;
       case 3:
          *bitUC = 0;        // Moving to a new uint8_t
          ++(*twoBitAryUC);  // Move to the next element
          return;
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveToNextElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: towBitAryMoveBackOneElm
   '  - Moves back one element in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move the bit above 
   switch(*bitUC)
   { // Switch; check if need to move to a new element
       case 0:
          *bitUC = 3;        // One last bit in the previous unit8_t
          --(*twoBitAryUC);  // Move to previous uint8_t in the array
          return;
       case 1:
       case 2:
       case 3:
          --(*bitUC);        // Move to the previous two bits in uint8_t
          return;
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveToNextElm

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: towBitAryMoveBackXElm
   '  - Moves back X elements in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Move back by the indes
   (*twoBitAryUC) -= (shiftByI >> 2);

   // Move the bit above 
   switch(*bitUC)
   { // Switch; check if need to move to a new element
       case 0:
       // Case 0: Current position is on the first bit of the array
           switch(shiftByI & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   *bitUC = 3; // Just starting on the next array
                   --(*twoBitAryUC); // Acount for off bit
                   return;
               case 2:
                   *bitUC = 2;
                   --(*twoBitAryUC); // Acount for off bit
                   return;
               case 3:
                   *bitUC = 1;
                   --(*twoBitAryUC); // Acount for off bit
                   return;
           } // Switch; find out how many bits to adjust
       // Case 0: Current position is on the first bit of the array

       case 1:
       // Case 1: Have one bit to move back by before starting new array
           switch(shiftByI & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   *bitUC = 0; // Had an extra bit to spare
                   return;
               case 2:
                   *bitUC = 3; // On last element of next uint8_t
                   --(*twoBitAryUC); // Acount for off bit
                   return;
               case 3:
                   *bitUC = 2;
                   --(*twoBitAryUC); // Acount for off bit
                   return;
           } // Switch; find out how many bits to adjust
       // Case 1: Have one bit to move back by before starting new array

       case 2:
       // Case 2: Have two bits to move back by before starting new uint
           switch(shiftByI & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   *bitUC = 1; // Had an extra bit to spare
                   return;
               case 2:
                   *bitUC = 0; // On last element of shift uint8_t
                   return;
               case 3:
                   *bitUC = 3;
                   --(*twoBitAryUC); // Acount for off bit
                   return;
           } // Switch; find out how many bits to adjust
       // Case 2: Have two bits to move back by before starting new uint

       case 3:
       // Case 3: was on the last element of the uint8
           switch(shiftByI & (1 | 2))
           { // Switch; find out how many bits to adjust
               case 0: return; // Finished
               case 1:
                   *bitUC = 2;
                   return;
               case 2:
                   *bitUC = 1;
                   return;
               case 3:
                   *bitUC = 0;
                   return;
           } // Switch; find out how many bits to adjust
       // Case 3: was on the last element of the uint8
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveBackXElm

/*---------------------------------------------------------------------\
| Output: Returns: Initalized alnSet structer or 0 for memory error
| Note: I prefer using stack, so this is more here for someone else
\---------------------------------------------------------------------*/
struct alnSet * makeAlnSetST(
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: makeAlnSet
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
int8_t getBasePairScore(
    const char *queryBaseC, // Query base of pair to get score for
    const char *refBaseC,   // Reference base of pair to get score for
    struct alnSet *alnSetST // structure with scoring matrix to change
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-12 TOC: Sec-1 Sub-1: getBasePairScore
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
   ' Fun-13 TOC: Sec-1 Sub-1: initAlnSet
   '  - Set values in altSet (alingment settings) structure to defaults
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   // Kmer mapping settings
   alnSetST->percOverlapD = defPercOverlap;
   alnSetST->maxChunksI = defMaxChunks;
   alnSetST->minChunkLenI = defMinChunkLen;
   alnSetST->minPercKmerD = defMinPercKmer;

   // Aligment Needleman-Wunsch settings
   alnSetST->gapStartPenaltyI = defGapStartPenalty;
   alnSetST->gapExtendPenaltyI = defGapExtendPenalty;
   alnSetST->minScoreUI = defMinScore;

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
    char newScoreC,         // New value for [query][ref] combination
    struct alnSet *alnSetST   // structure with scoring matrix to change
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-14 TOC: Sec-1 Sub-1: setBasePairScore
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
   ' Fun-15 TOC: Sec-1 Sub-1: freeAlnSet
   '  - Frees and alnSet (alignment settings) structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   if(stackBl & 1) return; // Nothing to do
   if(alnSetST != 0) free(alnSetST); // No heap variables
   return; 
} // freeAlnSet
    
/*---------------------------------------------------------------------\
| Output: 1: if bases were a match; 0 if not
\---------------------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC, // Query base to see if same as reference
    char *refBaseC    // Reference base to see if same as query
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-16 TOC: Sec-1 Sub-1: checkIfBasesMatch
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
| Output: Modifies: alnErrUCAry to have letters instead of codes
\---------------------------------------------------------------------*/
void cnvtAlnErrAryToLetter(
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
                break;
            case defInsFlag:                    // insertion
                *alnErrUCAry = 'I';
                break;
            case defMatchFlag:                    // match
                *alnErrUCAry = '=';
                break;
            case defBaseFlag:                   // snp
                *alnErrUCAry = 'X';
                break;
        } // Switch; check the error type

        ++alnErrUCAry;
    } // While have codes to convert

    return;
} // cnvtAlnErryAryToLetter

/*---------------------------------------------------------------------\
| Output:
|  - Returns:
|    o array with flags for snp/match, insertion, and deletions at each
|      position. (1 = snp/match, 2 = insertion, 4 = deletion)
|    o 0 for memory allocation errors
|  - Modifies:
|    o lenErrAryUI to hold the length of the returned array
\---------------------------------------------------------------------*/
uint8_t * NeedleManWunschAlnCostly(
    char *queryCStr,        // Full query sequence as c-string
    int32_t queryStartI,    // Starting query coordinate for alignment
    int32_t queryEndI,      // Ending query coordinate for alignment
    char *refCStr,          // Full reference sequence as c-string
    int32_t refStartI,      // Starting reference coordinate for aln
    int32_t refEndI,        // Ending reference coordinate for alignment
    struct alnSet *settings,// Settings for the alignment
    uint32_t *lenErrAryUI,  // Will hold the return arrays length
    int32_t *scoreI        // Score for the alignment (bottom right)
    // *startI and *endI paramaters should be index 1
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-0? TOC: NeedleManWunschAlnCostly
   '  - Perform a Needleman-Wunsch alignment on input sequences
   '    this variation uses a scoring array and direction array
   '  o fun-0? sec-1: Variable declerations
   '  o fun-0? sec-2: Allocate memory for alignment
   '  o fun-0? sec-3: Fill in the initial negatives for the reference
   '  o fun-0? sec-4: Fill the matrix with scores
   '  o fun-0? sec-5: Find the best path
   '  o fun-0? sec-6: Clean up and invert error array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-1: Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *tmpRefCStr = refCStr;
   char *tmpQueryCStr = queryCStr;

   // Find the length of the reference and query
   int32_t lenQueryI = queryEndI - queryStartI + 1;
       // +1 to convert back to 1 index (subtraction makes 0)
   int32_t lenRefI = refEndI - refStartI + 1;
       // +1 to convert back to 1 index (subtraction makes 0)

   int32_t snpScoreI = 0;     // Score for single base pair
   int32_t scoreTopI = 0;     // Score when using the top cell
   int32_t scoreDiagnolI = 0; // Score when using the diagnol cell
   int32_t scoreLeftI = 0;    // Score when using the left cell

   int32_t *scoreMatrixI = 0; // matrix to use in alignment
   int32_t *scoreOnIPtr = 0;  // Score I am currently working on
   int32_t *lastBaseIPtr = 0; // Pointer to cell with last base
   int32_t shiftByI = 0;      // How much to shift back when aliging

   uint8_t *dirMatrixUC = 0;  // Directions for each score cell
   uint8_t bitElmUC = 0;      // The value of a single bit element
   uint8_t *dirOnUCPtr = 0;   // Score working on
   uint8_t bitUC = 0;  // Keep track of if need to change direction elm

   uint8_t *topDirUCPtr = 0; // Direction of last score in the matrix
   uint8_t topBitUC = 0;     // Element on for last direction

   uint8_t *leftDirUCPtr = 0;
   uint8_t leftBitUC = 0;
     // Every 2 bits tells were to move next 11=top,10=diagnol,01=left
     // Find next cell: top: scoreOnIPtr - lenRefI;
     // Find next cell: diagnol: scoreOnIPtr - lenRefI - 1;
     // Find next cell: left: --scoreOnIPtr;

   uint8_t *alnErrAryUC = 0;   // Error type array (match/snp, ins, del)
   uint32_t numErrUI = 0;      // Number of error detected
   uint8_t *startUCPtr = 0;    // For inverting the alignment array
   uint8_t *endUCPtr = 0;     // For inverting the alignment array
   uint8_t swapUC = 0;         // For swaping elements in the array

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-2: Allocate memory for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   scoreMatrixI =
       malloc(sizeof(int32_t) * ((lenRefI + 1) * (lenQueryI + 1)));
       // lenRefI + 1 is to account for the insertion  reference row
       // lenQeurI + 1 is to account for insertion zero query column
   if(scoreMatrixI == 0) return 0;

   // This gives me an array of 0's (so no extra clearing)
   dirMatrixUC =
       calloc((
           1 + (((lenRefI + 1) * (lenQueryI + 1)) >> 2)),
           sizeof(uint8_t)
   ); // Make the direction array for the scoring array
       // 1 + is to make sure have enough chars
       // lenRefI + 1 is to account for the insertion  reference row
       // lenQeurI + 1 is to account for insertion zero query column
       // x >> 2 = x / 4 & accounts for taking 2 bits (char = 8 bits)
          // per cell

   if(dirMatrixUC == 0)
   { // If I do not have a direction matrix for each cell
       free(scoreMatrixI);
       return 0;
   } // If I do not have a direction matrix for each cell

   // Longest path in an n^2 matrix is query + reference
   if(*lenErrAryUI == 0) *lenErrAryUI = lenQueryI + lenRefI;
   alnErrAryUC = malloc(sizeof(int8_t) * *lenErrAryUI);

   if(alnErrAryUC == 0) 
   { /*If I had a memory allocation error*/
       free(scoreMatrixI);
       free(dirMatrixUC);
       return 0;
   } /*If I had a memory allocation error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-3: Fill in the initial negatives for the reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // Build up the indels for the reference row
   scoreOnIPtr = scoreMatrixI;
   *scoreOnIPtr = 0;         // Top left cell
   ++scoreOnIPtr;            // Move to first reference base cell
   *scoreOnIPtr = settings->gapStartPenaltyI;
       // Second column of the first row holds the first indel
   ++scoreOnIPtr;            // Move to first reference base cell

   // Get direction matrix in sync with scoring matrix
   dirOnUCPtr = dirMatrixUC; // Move left
   *dirOnUCPtr |= defMoveLeft;
   *dirOnUCPtr = *dirOnUCPtr << 2;

   *dirOnUCPtr |= defMoveLeft;
   *dirOnUCPtr = *dirOnUCPtr << 2;
   bitUC = 2;                 // Have one value in the two bit array

   // <= to deal with extra blank column with indel penalty
   for(int32_t iCol = 2; iCol <= lenRefI; ++iCol)
   { // loop; till have initalized the first row
       *scoreOnIPtr = *(scoreOnIPtr - 1) +settings->gapExtendPenaltyI;
       *dirOnUCPtr |= defMoveLeft;

       ++scoreOnIPtr; // Move to the next element
       // Move past elements in the bit array
       // Not worried about specifiying left shift (0), since calloc
       // already set everything to 0.

       twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
   } // loop; till have initalized the first row

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-4: Fill the matrix with scores
   ^  o fun-0? sec-4 sub-1: Fill in the indel column
   ^  o fun-0? sec-4 sub-2: Get scores for insertion, deletion, match
   ^  o fun-0? sec-4 sub-3: Find direction from scores (best score)
   ^  o fun-0? sec-4 sub-4: Move to the next refernce/query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-0? Sec-4 Sub-1: Fill in the indel column
   \*******************************************************************/

   lastBaseIPtr = scoreMatrixI; // Starting on the 0 cell

   // Direction of the upper cell
   topDirUCPtr = dirMatrixUC;
   topBitUC = 0;

   // Starting on the first sequence row
   for(int32_t iRow = 0; iRow < lenQueryI; ++iRow)
   { // loop; compare one query base against all reference bases

       // Set up the indel column
       if(lastBaseIPtr != scoreMatrixI)
           *scoreOnIPtr = *lastBaseIPtr + settings->gapExtendPenaltyI;
       else
           *scoreOnIPtr = settings->gapStartPenaltyI;
           // Else is the first indel for the indel column
 
       *dirOnUCPtr |= defMoveUp; // Set to move to top
       leftDirUCPtr = dirOnUCPtr; // Previous direction
       leftBitUC = bitUC;         // bit of previous direction

       // Move to the first base comparison
       ++scoreOnIPtr;  // Get of negative column for the new query base
       ++lastBaseIPtr; // Get of negative column for last query base

       twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
       twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);

       tmpRefCStr = refCStr; // Start at the beging of the reference

       /***************************************************************\
       * Fun-0? Sec-4 Sub-2: Get scores for insertion, deletion, match
       \***************************************************************/

       // First reference bases column (indel column already handled)
       for(int32_t iCol = 0; iCol < lenRefI; ++iCol)
       { // loop; compare one query to one reference base

           snpScoreI =
              getBasePairScore(
                  tmpQueryCStr,
                  tmpRefCStr,
                  settings
           ); // Find the score for the two base pairs

           // Find the score for the diagnol cell (snp/match)
           scoreDiagnolI = *(lastBaseIPtr - 1) + snpScoreI;

           // Check if I need to find the other scores
           if(checkIfBasesMatch(tmpQueryCStr, tmpRefCStr) & 1)
           { //If I have a match
               *dirOnUCPtr |= defMoveMatch;
               *scoreOnIPtr = scoreDiagnolI;
           } //If I have a match

           else
           { // else the bases do not match

               // Find the score for the top cell (insertion)
               bitElmUC = getTwoBitAryElm(topDirUCPtr, &topBitUC);

               if(bitElmUC & 1)     // If is the first indel
                 scoreTopI = *lastBaseIPtr + settings->gapStartPenaltyI;
               else                 // Else is part of a larger indel
                 scoreTopI = *lastBaseIPtr +settings->gapExtendPenaltyI;

               // find the score for the left cell (deletion)
               if(bitUC > 0) bitElmUC = (*dirOnUCPtr & (4 | 8)) >> 2;
               else bitElmUC = getTwoBitAryElm(leftDirUCPtr,&leftBitUC);
                  // If handles when have an incomplete element

               if(bitElmUC & 1)     // If is the first indel
                 scoreLeftI =
                     *(scoreOnIPtr - 1) + settings->gapStartPenaltyI;
               else                 // Else is part of a larger indel
                 scoreLeftI =
                     *(scoreOnIPtr - 1) + settings->gapExtendPenaltyI;

               /********************************************************\
               * Fun-0? Sec-4 Sub-3: Find direction (best score)
               \********************************************************/
               
               // The logic here is that I decide the best path as I
               // score, since the best of equal alternatives will
               // always follow an arbitary decision or report all
               // possible paths (I only care about 1). This allows me
               // to use a smaller direction matrix, but does add in
               // some more time.
               if(scoreTopI > scoreDiagnolI)
               { // If top score is better than the diagnol (or as good)
                   if(scoreTopI >= scoreLeftI)
                   { // If have a better top score
                       *dirOnUCPtr |= defMoveUp;
                       *scoreOnIPtr = scoreTopI;
                   } // If have a better top score

                   else
                   { // Else if have a better left score
                       *dirOnUCPtr |= defMoveLeft;
                       *scoreOnIPtr = scoreLeftI;
                   } // Else if have a better left score
               } // If top score is better than the diagnol (or as good)

               else if(scoreLeftI <= scoreDiagnolI)
               { //else if have a better or as good diagnol score
                   *dirOnUCPtr |= defMoveDiagnol;
                   *scoreOnIPtr = scoreDiagnolI;
               } //else if have a better or as good diagnol score

               else
               { //Else have a better left score
                   *dirOnUCPtr |= defMoveLeft;
                   *scoreOnIPtr = scoreLeftI;
               } //Else have a better left score
           } // If the bases do not match

           /************************************************************\
           * Fun-0? Sec-4 Sub-4: Move to the next refernce/query base
           \************************************************************/
       
           // Move to the next cell to score
           ++scoreOnIPtr; // Move to next comparison for this query base
           ++lastBaseIPtr; // Move to next element
           ++tmpRefCStr;   // Move to the next reference base

           leftDirUCPtr = dirOnUCPtr; // Previous direction
           leftBitUC = bitUC;         // bit of previous direction

           twoBitAryShiftBitsForNewElm(&dirOnUCPtr, &bitUC);
           twoBitAryMoveToNextElm(&topDirUCPtr, &topBitUC);
       } // loop; compare one query to one reference base

       ++tmpQueryCStr; // Move to the next query base
   } // loop; compare one query base against all reference bases

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-5: Find the best path
   ^   o fun-0? sec-5 sub-1: Get to the very last score (bottom right)
   ^   o fun-0? sec-5 sub-2: Find the best path
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*******************************************************************\
   * Fun-0? Sec-5 Sub-1: Get to the very last score (bottom right)
   \*******************************************************************/

   // Move to the bottom right base
   --scoreOnIPtr;
   *scoreI = *scoreOnIPtr;           // Get the socre for the alignment

   switch(bitUC)
   { // Switch; check if need to move to new direction element
       case 0: // moved to new char, so last score is 1st & 2nd bits
          --dirOnUCPtr;  // Move back to the char for the last score
          bitElmUC = *dirOnUCPtr & (1 | 2);
          break;
       case 1:           // Is the last bit 
       case 2:
       case 3:                                // The first two bitsk
          bitElmUC = *dirOnUCPtr >> 2; //Get off score after final score
          bitElmUC = *dirOnUCPtr & (1 | 2);
          break;
   } // Switch; check if need to move to new direction element

   /*******************************************************************\
   * Fun-0? Sec-5 Sub-2: Find the best path
   \*******************************************************************/

   while(scoreOnIPtr > scoreMatrixI)
   { // While I have more bases to add to the path
       switch(bitElmUC)
       { // Switch: check what the next base is in the sequence
           case defMoveUp:                    // Move to top (insertion)
               shiftByI = (scoreOnIPtr - scoreMatrixI) - lenRefI - 1;
               *(alnErrAryUC + numErrUI) = defInsFlag;  // flag for ins
               break;

           case defMoveMatch:
               *(alnErrAryUC + numErrUI) = defMatchFlag; // match
               shiftByI = (scoreOnIPtr - scoreMatrixI) - lenRefI - 2;
               break;

           case defMoveDiagnol:           // Move to diagnol (match/snp)
               *(alnErrAryUC + numErrUI) = defBaseFlag; // snp
               shiftByI = (scoreOnIPtr - scoreMatrixI) - lenRefI - 2;
               break;

           case defMoveLeft:              // Move to left (deletion)
               shiftByI = (scoreOnIPtr - scoreMatrixI) - 1;
               *(alnErrAryUC + numErrUI) = defDelFlag;  // flag for del
               break;
       } // Switch: check what the next base is in the sequence

       scoreOnIPtr = scoreMatrixI + shiftByI;

       // Get the next direction
       dirOnUCPtr = dirMatrixUC + (shiftByI >> 2);
         // >> 2 accounts for having 4 scores per char
       bitElmUC = shiftByI & (1 | 2); // Get the bit on
       bitElmUC = getTwoBitAryElm(dirOnUCPtr, &bitElmUC);
           // Gives the new direction

       ++numErrUI; // Account for the added error
   } // While I have more bases to add to the path

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-0? Sec-6: Clean up and invert error array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(scoreMatrixI);
   free(dirMatrixUC);

   *(alnErrAryUC + numErrUI) = 0;    // flag the end of the error array

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

