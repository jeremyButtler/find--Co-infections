/*######################################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Includes:
# C Standard Includes:
#   - <stdint.h>
######################################################################*/

#include "twoBitArrays.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  - fun-01 getElmFromToBitUCAry:
'     o Get an element from a two bit array
'  - fun-02 twoBitAryShiftBytsForNewElm:
'     o Make room in a unit8_t for two more two-bit elements
'  - fun-03 twoBitAryMoveToNextElm:
'     o Moves to the next element in a two-bit array
'  - fun-04 twoBitAryMoveForXElm:
'     o Moves forward in a two bit array by x two bit elements
'  - fun-05 twoBitAryMoveBackOneElm:
'     o Moves back one element in a 2-bit array
'  - fun-06 twoBitAryMoveBackXElm:
'     o Moves back X elements in a 2-bit array
'  - fun-07 finishTowBitAry:
'     o Makes sure the last element in a two bit array is buffered. This
'       aviods acess errors when trying to get at an element
'  - fun-08 changeTwoBitElm:
'     o Changes a single two bit value in a two bit array.
'  - fun-09 blankLimb:
'     o Sets a uint8_t (a limb) in a two bit array to 0. Each limb holds
'       four two bit elments.
'  - fun-10 moveToNextLimb:
'     o Moves to the start of the next limb (uint8_t) in the two bit 
'       array
'  - fun-11 moveToLastLimb:
'     o Moves to the start of the previous limb (uint8_t) in the two bit 
'       array
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Two bits of interest from the two bit array
\---------------------------------------------------------------------*/
uint8_t getTwoBitAryElm(
    uint8_t *twoBitAryUC,// Array of bytes that contain 4 2 bit elements
    uint8_t *charElmOnUC // The two bits on in the unit8_t
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
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
   ' Fun-02 TOC: Sec-1 Sub-1: twoBitAryShiftBytsForNewElm
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
   ' Fun-03 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
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
|    o twoBitAryUC to point to the unit8_t that is shiftByI two bit
|      elements ahead.  (one uint8_t holds four 2-bit elements)
|    o Changes bitUC to be on the target 2-bit element in the uint8_t
|      element twoBitAryUC will point to
\---------------------------------------------------------------------*/
void twoBitAryMoveForXElm(
    uint8_t **twoBitAryUC, // To bit array to move to next element in
    uint8_t *bitUC,        // Which bit am I on in the current uint8_t
    int32_t shiftByI       // How many elements to shift back by
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: twoBitAryMoveForXElm
   '  - Moves forward in a two bit array by x two bit elements
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   *twoBitAryUC += (shiftByI >> 2); // Get whole shifts to perform

   // Move the bit above 
   switch(*bitUC + (shiftByI & (1 | 2)))
   { // Switch; check if need to move to a new element
       case 0:
           *bitUC = 0;
           return;
       case 1:
           *bitUC = 1;
           return;
       case 2:
          *bitUC = 2;
          return;
       case 3:
          *bitUC = 3;   // On the last element
          return;
       case 4:
          ++(*twoBitAryUC);
          *bitUC = 0;   // Starting a new element
          return;
       case 5:
          ++(*twoBitAryUC);
          *bitUC = 1;   // On the second element in the next uint8_t
          return;
       case 6:
          ++(*twoBitAryUC);
          *bitUC = 2;   // On the third element in the next uint8_t
          return;
   } // Switch; check if need to move to a new element

   return;
} // twoBitAryMoveForXElm

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
   ' Fun-05 TOC: Sec-1 Sub-1: twoBitAryMoveBackOneElm
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
   ' Fun-06 TOC: Sec-1 Sub-1: towBitAryMoveBackXElm
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
| Output: Modifies end of two bit array to be shifted correctly
\---------------------------------------------------------------------*/
void finishTwoBitAry(
    uint8_t *endOfTwoBitUCArray,
    uint8_t lastBitUC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: finishTowBitAry
   '  - Makes sure the last element in a two bit array is buffered. This
   '    aviods acess errors when trying to get at an element
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   switch(lastBitUC)
   { // Switch; check if need to finish off the bit elment
       case 0: return; // On a new element, no need to shift
       case 1:         // Current value is on the last element
           *endOfTwoBitUCArray = *endOfTwoBitUCArray << 4;
           return;
       case 2:         // Current value is one element off
           *endOfTwoBitUCArray = *endOfTwoBitUCArray << 2;
           return;
       case 3: return;  // Finshed the element, no need to shift
   } // Switch; check if need to finish off the bit elment
} // finishTwoBitAry

/*---------------------------------------------------------------------\
| Output: Changes the taget elemtn in twoBitUCArray to newValueUC
| Note: This function is here to complete the two bit array functions,
|       but has not been tested. It will likely work.
\---------------------------------------------------------------------*/
void changeTwoBitElm(
    uint8_t *twoBitUCArray,
    uint8_t bitUC,
    uint8_t newValueUC // Two bit value to change to. Only the first two
                       // bits should be filled (otherwise will mess up
                       // array). This function will put the 1st 2 bits
                       // into the correct position, but will not clear
                       // any bits in newValueUC.
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: changeTwoBitElm
   '  - Changes a single two bit value in a two bit array.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    switch(bitUC)
    { // Switch: find the element on
       case 0:                               // Working on 7th/8th bits
           *twoBitUCArray &= (~(64 | 128));
           *twoBitUCArray |= (newValueUC  << 6);
           return;
       case 1:                               // Working on 5th/6th bits
           *twoBitUCArray &= (~(16 | 32));
           *twoBitUCArray |= (newValueUC << 4); 
           return;
       case 2:                               // Working on 5th/6th bits
           *twoBitUCArray &= (~(4 | 8));
           *twoBitUCArray |= (newValueUC << 2); 
           return;
       case 3: // The first two bits
           *twoBitUCArray &= (~(1 | 2));
           *twoBitUCArray |= newValueUC; 
           return;
    } // Switch: find the element on

    return;
} // changeTwoBitElm

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray and bitUC to 0
\---------------------------------------------------------------------*/
void blankLimb(
    uint8_t *twoBitUCArray, // Limb (uint8_t) on in the two bit array
    uint8_t *bitUC          // Element on in the limb (set to 0)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: blankLimb
   '  o Sets a uint8_t (a limb) in a two bit array to 0. Each limb holds
   '    four two bit elments.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   *twoBitUCArray = 0;
   *bitUC = 0;         // Limb has no data in it
} // blankLimb

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray to the next limb (uint8_t) and bitUC to 0
\---------------------------------------------------------------------*/
void moveToNextLimb(
    uint8_t *twoBitUCArray, // Limb to move from
    uint8_t *bitUC          // element on in the limb (set to 0)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: moveToNextLimb
   '  - Moves to the start of the next limb (uint8_t) in the two bit 
   '    array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   ++(*twoBitUCArray);
   *bitUC = 0;
} // moveToNextLimb

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray to the previous limb (uint8_t) & bitUC to 0
\---------------------------------------------------------------------*/
void moveToLastLimb(
    uint8_t *twoBitUCArray, // Limb to move from
    uint8_t *bitUC          // element on in the limb (set to 0)
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: moveToLastLimb
   '  - Moves to the start of the previous limb (uint8_t) in the two bit 
   '    array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   --(*twoBitUCArray);
   *bitUC = 0;
} // moveToLastLimb
