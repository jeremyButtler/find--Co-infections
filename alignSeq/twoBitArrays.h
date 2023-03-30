/*######################################################################
# Name: twoBitArrays
# Use:
#   o Holds functions to handle two bit arrays
# Includes:
# C Standard Includes:
#   - <stdint.h>
######################################################################*/

#ifndef TWOBITARRAYS_H
#define TWOBITARRAYS_H

#include <stdint.h>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
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
| Interfacing: A quick guid to setting up your two bit array
|   - This is here to describe the order you might want to use these
|     functions.
|   o Make an uint8_t array for the twobit array
|     - Size should be: (maximum elemnents / 4) + 1
|     - This is an array of limbs, were each limb holds four elements
|   o make a uint8_t bit counter set to 0 (uint8_t bitUC = 0)
|     - This is to keep track of the element on in the limb
|   o TwoBitAryShiftBytsForNewElm is here to allow you to manipulate two
|     bit array elements directly. You can set new elements by
|     *twoBitArray |= newTwoBitElement and then calling
|     TwoBitAryShiftBytsForNewElm to shift the elements.
|     - This will result in the very end elements of the array not being
|       shifted correctly, so you will want to call finishTowBitAry at
|       the end.
|   o An alternative way to add elements would be call changeTwoBitElm,
|     which replaces the old elements with the target elements.
|   o After adding an element you will want to move to the next element
|     with twoBitAryMoveToNextElm
|     - Use twoBitAryMoveBackOneElm to move backwards
|   o Your can check an element using getElmFromToBitUCAry
\---------------------------------------------------------------------*/

/*---------------------------------------------------------------------\
| Output: Two bits of interest from the two bit array
\---------------------------------------------------------------------*/
uint8_t getTwoBitAryElm(
    uint8_t *twoBitAryUC,// Array of bytes that contain 4 2 bit elements
    uint8_t *charElmOnUC // The two bits on in the unit8_t
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Sec-1 Sub-1: getElmFromToBitUCAry
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
   ' Fun-02 TOC: Sec-1 Sub-1: twoBitAryShiftBytsForNewElm
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
   ' Fun-03 TOC: Sec-1 Sub-1: twoBitAryMoveToNextElm
   '  - Moves to the next element in a two-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: twoBitAryMoveForXElm
   '  - Moves forward in a two bit array by x two bit elements
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
   ' Fun-05 TOC: Sec-1 Sub-1: towBitAryMoveBackOneElm
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
   ' Fun-06 TOC: Sec-1 Sub-1: towBitAryMoveBackXElm
   '  - Moves back X elements in a 2-bit array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies end of two bit array to be shifted correctly
\---------------------------------------------------------------------*/
void finishTwoBitAry(
    uint8_t *endOfTwoBitUCArray,
    uint8_t lastBitUC
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: Sec-1 Sub-1: finishTowBitAry
   '  - Makes sure the last element in a two bit array is buffered. This
   '    aviods acess errors when trying to get at an element
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: Sec-1 Sub-1: changeTwoBitElm
   '  - Changes a single two bit value in a two bit array.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray and bitUC to 0
\---------------------------------------------------------------------*/
void blankLimb(
    uint8_t *twoBitUCArray, // Limb (uint8_t) on in the two bit array
    uint8_t *bitUC          // Element on in the limb (set to 0)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: Sec-1 Sub-1: blankLimb
   '  o Sets a uint8_t (a limb) in a two bit array to 0. Each limb holds
   '    four two bit elments.
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray to the next limb (uint8_t) and bitUC to 0
\---------------------------------------------------------------------*/
void moveToNextLimb(
    uint8_t *twoBitUCArray, // Limb to move from
    uint8_t *bitUC          // element on in the limb (set to 0)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: Sec-1 Sub-1: moveToNextLimb
   '  - Moves to the start of the next limb (uint8_t) in the two bit 
   '    array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Sets twoBitUCArray to the previous limb (uint8_t) & bitUC to 0
\---------------------------------------------------------------------*/
void moveToLastLimb(
    uint8_t *twoBitUCArray, // Limb to move from
    uint8_t *bitUC          // element on in the limb (set to 0)
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: Sec-1 Sub-1: moveToLastLimb
   '  - Moves to the start of the previous limb (uint8_t) in the two bit 
   '    array
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif

