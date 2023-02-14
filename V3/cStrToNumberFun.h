/*######################################################################
# Name: cStrToNumberFun.c
# Use: My numeric functions to convert numbers to specific data types.
#      This allows me to avoid issues of overflows form the base 
#      c functions, wich always convert to unsigned longs or longs.
# Note: This uses set sizes for each data type & so is a bit less
#       universall than the base functions, but makes these a bit
#       faster to. I also have unrolled loops so they will be faster,
#       however, their are still some good performance gains to be had
#       with O3 compiling.
# Includes:
#    - <stdint.h>
######################################################################*/

#ifndef CSTRTONUMBERFUN_H
#define CSTRTONUMBERFUN_H

#include <stdint.h> /*uintx_t variables*/

/*######################################################################
# Output:
#    Returns: pionter to character after last converted number
#    Modifies: retUInt to hold the uint32_t integer
# Note:
#    This function will only convert the first number of digits in an
#    uint32_t unsigned integer or till non-numeric character
######################################################################*/
char * cStrToUInt(
    char *charUCStr, /*C-string to convert to number*/
    uint32_t *retUInt   /*Holds converted number*/
); /*converst a c-string into an uint32_teger*/

/*######################################################################
# Output:
#    Returns: pionter to character after last converted number
#    Modifies: retUInt to hold the uint16_t unsigned short
# Note:
#    This function will only convert the first number of digits in an
#    uint16_t unsigned short or till non-numeric character
######################################################################*/
char * cStrToUSht(
    char *charUCStr, /*C-string to convert to number*/
    uint16_t *retUSht   /*Holds converted number*/
); /*converst a c-string into an uint16_t integer*/

/*######################################################################
# Output:
#    Returns: pionter to character after last converted number
#    Modifies: retUInt to hold the char unsigned character
# Note:
#    This function will stop converting at a buffer overlflow or till a
#    non-numeric character
######################################################################*/
char * cStrToUChar(
    char *charUCStr, /*C-string to convert to number*/
    unsigned char *retUChar   /*Holds converted number*/
); /*converst a c-string into an uint32_teger*/

/*######################################################################
# Output:
#    Returns: pionter to character after last converted number
#    Modifies: retUInt to hold the uint32_t integer
# Note:
#    This function will only convert the first number of digits in an
#    uint32_t unsigned integer or till non-numeric character
######################################################################*/
char * backwarsCStrToUInt(
    char *charUCStr, /*C-string to convert to number*/
    uint32_t *retUInt   /*Holds converted number*/
); /*converst a backwards c-string into an uint32_t integer*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - buffCStr to be c-string with the converted number
|    Returns:
|        - pointer to end of bufferCStr (will point to '\0')
\----------------------------------------------------------------------*/
char * uCharToCStr(
    char *buffCStr,  /*Buffer to hold output c-string (4 elements)*/
    char uCharToCnvt /*Chacter to convert to c-string*/
); /*converts unsigned character value to a numeric c-string*/

#endif
