/*##############################################################################
# Name: scoreReadsConversion
# Use: Holds functions to convert data types (ex: c-string to number)
##############################################################################*/

#include "scoreReadsConversion.h"

/*##############################################################################
#    char checkIfNumChar: Convert number to char (if by 256) [fail: -1]
#       fun-1 sec-1: variable declarations
#       fun-1 sec-2: Check if string is number (returns 0 if not)
#       fun-1 sec-3: Convert the string to an char
#    char checkIfNumChar: Convert number to char (if by 256) [fail: -1]
#       fun-2 sec-1: variable declarations
#       fun-2 sec-3: Convert the string to an char
##############################################################################*/

char NUM_DIGS_UINT = log10(UINT_MAX) + 1; /*number digits in unisigned int*/
char NUM_DIGS_CHAR = log10(CHAR_MAX) + 1; /*number digits in character*/

/*##############################################################################
# Name: cStrToChar
# Use: Checks if the input is a numeric number roughly the size of character
# Input:
#    charCStr: Numeric c-string with to convert to char with number
#    delminChar: Character to stop reading line on
# Output:
#    returns char holding the input number
#    returns -1 if failed (Assuming user will never input -1)
##############################################################################*/
char cStrToChar(char *charCStr, char delimChar)
{ /*cStrToChar*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: variable declarations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = charCStr, numDigChar = 0;
    unsigned long tmpULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Check if string is number (returns 0 if not)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpCStr != delimChar)
    { /*loop through c-string and check if char is a number*/
        numDigChar++; /*keep track of the number of digits*/

        if(numDigChar > NUM_DIGS_CHAR)
            return -1; /*To long of a number for an unsigned int*/
        if(*tmpCStr < 48 || *tmpCStr > 57)
            return -1; /*is not a number*/

        tmpULng = tmpULng * 10 + (*tmpCStr - 48); /*Add the next digit in*/
        tmpCStr++; /*move to the next character*/
    } /*loop through c-string and check if char is a number*/

    if(tmpULng > CHAR_MAX)
        return -1; /*if number is longer than max char (~256)*/

    return tmpULng; /*Return success*/
} /*cStrToChar*/

/*##############################################################################
# Name: cStrToUInt
# Use: Checks if the input is a numeric number roughly the size of unsigned int
# Input:
#    charCStr: Numeric c-string with to convert to unsigned int (c-string)
#    retUInt: Modified to hold the conversion, if successful (unsigned int)
#    delminChar: Character to stop reading line on
# Output:
#    returns 1 if successful, 0 if failed
#    Modifies: retUInt to hold the unsigned integer
##############################################################################*/
char cStrToUInt(char *charCStr, unsigned int *retUInt, char delimChar)
{ /*cStrToUInt*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: variable declarations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = charCStr, /*temporary pointer to c-string*/
         numDigChar = 0;      /*holds the number of digits added*/
    unsigned long tmpULng = 0; /*holds the conversion, till confirmed*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Convert and check the number
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*tmpCStr != delimChar)
    { /*loop through c-string and check if char is a number*/
        numDigChar++; /*keep track of the number of digits*/

        if(numDigChar > NUM_DIGS_UINT)
            return 0; /*To long of a number for an unsigned int*/
        if(*tmpCStr < 48 || *tmpCStr > 57)
            return 0; /*is not a number*/

        tmpULng = tmpULng * 10 + (*tmpCStr - 48); /*Add the next digit in*/
        tmpCStr++; /*move to the next character*/
    } /*loop through c-string and check if char is a number*/

    if(tmpULng > UINT_MAX)
        return -1; /*if number is longer than max char (~256)*/

    *retUInt = tmpULng; /*Send back the correct number*/
    return 1; /*Return success*/
} /*cStrToUInt*/
