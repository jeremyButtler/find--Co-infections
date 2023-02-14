/*##############################################################################
# Name: scoreReadsConversion
# Use: Holds functions to convert data types (ex: c-string to number)
##############################################################################*/

#ifndef SCOREREADSCONVERSION_H
#define SCOREREADSCONVERSION_H

#include <limits.h> /*get max data type sizes*/
#include <math.h>   /*For log10*/

/*Convert c-string to numeric character*/
char cStrToChar(char *charCStr, char delimChar);

/*Convert c-string to unsinged integer*/
char cStrToUInt(char *charCStr, unsigned int *retUInt, char delimChar);
#endif

