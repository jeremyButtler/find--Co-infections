/*######################################################################
# Name: printErrors.c
# Use:
#     - These are my functions for common error messages I might print
#       out. Right know this is only the memory allocation error.
# Includes: 
#    - <stdio.h>
######################################################################*/

#ifndef PRINTERRORS_H
#define PRINTERRORS_H

#include <stdio.h>
#include <stdint.h>

/*######################################################################
# Output:
#    prints: Memory allocation error message to stdout
######################################################################*/
void printMemAlocErr(
    char * fileUCStr,          /*name of file error was in*/
    char * funNameUCStr,       /*name of function with error*/
    char funNumUChar,          /*Number of function error out*/
    unsigned long lineNumULng     /*Line number error out at*/
); /*Prints out error message for memory allocation failure*/

#endif
