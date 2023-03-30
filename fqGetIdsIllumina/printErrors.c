/*######################################################################
# Name: printErrors.c
# Use:
#     - These are my functions for common error messages I might print
#       out. Right know this is only the memory allocation error.
# Includes: 
#    - <stdio.h>
######################################################################*/

#include "printErrors.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    fun-1 printMemAlocErr: Prints memorry allocation error message
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*######################################################################
# Output:
#    prints: Memory allocation error message to stdout
######################################################################*/
void printMemAlocErr(
    char * fileUCStr,          /*name of file error was in*/
    char * funNameUCStr,       /*name of function with error*/
    char funNumUChar,          /*Number of function error out*/
    unsigned long lineNumULng     /*Line number error out at*/
)/*Prints out error message for memory allocation failure*/
{ /*Prints error message about memmory allocation failure*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-13 Sec-1 Sub-1 TOC: printMemAlocErr: print memory error
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fprintf(
        stdout,
        "Failed memory allocation: %s %s (Fun-%u) line %lu\n",
        fileUCStr,
        funNameUCStr,
        funNumUChar,
        lineNumULng
   ); /*Print out error*/

   return; /*Done*/
} /*Prints error message about memmory allocation failure*/
