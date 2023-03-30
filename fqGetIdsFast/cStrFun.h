/*######################################################################
# Use:
#   o Holds functions for copying or manipualting c-strings
######################################################################*/

#ifndef CSTRFUN_H
#define CSTRFUN_H

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold the copied C-string
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cStrCpInvsDelm(
    char *cpToCStr,  /*C-string to copy values to*/
    char *cpFromCStr /*C-string to copy*/
); /*Copy one c-string till an tab, newline, or '\0' (keeps spaces)*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - cpToCStr to hold space, parameter, space, & argument
|    Returns:
|        - pointer to null at end of cpToCStr
\---------------------------------------------------------------------*/
char * cpParmAndArg(
    char *cpToCStr,   /*Holds copied parameter and argement*/
    char *cpParmCStr, /*Paramater to copy*/
    char *cpArgCStr   /*Argument to copy*/
); /*Copies adds a space and copies a paramater and a agrugment*/

#endif
