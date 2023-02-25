/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cStrFun SOF:
'   fun-1 cStrCpInvsDelm:
'     o Copy one c-string till an tab, newline, or '\0' (keeps spaces)
'   fun-2 cpParmAndArg:
'     o Copies adds a space and copies a paramater and a agrugment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
) /*Copy one c-string till an tab, newline, or '\0' (keeps spaces)*/
{ /*cStrCpInvsDelm*/


    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-20 TOC: Sec-1 Sub-1: cStrCpInvsDelim
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    while(*cpFromCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpFromCStr;
        ++cpToCStr;
        ++cpFromCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cStrCpInvsDelim*/

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
) /*Copies adds a space and copies a paramater and a agrugment*/
{ /*cpParmAndArg*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-21 TOC: Sec-1 Sub-1: cpSpaceCStr
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    *cpToCStr = ' ';
    ++cpToCStr;
    
    while(*cpParmCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpParmCStr;
        ++cpToCStr;
        ++cpParmCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = ' ';
    ++cpToCStr;

    while(*cpArgCStr > 31) /*Copy spaces & ever visible character*/
    { /*While have a  c-string to copy*/
        *cpToCStr = *cpArgCStr;
        ++cpToCStr;
        ++cpArgCStr;
    } /*While have a  c-string to copy*/

    *cpToCStr = '\0';
    return cpToCStr;
} /*cpParmAndArg*/
