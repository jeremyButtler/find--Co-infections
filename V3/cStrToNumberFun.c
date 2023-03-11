/*######################################################################
# Name: cStrToNumberFun.c
# Use: My numeric functions to convert numbers to specific data types.
#      This allows me to avoid issues of overflows form the base 
#      c functions, wich always convert to unsigned longs or longs.
# Note:
#    - These functions are only for base 10 conversions of their. This
#      makes these functions very specific, but also provides a speed
#      boost compared to strtoul, which is generalized for any base
#      conversion.
# Includes:
#    - <stdint.h>
######################################################################*/

#include "cStrToNumberFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' TOC:
'    fun-1 cStrToUInt:
'        - Converts c-string to uint32_t (unsigned integer)
'    fun-2 cStrToUSht:
'        - Converts c-string to uint16_t (unsigned short)
'    fun-3 cStrToUChar:
'        - Converts c-string to char (unsigned char)
'    fun-4 backwarsCStrToUInt:
'        - Converts backwards numeric c-string to uint32_t
'    fun-5 uCharToCStr:
'        - Converts number in char to c-string
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
) /*converst a c-string into an uint32_teger*/
{ /*cStrToUInt*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Sec-1 Sub-1: cStrToUInt
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*1st digit*/
        *retUInt = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUInt = 0;
        return charUCStr;
    } /*Else non-numeric*/

    ++charUCStr;

    for(uint32_t intCnt = 0; intCnt < 7; ++intCnt)
    { /*Convert digits with no overflow concern*/
        if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)
            *retUInt = *retUInt * 10 + *charUCStr - 48;
        else
            return charUCStr;
        ++charUCStr;
    } /*Convert digits with no overflow concern*/

    /*Convert last two digits, which could overflow*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*9th digit*/
    { /*If have one or tow more numbers*/
        if(10 * *retUInt + *charUCStr - 48 > 10000)
        { /*If can fit in one more number*/
            *retUInt = *retUInt * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or tow more numbers*/
    else
        return charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*10th digit*/
    { /*If have one more number*/
        if(10 * *retUInt + *charUCStr - 48 > 10000)
        { /*If would not cause an overflow*/
            /*Need to use > 1000 because will keep value as int*/
            *retUInt = *retUInt * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*cStrToUInt*/

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
) /*converst a c-string into an uint32_teger*/
{ /*cStrToUSht*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: Sec-1 Sub-1: cStrToUSht
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*1st digit*/
        *retUSht = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUSht = 0;
        return charUCStr;
    } /*Else non-numeric*/

    ++charUCStr;

    /*Convert digits were their is no risk of overlfow (2-8th)*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*2nd digit*/
        *retUSht = *retUSht * 10 + *charUCStr - 48;
    else
        return charUCStr;
    ++charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*3rd digit*/
        *retUSht = *retUSht * 10 + *charUCStr - 48;
    else
        return charUCStr;
    ++charUCStr;

    /*Convert last two digits, which could overflow*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*4th digit*/
    { /*If have one or tow more numbers*/
        if(10 * *retUSht + *charUCStr - 48 < 65536)
        { /*If can fit in one more number*/
            *retUSht = *retUSht * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or tow more numbers*/
    else
        return charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*5th digit*/
    { /*If have one more number*/
        if(10 * *retUSht + *charUCStr - 48 > 100)
        { /*If would not cause an overflow*/
            /*Need to use > 1000 because will keep value as int*/
            *retUSht = *retUSht * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*cStrToUSht*/

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
) /*converst a c-string into an uint32_teger*/
{ /*cStrToUChar*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-3 TOC: Sec-1 Sub-1: cStrToUChar
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    if(*charUCStr < 58 && *charUCStr > 47)
        *retUChar = *charUCStr - 48;
        /*& 15 sets 16 and 32 flags to 0, which converts char to num*/
    else
    { /*Else is a non-valid character*/
        *retUChar = 0;
        return charUCStr;
    } /*Else is a non-valid character*/

    ++charUCStr;

    if(*charUCStr < 58 && *charUCStr > 47)
        *retUChar = *retUChar * 10 + *charUCStr - 48;
    else
        return charUCStr;

    ++charUCStr;

    if(*charUCStr < 58 && *charUCStr > 47)
        if(*retUChar * 10 + *charUCStr - 48 < 256)
        { /*If we can add in one more digit*/
            *retUChar = *retUChar * 10 + *charUCStr - 48;
            ++charUCStr;
        } /*If we can add in one more digit*/

    return charUCStr;
} /*cStrToUInt*/

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
) /*converst a backwards c-string into an uint32_t integer*/
{ /*backwardsCStrToUInt*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-4 TOC: Sec-1 Sub-1: backwardsCStrToUInt
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*Convert first number*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*1st digit*/
        *retUInt = *charUCStr - 48;
    else
    { /*Else non-numeric*/
        *retUInt = 0;
        return charUCStr;
    } /*Else non-numeric*/

    --charUCStr;

    /*Convert digits were their is no risk of overlfow (2-8th)*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*2nd digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 10;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*3rd digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 100;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*4th digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 1000;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*5th digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 10000;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*6th digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 100000;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*7th digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 1000000;
    else
        return charUCStr;
    --charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*8th digit*/
        *retUInt = *retUInt + (*charUCStr - 48) * 10000000;
    else
        return charUCStr;
    --charUCStr;

    /*Convert last two digits, which could overflow*/
    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*9th digit*/
    { /*If have one or tow more numbers*/
        if(*retUInt + (*charUCStr - 48) * 10000000 > 10000)
        { /*If can fit in one more number*/
            *retUInt = *retUInt + (*charUCStr - 48) * 10000000;
            --charUCStr;
        } /*If can fit in one more number*/

        else
            return charUCStr;
    } /*If have one or tow more numbers*/
    else
        return charUCStr;

    if(*charUCStr - 58 < 0 && *charUCStr - 47 > 0)    /*10th digit*/
    { /*If have one more number*/
        if(*retUInt + (*charUCStr - 48) * 100000000 > 1000)
        { /*If would not cause an overflow*/
            /*Need to use > 1000 because will keep value as int*/
            *retUInt = *retUInt + (*charUCStr - 48) * 100000000;
            --charUCStr;
        } /*If would not cause an overflow*/
    } /*If have one more number*/

    return charUCStr;
} /*backwardsCStrToUInt*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|        - buffCStr to be c-string with the converted number
|    Returns:
|        - pointer to end of bufferCStr (will point to '\0')
\----------------------------------------------------------------------*/
char * uCharToCStr(
    char *buffCStr,  /*Buffer to hold output c-string (4 elements)*/
    unsigned char uCharToCnvt /*Chacter to convert to numeric c-string*/
) /*converts unsigned character value to a numeric c-string*/
{ /*uCharToCStr*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-5 TOC: Sec-1 Sub-1: uCharToCStr
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    char onesUC = 0;
    char tensUC = 0;

    onesUC = (uCharToCnvt % 10) + 48;
    uCharToCnvt = uCharToCnvt / 10; /*Move to next digit in character*/

    if(uCharToCnvt > 0)
        tensUC = (uCharToCnvt % 10) + 48;
    else
    { /*Else their is only a ones position*/
       *buffCStr = onesUC;
       ++buffCStr;

       *buffCStr = '\0';
       return buffCStr;
    } /*Else their is only a ones position*/

    uCharToCnvt = uCharToCnvt / 10; /*Move to next digit in character*/

    if(uCharToCnvt > 0)
    { /*If looking at the 100s position*/
        *buffCStr = (uCharToCnvt % 10) + 48;
        ++buffCStr;
    } /*If looking at the 100s position*/

    *buffCStr = tensUC;
    ++buffCStr;

    *buffCStr = onesUC;
    ++buffCStr;

    *buffCStr = '\0';
    return buffCStr;
} /*uCharToCStr*/
