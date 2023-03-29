/*######################################################################
# Name: alignSeq
# Use: Runs a Needlman Wunsch alignment on a pair of fasta files
# Includes:
#  - "alignmentsFun.h"
#  - "sequenceFun.h"
#  o "twoBitArrays.h"
# C standard libraries:
#  - <string.h>
#  o <stdlib.h>
#  o <stdio.h>
#  o <stdint.h>
######################################################################*/

#include "alignmentsFun.h"
#include "sequenceFun.h"
#include <string.h>

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOP: Start Of Program
'  - main: Run this program
'  - fun-01 checkInput;
'    o Gets user input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to -score-matrix if the scoring matrix is invalid
|  - Prints to stdout when the scoring file is invalid
\---------------------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,               /*Number arguments user input*/
    char *argsCStr[],              /*Array with user arguments*/
    char **refFileCStr,            /*file name of the reference file*/
    char **queryFileCStr,          /*File name of the query file*/
    char **outFileCStr,            // Name of the output file
    struct alnSet *alnSetST        // Aligment settings
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


int main(
    int lenArgsInt,
    char *argsCStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '  - Currently dummy wrapper for testing majCon, but will eventually
   '    handle user input
   '  o main sec-1: Variable declerations
   '  o main sec-2: Read in user input and check input
   '  o main sec-3: read in the reference sequence
   '  o main sec-4: read in the query sequence
   '  o main sec-5: Do the alingment
   '  o main sec-6: Make aligned query sequence
   '  o main sec-7: Make aligned reference sequence
   '  o main sec-8: Print out the alignment
   '  o main sec-9: Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-1: Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   // User input
   char *refFileCStr = 0;
   char *queryFileCStr = 0;
   char *outFileCStr = 0;
   char *inputCStr = 0;

   // Holds the reference sequence
   char *refCStr = 0;
   char *refHeadCStr = 0;
   uint32_t lenRefBuffUI = 0;
   uint32_t lenRefHeadBuffUI = 0;
   uint32_t lenRefSeqUI = 0;

   // Holds the query sequence
   char *queryCStr = 0;
   char *queryHeadCStr = 0;
   uint32_t lenQueryBuffUI = 0;
   uint32_t lenQueryHeadBuffUI = 0;
   uint32_t lenQuerySeqUI = 0;

   // Output aligned sequences
   char *queryAlnCStr = 0;
   char *refAlnCStr = 0;

   // Holds the alignment details
   uint8_t *alnErrUCAry = 0;
   uint32_t lenErrAryUI = 0;
   uint8_t *tmpUCPtr = 0;    // For iterating through alnErrUCAry
   long alnScoreL = 0;

   char *refPtrCStr = 0;
   char *queryPtrCStr = 0;
   uint8_t charOnUC = 0;

   unsigned char errUC = 0; // Hold error output from functions
    
   struct alnSet alnSetST;

   FILE *faFILE = 0;
   FILE *outFILE = 0;

   char *helpCStr = "\
       \n Use: Runs a Needleman-Wunsch alignment on input sequences\
       \n Run: alingSeq -query query.fasta -ref reference.fasta [...]\
       \n Input:\
       \n   -query: [Required]\
       \n     o Fasta sequence to align.\
       \n   -ref: [Required]\
       \n     o Second fasta sequence to align.\
       \n   -out: [stdout]\
       \n     o File to output alignment to (default to screen)\
       \n   -gapopen: [-1]\
       \n     o Cost of starting an indel (as integer).\
       \n     o A negative value is a penalty, while postives values\
       \n       are correct (+ favors gaps, - disfavors).\
       \n   -gapextend: [-4]\
       \n     o Cost of extending an indel.\
       \n     o A negative value is a penalty, while postives values\
       \n       are correct (+ favors gaps, - disfavors).\
       \n   -score-matrix: [ENDNAFULL]\
       \n     o File with matrix to use. It should have one line for\
       \n       each possible score. Scores not supplied will be set to\
       \n       defaults. (see scoring-matrix.txt file in source).\
       \n     o Format:\
       \n       a t -4\
       \n       a a 5 \
       \n       \\\\ This is a comment\
       \n  -use-water: [Needleman Wunsch]\
       \n     o Use a Waterman Smith alignment instead of the Needleman\
       \n       Wunsch (Waterman is a local, Needle is an global).\
       \n     o This only tracks one best alignment, not all best\
       \n       (or all) alignments.\
       \n  -enable-match-priority: [No]\
       \n     o Always go with matches for best score, even when an\
       \n       indel would give a better score.\
       \n     o This does not apply to SNPs.\
       \n  -match-ins-del: [ins-match-del]\
       \n     o For equal scores choose matches/SNPs over insertions.\
       \n       Choose insertions over deletions.\
       \n  -match-del-ins: [ins-match-del]\
       \n     o For equal scores choose matches/SNPs over deletions.\
       \n       Choose deletions over insertions.\
       \n  -ins-match-del: [ins-match-del]\
       \n     o For equal scores choose insertions over matches/SNPs.\
       \n       Choose matchs/SNPs over deletions.\
       \n  -del-match-ins: [ins-match-del]\
       \n     o For equal scores choose deletions over matches/SNPs.\
       \n       Choose matchs/SNPs over insertions.\
       \n  -ins-del-match: [ins-match-del]\
       \n     o For equal scores choose insertions over deletions.\
       \n       Choose deletions over matchs/SNPs.\
       \n  -del-ins-match: [ins-match-del]\
       \n     o For equal scores choose deletions over insertions.\
       \n       Choose insertions over matches/SNPs.\
       \n     o When scores are equal choose deletions then insertions.\
       \n Output:\
       \n  - Alignment with query on top, reference, and error line to\
       \n    stdout or specified file.\
       \n  - Error line format:\
       \n    o I = Insertion\
       \n    o D = Deletion\
       \n    o X = mismatch\
       \n    o = = match\
       ";

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-2: Read in user input and check input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   initAlnSet(&alnSetST);

   inputCStr =
       checkInput(
           &lenArgsInt,         /*Number arguments user input*/
           argsCStr,        /*Array with user arguments*/
           &refFileCStr,            /*file name of the reference file*/
           &queryFileCStr,          /*File name of the query file*/
           &outFileCStr,            // Name of the output file
           &alnSetST        // Aligment settings
   );

   if(inputCStr != 0)
   { // If had problematic input
        if(strcmp(inputCStr, "-h") == 0 ||
           strcmp(inputCStr, "--h") == 0 ||
           strcmp(inputCStr, "-help") == 0 ||
           strcmp(inputCStr, "--help") == 0 ||
           strcmp(inputCStr, "help") == 0
        ) { /*If user wanted the help message*/
            fprintf(stdout, "%s\n", helpCStr);
            exit(0);
        } /*If user wanted the help message*/

        if(strcmp(inputCStr, "-V") == 0 ||
           strcmp(inputCStr, "-v") == 0 ||
           strcmp(inputCStr, "--V") == 0 ||
           strcmp(inputCStr, "--v") == 0 ||
           strcmp(inputCStr, "--version") == 0 ||
           strcmp(inputCStr, "--Version") == 0 ||
           strcmp(inputCStr, "-version") == 0 ||
           strcmp(inputCStr, "-Version") == 0
        ) { /*if the user wanted the version number*/
            fprintf(
                stdout,
                "alignSeq from findCoInft version: %.8f\n",
                defVersion
            ); /*Print out the closest thing to a version*/
            exit(0);
        } /*Else if the user wanted the version number*/

        // If their was an invalid scoring file (error already printed)
        else if(strcmp(inputCStr, "-score-matrix") == 0)
            exit(1);

        else if(inputCStr != 0)
        { /*If user had invalid input*/
            fprintf(
                stderr,
                "%s\n%s is invalid\n",
                helpCStr,
                inputCStr
            ); /*Print out the problem*/
            exit(1); /*Let user know their was an error*/
        } /*If user had invalid input*/
   } // If had problematic input

   if(outFileCStr != 0)
   { // If printing output to a file
        outFILE = fopen(outFileCStr, "w");

        if(outFILE == 0)
        { // If an invalid output file
            printf("Output (-out %s) file is invalid.\n", outFileCStr);
            exit(-1);
        } // If an invalid output file

        fclose(outFILE);
        outFILE = 0;
   } // If printing output to a file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-3: read in the reference sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   faFILE = fopen(refFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       printf("Reference (-ref %s) could not be opend\n", refFileCStr);
       exit(-1);
   } // If reference file could not be opened

   errUC = 
       readFaSeq(
           faFILE,
           &refHeadCStr,
           &lenRefHeadBuffUI,
           &refCStr,
           &lenRefBuffUI,
           &lenRefSeqUI
   ); // Read in the reference sequence

   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       if(refCStr != 0) free(refCStr);
       if(refHeadCStr != 0) free(refHeadCStr);
       printf("Reference (-ref %s) is not valid\n", refFileCStr);
       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
       if(refCStr != 0) free(refCStr);
       if(refHeadCStr != 0) free(refHeadCStr);
       printf("Memory allocation error (not enough memory for ref)\n");
       exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-4: read in the query sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   faFILE = fopen(queryFileCStr, "r");

   if(faFILE == 0) 
   { // If reference file could not be opened
       free(refCStr);
       free(refHeadCStr);
       printf("Query (-query %s) could not be opend\n", queryFileCStr);
       exit(-1);
   } // If reference file could not be opened

   errUC = 
       readFaSeq(
           faFILE,
           &queryHeadCStr,
           &lenQueryHeadBuffUI,
           &queryCStr,
           &lenQueryBuffUI,
           &lenQuerySeqUI
   ); // Read in the reference sequence

   fclose(faFILE);
   faFILE = 0;

   if(errUC & 2)
   { // Invalid fasta file
       free(refCStr);
       free(refHeadCStr);
       if(queryHeadCStr != 0) free(queryHeadCStr);
       if(queryCStr != 0) free(queryCStr);
 
       printf("Query (-query %s) is not valid\n", refFileCStr);
       exit(-1);
   } // Invalid fasta file

   if(errUC & 64)
   { // Invalid fasta file
      free(refCStr);
      free(refHeadCStr);
      if(queryHeadCStr != 0) free(queryHeadCStr);
      if(queryCStr != 0) free(queryCStr);
      printf("Memory allocation error (not enough memory for query)\n");
      exit(-1);
   } // Invalid fasta file

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-5: Do the alingment
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(alnSetST.useNeedleBl != 0)
   { // If doing a Needleman Wunsch alignment
       alnErrUCAry =
           NeedleManWunschAln(
               queryCStr,
               1,
               lenQuerySeqUI,
               refCStr,
               1,
               lenRefSeqUI,
               &alnSetST,
               &lenErrAryUI,
               &alnScoreL
       ); // Get the alignment
   } // If doing a Needleman Wunsch alignment

   else
   { // Else doing a Waterman Smith alignment
       alnErrUCAry =
           WatermanSmithAln(
               queryCStr,
               1,
               lenQuerySeqUI,
               refCStr,
               1,
               lenRefSeqUI,
               &alnSetST,
               &lenErrAryUI,
               &alnScoreL
       ); // Get the alignment
   } // Else doing a Waterman Smith alignment

   if(alnErrUCAry == 0)
   { // If did not have enough memory
       free(refCStr);
       free(refHeadCStr);
       free(queryHeadCStr);
       free(queryCStr);

      fprintf(stderr, "Memory allocation error in aligment step\n");
       exit(-1);
   } // If did not have enough memory

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-6: Make aligned query sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(outFileCStr != 0) outFILE = fopen(outFileCStr, "w");
   else outFILE = stdout;

   queryAlnCStr =
       cnvtAlnErrToSeq(queryCStr, 1, 1, alnErrUCAry, lenErrAryUI);

   if(queryAlnCStr == 0)
   { // If had a memory allocation error
       cnvtAlnErrAryToLetter(refCStr, queryCStr, alnErrUCAry);
       fprintf(outFILE, "%s\n", alnErrUCAry);

       fprintf(stderr, "Memory error, while outputing sequences\n");
       fprintf(stderr, "Printed error array\n");

       free(alnErrUCAry);
       free(refCStr);
       free(refHeadCStr);
       free(queryHeadCStr);
       free(queryCStr);
       fclose(outFILE);
       outFILE = 0;

       exit(-1);
   } // If had a memory allocation error

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-7: Make aligned reference sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   refAlnCStr =
       cnvtAlnErrToSeq(refCStr, 1, 0, alnErrUCAry, lenErrAryUI);

   if(refAlnCStr == 0)
   { // If had a memory allocation error
       cnvtAlnErrAryToLetter(refCStr, queryCStr, alnErrUCAry);
       fprintf(outFILE, "%s\n%s\n", queryAlnCStr, alnErrUCAry);

       fprintf(
            stderr,
            "Memory error, while outputing sequences\n"
            "Printed error array and query array\n"
       );

       free(alnErrUCAry);
       free(queryAlnCStr);
       free(refCStr);
       free(refHeadCStr);
       free(queryHeadCStr);
       free(queryCStr);
       fclose(outFILE);
       outFILE = 0;

       exit(-1);
   } // If had a memory allocation error

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-8: Print out the alignment
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   cnvtAlnErrAryToLetter(refCStr, queryCStr, alnErrUCAry);

   fprintf(outFILE, "# Query = %s", queryHeadCStr);
   fprintf(outFILE, "# Ref = %s", refHeadCStr);
   fprintf(outFILE, "# Eqx = Error line\n");
   fprintf(outFILE, "#   - = is match\n");
   fprintf(outFILE, "#   - X is mismatch\n");
   fprintf(outFILE, "#   - I is insertion\n");
   fprintf(outFILE, "#   - D is deletion\n");
   fprintf(outFILE, "#   - S is soft mask on query and reference\n");
   fprintf(outFILE, "#   - s is soft mask on query only\n");
   fprintf(outFILE, "#   - P is soft mask on reference only\n\n");

   queryPtrCStr = queryAlnCStr;
   refPtrCStr = refAlnCStr;
   tmpUCPtr = alnErrUCAry;

   // This is not very elegent or efficent, but it works for now
   while(*queryPtrCStr != 0)
   { // While I have bases to print out
       charOnUC = 0;
       fprintf(outFILE, "Ref:                 ");
       while(charOnUC < 49)
       { // While I have query bases to print out
           fprintf(outFILE, "%c", *refPtrCStr);
           ++charOnUC;
           ++refPtrCStr;
           if(*refPtrCStr == 0) break;
       } // While I have query bases to print out
       fprintf(outFILE, "\n");

       charOnUC = 0;
       fprintf(outFILE, "Query:               ");
       while(charOnUC < 49)
       { // While I have query bases to print out
           fprintf(outFILE, "%c", *queryPtrCStr);
           ++charOnUC;
           ++queryPtrCStr;
           if(*queryPtrCStr == 0) break;
       } // While I have query bases to print out

       fprintf(outFILE, "\n");

       charOnUC = 0;
       fprintf(outFILE, "Eqx:                 ");
       while(charOnUC < 49)
       { // While I have query bases to print out
           fprintf(outFILE, "%c", *tmpUCPtr);
           ++charOnUC;
           ++tmpUCPtr;
           if(*tmpUCPtr == 0) break;
       } // While I have query bases to print out

       fprintf(outFILE, "\n\n");
   } // While I have bases to print out

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec-9: Clean up
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   fclose(outFILE);
   outFILE = 0;

   free(queryAlnCStr);
   free(refAlnCStr);
   free(alnErrUCAry);
   free(refCStr);
   free(refHeadCStr);
   free(queryHeadCStr);
   free(queryCStr);

   exit(0);
} // main

/*---------------------------------------------------------------------\
| Output: Modifies: Each input variable to hold user input
|  - Returns:
|     o 0: If no errors
|     o Pointer to parameter that was an issue
|     o Pointer to -score-matrix if the scoring matrix is invalid
|  - Prints to stdout when the scoring file is invalid
\---------------------------------------------------------------------*/
char * checkInput(
    int *lenArgsInt,               /*Number arguments user input*/
    char *argsCStr[],              /*Array with user arguments*/
    char **refFileCStr,            /*file name of the reference file*/
    char **queryFileCStr,          /*File name of the query file*/
    char **outFileCStr,            // Name of the output file
    struct alnSet *alnSetST        // Aligment settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: checkInput
   '  - Checks & extracts user input
   '    fun-01 sec-1: Variable declerations
   '    fun-01 sec-2: Look through user input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpCStr = 0;
    char *singleArgCStr = 0;
    unsigned long scoreFileErrUL = 0;
    FILE *scoreFILE = 0;  // For loading the scoring matrix

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-01 Sec-2: Look through user input
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    for(int intArg = 1; intArg < *lenArgsInt; intArg++)
    { /*loop through all user input arguments*/    /*0 is program name*/
        singleArgCStr = *(argsCStr +intArg + 1);   /*supplied argument*/
        tmpCStr = *(argsCStr + intArg);            /*Paramter*/

        if(strcmp(tmpCStr, "-ref") == 0)
            *refFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-query") == 0)
            *queryFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-out") == 0)
            *outFileCStr = singleArgCStr;

        else if(strcmp(tmpCStr, "-gapopen") == 0)
            alnSetST->gapStartPenaltyI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-gapextend") == 0)
            alnSetST->gapExtendPenaltyI =
                strtol(singleArgCStr, &tmpCStr, 10);

        else if(strcmp(tmpCStr, "-use-water") == 0)
        { // Else if disabling match priority
            alnSetST->useNeedleBl = !defUseNeedle;
            --intArg;
        } // Else if disabling match priority

        else if(strcmp(tmpCStr, "-enable-match-priority") == 0)
        { // Else if disabling match priority
            alnSetST->matchPriorityBl = !defMatchPriority;
            --intArg;
        } // Else if disabling match priority

        else if(strcmp(tmpCStr, "-match-ins-del") == 0)
        { // Else if user wants matches->insertions->deletions
            alnSetST->diagnolPriorityC = 0;
            alnSetST->topPriorityC = 1;
            alnSetST->leftPriorityC = 2;
            --intArg;
        } // Else if user wants matches->insertions->deletions

        else if(strcmp(tmpCStr, "-match-del-ins") == 0)
        { // Else if user wants matches->deletions->insertions
            alnSetST->diagnolPriorityC = 0;
            alnSetST->topPriorityC = 2;
            alnSetST->leftPriorityC = 1;
            --intArg;
        } // Else if user wants matches->deletions->insertions

        else if(strcmp(tmpCStr, "-ins-match-del") == 0)
        { // Else if user wants insertions->matches->deletions
            alnSetST->diagnolPriorityC = 1;
            alnSetST->topPriorityC = 0;
            alnSetST->leftPriorityC = 2;
            --intArg;
        } // Else if user wants insertions->matches->deletions

        else if(strcmp(tmpCStr, "-del-match-ins") == 0)
        { // Else if user wants deletions->matches->insertions
            alnSetST->diagnolPriorityC = 1;
            alnSetST->topPriorityC = 2;
            alnSetST->leftPriorityC = 0;
            --intArg;
        } // Else if user wants deletions->matches->insertions

        else if(strcmp(tmpCStr, "-ins-del-match") == 0)
        { // Else if user wants insertions->deletions->matches
            alnSetST->diagnolPriorityC = 2;
            alnSetST->topPriorityC = 0;
            alnSetST->leftPriorityC = 1;
            --intArg;
        } // Else if user wants insertions->deletions->matches

        else if(strcmp(tmpCStr, "-del-ins-match") == 0)
        { // Else if user wants deletions->insertions->matches
            alnSetST->diagnolPriorityC = 2;
            alnSetST->topPriorityC = 1;
            alnSetST->leftPriorityC = 0;
            --intArg;
        } // Else if user wants deletions->insertions->matches

        else if(strcmp(tmpCStr, "-score-matrix") == 0)
        { // else if the user supplied a scoring matrix
            scoreFILE = fopen(singleArgCStr, "r");

            if(scoreFILE == 0)
            { // If I could not open the scoring file
                return tmpCStr; // So user knows invalid

                fprintf(stderr,
                    "-score-matrix %s is an invalid scoring file\n",
                    singleArgCStr
                ); /*Print out the problem*/
            } // If I could not open the scoring file

            scoreFileErrUL = readInScoreFile(alnSetST, scoreFILE);

            if(scoreFileErrUL != 0)
            { // If the scoring file had an error
                fprintf(stderr,
                    "Line: %lu in -score-matrix %s is invalid\n",
                    scoreFileErrUL,
                    singleArgCStr
                ); /*Print out the problem*/

                return tmpCStr;    // invalid file
            } // If the scoring file had an error
        } // else if the user supplied a scoring matrix
            
        else return tmpCStr; // Invalid parameter

        intArg++; /*Move to the parameter, so next input is a flag*/
    } /*loop through all user input arguments*/

    return 0; /*input is valid*/
} /*checkInput*/
