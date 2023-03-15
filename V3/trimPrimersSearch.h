/*######################################################################
# Name: trimPrimersSearch
# Use:
#   o Sets up hash table of reads with primer mappings, searchs for
#     the reads in the hash table, & trims off the primers from the 
#     target reads.
# Output
#   o Trimmed reads are output to a fastq file
# Includes:
#   - "trimPrimersHash.h"
#   o "fqGetIdsHash.h"
#   o "trimPrimersAVLTree.h"
#   o "defaultSettings.h"
#   o "cStrFun.h"
#   o "fqGetIdsFqFun.h"
#   o "trimPrimersStructs.h"
#   o "fqGetIdsStructs.h"
#   o "cStrToNumberFun.h"
#   o "fqAndFaFun.h"
#   o "FCIStatsFun.h"     (fqAndFqFun.h)
#   o "minAlnStats.h"     (fqAndFaFun.h)
#   o "samEntryStruct.h"  (fqAndFaFun.h->FCIstatsFun.h)
#   o "printError.h"      (fqAndFaFun.h)
# C standard includes:
#   o <stdlib.h>
#   o <stdio.h>
#   o <stdint.h>
######################################################################*/

#ifndef TRIMPRIMERSSEARCH_H
#define TRIMPRIMERSSEARCH_H

#include "trimPrimersHash.h" /*includes fqGetIdsStructs.h*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimPrimersSerach SOH: Start Of Header
'  - fun-1 trimPrimers:
'    o Wrapper functions for a series of functions that map the reads
'      to primers, builds a hash table for mappings with primers,
'      & extracts & trims reads to their primer mappings.
'   - fun-2 primReadsExtReads:
'     o Extract target reads from fastq file with hash table or tree
'   - fun-3 trimAndPrintRead
'     o Trims off primer regions and prints out the untrimmed region.
'       Multiple fastq entrries are printed out for reads with multiple
'       primer targets. 
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 2 if could not open the fasta file
|     o 4 if could not open the fastq file
|     o 8 if could not open the output file
|     o 16 if both primer fasta and fastq file coming from stdin
|     o 32 if was not a valid fastq file
|     o 64 for memory allocation errors
\---------------------------------------------------------------------*/
unsigned char trimPrimers(
    char *faPathCStr,  /*Path to fasta file with primers to map*/
    char *pafFileCStr,     /*Paf with read mappings; skips minimap2*/
    char stdinPafBl,       /*Paf with read mappings from stdin*/
    char *fqPathCStr,      /*Path to fastq file with reads to trim*/
    char *outPathCStr,     /*Path to fastq file to write trimmed reads*/
    char *threadsCStr,     /*Number of threads to use with minimap2*/
    char hashSearchBl      /*1: do hash search, 0: do Tree search*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-1 TOC: trimPrimers
   '   - Wrapper functions for a series of functions that map the reads
   '     to primers, builds a hash table for mappings with primers,
   '     & extracts & trims reads to their primer mappings.
   '   o fun-1 sec-1: Variable declerations
   '   o fun-1 sec-2: Check if files exist
   '   o fun-1 sec-3: Build the tree or hash table for reads
   '   o fun-1 sec-4: extract/trim reads, handle errors, clean up,& exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|    stdout: Prints out reads in hash table
|    Returns:
|      o 1 for success
|      o 32 if was not a valid fastq file
|      o 64 for memory allocation errors
\---------------------------------------------------------------------*/
unsigned char extracAndTrimReads(
    FILE *fqFILE,               /*fastq file to search through*/
    FILE *outFILE,              /*File to write extracted reads to*/
    struct readPrimHash *hashST /*Holds tree/hash table variables*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-2 TOC: primReadsExtReads
   '  - Extract target reads from fastq file with hash table or tree
   '  o fun-2 sec-1: Variable declerations
   '  o fun-2 sec-2: Extract target reads from fastq file
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output: Prints trimmed reads to outFILE
\---------------------------------------------------------------------*/
void trimAndPrintRead(
    struct samEntry *samST,  /*Has buffer and sequence to process*/
    struct readPrim *readIn, /*Has sorted primer coordinates to cut at*/
    FILE *outFILE            /*File to output everything to*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-3 TOC: trimAndPrintRead
   '  - Trims off primer regions and prints out the untrimmed region.
   '    Multiple fastq entrries are printed out for reads with multiple
   '    primer targets. 
   '  o fun-3 sec-1: Variable declerations
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
