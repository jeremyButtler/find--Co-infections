#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Trim & bin reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: trimReads.sh
# Use:
#   Trims off start & end of a sequence that do not map to a reference genome.
# Input:
#   -f: Fasta or fastq File with sequences to trim       [Required]
#   -ref: Reference genome(s) to trim the sequences with [Required]
#   -t: Number of threads to use with minimap2           [3]
#   -p: Prefix for output files                          [out]
# Output:
#   File: fastx file with trimed sequences [prefix--reference.fastx]
#     - Creates a separate file for each refernce with mapped sequences
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#   sec-1 sub-1: User input & script variables
#   sec-1 sub-2: Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 sub-1: User input & script variables
#*******************************************************************************

fastxFileStr="";
refFastaStr="";
threadsInt=3;
prefixStr="out";

scriptDirStr="$(dirname "$0")";

helpStr="$(basename "$0")
  Use:
    Trims off start & end of a sequence that do not map to a reference genome.
  Input:
    -f: Fasta or fastq File with sequences to trim       [Required]
    -ref: Reference genome(s) to trim the sequences with [Required]
    -t: Number of threads to use with minimap2           [3]
    -p: Prefix for output files                          [out]
  Output:
    File: fastx file with trimed sequences [prefix--reference.fastx]
      - Creates a separate file for each refernce with mapped sequences
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get and check user input
#   sec-2 sub-1: Get user input
#   sec-2 sub-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Get user input
#*******************************************************************************

while [ $# -gt 0 ]; do                     # whild there is user input to look at
# whild have cmd arguments

  if [[ "$2" == "" ]]; then             # parameter missing argument
    printf "%s\n%s has no argument\n" "$helpStr" "$1";
    exit;
  fi;
        
  case $1 in                                # $1 is the parameter
    -f) fastxFileStr="$2";;                 # file with sequences to trim
    -t) threadsInt="$2";;                   # number of threads to use
    -p) prefixStr="$2";;                    # prefix for output file names
    -ref) refFastaStr="$2";;                # file with reference to trim with
    -h) printf "%s\n" "$helpStr"; exit;; # print help message
    ?) printf "%s\n-%s is not valid\n" "$helpStr" "$1"; exit;;
  esac

  shift;                              # move to argument
  shift;                              # move to next parameter
done # whild have cmd arguments

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastxFileStr" ]]; then
    printf "%s is not a file.\n Input fastx file with sequences using -f\n" \
      "$fastxFileStr";
    exit;
fi # Check if input file was valid

if [[ ! -f "$refFastaStr" ]]; then
    printf "%s is not a file.\n Input fasta file with reference using -ref\n" \
      "$refFastaStr";
    exit;
fi

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Trim & bin reads
#   sec-3 sub-1: Trim & bin reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 sub-1: Trim & bin reads
#*******************************************************************************

# --sam-hit-only: Only output the mapped reads
# --secondary=no: Only output the best alignment

minimap2 \
    --sam-hit-only \
    --secondary=no \
    --eqx \
    -ax map-ont \
    -t "$threadsInt" \
    "$refFastaStr" \
    "$fastxFileStr" |
  awk \
    -v prefStr="$prefixStr" \
    -f "$scriptDirStr/trimMaskedStartEnd.awk";

exit;
