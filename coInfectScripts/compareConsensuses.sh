#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Compare consensues
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: compareConsensuses.sh
# Use:
#   Determines if there is a consensus in a set of consensus that is similar
#    to the input reference
# Input:
#   -f: fasta file with consensuses to compare reference to           [Required]
#   -ref: fasta file with single reference compare to consensuses     [Required]
#   -t: number threads to use with minimap2                           [2]
#   -min-diff: Min % difference for reference & consensus to be differnt [0.01 (1%)]
#   -min-mismatches: Min % mismatchs for reference to be different    [0.003]
# Output:
#   Prints name of consensus that matches reference to stdout
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#   sec-1 sub-1: user input
#   sec-1 sub-2: help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: User input
#*******************************************************************************

fastaFileStr="";            # fasta file with consensuses to compare refernce to
refFastaFileStr="";         # fasta with reference to compare to consensuses
threadsInt=3;               # number threads to use with minimap2
minDiffDbl=0.01;            # consensuses under 1% difference are the same
minMisDbl=0.003;            # consensuses with < 0.3% mismatches are the same
scriptDirStr="$(dirname "$0")";

#*******************************************************************************
# Sec-1 Sub-2: help message
#*******************************************************************************

helpStr="$(basename "$0") -consensuses file.fasta -ref file.fasta
  Use:
    Determines if there is a consensus in a set of consensus that is similar
     to the input reference
  Input:
    -f: fasta file with consensuses to compare reference to           [Required]
    -ref: fasta file with single reference compare to consensuses     [Required]
    -t: number threads to use with minimap2                           [2]
    -min-diff: Min % difference for reference & consensus to be differnt [0.01 (1%)]
    -min-mismatches: Min % mismatchs for reference to be different    [0.003]
  Output:
    Prints name of consensus that matches reference to stdout
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get and check user input
#   sec-2 sub-1: Get user input
#   sec-2 sub-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Get user input
#*******************************************************************************

while [ $# -gt 0 ]; do
# while have user input arguments

  if [[ "$2" == "" ]]; then
    printf "%s\n%s has no argument\n" "$helpStr" "$1";
    exit;
  fi

  case $1 in
    -f) fastaFileStr="$2";;           # fasta file with consensuses
    -ref) refFastaFileStr="$2";;      # fasta with reference to compare
    -t) threadsInt="$2";;             # number threads to use with minimap2
    -min-diff) minDiffDbl="$2";;      # min % difference between ref & consensus
    -min-mismatches) minMisDbl="$2";; # min % mismatche difference
    -h)printf "%s\n" "$helpStr"; exti;;  # print help message
    ?) printf "%s\n%s is not valid\n" "$helpStr" "$1"; exit;;
  esac

  # move to next arugments
  shift;
  shift;
done # while have user input arguments

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastaFileStr" ]]; then

  { # print out error
    printf "%s is not a valid fasta file\n" \
      "$fastaFileStr";
    printf "Provided a fasta with a set of consensuses with -consensuses\n"
  } # print out error

  exit;
fi # if no valid consensus supplied

if [[ ! -f "$refFastaFileStr" ]]; then

  { # print out error
    printf "%s is not a valid fasta file\n" \
      "$refFastaFileStr";
    printf "Provided a fasta with a reference with -ref\n"
  } # print out error

  exit;
fi # if no valid consensus supplied

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Compare consensues
#   sec-3 sub-1: Compare consensues
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Compare consensues
#*******************************************************************************

printf "%s\n" \
  "$(
    minimap2 \
      --eqx \
      -t "$threadsInt" \
      -a \
      "$refFastaFileStr" \
      "$fastaFileStr" |
    awk -f "$scriptDirStr/printConsensusScores.awk" |
    awk -v minDiffDbl="$minDiffDbl" \
      -v minMisDbl="$minMisDbl" \
      '{ # Main
         if($3 < minMisDbl || $4 < minDiffDbl)
         { # If very similar to a previously build consenuses
           sub(/>[ \t]*/, "", $0);   # remove header marker
           print $1;
           exit;
         } # If very similar to a previously build consenuses
    }' \
)"; # print matching consensuses name to stdout (if similar enough)

exit;
