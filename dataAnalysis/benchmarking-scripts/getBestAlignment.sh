#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: variable declerations
#   sec-2: Get the best aligment
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: getBestAligment.sh
# Use:
#   Alignes a consensus to a reference and keeps alignment with fewest errors
# Input:
#   $1: consensus to get errors for
#   $2: Reference(s) to aligne consensus to
# Output:
#   stdout: reference-name\tnumber-mismatches\tnumber-indels
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#    sec-1 sub-1: varaibles holding user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: Variables holding user input
#*******************************************************************************

fastaFileStr="$1";
refFastaStr="$2";               # references to map reads to
scriptDirStr="$(dirname "$0")"; # this script's directory

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Check user input
#   sec-2 sub-1: check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Check user input
#*******************************************************************************

if [[ ! -f "$fastaFileStr" ]]; then
   printf "%s is not a fasta file\nProvide fasta file \$1\n" "$fastaFileStr";
   exit;
fi # if fasta file not valid

if [[ ! -f "$refFastaStr" ]]; then
  printf "%s is not a fasta file with referenceswith \n provide refernces \$2\n" \
      "$refFastaStr";
  exit;
fi # if references to map to are not valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Align consensus to references & find top scoring alignment
#   sec-3 sub-1: Align consensus to references & find top scoring alignment
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Align consensus to references & find top scoring alignemnt
#*******************************************************************************

minimap2 \
    --eqx \
    -ax map-ont \
    "$refFastaStr" \
    "$fastaFileStr" |
  awk \
    -f "$scriptDirStr/getMisIndelCnt.awk" |
  awk \
    '{ # MAIN
     if(NR == 1)
     { # if on the first row
       errorInt = $2 + $3; # mis + indel
       lineStr = $0; # save the alignment
       next;
     } # if on the first row

     if($2 + $3 < errorInt)
     { # if less errors on this row
       errorInt = $2 + $3; # mis + indel
       lineStr = $0; # better alignemnt
       next;
     } # if this alignment had less errors
     } # MAIN

    END{print lineStr}
'
