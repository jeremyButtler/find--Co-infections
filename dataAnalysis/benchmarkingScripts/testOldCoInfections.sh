#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: variable declerations
#   sec-2: Get and check user input
#   sec-3: Build consensus
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: allReadConsensus.sh
# Use:
#   builds consensus genome from input reads
# Input:
#   -fastq: Fastq file to build consensus from
#   -ref: Reference used to filter out unmapped reads
#   -p: prefix to name output files                           [out]
# Output:
#   file: fasta file with consensus genome
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#    sec-1 sub-1: varaibles holding user input
#    sec-1 sub-2: script variables
#    sec-1 sub-3: Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: Variables holding user input
#*******************************************************************************

fastqFileStr="";              # fasta file with references
refFastaStr="";               # references to map reads to
prefixStr="out";

#*******************************************************************************
# Sec-1 Sub-2: script varaibles
#*******************************************************************************

scriptDirStr="$(dirname "$0")"; # this script's directory

#*******************************************************************************
# Sec-1 Sub-3: Help message
#*******************************************************************************

helpStr="$(basename "$0") -fasta file.fasta -ref references.fasta [options]
  Use:
    builds consensus genome from input reads
  Input:
    -fastq: Fastq file to build consensus from
    -ref: Reference used to filter out unmapped reads
    -p: prefix to name output files                           [out]
  Output:
    file: fasta file with consensus genome
"; # help message

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get user input
#   sec-2 sub-1: get user input
#   sec-2 sub-2: check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: get user input
#*******************************************************************************

while [ $# -gt 0 ]; do
# While there is input to read

  if [[ "$2" == "" ]]; then
    printf "%s\n%s has no arguments\n" \
        "$helpStr" \
        "$1";
    exit;
  fi # if argument is blank

  case $1 in
    -fastq) fastqFileStr="$2";;
    -ref) refFastaStr="$2";;
    -p) prefixStr="$2";;
    -h) printf "%s\n" \
            "$helpStr" &&
          exit;;
   ?) printf "%s\n%s is not valid\n" \
          "$helpStr" \
          "$1" &&
        exit;;
  esac

  shift;  # move to parameter
  shift;  # move to next argument
done # while their are user arguemnts to check

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastqFileStr" ]]; then
   printf "%s is not a fastq file\nProvide fastq file as -fastq" "$fastqFileStr";
   exit;
fi # if fasta file not valid

if [[ ! -f "$refFastaStr" ]]; then
  printf "%s is not a fasta file with references\n provide refernces as -ref\n" \
      "$refFastaStr";
  exit;
fi # if references to map to are not valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Build consensus
#   sec-3 sub-1: Build consensus with old findCoInfections pipeline
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Run old co-infection pipeline
#*******************************************************************************
bash "$scriptDirStr/../oldFindCoInfections/binAndBuildConsensus.sh" \
    -i "$fastqFileStr" \
    -r "$refFastaStr" \
    -p "$prefixStr" \
    -t 3;

exit;

