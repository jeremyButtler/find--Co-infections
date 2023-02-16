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
mapqDbl=20;

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
    -Q: min mapq to keep read  [20]
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
    -Q) mapqDbl="$2";;
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
#   sec-3 sub-1: filter out low quality & unmapped reads
#   sec-3 sub-2: Find the longest read
#   sec-3 sub-3: Build the consensus genome
#   sec-3 sub-4: Clean up
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: filter out low qualtiy & unmapped reads
#*******************************************************************************

filtlong \
    --min_mean_q 13 \
    "$fastqFileStr" \
  > "$prefixStr--tmp.fastq";

# remove reads that did not map to references
minimap2 \
    --eqx \
    -ax map-ont \
    "$refFastaStr" \
    "$prefixStr--tmp.fastq" |
  samtools view \
    -F 0x04 \
    -b \
    --min-qlen 600 \
    -q "$mapqDbl" |
  samtools fastq \
  > "$prefixStr--tmp--minimap2.fastq";

rm "$prefixStr--tmp.fastq"; # no longer need

#*******************************************************************************
# Sec-3 Sub-2: Find the longest read
#*******************************************************************************

awk \
    -f "$scriptDirStr/../findCoInfections/scripts/extractLongestRead.awk" \
    -v filePrefixStr="$prefixStr--tmp" \
    < "$prefixStr--tmp--minimap2.fastq";

#*******************************************************************************
# Sec-3 Sub-3 Build the consenus genome
#*******************************************************************************

bash "$scriptDirStr/../findCoInfections/scripts/buildConsensus.sh" \
    -i "$prefixStr--tmp--other-reads.fastq" \
    -I "$prefixStr--tmp--longest-read.fasta" \
    -p "$prefixStr" \
    -t 3;
exit;
#*******************************************************************************
# Sec-3 Sub-4: Clean up
#*******************************************************************************

rm \
    "$prefixStr--tmp--longest-read.fasta" \
    "$prefixStr--tmp--other-reads.fastq";
