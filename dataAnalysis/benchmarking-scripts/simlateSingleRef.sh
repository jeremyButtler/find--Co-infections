#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: variable declerations
#   sec-2: Get and check user input
#   sec-3: simulate reads & run scripts
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name:
# Use:
#   simlulates read for input pair of references with badread & runs
#   findCoInfections scripts on simulated data
# Input:
#   -sim-fasta: Fasta file with two references to simluate reads for [Required]
#   -map-ref: Fasta file to use in binning reads
#   -num-reads: Number reads to simulate                      [20000]
#   -p: prefix to name output files                           [out]
#   -depth: comma deliminated strings with depths to check       ["50,25,10,5,1"]
#   -seed: Seed to build reads with                              [1026]
# Output:
#   stdout: tsv with stats for the simulated run
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

fastaFileStr="";              # fasta file with references
refFastaStr="";               # references to map reads to
numReadsInt=20000;            # number of reads to simulate with badread
prefixStr="out";
depthsStr="50,25,10,5,1";                 # comma deliminated strings with depths to check
seedInt=1026;                   # seed for badread

#*******************************************************************************
# Sec-1 Sub-2: script varaibles
#*******************************************************************************

pIdDbl=0; # percent idenity between references
majVarStr=""; # major variant
minVarStr=""; # minor variant
numMajVarReadsInt=0; # number of major variant reads
numMinVarReadsInt=0; # number of reads from minor variant
pLineStr="";         # Line to print out
scriptDirStr="$(dirname "$0")"; # this script's directory

#*******************************************************************************
# Sec-1 Sub-3: Help message
#*******************************************************************************

helpStr="$(basename "$0") -fasta file.fasta -map-ref references.fasta [options]
  Use:
    simlulates read for input pair of references with badread & runs
    findCoInfections scripts on simulated data
  Input:
    -sim-fasta: Fasta file with two references to simluate reads for [Required]
    -map-ref: Fasta file to use in binning reads                     [Required]
    -num-reads: Number of reads to simulate                          [20000]
    -p: prefix to name output files                                  [out]
    -depth: comma deliminated strings with depths to check      ["50,25,10,5,1"]
    -seed: Seed to build reads with                                  [1026]
  Output:
    stdout: tsv with stats for the simulated run
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
    -sim-fasta) fastaFileStr="$2";;
    -map-ref) refFastaStr="$2";;
    -num-reads) numReadsInt="$2";;
    -p) prefixStr="$2";;
    -depth) depthStr="$2";;
    -seed) seedInt="$2";;
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

if [[ ! -f "$fastaFileStr" ]]; then
   printf "%s is not a fasta file\nProvide fasta file as -sim-fasta" "$fastaFileStr";
   exit;
fi # if fasta file not valid

if [[ ! -f "$refFastaStr" ]]; then
  printf "%s is not a fasta file with references\n provide refernces as -map-ref\n" \
      "$refFastaStr";
  exit;
fi # if references to map to are not valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Run badread & find co-infections
#   sec-3 sub-1: put read depth ratios in temporary file
#   sec-3 sub-2: Simulate reads with badread
#   sec-3 sub-3: build consensuses from all reads
#   sec-3 sub-3: build consensues from 300 reads, ignoring co-infections
#   sec-3 sub-4: find co-infections & build consensuses using old pipeline
#   sec-3 sub-5: find co-infections & build consensuses with new pipepline
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: put read depth ratios in temporary file
#*******************************************************************************

for depthDbl in "$(printf "$depthStr" | awk 'BEGIN{FS=","}; {print $0}')"; do
# for all depths run badread

  awk \
      -v depthDbl="$depthDbl" \
      'BEGIN{
         print $0, " depth=" depthDbl; # print out depth for badread
         next;                          # move to sequence line
       } # BEGIN block

       { # MAIN block
         if($0 ~ /^>/)
         { # if on the next sequences header
            print $0, " depth=" 100 - depthDbl; # print out minor var depth
            next;                           # move to the next sequence
         } # if on the next sequences header

         print $0;         # on sequence
        }' \
    > tmp.fasta


  #*****************************************************************************
  # Sec-3 Sub-2: Run badread
  #*****************************************************************************

  bash "$scriptDirStr/twoRefRunBadread.sh" \
      -fasta tmp.fasta \
      -num-reads "$numReadsInt" \
      -seed "$seedInt" \
      -p "$prefixStr--badread";

  #*****************************************************************************
  # Sec-3 Sub-3: build consensus from all reads (ignore co-infections)
  #*****************************************************************************

  filtlong \
      --min_mean_q 13 \
      "$prefixStr--badread"*".fastq"  \
    > "$prefixStr--tmp.fastq";

  # remove reads that did not map to references
  minimap2 \
      --eqx \
      -ax map-ont \
      "$refFastaStr" \
      "$prefixStr--tmp.fastq" |
    samtools view -F 0x04 |
    samtools fastq \
    > "$prefixStr--tmp--minimap2.fastq";

  rm "$prefixStr--tmp.fastq"; # no longer need

  # find the longest read
  awk \
      -f "$scirptDirStr/findCoInfections/scripts/extractLongestRead.awk" \
      -v filePrefixStr="$prefixStr--badread" \
      < "$prefixStr--tmp--minimap2.fastq";

  "$scriptDirStr/findCoInfections/scripts/buildConsensus.sh \
      -i "$prefixStr--badread--other-reads.fastq" \
      -I "$prefixStr--badread--longest-read.fasta \
      -p "$prefixStr--badread" \
      -t 3;
  # clean up
  rm \
      "$prefixStr--badread--longest-read.fasta" \
      "$prefixStr--tmp--other-reads.fastq" \
      "$prefixStr--tmp--minimap2.fastq";

  # get scores for consensuses
  pLineStr="$(
    minimap2 \
      --eqx \
      -a \
      -t 3 \
      "$refFastaStr" \
      "$prefixStr--badread"*"consensus.fasta" |
    awk \
      -f "$scriptDirStr/getMisIndelCnt.awk"
  )";
