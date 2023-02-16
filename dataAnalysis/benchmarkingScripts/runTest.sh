#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   Sec-1: variable declerations
#   Sec-1: run single test
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Use:
#   Quick testing script for my pipeline

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: variable declerations
#   sec-1 sub-1: variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: variable declerations
#*******************************************************************************

fastaDirStr="$1";
refFastaStr="$2";
prefStr="$3";
scriptDirStr="$(dirname "$0")";
loopBool=0; # 0 is first loop
mapRefQualDbl="$4";
maxErrorRateDbl="$5";
maxClustErrorRateDbl="$7";
readDepthInt=20000;  # depth of reads to simulate

if [[ "$mapRefQualDbl" == "" ]]; then
    mapRefQualDbl=20; # min mapq to keep read when mapped to ref
fi

if [[ "$maxErrorRateDbl" == "" ]]; then
    maxErrorRateDbl=0.07; # min mapq to keep read when mapped to ref
fi

if [[ "$maxClustErrorRateDbl" == "" ]]; then
    maxClustErrorRateDbl=0.07; # min mapq to keep read when mapped to ref
fi

if [[ "$prefStr" == "" ]]; then
    prefStr="out";
fi # if no prefix provided

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1:
#   sec-2 sub-1: Run a single test
#   sec-2 sub-2: Remove the extra headers in --stat.tsv file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Run a single test
#*******************************************************************************

# make log file
printf \
    "fasta\trefs\tminMapq\tminDiff\tminClustDiff\treadDepth\n" \
  >> "$prefStr--log.tsv";

for strFasta in ./"$fastaDirStr/"*.fasta; do
# for all fasta files to simulate reads for

  if [[ ! -f "$strFasta" ]]; then
    continue;
  fi # if is the null case

  bash "$scriptDirStr/runSingleTest.sh" \
      -fasta "$strFasta" \
      -map-ref "$refFastaStr" \
      -ref-mapq "$mapRefQualDbl" \
      -max-dif "$maxErrorRateDbl" \
      -max-clust-read-dif "$maxClustErrorRateDbl" \
      -p "out" \
      -num-reads "$readDepthInt";

  cat \
      < "out--stats.tsv" \
      >> "$prefStr--stats.tsv";

  # add entry to log file
  printf \
      "%s\t%s\t%s\t%s\t%s\t%s\n";
      "$strFasta" \
      "$refFastaStr" \
      "$mapRefQualDbl" \
      "$maxErrorRateDbl" \
      "$maxClustErrorRateDbl" \
      "$readDepthInt" \
    >> "$prefStr--log.tsv";

  rm -r out*;  # remove all the previous rounds files
done # for all fasta files to simulate reads for

#*******************************************************************************
# Sec-2 Sub-2: Remove the extra headers in --stat.tsv file
#*******************************************************************************

# remove extra headers in the statas file
awk \
    'BEGIN{
       OFS="\t";
       getline; # get the first line
       print $0; # print the header
     } # BEGIN block

     { # MAIN
       if($1 != "pipeline")
         print $0;
     } # MAIN
  ' < "$prefStr--stats.tsv" \
  > "$prefStr--tmp--stats.tsv";

# remane tmp--stats to be the orignal stats file
mv \
    "$prefStr--tmp--stats.tsv" \
    "$prefStr--stats.tsv";

exit;
