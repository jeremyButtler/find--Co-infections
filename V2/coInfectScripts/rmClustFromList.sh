#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: rmClustFromList.sh
# Use:
#   Removes a percentage of clustered reads from a list of reads
# Input:
#   $1: file with list of reads to remove cluster from              [Required]
#   $2: tsv with reference counts (ref-name\tnum-alignments)        [Required]
#   $3: File with all sequences clustered                           [Required\
#   $4: Percent of reads to remove                                  [0.5 = 50%]
# Output:
#   Modifies: read list (-f) to have the reads from cluster tsv (-clust-tsv)
#             removed
#   Modifies: $2 to have all reads in file (reads not clustered assigned 1)
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#   sec-1 sub-1: variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: variable declerations
#*******************************************************************************

readListStr="$1";     # list with reads to remove
refCountFileStr="$2";   # tsv with refInCluster\tnumReads
clusteredReadsStr="$3"; # list of all reads clustered
percClustRmDbl="$4";  # precent of reads to remove
numClustToRm=0;       # number of clusters to remove (found later)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Check user input
#   sec-2 sub-1: check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: check user input
#*******************************************************************************

if [[ ! -f "$readListStr" ]]; then
  printf \
      "%s is not a file. Provide a list of reads as 1st arugment\n" \
      "$readListStr";
  exit;
fi # if user did not provide a percentage of reads to remove

if [[ ! -f "$refCountFileStr" ]]; then
  printf \
      "%s is not a file. Provide tsv with read\tcluster 2nd argument\n" \
      "$refCountFileStr";
  exit;
fi # if user did not provide a percentage of reads to remove

if [[ ! -f "$clusteredReadsStr" ]]; then
  printf \
    "%s is not a file. Provide txt file with list of reads as 3rd argument\n" \
    "$clusteredReadsStr";
  exit;
fi # if user did nto provide file with list of reads

if [[ "$percClustRmDbl" == "" ]]; then
  percClustRmDbl=0.5; # set to default
fi # if user did not provide a percentage of reads to remove

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Remove reads from read list
#   sec-3 sub-1: Find the number of clusters to discard
#   sec-3 sub-2: Find the clusters & reads to remove
#   sec-3 sub-3: Remove reads from read list
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Find the number of clusters to discard
#*******************************************************************************

# Get the reads that were not kept after scoreReads
awk \
    '{print $1}' \
    < "$refCountFileStr" \
  > "tmp-filter.grep";
grep \
    -v \
    -f "tmp-filter.grep" \
    < "$clusteredReadsStr" |
  sed 's/$/	1/' \
  >> "$refCountFileStr";

sort \
    -n \
    -k 2 \
    < "$refCountFileStr" \
  > "tmp-ref-count.tsv";

rm "tmp-filter.grep";

mv \
    "tmp-ref-count.tsv" \
    "$refCountFileStr";

numClustToRm="$(
  wc \
      -l \
      < "$refCountFileStr" |
    awk \
      -v percToRmDbl="$percClustRmDbl" \
      '{printf "%i", $1 * percToRmDbl}';
)"; # get the number of reads to remove

#*******************************************************************************
# Sec-3 Sub-2: Find the reads to remove
#*******************************************************************************

head \
    -n+"$numClustToRm" \
    < "$refCountFileStr" |
  awk '{print $1}' \
  > "tmp--rm-reads.txt"; # is the half of reads with least amount of alignments

#*******************************************************************************
# Sec-3 Sub-3: Remove reads from read list
#*******************************************************************************

# make temporary file without extracted reads
grep \
    -v \
    -f "tmp--rm-reads.txt" \
    < "$readListStr" \
  > "tmp--read-list.tsv";

rm "tmp--rm-reads.txt"; # no longer need

# rename file
mv "tmp--read-list.tsv" "$readListStr";

exit;
