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
#   -fasta: Fasta file with two references to simluate reads for     [Required]
#   -map-ref: Fasta file to use in binning reads
#   -num-reads: Number reads to simulate                             [20000]
#   -p: prefix to name output files                                  [out]
#   -seed: Seed to build reads with                                  [1026]
# Output:
#   File: with stats for each pipeline                 ["$prefixStr--stats.tsv"]
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
seedInt=1026;                   # seed for badread
mapRefQualDbl=20;             # min mapq to keep read when mapped to ref
maxErrorRateDbl=0.07;         # max difference between read and reference (7%)
maxClustErrorRateDbl=0.07;    # max diff between reads for clustering (7%)

#*******************************************************************************
# Sec-1 Sub-2: script varaibles
#*******************************************************************************

majVarStr=""; # major variant
minVarStr=""; # minor variant
numMajVarReadsInt=0; # number of major variant reads
numMinVarReadsInt=0; # number of reads from minor variant
numJunkInt=0;        # number of junk sequences
numRandInt=0;        # number of random sequences
pLineStr="";         # Line to print out
refStatStr="";       # holds number of errors in reference
mapRefStr="";        # reference consensus mapped to in reference library
simFastqFileStr="";  # points to simulated fastq file
scriptDirStr="$(dirname "$0")"; # this script's directory

#*******************************************************************************
# Sec-1 Sub-3: Help message
#*******************************************************************************

helpStr="$(basename "$0") -fasta file.fasta -map-ref references.fasta [options]
  Use:
    simlulates read for input pair of references with badread & runs
    findCoInfections scripts on simulated data
  Input:
    -fasta: Fasta file with two references to simluate reads for [Required]
    -map-ref: Fasta file to use in binning reads                     [Required]
    -num-reads: Number of reads to simulate                          [20000]
    -p: prefix to name output files                                  [out]
    -seed: Seed to build reads with                                  [1026]
  Output:
    File: with stats for each pipeline                 [prefix--stats.tsv]

  -ref-mapq: Min mapping quality to keep read in reference mapping [20]
  -max-dif: Max difference between reads in reference alignment    [0.07 (7%)]
  -max-clust-read-dif: max diference beteween reads in clustering  [0.07 (7%)]
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
    -fasta) fastaFileStr="$2";;
    -map-ref) refFastaStr="$2";;
    -num-reads) numReadsInt="$2";;
    -p) prefixStr="$2";;
    -seed) seedInt="$2";;
    -ref-mapq) mapRefQualDbl="$2";;
    -max-dif) maxErrorRateDbl="$2";;  # max difference between query & ref
    -max-clust-read-dif) maxClustErrorRateDbl="$2";; # max diference beteween reads in clustering
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
   printf "%s is not a fasta file\nProvide fasta file as -fasta" "$fastaFileStr";
   exit;
fi # if fasta file not valid

if [[ ! -f "$refFastaStr" ]]; then
  printf "%s is not a fasta file with references\n provide refernces as -map-ref\n" \
      "$refFastaStr";
  exit;
fi # if references to map to are not valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Run badread & find co-infections
#   sec-3 sub-1: Simulate reads with badread
#   sec-3 sub-2: get metadata & major + minor variant ref names
#   sec-3 sub-3: print out the stats file header
#   sec-3 sub-4: build consensuses from all reads
#   sec-3 sub-6: Build consensuses using old pipeline
#   sec-3 sub-7: get stats for old pipeline references
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*****************************************************************************
# Sec-3 Sub-1: Run badread
#*****************************************************************************

bash "$scriptDirStr/twoRefRunBadread.sh" \
    -fasta "$fastaFileStr" \
    -num-reads "$numReadsInt" \
    -seed "$seedInt" \
    -p "$prefixStr";

simFastqFileStr="$(
  printf "%s" ./"$prefixStr"*".fastq" |
  sed 's/\.\///;' \
)"; # get fastq file name

#*******************************************************************************
# Sec-3 Sub-2: Get meta data & reference names
#*******************************************************************************

pLineStr="$(
  awk \
    'BEGIN{FS="_"; OFS="\t"};
     {print $2, $3, $4, $5, $6}
  ' <(printf "%s" "$simFastqFileStr") \
)"; # get meta data from name (id, major ref, minor ref, read maj, reads min)

majVarStr="$(
  printf "%s"  "$pLineStr" |
  awk '{print $2}' \
)"; # get the major references name

minVarStr="$(
  printf "%s" "$pLineStr" |
  awk '{print $3}' \
)"; # get the minor variants name

#*******************************************************************************
# Sec-3 Sub-3: Print out stats file header (makes new stat file)
#*******************************************************************************

{ # stats header print block
  printf "pipeline\tseed\tpercId\tmajorRef\tminorReft\tnoMajorRef\tnoMinorRef";
  printf "\tnoBinMajorRef\tnoBinMinorRef\tnoBinJunk\tnoBinRand\ttime\tmemory";
  printf "\tmappedRef\tacutalRef\tmismatches\tindels\n"
} > "$prefixStr--stats.tsv" # stats header print block

#*****************************************************************************
# Sec-3 Sub-4: build consensus from all reads (ignore co-infections)
#*****************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "$prefixStr--time.tsv" \
  bash "$scriptDirStr/allReadConsensus.sh" \
    -fastq "$simFastqFileStr" \
    -ref "$refFastaStr" \
    -Q "$mapRefQualDbl" \
    -p "$prefixStr--all";

#*****************************************************************************
# Sec-3 Sub-5: get consensus stats
#*****************************************************************************

# Get fastq file with metadata
sed \
    -n \
    'p;n;n;n;' \
    < "$prefixStr--all--tmp--minimap2.fastq" |
  "$scriptDirStr/../findCoInfections/scripts/fastqGrep" \
    -stdin-filt \
    -fastq "$simFastqFileStr" \
  > "$prefixStr--all--bin.fastq";

rm "$prefixStr--all--tmp--minimap2.fastq"; # no longer need

numMajVarReadsInt="$(
  awk \
     -f "$scriptDirStr/getNumRefReads.awk" \
     -v refStr="$majVarStr" \
   < "$prefixStr--all--bin.fastq" \
)"; # find number of major variant reads

numMinVarReadsInt="$(
  awk \
    -f "$scriptDirStr/getNumRefReads.awk" \
    -v refStr="$minVarStr" \
    < "$prefixStr--all--bin.fastq" \
)"; # Find the number of reads in bin from minor variant

numJunkInt="$(
  awk \
    -f "$scriptDirStr/getNumRefReads.awk" \
    -v refStr="junk_seq" \
    < "$prefixStr--all--bin.fastq" \
)"; # Find the number of reads in bin from minor variant

numRandInt="$(
  awk \
    -f "$scriptDirStr/getNumRefReads.awk" \
    -v refStr="random_seq" \
    < "$prefixStr--all--bin.fastq" \
)"; # Find the number of reads in bin from minor variant

refStatStr="$(
  bash "$scriptDirStr/getBestAlignment.sh" \
      "$prefixStr--all--consensus.fasta" \
      "$fastaFileStr" \
)"; # get the errors for the consensus

mapRefStr="$(
  bash "$scriptDirStr/getBestAlignment.sh" \
      "$prefixStr--all--consensus.fasta" \
      "$refFastaStr" |
  awk '{print $1}' \
)"; # Get the closest reference in reference database

# print out data for this consensus
printf "allReads\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$seedInt" \
    "$pLineStr" \
    "$numMajVarReadsInt" \
    "$numMinVarReadsInt" \
    "$numJunkInt" \
    "$numRandInt" \
    "$(cat "$prefixStr--time.tsv")" \
    "$mapRefStr" \
    "$refStatStr" \
  >> "$prefixStr--stats.tsv";

rm "$prefixStr--time.tsv"; # Just to be safe
rm -r "$prefixStr--all"*;  # remove all files made by all pipeline

#*****************************************************************************
# Sec-3 Sub-6: Build consensuses using old pipeline
#*****************************************************************************

mkdir "$prefixStr--tmp";
mv "$simFastqFileStr" "$prefixStr--tmp"; # protect my fastq file


/usr/bin/time \
    -f "%e\t%M" \
    -o "$prefixStr--time.tsv" \
  bash "$scriptDirStr/testOldCoInfections.sh" \
      -fastq "$prefixStr--tmp/$simFastqFileStr" \
      -ref "$refFastaStr" \
      -Q "$mapRefQualDbl" \
      -p "$prefixStr--old";

mv "$prefixStr--tmp/$simFastqFileStr" \
   "$simFastqFileStr";   # move fastq file back to working directory

rm -r "$prefixStr--tmp"; # no longer need

#*******************************************************************************
# Sec-3 Sub-7: Get stats for consensuses build using old coinfection pipeline
#*******************************************************************************

for strCon in ./"$prefixStr--old--"*"consensus.fasta"; do
# for all kept consensuses, get stats

  if [[ ! -f "$strCon" ]]; then
    continue; # null case
  fi # if the null case

  conFastqStr="$(
    printf "%s" "$strCon" |
    sed "
      s/\.\///;
      s/--consensus.fasta//;
      s/^/$prefixStr--old--bins\//;
    " \
  ).fastq"; # find the fastq file used to build ref

  # Get fastq file with metadata
  sed \
      -n \
      'p;n;n;n;' \
      < "$conFastqStr" |
    "$scriptDirStr/../findCoInfections/scripts/fastqGrep" \
      -stdin-filt \
      -fastq "$simFastqFileStr" \
    > "$prefixStr--old--bin.fastq";

  rm "$conFastqStr"; # no longer need

  numMajVarReadsInt="$(
    awk \
       -f "$scriptDirStr/getNumRefReads.awk" \
       -v refStr="$majVarStr" \
       < "$prefixStr--old--bin.fastq" \
  )"; # find number of major variant reads

  numMinVarReadsInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="$minVarStr" \
      < "$prefixStr--old--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  numJunkInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="junk_seq" \
      < "$prefixStr--old--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant
  
  numRandInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="random_seq" \
      < "$prefixStr--old--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  refStatStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$fastaFileStr" \
  )"; # get the errors for the consensus

  mapRefStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$refFastaStr" |
    awk '{print $1}' \
  )"; # Get the closest reference in reference database

  # print out data for this consensus
  printf "oldPipe\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$seedInt" \
      "$pLineStr" \
      "$numMajVarReadsInt" \
      "$numMinVarReadsInt" \
      "$numJunkInt" \
      "$numRandInt" \
      "$(cat "$prefixStr--time.tsv")" \
      "$mapRefStr" \
      "$refStatStr" \
    >> "$prefixStr--stats.tsv";

  rm "$prefixStr--old--bin.fastq"; # make sure removed
done # for all kept consensuses, get stats

rm -r "$prefixStr--old"*; # remove all made files
rm "$prefixStr--time.tsv"; # Just to be safe

#*******************************************************************************
# Sec-3 Sub-8: Test my new pipeline whithout clustering
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "$prefixStr--time.tsv" \
  bash "$scriptDirStr/../findCoInfections/findCoInfections.sh" \
    -f "$simFastqFileStr" \
    -ref "$refFastaStr" \
    -p "$prefixStr--new" \
    -cluster "no" \
    -ref-mapq "$mapRefQualDbl" \
    -max-dif "$maxErrorRateDbl" \
    -max-clust-read-dif "$maxClustErrorRateDbl";
  # lowering min-con-diff to help with 1% case

#*******************************************************************************
# Sec-3 Sub-9: Get stats for new pipeline without clustering
#*******************************************************************************

for strCon in ./"$prefixStr--new--"*"consensus.fasta"; do
# for all kept consensuses, get stats

  if [[ ! -f "$strCon" ]]; then
    continue; # null case
  fi # if the null case

  conFastqStr="$(
    printf "%s" "$strCon" |
    sed "
      s/\.\///;
      s/--consensus.fasta//;
    " \
  ).fastq"; # find the fastq file used to build ref

  # Get fastq file with metadata
  sed \
      -n \
      'p;n;n;n;' \
      < "$conFastqStr" |
    "$scriptDirStr/../findCoInfections/scripts/fastqGrep" \
      -stdin-filt \
      -fastq "$simFastqFileStr" \
    > "$prefixStr--new--bin.fastq";

  rm "$conFastqStr"; # no longer need

  numMajVarReadsInt="$(
    awk \
       -f "$scriptDirStr/getNumRefReads.awk" \
       -v refStr="$majVarStr" \
       < "$prefixStr--new--bin.fastq" \
  )"; # find number of major variant reads

  numMinVarReadsInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="$minVarStr" \
      < "$prefixStr--new--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  numJunkInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="junk_seq" \
      < "$prefixStr--new--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant
  
  numRandInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="random_seq" \
      < "$prefixStr--new--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  refStatStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$fastaFileStr" \
  )"; # get the errors for the consensus

  mapRefStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$refFastaStr" |
    awk '{print $1}' \
  )"; # Get the closest reference in reference database

  # print out data for this consensus
  printf "newPipe\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$seedInt" \
      "$pLineStr" \
      "$numMajVarReadsInt" \
      "$numMinVarReadsInt" \
      "$numJunkInt" \
      "$numRandInt" \
      "$(cat "$prefixStr--time.tsv")" \
      "$mapRefStr" \
      "$refStatStr" \
    >> "$prefixStr--stats.tsv";

  rm "$prefixStr--new--bin.fastq"; # make sure removed
done # for all kept consensuses, get stats

rm -r "$prefixStr--new"*; # remove all files
rm "$prefixStr--time.tsv";

#*******************************************************************************
# Sec-3 Sub-10: Test my new pipeline with clustering
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "$prefixStr--time.tsv" \
  bash "$scriptDirStr/../findCoInfections/findCoInfections.sh" \
    -f "$simFastqFileStr" \
    -ref "$refFastaStr" \
    -p "$prefixStr--new-clust" \
    -ref-mapq "$mapRefQualDbl" \
    -max-dif "$maxErrorRateDbl" \
    -max-clust-read-dif "$maxClustErrorRateDbl";
  # lowering min-con-diff to help with 1% case

#*******************************************************************************
# Sec-3 Sub-11: Get stats for new pipeline with clustering
#*******************************************************************************

for strCon in ./"$prefixStr--new-clust--"*"consensus.fasta"; do
# for all kept consensuses, get stats

  if [[ ! -f "$strCon" ]]; then
    continue; # null case
  fi # if the null case

  conFastqStr="$(
    printf "%s" "$strCon" |
    sed "
      s/\.\///;
      s/--consensus.fasta//;
    " \
  ).fastq"; # find the fastq file used to build ref

  # Get fastq file with metadata
  sed \
      -n \
      'p;n;n;n;' \
      < "$conFastqStr" |
    "$scriptDirStr/../findCoInfections/scripts/fastqGrep" \
      -stdin-filt \
      -fastq "$simFastqFileStr" \
    > "$prefixStr--new-clust--bin.fastq";


  rm "$conFastqStr"; # no longer need

  numMajVarReadsInt="$(
    awk \
       -f "$scriptDirStr/getNumRefReads.awk" \
       -v refStr="$majVarStr" \
       < "$prefixStr--new-clust--bin.fastq" \
  )"; # find number of major variant reads

  numMinVarReadsInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="$minVarStr" \
      < "$prefixStr--new-clust--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  numJunkInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="junk_seq" \
      < "$prefixStr--new-clust--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant
  
  numRandInt="$(
    awk \
      -f "$scriptDirStr/getNumRefReads.awk" \
      -v refStr="random_seq" \
      < "$prefixStr--new-clust--bin.fastq" \
  )"; # Find the number of reads in bin from minor variant

  refStatStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$fastaFileStr" \
  )"; # get the errors for the consensus

  mapRefStr="$(
    bash "$scriptDirStr/getBestAlignment.sh" \
        "$strCon" \
        "$refFastaStr" |
    awk '{print $1}' \
  )"; # Get the closest reference in reference database


  # print out data for this consensus
  printf "newClustPipe\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$seedInt" \
      "$pLineStr" \
      "$numMajVarReadsInt" \
      "$numMinVarReadsInt" \
      "$numJunkInt" \
      "$numRandInt" \
      "$(cat "$prefixStr--time.tsv")" \
      "$mapRefStr" \
      "$refStatStr" \
    >> "$prefixStr--stats.tsv";

  rm "$prefixStr--new-clust--bin.fastq"; # make sure removed
done # for all kept consensuses, get stats

rm -r "$prefixStr--new-clust"*; # remove all files
rm "$prefixStr--time.tsv";
