#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: variable declerations
#   sec-2: read user input
#   sec-3: run badread
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: twoRefRunBadread.sh
# Use:
#   Runs badread on fasta file with two rerferences
# Input:
#   -fasta: fasta file name to run badread on [Required]
#   -num-reads: Number reads to simulate      [20000]
#   -seed: seed to use with badread           [1026]
#   -p: prefix to name file                   [out]
# Output:
#    File: fastq file with simulated read
#      - percent-identity_major-variant_minor-variant_\
#        number-major-variant-reads_number-minor-var-reads_\
#        reads.fastq
################################################################################


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#   sec-1 sub-1: user input
#   sec-1 sub-2: script variables
#   sec-1 sub-3: help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: Variables holding user input
#*******************************************************************************

fastaFileStr="";        # fasta file with two references
prefixStr="out";        # prefix to name everything
numReadsInt=20000;      # number reads to simulate from refs
seedInt=1026;           # seed for badread

#*******************************************************************************
# Sec-1 Sub-2: Script variables
#*******************************************************************************

percIdDbl=0;
majVarStr="";        # holds major variant name
minVarStr="";        # holds minor variant name
numMajVarInt=0;      # Holds number of reads from major variant
numMinVarInt=0;      # Holds number reads from minor var
scriptDirStr="$(dirname "$0")"; # this script's directory

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -fasta file.fasta [options ...]
  Use:
    Runs badread on fasta file with two rerferences
  Input:
    -fasta: fasta file name to run badread on [Required]
    -num-reads: Number reads to simulate      [20000]
    -seed: seed to use with badread           [1026]
    -p: prefix to name file                   [out]
  Output:
     File: fastq file with simulated read
       - Name: percent-identity_major-variant_minor-variant_\
         number-major-variant-reads_number-minor-var-reads_\
         reads.fastq
"; # help message

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get user input
#   sec-2 sub-1: get user input
#   sec-2 sub-2: check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Get user input
#*******************************************************************************

while [ $# -gt 0 ]; do
# while there is user input to check

  if [[ "$2" == "" ]]; then
    printf "%s\n%s has no argument\n" \
        "$helpStr" \
        "$2";
    exit;
  fi # if parameter is missing argument

  case $1 in
    -fasta) fastaFileStr="$2";;    # references to simulate reads for
    -num-reads) numReadsInt="$2";; # number of reads to simulate
    -seed) seedInt="$2";;          # seed for badread
    -p)prefixStr="$2";;            # prefix to name everything
    -h) printf "%s\n" \
            "$helpStr" &&
          exit;;
    ?) printf "%s\n%s is not valid\n" \
           "$helpStr" \
           "$fastaFileStr" &&
         exit;;
  esac

  shift; # move to argument
  shift; # move to next parameter
done # while there is user input to check

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastaFileStr" ]]; then
  printf "%s is not a fasta file.\nProvide references with -fasta\n" \
      "$fastaFileStr";
  exit;
fi # if references do not exist

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Run tests
#   sec-3 sub-1: Get minor & major variant names
#   sec-3 sub-2: Extract the major variant
#   sec-3 sub-3: Extract the minor variant
#   sec-3 sub-4: Find percent idenitity between variants
#   sec-3 sub-5: run badread
#   sec-3 sub-6: Get number reads for major & minor variant
#   sec-3 sub-7: Rename fastq file to hold metadata
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Get major & minor variant names
#*******************************************************************************

# Find reference names
majVarStr="$(
  sed -n '
    1{
       s/^>[ \t]*//;
       s/[ \t].*//;
       p;
       q;
    }' < "$fastaFileStr"
)"; # Get major variant name

minVarStr="$( \
  awk \
    'BEGIN{getline;}; # move of first sequence

     { # MAIN block
      if($0 ~ /^>/)
      { # if found the second reference
        sub(/^>/, "", $1);   # make sure > not in name
        print $1;
        exit;
      } # if found the second reference
     } # MAIN block
  ' < "$fastaFileStr" \
)"; # get minor varaint name


#*******************************************************************************
# Sec-3 Sub-2: Extract the major variant
#*******************************************************************************

# Make major variant file
awk \
    'BEGIN{
       getline;
       print $0;
     }; # BEGIN block

     { # MAIN block
       if($1 ~ /^>/){exit;};
       print $0;
     } # MAIN block
  ' < "$fastaFileStr" \
 > "$prefixStr--tmp--major.fasta";

#*******************************************************************************
# Sec-3 Sub-3: Extract minor variant sequence
#*******************************************************************************

# Make minor variant file
awk \
    'BEGIN{
      getline; # move onto header (first line)
      getline; # move off header
     } # BEGIN block

     { # MAIN block
       if($1 ~ /^>/){pBool = 1;};
       if(pBool == 1){print $0}
     } # MAIN block
  ' < "$fastaFileStr" \
  > "$prefixStr--tmp--minor.fasta";

#*******************************************************************************
# Sec-3 Sub-4: Find percent idenitity between variants
#*******************************************************************************

# find percent idenity
percIdDbl="$( \
  minimap2 \
      -a \
      "$prefixStr--tmp--major.fasta" \
      "$prefixStr--tmp--minor.fasta" |
    awk '{
      if($1 ~ /^@/){next;};
      sub(/NM:i:/, "", $12); # remove uneeded characters
      printf "%.3f", $12 / length($10); # $12 is indeles + mis
      exit;                    # only one match
    }' |
  sed '
    s/0\.//;
    s/\([0-9][0-9]\)/\1./;
  ' \
)"; # get percent id between both refs

rm "$prefixStr--tmp--major.fasta" \
   "$prefixStr--tmp--minor.fasta";


#*******************************************************************************
# Sec-3 Sub-5: Run bad read
#*******************************************************************************

badread \
    simulate \
    --reference "$fastaFileStr" \
    --quantity "$numReadsInt"x \
    --seed "$seedInt" \
  > "$prefixStr--tmp.fastq";

#*******************************************************************************
# Sec-3 Sub-6: Get number minor & major varaint reads
#*******************************************************************************

numMajVarInt="$( \
  awk \
    -f "$scriptDirStr/getNumRefReads.awk" \
    -v refStr="$majVarStr" \
  < "$prefixStr--tmp.fastq" \
)"; # get the number of major variant reads

numMinVarInt="$( \
  awk \
    -f "$scriptDirStr/getNumRefReads.awk" \
    -v refStr="$minVarStr" \
  < "$prefixStr--tmp.fastq" \
 )"; # get the number of minor variant reads

#*******************************************************************************
# Sec-3 Sub-7: Rename file
#*******************************************************************************

mv "$prefixStr--tmp.fastq" \
  "$prefixStr""_$percIdDbl""_$majVarStr""_$minVarStr""_$numMajVarInt""_$numMinVarInt""_reads.fastq";

exit;
