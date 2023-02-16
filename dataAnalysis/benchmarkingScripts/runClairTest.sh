#!/usr/bin/bash

################################################################################
# Name: runClairTest.sh
# Use:
#    Does clair varaint calling on a single set of simulated reads.
#    Does variant calling using major and minor variant refs and consensuses.
# Input:
#    $1: 
#        References to simulate reads with using badread [Required: fasta]
#    $2:
#        non-variant Reference to call variants with     [Required: fasta]
#    $3:
#        - Number of reads to simulate with badread      [Default: 15000]
# Output:
#    Appends: stats to test--stats.tsv
# Requires:
#    minimap2
#    samtools
#    /usr/bin/time
#    twoRefRunBadread.sh: badread
#    runClair3: Clair3 [hkubal/clair3:latest] (by docker)
#    getBestAlignment.sh: (has an awk script required)
#    ../findCoInfections/scripts/extractLongestRead.awk (for consensus building)
#    ../findCoInfections/scripts/buildConsensus.sh   (for consensus building)
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: 
#    sec-1: Variable declerations
#    sec-2: Check user input
#    sec-3: Simulate reads
#    sec-4: Detect varianst with major ref
#    sec-5: Detect variants with minor ref
#    sec-6: Detect varianst with consensus ref
#    sec-7: Detect variants with more distant ref provided by user
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

simRefStr="$1";   # reference to simulate reads for
compRefStr="$2";  # reference to be used comparing long shot to distant ref
statFileStr="test-stats"; # file to hold data

numReadsInt="$3";  # number of reads to simulate
seedInt=1026;      # seed for badread
varNumInt=0;       # number of variants found
majVarStr="";      # holds name of major variant
minVarStr="";      # holds name of minor variant

scriptDirStr="$(dirname "$0")";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ ! -f "$simRefStr" ]]; then
    printf \
        "%s is not a file, provided a fasta file to simulate reads with\n" \
        "$simRefStr"; 
    exit; 
fi # if no refereence provided to simulate reads with

if [[ "$numReadsInt" -le 1 ]]; then
    printf \
        "Read depth not provided for simulation, simulating 15000 reads\n"
    numReadsInt=15000
fi # if need to use default read read depth

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Simulate reads
#     sec-3 sub-1: simulate reads
#     sec-3 sub-2: get stats from simulated reads
#     sec-3 sub-3: Print out stats file header (makes new stat file)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: simulate reads
#*******************************************************************************

bash "$scriptDirStr/twoRefRunBadread.sh" \
    -fasta "$simRefStr" \
    -num-reads "$numReadsInt" \
    -seed "$seedInt" \
    -p "clair-sim-reads";

simFastqFileStr="$(
  printf "%s" ./"clair-sim-reads"*".fastq" |
  sed 's/\.\///;' \
)"; # get fastq file name

#*******************************************************************************
# Sec-3 Sub-2: get stats from simulated reads
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

if [[ ! -f "$statFileStr.tsv" ]]; then 
    { # stats header print block
      printf "pipeline\tseed\tpercId\tmajorRef\tminorReft\tnoMajorRef";
      printf "\tnoMinorRef\tnoBinMajorRef\tnoBinMinorRef\tnoBinJunk\tnoBinRand";
      printf "\ttime\tmemory\t\tmappedRef\tacutalRef\tmismatches\tindels\n"
    } >> "$statFileStr.tsv" # stats header print block
fi # if need to make the stats file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Detect varianst with major ref
#    sec-4 sub-1: run clair with major ref
#    sec-4 sub-2: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

sed \
    -n \
    "/$majVarStr/{
        p;         # print the first line with major reference (header)
        n;         # move to the sequence line
        p;         # print the seqeuence line of the major reference
        q;         # quit sed (no longer need to run)
     }
    " \
    < "$simRefStr" \
    > "clair--majorRef.fasta"; # put major reference in fasta for clair

/usr/bin/time \
    -f "%e\t%M" \
    -o "clair--time.tsv" \
  bash "$scriptDirStr/runClair3.sh" \
    "$simFastqFileStr" \
    "clair--majorRef.fasta"; # call variants

rm "clair--majorRef.fasta"; # no longer need reference

#*******************************************************************************
# Sec-4 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants clair detect, get number of errors
    if [[ ! -f "$strVar" ]]; then
        continue;
    fi # if null case

    # get the errors in the reference
    refStatStr="$(
      bash "$scriptDirStr/getBestAlignment.sh" \
          "$strVar" \
          "$simRefStr" \
    )"; # get the errors for the consensus

    varNumInt="$( \
        printf \
            "%s" \
            "$strVar" |
          sed '
            s/^\.\/variant-//; # remove file prefix, leaves number.fasta
            s/\.fasta$//;      # remove .fasta, leaving number
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "clairMajorRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "clair--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants clair detect, get number of errors

rm "clair--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Detect variants with minor ref
#    sec-5 sub-1: run clair with minor reference
#    sec-5 sub-2: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

sed \
    -n \
    "/$minVarStr/{
        p;         # print the first line with minor reference (header)
        n;         # move to the sequence line
        p;         # print the seqeuence line of the minor reference
        q;         # quit sed (no longer need to run)
     }
    " \
    < "$simRefStr" \
    > "clair--minorRef.fasta";

/usr/bin/time \
    -f "%e\t%M" \
    -o "clair--time.tsv" \
  bash "$scriptDirStr/runClair3.sh" \
    "$simFastqFileStr" \
    "clair--minorRef.fasta"; # call variants using minor variant as ref

rm "clair--minorRef.fasta"; # no longer need

#*******************************************************************************
# Sec-5 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants clair detect, get number of errors
    if [[ ! -f "$strVar" ]]; then
        continue;
    fi # if null case

    # get the errors in the reference
    refStatStr="$(
      bash "$scriptDirStr/getBestAlignment.sh" \
          "$strVar" \
          "$simRefStr" \
    )"; # get the errors for the consensus

    varNumInt="$( \
        printf \
            "%s" \
            "$strVar" |
          sed '
            s/^\.\/variant-//; # remove file prefix, leaves number.fasta
            s/\.fasta$//;      # remove .fasta, leaving number
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "clairMinorRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "clair--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants clair detect, get number of errors

rm "clair--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Detect varianst with consensus ref
#    sec-6 sub-1: build consensus to use as reference
#    sec-6 sub-2: prep consensus for long shot
#    sec-6 sub-3: run clair
#    sec-6 sub-4: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-6 Sub-1: build consensus to use as reference
#*******************************************************************************

# Filter reads
minimap2 \
    -ax map-ont \
    "$simRefStr" \
    "$simFastqFileStr" |
  samtools \
    view \
    -b \
    --min-qlen 600 \
    -q 20 \
    -F 0x04 |
  samtools \
    fastq \
  > "buldingCon.fastq";

# extract the longest read
awk \
    -f "$scriptDirStr/../findCoInfections/scripts/extractLongestRead.awk" \
    -v filePrefixStr="clair--tmp" \
    < "buldingCon.fastq";

# extract first 300 reads
head \
    -n 1200 \
    "clair--tmp--other-reads.fastq" \
    > "clair--tmp--other-reads2.fastq";

mv \
    "clair--tmp--other-reads2.fastq" \
    "clair--tmp--other-reads.fastq";

# Build the consensuses to use as a reference
bash "$scriptDirStr/../findCoInfections/scripts/buildConsensus.sh" \
    -i "clair--tmp--other-reads.fastq" \
    -I "clair--tmp--longest-read.fasta" \
    -p "clair" \
    -t 3;

rm \
    "clair--tmp--other-reads.fastq" \
    "clair--tmp--longest-read.fasta" \
    "clair--racon-round-1.fasta.fai" \
    "buldingCon.fastq";

#*******************************************************************************
# Sec-6 Sub-2: prep consensus for long shot
#*******************************************************************************

# trim the N's off the cosensuses
sed '
     s/^N*//; # trim anonymous (masked) bases at start of sequence
     s/N*$//; # trim masked bases at end of the sequence
    ' \
    < "clair--consensus.fasta" \
  > "clair--consensus2.fasta";

mv \
    "clair--consensus2.fasta" \
    "clair--consensus.fasta";

#*******************************************************************************
# Sec-6 Sub-3: run clair
#*******************************************************************************

# call variants
/usr/bin/time \
    -f "%e\t%M" \
    -o "clair--time.tsv" \
  bash "$scriptDirStr/runClair3.sh" \
    "$simFastqFileStr" \
    "clair--consensus.fasta";

rm \
   -r \
   "clair--consensus.fasta" \
   "clair--read-depth.txt" \
   "clair--medaka" \
   "clair--racon-round-1.fasta" \
   "clair--racon-round-1.fasta.map-ont.mmi";

#*******************************************************************************
# Sec-6 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants clair detect, get number of errors
    if [[ ! -f "$strVar" ]]; then
        continue;
    fi # if null case

    # get the errors in the reference
    refStatStr="$(
      bash "$scriptDirStr/getBestAlignment.sh" \
          "$strVar" \
          "$simRefStr" \
    )"; # get the errors for the consensus

    varNumInt="$( \
        printf \
            "%s" \
            "$strVar" |
          sed '
            s/^\.\/variant-//; # remove file prefix, leaves number.fasta
            s/\.fasta$//;      # remove .fasta, leaving number
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "clairConsensus\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "clair--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants clair detect, get number of errors

rm "clair--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Detect variants with more distant ref provided by user
#    sec-7 sub-1: detect variants using clair & distance refs
#    sec-7 sub-2: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-7 Sub-1: detect variants using clair & distance refs
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "clair--time.tsv" \
  bash "$scriptDirStr/runClair3.sh" \
    "$simFastqFileStr" \
    "$compRefStr";

#*******************************************************************************
# Sec-7 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants clair detect, get number of errors
    if [[ ! -f "$strVar" ]]; then
        continue;
    fi # if null case

    # get the errors in the reference
    refStatStr="$(
      bash "$scriptDirStr/getBestAlignment.sh" \
          "$strVar" \
          "$simRefStr" \
    )"; # get the errors for the consensus

    varNumInt="$( \
        printf \
            "%s" \
            "$strVar" |
          sed '
            s/^\.\/variant-//; # remove file prefix, leaves number.fasta
            s/\.fasta$//;      # remove .fasta, leaving number
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "clairDistRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "clair--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants clair detect, get number of errors

rm "clair--time.tsv"; # no longer need the time file
