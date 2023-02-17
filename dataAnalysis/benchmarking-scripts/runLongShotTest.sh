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
    -p "longshot-sim-reads";

simFastqFileStr="$(
  printf "%s" ./"longshot-sim-reads"*".fastq" |
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
#    sec-4 sub-1: run longshot with major ref
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
    > "longshot--majorRef.fasta"; # put major reference in fasta for longshot

/usr/bin/time \
    -f "%e\t%M" \
    -o "longshot--time.tsv" \
  bash "$scriptDirStr/runLongShot.sh" \
    "$simFastqFileStr" \
    "longshot--majorRef.fasta"; # call variants

rm "longshot--majorRef.fasta"; # no longer need reference

#*******************************************************************************
# Sec-4 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants longshot detect, get number of errors
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
        "longshotMajorRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "longshot--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants longshot detect, get number of errors

rm "longshot--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Detect variants with minor ref
#    sec-5 sub-1: run longshot with minor reference
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
    > "longshot--minorRef.fasta";

/usr/bin/time \
    -f "%e\t%M" \
    -o "longshot--time.tsv" \
  bash "$scriptDirStr/runLongShot.sh" \
    "$simFastqFileStr" \
    "longshot--minorRef.fasta"; # call variants using minor variant as ref

rm "longshot--minorRef.fasta"; # no longer need

#*******************************************************************************
# Sec-5 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants longshot detect, get number of errors
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
        "longshotMinorRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "longshot--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants longshot detect, get number of errors

rm "longshot--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Detect varianst with consensus ref
#    sec-6 sub-1: build consensus to use as reference
#    sec-6 sub-2: prep consensus for long shot
#    sec-6 sub-3: run longshot
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
    -v filePrefixStr="longshot--tmp" \
    < "buldingCon.fastq";

# extract first 300 reads
head \
    -n 1200 \
    "longshot--tmp--other-reads.fastq" \
    > "longshot--tmp--other-reads2.fastq";

mv \
    "longshot--tmp--other-reads2.fastq" \
    "longshot--tmp--other-reads.fastq";

# Build the consensuses to use as a reference
bash "$scriptDirStr/../findCoInfections/scripts/buildConsensus.sh" \
    -i "longshot--tmp--other-reads.fastq" \
    -I "longshot--tmp--longest-read.fasta" \
    -p "longshot" \
    -t 3;

rm \
    "longshot--tmp--other-reads.fastq" \
    "longshot--tmp--longest-read.fasta" \
    "longshot--racon-round-1.fasta.fai" \
    "buldingCon.fastq";

#*******************************************************************************
# Sec-6 Sub-2: prep consensus for long shot
#*******************************************************************************

# trim the N's off the cosensuses
sed '
     s/^N*//; # trim anonymous (masked) bases at start of sequence
     s/N*$//; # trim masked bases at end of the sequence
    ' \
    < "longshot--consensus.fasta" \
  > "longshot--consensus2.fasta";

mv \
    "longshot--consensus2.fasta" \
    "longshot--consensus.fasta";

#*******************************************************************************
# Sec-6 Sub-3: run longshot
#*******************************************************************************

# call variants
/usr/bin/time \
    -f "%e\t%M" \
    -o "longshot--time.tsv" \
  bash "$scriptDirStr/runLongShot.sh" \
    "$simFastqFileStr" \
    "longshot--consensus.fasta";

rm \
   -r \
   "longshot--consensus.fasta" \
   "longshot--read-depth.txt" \
   "longshot--medaka" \
   "longshot--racon-round-1.fasta" \
   "longshot--racon-round-1.fasta.map-ont.mmi";

#*******************************************************************************
# Sec-6 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants longshot detect, get number of errors
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
        "longshotConsensus\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "longshot--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants longshot detect, get number of errors

rm "longshot--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Detect variants with more distant ref provided by user
#    sec-7 sub-1: detect variants using longshot & distance refs
#    sec-7 sub-2: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-7 Sub-1: detect variants using longshot & distance refs
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "longshot--time.tsv" \
  bash "$scriptDirStr/runLongShot.sh" \
    "$simFastqFileStr" \
    "$compRefStr";

#*******************************************************************************
# Sec-7 Sub-2: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"variant-"*".fasta"; do
# for all variants longshot detect, get number of errors
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
        "longshotDistRef\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "longshot--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

    rm "$strVar";  # No longer need the variant file (already have data)
done # for all variants longshot detect, get number of errors

rm "longshot--time.tsv"; # no longer need the time file
