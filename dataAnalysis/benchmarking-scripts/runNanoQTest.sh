#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

simRefStr="$1";   # reference to simulate reads for
refDbStr="$2";    # database of references to run NanoQ with
statFileStr="Nano-Q--stats"; # file to hold data

numReadsInt="$3";  # number of reads to simulate
minMapqDbl=20;      # mapping quality to use
medHamDistInt=170;     # hamming distance for Nano-Q (80-95 < 5%, 160-180 5-15%)
minClustDepthInt=100; # at least 100 reads per cluster
seedInt=1026;      # seed for badread
trimReadsToInt=1100; # length for Nano-Q to trim reads to
jumpInt=500;       # have Nano-Q process 500 reads per thread
varNumInt=0;       # number of variants found
majVarStr="";      # holds name of major variant
minVarStr="";      # holds name of minor variant
prefixStr="Nano-Q"; # name of directory NanoQ outputs
threadsInt=4;      # number threads to use

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
    -p "Nano-Q-sim-reads";

simFastqFileStr="$(
  printf "%s" ./"Nano-Q-sim-reads"*".fastq" |
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
#    sec-4 sub-1: run Nano-Q with major ref
#    sec-4 sub-2: get percent similarity to consensus & print out data
#    sec-4 sub-3: Run major consensus with median (190) hamming distance
#    sec-4 sub-4: get percent similarity to consensus & print out data
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
    > "Nano-Q--majorRef.fasta"; # put major reference in fasta for Nano-Q

#/usr/bin/time \
#    -f "%e\t%M" \
#    -o "Nano-Q--time.tsv" \
#  bash "$scriptDirStr/runNanoQ.sh" \
#    -ref "Nano-Q--majorRef.fasta" \
#    -fastq "$simFastqFileStr" \
#    -min-q "$minMapqDbl" \
#    -p "$prefixStr" \
#    -j "$jumpInt" \
#    -ht "$closeHamDistInt" \
#    -t "$threadsInt" \
#    -mt "$minClustDepthInt" \
#    -racon \
#    -l "$trimReadsToInt";
#
##*******************************************************************************
## Sec-4 Sub-2: get percent similarity to consensus & print out data
##*******************************************************************************
#
#for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
## for all variants Nano-Q detect, get number of errors
#    printf "HI %s\n" "$strVar";
#
#    if [[ ! -f "$strVar" ]]; then
#        continue;
#    fi # if null case
#
#    # get the errors in the reference
#    refStatStr="$(
#      bash "$scriptDirStr/getBestAlignment.sh" \
#          "$strVar" \
#          "$simRefStr" \
#    )"; # get the errors for the consensus
#
#    varNumInt="$( \
#        printf \
#            "%s" \
#            "$strVar" |
#          sed '
#            s/^.*\///; # remove directory info
#            s/Cluster//; # remove Cluster leavning #_Consensus.fa
#            s/_Con.*$//; # remove Ending leaving ref + cluster
#          '\
#    )"; # get the variant number
#
#    # print out data for this consensus
#    printf \
#        "Nano-QMajorRef-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
#        "$closeHamDistInt" \
#        "$seedInt" \
#        "$pLineStr" \
#        "$(cat "Nano-Q--time.tsv")" \
#        "$varNumInt" \
#        "$refStatStr" \
#      >> "$statFileStr.tsv";
#done # for all variants Nano-Q detect, get number of errors
#
#rm \
#    -r \
#    "$prefixStr" \
#    "Nano-Q--time.tsv"; # no longer need the time file
#
#*******************************************************************************
# Sec-4 Sub-3: Run major consensus with median (190) hamming distance
#*******************************************************************************


/usr/bin/time \
    -f "%e\t%M" \
    -o "Nano-Q--time.tsv" \
  bash "$scriptDirStr/runNanoQ.sh" \
    -ref "Nano-Q--majorRef.fasta" \
    -fastq "$simFastqFileStr" \
    -min-q "$minMapqDbl" \
    -p "$prefixStr" \
    -j "$jumpInt" \
    -ht "$medHamDistInt" \
    -t "$threadsInt" \
    -mt "$minClustDepthInt" \
    -racon \
    -l "$trimReadsToInt";

#*******************************************************************************
# Sec-4 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
# for all variants Nano-Q detect, get number of errors

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
            s/^.*\///; # remove directory info
            s/Cluster//; # remove Cluster leavning #_Consensus.fa
            s/_Con.*$//; # remove Ending leaving ref + cluster
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "Nano-QMajorRef-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$medHamDistInt" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "Nano-Q--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";
done # for all variants Nano-Q detect, get number of errors

rm \
    -r \
    "$prefixStr" \
    "Nano-Q--time.tsv"; # no longer need the time file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Detect variants with minor ref
#    sec-5 sub-1: run Nano-Q with minor reference
#    sec-5 sub-2: get percent similarity to consensus & print out data
#    sec-5 sub-3: Find co-infections with Nano-Q median distance
#    sec-5 sub-4: get percent similarity to consensus & print out data
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
    > "Nano-Q--minorRef.fasta";

#/usr/bin/time \
#    -f "%e\t%M" \
#    -o "Nano-Q--time.tsv" \
#  bash "$scriptDirStr/runNanoQ.sh" \
#    -ref "Nano-Q--minorRef.fasta" \
#    -fastq "$simFastqFileStr" \
#    -min-q "$minMapqDbl" \
#    -p "$prefixStr" \
#    -j "$jumpInt" \
#    -ht "$closeHamDistInt" \
#    -t "$threadsInt" \
#    -mt "$minClustDepthInt" \
#    -racon \
#    -l "$trimReadsToInt";
#
#
##*******************************************************************************
## Sec-5 Sub-2: get percent similarity to consensus & print out data
##*******************************************************************************
#
#for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
## for all variants Nano-Q detect, get number of errors
#    if [[ ! -f "$strVar" ]]; then
#        continue;
#    fi # if null case
#
#    # get the errors in the reference
#    refStatStr="$(
#      bash "$scriptDirStr/getBestAlignment.sh" \
#          "$strVar" \
#          "$simRefStr" \
#    )"; # get the errors for the consensus
#
#    varNumInt="$( \
#        printf \
#            "%s" \
#            "$strVar" |
#          sed '
#            s/^.*\///; # remove directory info
#            s/Cluster//; # remove Cluster leavning #_Consensus.fa
#            s/_Con.*$//; # remove Ending leaving ref + cluster
#          '\
#    )"; # get the variant number
#
#    # print out data for this consensus
#    printf \
#        "Nano-QMinorRef-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
#        "$closeHamDistInt" \
#        "$seedInt" \
#        "$pLineStr" \
#        "$(cat "Nano-Q--time.tsv")" \
#        "$varNumInt" \
#        "$refStatStr" \
#      >> "$statFileStr.tsv";
#done # for all variants Nano-Q detect, get number of errors
#
#rm \
#    -r \
#    "$prefixStr" \
#    "Nano-Q--time.tsv"; # no longer need the time file

#*******************************************************************************
# Sec-5 Sub-3: Find co-infections with Nano-Q median distance
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "Nano-Q--time.tsv" \
  bash "$scriptDirStr/runNanoQ.sh" \
    -ref "Nano-Q--minorRef.fasta" \
    -fastq "$simFastqFileStr" \
    -min-q "$minMapqDbl" \
    -p "$prefixStr" \
    -j "$jumpInt" \
    -ht "$medHamDistInt" \
    -t "$threadsInt" \
    -mt "$minClustDepthInt" \
    -racon \
    -l "$trimReadsToInt";


#*******************************************************************************
# Sec-5 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
# for all variants Nano-Q detect, get number of errors
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
            s/^.*\///; # remove directory info
            s/Cluster//; # remove Cluster leavning #_Consensus.fa
            s/_Con.*$//; # remove Ending leaving ref + cluster
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "Nano-QMinorRef-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$medHamDistInt" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "Nano-Q--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";
done # for all variants Nano-Q detect, get number of errors

rm \
    -r \
    "$prefixStr" \
    "Nano-Q--time.tsv"; # no longer need the time file


rm "Nano-Q--minorRef.fasta"; # no longer need

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Detect varianst with consensus ref
#    sec-6 sub-1: build consensus to use as reference
#    sec-6 sub-2: prep consensus for long shot
#    sec-6 sub-3: run Nano-Q
#    sec-6 sub-4: get percent similarity to consensus & print out data
#    sec-6 sub-5: run Nano-Q for median hamming distance
#    sec-6 sub-6: get percent similarity to consensus & print out data
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
    -v filePrefixStr="Nano-Q--tmp" \
    < "buldingCon.fastq";

# extract first 300 reads
head \
    -n 1200 \
    "Nano-Q--tmp--other-reads.fastq" \
    > "Nano-Q--tmp--other-reads2.fastq";

mv \
    "Nano-Q--tmp--other-reads2.fastq" \
    "Nano-Q--tmp--other-reads.fastq";

# Build the consensuses to use as a reference
bash "$scriptDirStr/../findCoInfections/scripts/buildConsensus.sh" \
    -i "Nano-Q--tmp--other-reads.fastq" \
    -I "Nano-Q--tmp--longest-read.fasta" \
    -p "Nano-Q" \
    -t 3;

rm \
    "Nano-Q--tmp--other-reads.fastq" \
    "Nano-Q--tmp--longest-read.fasta" \
    "Nano-Q--racon-round-1.fasta.fai" \
    "buldingCon.fastq";

#*******************************************************************************
# Sec-6 Sub-2: prep consensus for long shot
#*******************************************************************************

# trim the N's off the cosensuses
sed '
     s/^N*//; # trim anonymous (masked) bases at start of sequence
     s/N*$//; # trim masked bases at end of the sequence
    ' \
    < "Nano-Q--consensus.fasta" \
  > "Nano-Q--consensus2.fasta";

# Trim the consensues
minimap2 \
     -a \
     --eqx \
     --secondary=no \
     "Nano-Q--majorRef.fasta" \
     "Nano-Q--consensus2.fasta" |
    awk '
        { # MAIN
            if($0 ~ /^@/)
                next; # ignore the headers

            # split cigar into separate entries
            lenCigInt = split($6, cigAryStr, "=");

            lenSeqInt = length($10);              # length of sequence
            seqStr = $10;             # save sequence

            if(cigAryStr[lenCigInt] ~ /S$/)
            { # if need to trim the consensus
                sub(/S/, "", cigAryStr[lenCigInt]); # get number bases to trim
                trimInt = lenSeqInt - cigAryStr[lenCigInt]; # get stop of trim

                seqStr = substr(seqStr, 1, trimInt);  # trim bases from end
            } # if need to trim the consensus

            if(cigAryStr[1] ~ /S/)
            { # if need to trim the consensus at the start
                sub(/S/, "", cigAryStr[1]);    # get number of bases to trim
                seqStr = substr(seqStr, 1, cigAryStr[1] + 1); # trim bases
            } # if need to trim the consensus at the start

            printf ">Consensus\n%s\n", seqStr;
	    exit;
        } # MAIN
' > "Nano-Q--consensus.fasta";

# remove files that I no longer need
rm \
    "Nano-Q--majorRef.fasta" \
    "Nano-Q--consensus2.fasta";

##*******************************************************************************
## Sec-6 Sub-3: run Nano-Q
##*******************************************************************************
#
## call variants
#/usr/bin/time \
#    -f "%e\t%M" \
#    -o "Nano-Q--time.tsv" \
#  bash "$scriptDirStr/runNanoQ.sh" \
#    -ref "Nano-Q--consensus.fasta" \
#    -fastq "$simFastqFileStr" \
#    -min-q "$minMapqDbl" \
#    -p "$prefixStr" \
#    -j "$jumpInt" \
#    -ht "$closeHamDistInt" \
#    -t "$threadsInt" \
#    -mt "$minClustDepthInt" \
#    -racon \
#    -l "$trimReadsToInt";
#
#rm \
#   -r \
#   "Nano-Q--read-depth.txt" \
#   "Nano-Q--medaka" \
#   "Nano-Q--racon-round-1.fasta" \
#   "Nano-Q--racon-round-1.fasta.map-ont.mmi";
#
##*******************************************************************************
## Sec-6 Sub-4: get percent similarity to consensus & print out data
##*******************************************************************************
#
#for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
## for all variants Nano-Q detect, get number of errors
#    if [[ ! -f "$strVar" ]]; then
#        continue;
#    fi # if null case
#
#    # get the errors in the reference
#    refStatStr="$(
#      bash "$scriptDirStr/getBestAlignment.sh" \
#          "$strVar" \
#          "$simRefStr" \
#    )"; # get the errors for the consensus
#
#    varNumInt="$( \
#        printf \
#            "%s" \
#            "$strVar" |
#          sed '
#            s/^.*\///; # remove directory info
#            s/Cluster//; # remove Cluster leavning #_Consensus.fa
#            s/_Con.*$//; # remove Ending leaving ref + cluster
#          '\
#    )"; # get the variant number
#
#    # print out data for this consensus
#    printf \
#        "Nano-QConsensus-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
#        "$closeHamDistInt" \
#        "$seedInt" \
#        "$pLineStr" \
#        "$(cat "Nano-Q--time.tsv")" \
#        "$varNumInt" \
#        "$refStatStr" \
#      >> "$statFileStr.tsv";
#done # for all variants Nano-Q detect, get number of errors
#
#rm \
#    -r \
#    "$prefixStr" \
#    "Nano-Q--time.tsv"; # no longer need the time file
#
##*******************************************************************************
## Sec-6 Sub-5: run Nano-Q for median hamming distance
##*******************************************************************************

# call variants
/usr/bin/time \
    -f "%e\t%M" \
    -o "Nano-Q--time.tsv" \
  bash "$scriptDirStr/runNanoQ.sh" \
    -ref "Nano-Q--consensus.fasta" \
    -fastq "$simFastqFileStr" \
    -min-q "$minMapqDbl" \
    -p "$prefixStr" \
    -j "$jumpInt" \
    -ht "$medHamDistInt" \
    -t "$threadsInt" \
    -mt "$minClustDepthInt" \
    -racon \
    -l "$trimReadsToInt";

rm \
   -r \
   "Nano-Q--read-depth.txt" \
   "Nano-Q--medaka" \
   "Nano-Q--racon-round-1.fasta" \
   "Nano-Q--racon-round-1.fasta.map-ont.mmi";

#*******************************************************************************
# Sec-6 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
# for all variants Nano-Q detect, get number of errors
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
            s/^.*\///; # remove directory info
            s/Cluster//; # remove Cluster leavning #_Consensus.fa
            s/_Con.*$//; # remove Ending leaving ref + cluster
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "Nano-QConsensus-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$medHamDistInt" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "Nano-Q--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";
done # for all variants Nano-Q detect, get number of errors

rm \
    -r \
    "$prefixStr" \
    "Nano-Q--time.tsv"; # no longer need the time file


rm "Nano-Q--consensus.fasta";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-8: Detect variants with reference database
#    sec-8 sub-1: detect variants using Nano-Q & reference database
#    sec-8 sub-2: get percent similarity to consensus & print out data
#    sec-8 sub-3: detect variants using Nano-Q & database using median hamming dist
#    sec-8 sub-4: get percent similarity to consensus & print out data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

##*******************************************************************************
## Sec-8 Sub-1: detect variants using Nano-Q & database
##*******************************************************************************
#
#/usr/bin/time \
#    -f "%e\t%M" \
#    -o "Nano-Q--time.tsv" \
#  bash "$scriptDirStr/runNanoQ.sh" \
#    -ref "$refDbStr" \
#    -fastq "$simFastqFileStr" \
#    -p "$prefixStr" \
#    -j "$jumpInt" \
#    -ht "$closeHamDistInt" \
#    -t "$threadsInt" \
#    -mt "$minClustDepthInt" \
#    -racon \
#    -l "$trimReadsToInt";
#
##*******************************************************************************
## Sec-8 Sub-2: get percent similarity to consensus & print out data
##*******************************************************************************
#
#for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
## for all variants Nano-Q detect, get number of errors
#    if [[ ! -f "$strVar" ]]; then
#        continue;
#    fi # if null case
#
#    # get the errors in the reference
#    refStatStr="$(
#      bash "$scriptDirStr/getBestAlignment.sh" \
#          "$strVar" \
#          "$simRefStr" \
#    )"; # get the errors for the consensus
#
#    varNumInt="$( \
#        printf \
#            "%s" \
#            "$strVar" |
#          sed '
#            s/^.*\///; # remove directory info
#            s/Cluster//; # remove Cluster leavning #_Consensus.fa
#            s/_Con.*$//; # remove Ending leaving ref + cluster
#          '\
#    )"; # get the variant number
#
#    # print out data for this consensus
#    printf \
#        "Nano-QDb-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
#        "$closeHamDistInt" \
#        "$seedInt" \
#        "$pLineStr" \
#        "$(cat "Nano-Q--time.tsv")" \
#        "$varNumInt" \
#        "$refStatStr" \
#      >> "$statFileStr.tsv";
#
#done # for all variants Nano-Q detect, get number of errors
#
#rm \
#    -r \
#    "$prefixStr" \
#    "Nano-Q--time.tsv"; # no longer need the time file
#
#*******************************************************************************
# Sec-8 Sub-3: detect variants using Nano-Q & database using median hamming dist
#*******************************************************************************

/usr/bin/time \
    -f "%e\t%M" \
    -o "Nano-Q--time.tsv" \
  bash "$scriptDirStr/runNanoQ.sh" \
    -ref "$refDbStr" \
    -fastq "$simFastqFileStr" \
    -p "$prefixStr" \
    -j "$jumpInt" \
    -ht "$medHamDistInt" \
    -t "$threadsInt" \
    -mt "$minClustDepthInt" \
    -racon \
    -l "$trimReadsToInt";

#*******************************************************************************
# Sec-8 Sub-4: get percent similarity to consensus & print out data
#*******************************************************************************

for strVar in ./"$prefixStr/"*"_Consensus.fa"; do
# for all variants Nano-Q detect, get number of errors
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
            s/^.*\///; # remove directory info
            s/Cluster//; # remove Cluster leavning #_Consensus.fa
            s/_Con.*$//; # remove Ending leaving ref + cluster
          '\
    )"; # get the variant number

    # print out data for this consensus
    printf \
        "Nano-QDb-ham-%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s\t%s\t%s\n" \
        "$medHamDistInt" \
        "$seedInt" \
        "$pLineStr" \
        "$(cat "Nano-Q--time.tsv")" \
        "$varNumInt" \
        "$refStatStr" \
      >> "$statFileStr.tsv";

done # for all variants Nano-Q detect, get number of errors

rm \
    -r \
    "$prefixStr" \
    "Nano-Q--time.tsv"; # no longer need the time file

