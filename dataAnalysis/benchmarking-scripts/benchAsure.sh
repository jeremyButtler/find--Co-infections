#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SOP: Start of program
#   sec-1: Variable declerations
#   sec-2: Set up the stats file
#   sec-3: Start loop and trim reads
#   sec-4: Find the number of references with reads in fastq file
#   sec-5: Detect co-infections/alternate genes
#   sec-6: Find the number of unique references that were detected
#   sec-7: Get consensus stats and print out entries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

fqPathStr="";
refsStr="";
prefixStr="out";
altConCmd=""; # store consensus buildng commands
modelStr="r941_min_high_g351";
numPolishRndsI=2; # number of rounds to rebuild/poilsh the consensus
keepFilesBl="FALSE";    # 1: keep all files, 0 discard

# hardcode variables
majMinBaseQI="7";
majMinBaseSupportI="0.35";
majMinInsQ="5";
majMinInsSupportI="0.30";

# script variables
statsFileStr="";
suffixStr="";
numRefsInFqStr="";
numRefsDetected="";
scoresStr="";
timeStatStr="";   # stats from /usr/bin/time
numClustReadsI=0; # number of reads assigend to a cluster
clustFqStr="";    # fastq assigned to a cluster
mapRefStr="";
refLenI="";

# flags
majConBl="TRUE";
raconBl="FALSE";
medakaBl="FALSE";

# paths to programs
trimSamPathStr="../trimSamFile";
findCoInftPathStr="../findCoInft";
scoreReadsPathStr="../scoreReads";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Set up the stats file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

while [ $# -gt 0 ]; do
# while have user input to read in

    case $1 in
        -fastq) fqPathStr="$2"; shift;;
        -prefix) prefixStr="$2"; shift;;
        -refs) refsStr="$2"; shift;;
        -model) modelStr="$2"; shift;;
        -enable-racon) useRaconCmd="-enable-racon"; raconBl="TRUE";;
        -enable-medaka) useMedakaCmd="-enable-medaka"; medakaBl="TRUE";;
        -disable-maj) noMajConCmd="-disable-majority-consensus"; majConBl="FALSE";;
        -rnds-polish) numPolishRndsI="$2"; shift;;
        -keep-files) keepFilesBl="TRUE";;
    esac

    shift;
done # while have user input to read in

if [[ "$majConBl" == "FALSE" ]]; then
   altConCmd="-disable-majority-consensus";
fi

if [[ "$raconBl" == "TRUE" ]]; then
   altConCmd="$altConCmd -enable-racon";
fi

if [[ "$medakaBl" == "TRUE" ]]; then
   altConCmd="$altConCmd -enable-medaka";
fi

if [[ "$prefixStr" == "" ]]; then
    prefixStr = "out";
fi

statsFileStr="$prefixStr--stats.tsv";

if [[ ! -f "$prefixStr--stats.tsv" ]]; then
# If I need to make the stats file
    { # print out header to stats file
        printf "Prefix\tfqFile\tReferenceFile\tnumTimesRebuildCon";
        printf "\tMajorityConsensus\tRaconUsed\tMedakaUsed\tMedakaModel";
        printf "\treferenceLength\tConsensus\tReference\treadLength";
        printf "\tAlignedReadLength\tSNPs\tInsertions\tDeletions";
        printf "\tReadsInCluster\tDetectedReferences";
        printf "\tDectableReferences\tTotalReferences\telapsedTime";
        printf "\tuserTime\tkernalTime\tmaxResidentMemory";
        printf "\tcpuPercentage\n";
    } > "$statsFileStr" # print out the header to the stats file
fi # If I need to make the stats file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: trim reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ ! -f "$fqPathStr" ]]; then
    printf "%s is not a file" "$fqPathStr";
    exit;
fi # if the null case

# Get the fastq files name
suffixStr="$(printf "%s" "$fqPathStr" | sed 's/.*\///; s/\.fastq//;')";

# trim the fastq file
minimap2 \
    --eqx \
    -a \
    --secondary=no \
    "$refsStr" \
    "$fqPathStr" |
  "$trimSamPathStr" -stdin |
  awk '{
    if($1 ~ /^@/ || and($2, 2048) || and($2, 256)){next;};
    printf "@%s\n%s\n+\n%s\n", $1, $10, $11
}' > "$prefixStr--$suffixStr.fastq";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Find the number of references with reads in fastq file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# find number references in fastq & references with over 100x reads
numRefsInFqStr="$( \
    minimap2 \
        --eqx \
        -a \
        --secondary=no \
        "$refsStr" \
        "$prefixStr--$suffixStr.fastq" |
     awk '{
         # Remove headers & supplemental, 2ndary, & unmapped alignments
         if($1 ~ /^@/ || and($2, 2048) || and($2, 256) || and($2, 4))
             {next;};
         print $3;
       }' |
     sort |
     uniq -c |
     awk '{ # main
         if($1 > 99){++detectRefI};# if reference has at least 100 reads
         ++totalRefI;     # count total number of references in the file
       } # main
       END{printf "%s\t%s", detectRefI, totalRefI};'
)"; # Find the number of references represented by the fastq file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Detect co-infections/alternate genes for majority consensus
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$altConCmd" == "" ]]; then
    /usr/bin/time \
        -f "%e\t%U\t%S\t%M\t%P" \
        -o "$prefixStr--$suffixStr--time.tsv" \
        -a \
      "$findCoInftPathStr" \
        -skip-bin \
        -prefix "$prefixStr--$suffixStr" \
	-extra-consensus-steps "$numPolishRndsI" \
        -model "$modelStr" \
        -maj-con-min-base-q "$majMinBaseQI" \
        -maj-con-min-bases "$majMinBaseSupportI" \
        -maj-con-min-ins-q "$majMinInsQ" \
        -maj-con-min-ins "$majMinInsSupportI" \
        -fastq "$prefixStr--$suffixStr.fastq";
else
    /usr/bin/time \
        -f "%e\t%U\t%S\t%M\t%P" \
        -o "$prefixStr--$suffixStr--time.tsv" \
        -a \
      "$findCoInftPathStr" \
        -skip-bin \
	-extra-consensus-steps "$numPolishRndsI" \
        -prefix "$prefixStr--$suffixStr" \
        -model "$modelStr" \
        -maj-con-min-base-q "$majMinBaseQI" \
        -maj-con-min-bases "$majMinBaseSupportI" \
        -maj-con-min-ins-q "$majMinInsQ" \
        -maj-con-min-ins "$majMinInsSupportI" \
        -fastq "$prefixStr--$suffixStr.fastq" \
        $altConCmd; # removeing "" because I want it to process spaces
fi # check if using alternative consensus commands

if [[ "$medakaBl" == "FALSE" ]]; then
    modelStr="NA";
fi # if not using medaka, blank the model

# get the stats from the time command
timeStatStr="$(cat "$prefixStr--$suffixStr--time.tsv")";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Find the number of unique references that were detected
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cat "$prefixStr--$suffixStr"*con.fasta \
  > "$prefixStr--$suffixStr-tmp.fasta"; # merge consensus for mapping

numRefsDetected="$( \
   minimap2 \
       --eqx \
       -a \
       --secondary=no \
       "$refsStr" \
       "$prefixStr--$suffixStr-tmp.fasta" |
     awk '{
         # Remove headers & supplemental, 2ndary, & unmapped alignments
         if($1 ~ /^@/ || and($2, 2048) || and($2, 256) || and($2, 4))
             {next;};
         print $3;
       }' |
     sort |
     uniq -c |
     wc -l \
)"; # find the number of consensuses made

rm "$prefixStr--$suffixStr-tmp.fasta"; # no longer need

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Get consensus stats and print out entries
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

for strCon in ./"$prefixStr--$suffixStr"*con.fasta; do
# loop and score consensuses for all consensus made

    if [[ ! -f "$strCon" ]]; then
        continue;
    fi # null case

    # Get the clusters fastq file
    clustFqStr="$(printf "%s" "$strCon" | sed 's/--con.fasta/.fastq/')";

    # find number of reads assigned to the cluster
    numClustReadsI="$(sed -n 'p;n;n;n;' "$clustFqStr" | wc -l)";

    scoresStr="$(
        minimap2 \
            --eqx \
            -a \
            --secondary=no \
            "$refsStr" \
            "$strCon" |
          "$scoreReadsPathStr" \
            -stdin \
            -min-read-length 0 \
            -max-read-length 0 |
          awk '{
              getline; # move past the header
              printf "%s\t%s\t%s\t%s", $1, $2, $4, $5;
              printf "\t%s\t%s\t%s", $15, $16, $17;
              exit;    # done with the file
          }' \
    )";

    mapRefStr="$(printf "%s" "$scoresStr" | awk '{print $2}')";

    if [[ "$mapRefStr" != "" ]]; then
        refLenI="$( \
            sed \
                -n \
                "/$mapRefStr/{n;p;q;}" \
                "$refsStr" |
              awk '{print length($0)}'
        )"; # get the length of the reference
    else
        mapRefStr="NA";
        scoresStr="NA\tNA\tNA\tNA\tNA\tNA\tNA"; # blank score line
    fi # see if mapped to a reference

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$prefixStr" \
        "$fqPathStr" \
        "$refsStr" \
        "$numPolishRndsI" \
        "$majConBl" \
        "$raconBl" \
        "$medakaBl" \
        "$modelStr" \
        "$refLenI" \
        "$scoresStr" \
        "$numClustReadsI" \
        "$numRefsDetected" \
        "$numRefsInFqStr" \
        "$timeStatStr" \
      >> "$statsFileStr";
done

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-8: Clean up and exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$keepFilesBl" == "FALSE" ]]; then
    rm "$prefixStr--$suffixStr"*; # remove all files, except stats file
fi

exit;
