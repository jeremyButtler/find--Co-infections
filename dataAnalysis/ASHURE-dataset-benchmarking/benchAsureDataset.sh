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
primPathStr="NA"; # fasta file with primers
prefixStr="out";
altConCmd=""; # store consensus buildng commands
modelStr="r941_min_high_g303";
   # ASURE data set on 9.5, guppy version 3.2.2, wich is less than
   #   guppy version 3.0.3. Rules for model selction, use a model that
   #   is equal to or less then the guppy version (g303 means guppy
   #    version 3.0.3)
   # I am assuming they use a high accuracy model.
   # Their model used with guppy was dna_r9.5_450bps_1d2_raw
   # for fast is "r941_min_fast_g303"
   # for high accurace is "r941_min_high_g303"
numPolishRndsI=2; # number of rounds to rebuild/poilsh the consensus
keepFilesBl="FALSE";    # 1: keep all files, 0 discard

# hardcode variables
majMinBaseQI="7";
majMinBaseSupportI="0.35";
majMinInsQ="5";
majMinInsSupportI="0.30";
maxReadLen="750";        # 0 is no limit
minPercReadsFlt="0.003"; # min percentage of reads to keep a bin

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
trimPrimPathStr="../trimPrimers";
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
	-primers) primPathStr="$2"; shift;;
        -refs) refsStr="$2"; shift;;
        -max-read-len) maxReadLen="$2"; shift;;
        -min-perc-reads) minPercReadsFlt="$2"; shift;;
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

if [[ "$primPathStr" != "NA" ]]; then
    altConCmd="$altConCmd -primers $primPathStr";
fi # if need to add in the primers command

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
        printf "\tcpuPercentage\tprimers\tmaxReadLength";
	printf "\tminPercReadsToKeepClust\n";
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

if [[ "$primPathStr" == "NA" ]]; then
# if trimming by references
    # trim the fastq file
    minimap2 \
        --eqx \
        -a \
        --secondary=no \
        "$refsStr" \
        "$fqPathStr" |
      "$trimSamPathStr" -stdin -keep-unmapped-reads |
      awk '{
        if($1 ~ /^@/ || and($2, 2048) || and($2, 256)){next;};
        printf "@%s\n%s\n+\n%s\n", $1, $10, $11
    }' > "$prefixStr--$suffixStr.fastq";
    
    fqPathStr="$prefixStr--$suffixStr.fastq";
    primPathStr="byActualRef";
fi # if trimming by references

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
        "$fqPathStr" |
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
        -min-perc-reads "$minPercReadsFlt" \
	-max-read-length "$maxReadLen" \
        -fastq "$fqPathStr";
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
        -min-perc-reads "$minPercReadsFlt" \
	-max-read-length "$maxReadLen" \
        -fastq "$fqPathStr" \
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
        scoresStr="NA	NA	NA	NA	NA	NA	NA"; # blank score line
    fi # see if mapped to a reference

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
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
        "$primPathStr" \
        "$maxReadLen" \
        "$minPercReadsFlt" \
      >> "$statsFileStr";
done

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-8: Clean up and exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ "$keepFilesBl" == "FALSE" ]]; then
    rm "$prefixStr--$suffixStr"*; # remove all files, except stats file
fi

exit;
