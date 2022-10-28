#!/usr/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Do an all versus all read mapping
#   sec-4: clean up & exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#
#Sam file table for first 11 columns (all sam files have)
#| Col | Field |  Type  |        Brief description              |
#|:---:|:-----:|:------:|:-------------------------------------:|
#|  1  | QNAME | String |       Query template NAME             |
#|  2  | FLAG  |  Int   |          bitwise FLAG                 |
#|  3  | RNAME | String |     Reference sequence NAME           |
#|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
#|  5  | MAPQ  |  Int   |          MAPping Quality              |
#|  6  | CIGAR | String |            CIGAR string               |
#|  7  | RNEXT | String | Reference name of the mate/next read  |
#|  8  | PNEXT |  Int   |   Position of the mate/next read      |
#|  9  | TLEN  |  Int   |      observed Template LENgth         |
#| 10  |  SEQ  | String |          segment SEQuence             |
#| 11  | QUAL  | String | ASCII of Phred-scaled base QUALity+33 |
#

################################################################################
# Name: allVsAll.sh
# Use: Does all versus all alignment with minimap2 & extracts mapq
# Input:
#    -f: reads to align               (Required)
#    -q: min mapq to keep an aligment (Defaut: 10)
# Output: stdout: reference\tread\tmapq
# Note: Why I used minimap2 in an non-optimal way
#    - Minimap2 does not output mapq scores for ava-ont (all versus all).
#    - For a non-ava alignment using the same fastq file, minimap2 will detect
#        that each read has a reference & only give the mapping of the self
#        aligment (this is usless)
#    - For a multiref aligment, minimap2 will only give one mapq score per read.
#        I need a mapq for all aligments (more than one read)
#    - To get around these problems I extracted and removed a single reference
#        read from the fastq & then used minimap2 with ava-ont settings. This
#        was repeated for each read in the file. This is inefficant, but works.
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# user input variables
fastqFileStr="";                              # reads to get mapq's for
readListFileStr="";          # file with read names to work on
minMapqDbl="13";                               # min mapq to keep an aligment
threadsInt=3;                                  # number threads to use
prefixStr="out";                               # prefix to name everything
maxErrorRateDbl=0.08;                       # max diference to keep alignemnt
maxMedakaReadsInt=300;                      # number reads to polish

# rachet variables (clustering small amount of reads)
minDiffDbl=0.01;                             # consensus must be 1% different
minMisDbl=0.002;                     # consensues must differ by 0.2% mismatches
numReadsForCluster=300;                      # max reads to input into medaka
maxMedakaClustRatio=2; # tells how many reads for each rachet (max medaka * ratio)
discardFailedClustPercInt=1; # number reads to discard on failed cluster [100%]

# polishing
minMaskDepthInt=30; # min depth to not mask bases in conesnsus
roundsRaconInt=1; # number of rounds to run racon
medakaModelStr="r941_min_high_g351"; # fill in later (default super accuracy)
 
# cluster variables
minClustReadsInt="100";                     # min number reads to keep a cluster
minMapReadsInt="30";                          # min number reads to keep ref
minPercEdgesDbl="0.1";                        # min edges from % of reads
minEdgesInt=0;                                # found using minPercEdgesDbl

# scoreReads parameters
baseQInt=8;                                   # min Qread-list for each base
minReadLenInt=600;                            # min read length to keep aligment
maxReadLenInt=1000;                           # max read length to keep aligment
minMeanQDbl=13;                              # min mean Qread-list to keep aligment
minMedainQDbl=13;                         # min medain Qread-list to keep alignment

# max homopolymer sizes to keep inserts from
maxAHomoLenInt=1;
maxTHomoLenInt=1;
maxGHomoLenInt=1;
maxCHomoLenInt=1;

# script variables
numReadsInt=0;                                 # number reads in fastq file
scriptDirStr="$(dirname "$0")";                # directory with this script
intClust=0;                                    # cluster on in loop
totalClustInt=0;                               # total clusters in fastq
numClustInt=0;                                 # number of clusters
fileStr="";                                    # holds a file name
tmpInt=0;

# help message
helpStr="$(basename "$0") -f reads.fastq [options ...]
        -h: print this help message & exit
        -f: file to do all versus all mapping on [Required]
        -read-list: file with list of read names to work on
            - Will be deleted
            - Only first column is used (spaces or tabs deliminate columns)
        -t: number of threads to use with minimap2 [Default: 3]
        -p: Prefix to name everything [Default: out]
          - Before runing this there sould be a prefix-read-lists.tsv
            file that has all reads to use
          - The scores.tsv file will be deleted

        Racheting:
        -min-con-dif: Percent difference consensus must be [Default: 0.01 (1%)]
        -min-con-perc-mis: Percentage of mismatches consensuses must differ by 
            [Default: 0.002 (0.2%)]
        -reads-per-rachet: Number reads to input into medaka [Default: 300]
        -medaka-clust-ratio: Used to find number of reads to cluster [2]
            - Number reads found by: -reads-per-rachet * medaka-clust-ratio
        -discard: Number reads to discared on failed cluster         [1 = 100%]
    
        Polishing
        -mask-depth: Min read depth to not mask a base [30]
        -rounds-racon: Number of rounds to run racon [1]
        -medaka-model: Model to use with medaka [r941_min_high_g351]

        Clustering Input:
        -min-reads: min number reads to keep a cluster
            - Default: 100
        -min-map-reads: min number of reads a single read mus     map     o
            - Default: 30
        -min-perc-edges: How many reads needed to assing a read to a cluster
            - Is a percentage of number reads in the file
            - Default: 0.1 (10%)

        Scoring Input:
        -max-dif: max difference between query & reference to keep aigment
            - Default: 0.08%
        -q-base: min Qread-list to keep a mismatch or insertion
            - Default: 8
        -min-length: min read length to keep aligment
            - Default: 600
        -max-length: max read length to keep aligment
            - Default: 1000
        -mapq: min mapping quality to keep aligment
            - Default: 13
        -mean-q: min mean Qread-list to keep alignemnt
            - Default: 20
        -med-q: min median Qread-list to keep aligment
            - Default: 20
        -ins-A: max size of an A homopolymer size to keep an insertion
            - Default: 1
        -ins-T: max size of an A homopolymer size to keep an insertion
            - Default: 1
        -ins-G: max size of an A homopolymer size to keep an insertion
            - Default: 1
        -ins-C: max size of an A homopolymer size to keep an insertion
            - Default: 1
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get and check user input
#    sec-2 sub-1: Get user input
#    sec-2 sub-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Get user input
#*******************************************************************************

while [ $# -gt 0 ]; do                     # whild there is user input to look at
# whild have cmd arguments

    if [[ "$2" == "" ]]; then             # parameter missing argument
        printf "%s\n%s has no argument\n" "$helpStr" "$1";
        exit;
     fi;
        
    case $1 in                            # $1 is the parameter
        -f) fastqFileStr="$2";;          # file for all versus all mapping
        -read-list) readListFileStr="$2";;     # file with read names to work on
        -t) threadsInt="$2";;             # number of threads to use
        -p) prefixStr="$2";;              # prefix for output file names

        # racheting
        -min-con-dif) minDiffDbl="$2";;   # consensus must be x% different
        -min-con-perc-mis) minMisDbl="$2";; # consensuses must differ by x% mis
        -reads-per-rachet) numReadsForCluster="$2";; # number reads per batch
        -medaka-clust-ratio) maxMedakaClustRatio="$2";; # readsPerRachet * clustRatio = clusterSize
        -discard) discardFailedClustPercInt="$2";; # number reads to discard

        # clustering input
        -min-reads) minClustReadsInt="$2";; # min number reads to keep a cluster
        -min-map-reads) minMapReadsInt="$2";; # min number reads to keep ref
        -min-perc-edges) minPercEdgesDbl="$2";; # min edges from % of reads

        # polishing
        -mask-depth) minMaskDepthInt="$2";;
        -rounds-racon) roundsRaconInt="$2";;
        -medaka-model) medakaModelStr="$2";;

        # Scoring input
        -max-dif) maxErrorRateDbl="$2";;  # max difference between query & ref
        -q-base) baseQInt="$2";;          # min Qread-list to keep a bases
        -min-length) minReadLenInt="$2";; # min read length to keep aligment
        -max-length) maxReadLenInt="$2";; # max read length to keep aligment
        -mapq) minMapqDbl="$2";;     # min mapping quality to keep aligment
        -mean-q) minMeanQDbl="$2";;       # min mean Qread-list to keep alignemnt
        -med-q) minMedainQDbl="$2";;      # min median Qread-list to keep aligment
        -ins-A) maxAHomoLenInt="$2";;    # max A homopolymer size to keep insert
        -ins-T) maxTHomoLenInt="$2";;    # max T homopolymer size to keep insert
        -ins-G) maxGHomoLenInt="$2";;    # max G homopolymer size to keep insert
        -ins-C) maxCHomoLenInt="$2";;    # max C homopolymer size to keep insert
        -h) printf "%s\n" "$helpStr"; exit;; # print help message
        ?) printf "%s\n-%s is not valid\n" "$helpStr" "$1"; exit;;
    esac

    shift;                              # move to argument
    shift;                              # move to next parameter
done # whild have cmd arguments

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastqFileStr" ]]; then
    printf "%s is not a file, please povide file with -f\n" "$fastqFileStr";
    exit;
fi # Check if input file was valid


if [[ ! -f "$readListFileStr" ]]; then
    { # print out error
        printf "%s is no a file\n" "$readListFileStr";
        printf "Provided a file with read names using -read-list\n";
    } # print out error

    exit;
fi # if the read name file was invalid

if [[ $minMapqDbl -lt 0 ]]; then
    printf "Mapping quality must be > -1\n";
    exit;
fi # check if the mapping qaulity is a valid

#*******************************************************************************
# Sec-2 Sub-3: Set variables for script
#*******************************************************************************

# Find number of reads to use in the clustering sept
numReadsForCluster=$((maxMedakaClustRatio * maxMedakaReadsInt));

# find number of reads in file & use to find number of egdes
numReadsInt="$(wc -l "$fastqFileStr" | awk '{print $1}')";
numReadsInt="$((numReadsInt / 4))"; #find num of reads
  # /4 to account for 4 lines for each read, subtraction so clusting min # reads

# find the number of edges needed to assign a read to a cluster
if [[ "$numReadsInt" -le "$numReadsForCluster" ]]; then
    minEdgesInt="$(printf "%s\t%s" "$minPercEdgesDbl" "$numReadsInt" |
                        awk '{printf "%i", $1 * $2;}')";
else
    minEdgesInt="$(printf "%s\t%s" "$minPercEdgesDbl" "$numReadsForCluster" |
                        awk '{printf "%i", $1 * $2;}')";
fi # check if should use the total read count or max cluser size for num edges

if [[ "$minEdgesInt" -lt 3 ]]; then
    { # print out the warning message*/
        printf "-min-perc-edges (%f) resulted" "$minPercEdgesDbl";
        printf " in only %i edges being kept\n" "$numReadsInt";
        printf "This will result in all reads being put in one cluster\n";
        printf "Please try again with a better percentage\n";
    } # print out my error message
    exit
fi # if their are to few edges per cluster

printf "copy file with read names\n"
cat "$readListFileStr" > "$prefixStr--read-lists.tsv";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Cluster reads
#    sec-3 sub-1: Subsample reads
#    sec-3 sub-2: Check if have enough reads to build a cluster & get reads
#    sec-3 sub-3: Cluster reads
#    Sec-3 Sub-4: Find how many clusters I have to build consensus for
#    sec-3 sub-5: Check if cluster has enough reads to build consensus
#    sec-3 sub-6: Extract longest read & top other top reads
#    sec-3 sub-7: Determine if should subsample reads
#    sec-3 sub-8: Polish reads and remove uneeded files
#    sec-3 sub-9: Compare consensus & decide if keep consensus
#    sec-3 sub-10: Assign clustered reads to their bin
#    sec-3 sub-11: Clean up subsamples & prepare for next subsample
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: subsample reads
#*******************************************************************************

while [[ "$(wc -l "$prefixStr--read-lists.tsv" | awk '{print $1}')" -gt 0 ]]; do
# while there are still reads to look at

  # split files up
  mkdir "$prefixStr--split";
  cd "$prefixStr--split";
  split \
      -l "$numReadsForCluster" \
      -a 7 \
      < "../$prefixStr--read-lists.tsv";

  if [[  "$(wc -l ./* | wc -l)" -lt 2 ]]; then
    finalClustBool=1;
  fi # if on the final cluster

  cd ..; # move back to my working directory

  for strReads in ./"$prefixStr--split/"*; do
  # For all subsamples, cluster reads

    #***************************************************************************
    # Sec-3 Sub-2: Check if have enough reads to build a cluster & get reads
    #***************************************************************************

    tmpInt="$(
      wc \
        -l \
        "$strReads" |
      awk '{print $1}'\
    )" # get number of reads in file

    if [[ "$tmpInt" -lt "$minClustReadsInt" ]]; then
      continue;
    fi # skip loop if there are not enough reads to build a consensuses

    # Grab the reads to cluster from the fastq file
    "$scriptDirStr/fastqGrep" \
            -f "$strReads" \
            -fastq "$fastqFileStr" \
        > "$prefixStr--rachet-reads.fastq";

    #***************************************************************************
    # Sec-3 Sub-3: Cluster reads
    #***************************************************************************

    printf "Scoring and clustering reads\n";
    bash "$scriptDirStr/prepForCluster.sh" \
            -f "$prefixStr--rachet-reads.fastq" \
            -t "$threadsInt" \
            -max-dif "$maxErrorRateDbl" \
            -q-base "$baseQInt" \
            -min-length "$minReadLenInt" \
            -max-length "$maxReadLenInt" \
            -mapq "$minMapqDbl" \
            -mean-q "$minMeanQDbl" \
            -med-q "$minMedainQDbl" \
            -align-mean-q "$minMeanQDbl" \
            -align-med-q "$minMedainQDbl" \
            -ins-A "$maxAHomoLenInt" \
            -ins-T "$maxTHomoLenInt" \
            -ins-G "$maxGHomoLenInt" \
            -ins-C "$maxCHomoLenInt" |
        awk \
            -f "$scriptDirStr/refCount.awk" \
            -v fileStr="$prefixStr--ref-counts.tsv" |
        "$scriptDirStr/cluster" \
            -stdin \
            -min-mapq "$minMapqDbl" \
            -min-reads-per-cluster "$minClustReadsInt" \
            -min-num-mapped-reads "$minMapReadsInt" \
            -min-shared-edges "$minEdgesInt" \
        > "$prefixStr--clusters.tsv";
        # prepForCluster.sh:
        #    maps & removes alignments that I do not want to keep
        # awk -f refCount.awk:
        #     - Counts & puts number of alignments per ref in file
        #     - Also passes the score onto cluster
        # Cluster:
        #    does the clustering

    #***************************************************************************
    # Sec-3 Sub-4: Find how many clusters I have to build consensus for
    #***************************************************************************

    numClustInt="$(
      sort \
        -k 2 \
        -n \
        "$prefixStr--clusters.tsv" |
      uniq \
        -c \
        -f 1 |
      awk \
        -v minReadsInt="$minClustReadsInt" \
        '{if($1 > minReadsInt){print $2;};}' |
      wc -l |
      sed 's/[ \t]//g' \
    )"; # Find the min number of clusters with the min number of reads

    if [[ "$numClustInt" -lt 1 ]]; then
    # Discard a portion of the reads (since no reads formed a cluster)
      bash "$scriptDirStr/rmClustFromList.sh" \
          "$prefixStr--read-lists.tsv" \
          "$prefixStr--ref-counts.tsv" \
          "$strReads" \
          "$discardFailedClustPercInt";
      continue; # move to next loop
    fi # Discard a portion of the reads (since no reads formed a cluster)

    #***************************************************************************
    # Sec-3 Sub-5: Check if cluster has enough reads to build consensus
    #***************************************************************************

    for intClust in $( \
      awk \
        '{print $2}' \
        < "$prefixStr--clusters.tsv" |
      sort |
      uniq \
    ); do # for all clusters in the cluster file
        # get reads in the cluster
        awk -v clustNumInt="$intClust" \
                '{if($2 == clustNumInt){print $1;};}' \
                < "$prefixStr--clusters.tsv" \
            > "$prefixStr--filt.grep";

        # check if have enough reads to continue
        tmpInt="$( \
          wc \
            -l \
            "$prefixStr--filt.grep" |
          awk '{print $1; exit;}' \
        )"; # get numver of reads in cluster

        if [[ "$tmpInt" -lt "$minClustReadsInt" ]]; then
          continue;
        fi # if have to few reads to make the cluster

        #***********************************************************************
        # Sec-3 Sub-6: Extract longest read & top other top reads
        #***********************************************************************

        # Get the reads for the cluster 
        "$scriptDirStr/fastqGrep" \
                -f "$prefixStr--filt.grep" \
                -fastq "$prefixStr--rachet-reads.fastq" \
            > "$prefixStr--tmp-reads--cluster-$totalClustInt.fastq";

        # Find & extract longest read
        # outputs: prefixStr--longest-read.fastq & prefixStr--other-reads.fastq
        awk -v filePrefixStr="$prefixStr--cluster-$totalClustInt" \
                -f "$scriptDirStr/extractLongestRead.awk" \
                < "$prefixStr--tmp-reads--cluster-$totalClustInt.fastq";
        rm "$prefixStr--tmp-reads--cluster-$totalClustInt.fastq"; 

        #***********************************************************************
        # Sec-3 Sub-7: Determine if should subsample reads
        #***********************************************************************

        # Get number of reads assigend to the cluster
        tmpInt="$(wc -l "$prefixStr--filt.grep" |
                        awk '{print $1}')";

        if [[ "$minClustReadsInt" -lt "$tmpInt" ]]; then
        # If need to down sample
            head -n+"$((maxMedakaReadsInt * 4))" \
                < "$prefixStr--cluster-$totalClustInt--other-reads.fastq"\
              > "$prefixStr--tmp.fastq";
            mv "$prefixStr--tmp.fastq" \
                "$prefixStr--cluster-$totalClustInt--other-reads.fastq";
        fi # If need to down sample

        #***********************************************************************
        # Sec-3 Sub-8: Polish reads and remove uneeded files
        #***********************************************************************

        bash "$scriptDirStr/buildConsensus.sh" \
                -i "$prefixStr--cluster-$totalClustInt--other-reads.fastq" \
                -I "$prefixStr--cluster-$totalClustInt--longest-read.fasta" \
                -t "$threadsInt" \
                -p "$prefixStr--cluster-$totalClustInt" \
                -c "$minMaskDepthInt" \
                -r "$roundsRaconInt" \
                -m "$medakaModelStr";

        # Remove the uneeded files
        rm -r "$prefixStr--cluster-$totalClustInt--medaka" \
                "$prefixStr--cluster-$totalClustInt--racon"* \
                "$prefixStr--cluster-$totalClustInt--other-reads.fastq" \
                "$prefixStr--cluster-$totalClustInt--longest-read.fasta" \
                "$prefixStr--cluster-$totalClustInt--read-depth.txt";

        #***********************************************************************
        # Sec-3 Sub-9: compare consensuses & decide if keep consensus
        #***********************************************************************

        fileStr="";

        if [[ -f "$prefixStr--consensuses.fasta" ]]; then
           fileStr="$( \
             bash "$scriptDirStr/compareConsensuses.sh" \
               -ref "$prefixStr--cluster-$totalClustInt--consensus.fasta" \
               -f "$prefixStr--consensuses.fasta" \
               -t "$threadsInt" \
               -min-diff "$minDiffDbl" \
               -min-mismatches "$minMisDbl").fastq"; # bins are fastq files
        fi # Check if I need to compare consenuses (past consensuses exist)

        if [[ -f "$fileStr" ]]; then
            rm "$prefixStr--cluster-$totalClustInt--consensus.fasta";
        else
            fileStr="$prefixStr--cluster-$totalClustInt.fastq";
            # Add new consensus to list of references
            cat "$prefixStr--cluster-$totalClustInt--consensus.fasta" \
                >> "$prefixStr--consensuses.fasta";
            # remove the old consensus file, since only taking up space
            rm "$prefixStr--cluster-$totalClustInt--consensus.fasta";
        fi # check if need to remove the new consensus

        #***********************************************************************
        # Sec-3 Sub-10: Assign clustered reads to their bin
        #***********************************************************************

        # Add mapped reads to the consensuses bin
        "$scriptDirStr/fastqGrep" \
                -f "$prefixStr--filt.grep" \
                -fastq "$prefixStr--rachet-reads.fastq" \
            >> "$fileStr";

       # remove the moved entries from my score file
       grep -v -f "$prefixStr--filt.grep" \
                < "$prefixStr--read-lists.tsv" \
            > "$prefixStr--tmp.tsv";
       mv "$prefixStr--tmp.tsv" "$prefixStr--read-lists.tsv";

       totalClustInt="$((totalClustInt + 1))";
       intClust="$((intClust + 1))";               # move to next cluster
    done # While I have clusters to build consensus for
  done # For all subsamples, cluster reads

  #***********************************************************************
  # Sec-3 Sub-11: Clean up subsamples & prepare for next subsample
  #***********************************************************************

  rm -r "$prefixStr--split";   # Remove old splits

  if [[ "$finalClustBool" -gt 0 ]]; then
    break;  # prevents an infinite loop, do to only on read in file
  fi # if done clustering reads (only one file left, & just split in half inf)
done # while there are still reads to look at


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: clean up & exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm \
    "$prefixStr--filt.grep" \
    "$prefixStr--trimmed-reads.fastq" \
    "$prefixStr--read-lists.tsv" \
    "$prefixStr--cluster-"*"-other-reads.fastq" \
    "$prefixStr--cluster-"*"-longest-reads.fastq" \
    "$prefixStr--clusters.tsv" \
    "$prefixStr--cluster-"*"--consensus.fasta" \
    "$prefixStr--read-lists.tsv" \
    "$prefixStr--rachet-reads.fastq" \
    "$prefixStr--tmp-reads--cluster"*".fastq";
exit;
