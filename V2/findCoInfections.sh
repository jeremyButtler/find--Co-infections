#!/usr/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Use reference to map, trim, & score reads
#   sec-4: Cluster reads
#   sec-5: Build the final consensus genomes
#   sec-6: get read counts & clean up
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: findCoInfections.sh
# Use: Finds co-infections in a fastq file
# Key Input:
#   -f: Fastq file to find co-infections in                    [Required]
#   -t: number of threads to use                               [3]
#   -p: prefix for output file names                           [out]
#   -ref: Fasta with references genomes to map reads to        [Required]
#     - References should only cover the region of interest
#   -full-ref: fasta file with full length references          [fasta from -ref]
#   -cluster:
#       yes: Do clustering step                                [Default no]
#       no: Skip clustering step
# For other input See the help message (sec-1 sub-7).
# Output: 
#   file: holding built consensuses [prefix--final--consensuses.fasta]
#   file: holding the read counts for each bin [prefixStr--counts.tsv]
#   file: fastq files with reads for each bin
#        - [prefix--reference--cluster-#.fastq]
################################################################################

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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#    sec-1 sub-1: general user input
#    sec-1 sub-2: Variables for scoring reads or references
#    sec-1 sub-3: Variables for racheting steps (clustering subsamples)
#    sec-1 sub-4: Variables for the polshing steps
#    sec-1 sub-5: Variables for checking read or alignemnt quality
#    sec-1 sub-6: Variables not set by user (script variables)
#    sec-1 sub-7: help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: general user input
#*******************************************************************************

fastqFileStr="";                              # reads to get mapq's for
threadsInt=3;                                  # number threads to use
prefixStr="out";                               # prefix to name everything
refFastaStr="";                             # reference to map & trim reads with
fullRefFastaStr="";                # holds full length reference for final polish
clustStr="no";                                # do clustering step yes, no 

#*******************************************************************************
# Sec-1 Sub-2: Variables for scoring reads or references
#*******************************************************************************

lenWeightDbl=1.25;                          # Weight of read length in scoring
mapqWeightDbl=1;                          # Weight of mapping quality in scoring
maxErrorRateDbl=0.07;                       # max diference to keep alignemnt
refFiltMapqDbl=20;                        # min mapq in reference filtering step

#*******************************************************************************
# Sec-1 Sub-3: Variables for racheting steps (clustering subsamples)
#*******************************************************************************

minDiffDbl=0.01;                             # consensus must be 1% different
minMisDbl=0.003;                     # consensues must differ by 0.3% mismatches
numReadsForCluster=300;                      # max reads to input into medaka
minRachetMapqDbl=20;                         # min read mapq for racheting
maxMedakaClustRatio=2; # tells how many reads for each rachet (max medaka * ratio)
discardFailedClustPercInt=1; # number reads to discard on failed cluster [100%]

#*******************************************************************************
# Sec-1 Sub-4: Variables for the polshing steps
#*******************************************************************************

minMaskDepthInt=30; # min depth to not mask bases in conesnsus
roundsRaconInt=1; # number of rounds to run racon
medakaModelStr="r941_min_high_g351"; # fill in later (default super accuracy)
maxMedakaReadsInt=300;                      # number reads to polish
 
#*******************************************************************************
# Sec-1 Sub-4: Variables for the polshing steps
#*******************************************************************************

minMapqDbl="13";                               # min mapq to keep an aligment
minClustReadsInt="100";                     # min number reads to keep a cluster
minMapReadsInt="30";                          # min number reads to keep ref
minPercEdgesDbl="0.1";                        # min edges from % of reads

#*******************************************************************************
# Sec-1 Sub-5: Variables for checking read or alignemnt quality
#*******************************************************************************

baseQInt=8;                                   # min Q-score for each base
maxClustErrorRateDbl=0.07;                    # max error rate for clustering reads
minReadLenInt=600;                            # min read length to keep aligment
maxReadLenInt=0;                              # max read length to keep aligment
minMeanQDbl=13;                              # min mean Q-score to keep aligment
minMedainQDbl=13;                         # min medain Q-score to keep alignment
minAlignMeanQDbl=13;                         # min aligned mean Q-score
minAlignMedainQDbl=13;                       # min aligned median Q-score
maxAHomoLenInt=1;
maxTHomoLenInt=1;
maxGHomoLenInt=1;
maxCHomoLenInt=1;

#*******************************************************************************
# Sec-1 Sub-6: Variables not set by user
#*******************************************************************************

scriptDirStr="$(dirname "$0")/coInfectScripts"; # directory with other scripts
refStr="";                                     # holds name of reference
clusterStr="";                               # holds cluster--number of cluster
tmpInt=0;

#*******************************************************************************
# Sec-1 Sub-7: help message
#*******************************************************************************

keyHelpStr="$(basename "$0") -f reads.fastq -ref refences.fasta [options ...]
  Use: Finds co-infections in a fastq file
  Key Input:
    -f: Fastq file to find co-infections in                  [Required]
    -t: number of threads to use                             [3]
    -p: prefix for output file names                         [out]
    -ref: Fasta with references genomes to map reads to      [Required]
       - References should only cover the region of interest
    -full-ref: fasta file with full length references        [fasta from -ref]
    -cluster:
        yes: Do clustering step                              [Default no]
        no: Skip clustering step

    -min-length: min read length to keep aligment                   [600]
    -max-length: max read length to keep aligment                   [1000]

    -rounds-racon: Number of rounds to run racon                    [1]
    -medaka-model: Model to use with medaka             [r941_min_high_g351]

    -h: print this help message & exit
    -h-detail: print a more detailed help message
  Output: 
    file: holding built consensuses [prefix--final--consensuses.fasta]
    file: holding the read counts for each bin [prefixStr--counts.tsv]
    file: fastq files with reads for each bin 
      - [prefix--reference--cluster-#.fastq]
";

helpStr="$(basename "$0") -f reads.fastq -ref refences.fasta [options ...]
  Use: Finds co-infections in a fastq file
  Key Input:
    -f: Fastq file to find co-infections in                  [Required]
    -t: number of threads to use                             [3]
    -p: prefix for output file names                         [out]
    -ref: Fasta with references genomes to map reads to      [Required]
       - References should only cover the region of interest
    -full-ref: fasta file with full length references        [fasta from -ref]
    -cluster:
        yes: Do clustering step                              [Default]
        no: Skip clustering step
    -h: print simpler help message & exit
    -h-detail: print this help message & exit
  Output: 
    file: holding built consensuses [prefix--final--consensuses.fasta]
    file: holding the read counts for each bin [prefixStr--counts.tsv]
    file: fastq files with reads for each bin 
      - [prefix--reference--cluster-#.fastq]

  Other parameters:
     Read scoring:
       -len-weight: Weight of aligned read length in scoring            [1.25]
         - Score is square root of read length
       -mapq-weight: Weight of mapping quality in scoring               [1]
       -ref-mapq: Min mapping quality to keep read in reference mapping [20]

     Racheting (clustering on subsampe of reads):
       -min-con-dif: % difference between consensus                  [0.01 (1%)]
       -min-con-perc-mis: Min % mismatch difference between consensuses [0.003]
       -min-rachet-mapq: min mapq for mapping read to consensuses       [30]
       -reads-per-rachet: Max number of reads to input into medaka      [300]
       -medaka-clust-ratio: Used to find subsample depth in cluster step [2]
         - Number reads to cluster = -reads-per-rachet * -medaka-clust-ratio
       -discard: Number reads to discared on failed cluster         [1 = 100%]

     Polishing:
        -mask-depth: Min read depth to not mask a base                  [30]
        -rounds-racon: Number of rounds to run racon                    [1]
        -medaka-model: Model to use with medaka             [r941_min_high_g351]

     Clustering:
        -min-reads: min number reads to keep a cluster                  [100]
        -min-map-reads: min number of reads a single read must map to   [30]
	-max-clust-read-dif: max diference beteween reads in clustering [0.07 (7%)]
        -min-perc-edges: Used to find number of mappings needed to
                         assign read to a cluster                       [0.1 (10%)]
          - Min number edges = -min-perc-edges * number of reads

     Scoring:
        -max-dif: max difference between query & other read             [0.07 (7%)]
        -q-base: min Q-score to keep a mismatch or insertion            [8]
        -min-length: min read length to keep aligment                   [600]
        -max-length: max read length to keep aligment                   [1000]
        -mapq: min mapping quality to keep aligment                     [13]
        -mean-q: min mean Q-score to keep alignemnt                     [20]
        -med-q: min median Q-score to keep aligment                     [20]
        -ins-A: max size of an A homopolymer size to keep an insertion  [1]
        -ins-T: max size of an T homopolymer size to keep an insertion  [1]
        -ins-G: max size of an G homopolymer size to keep an insertion  [1]
        -ins-C: max size of an C homopolymer size to keep an insertion  [1]

  Output: 
    file: holding built consensuses [prefix--final--consensuses.fasta]
    file: holding the read counts for each bin [prefixStr--counts.tsv]
    file: fastq files with reads for each bin 
      - [prefix--reference--cluster-#.fastq]

";

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
        -t) threadsInt="$2";;             # number of threads to use
        -p) prefixStr="$2";;              # prefix for output file names
        -ref) refFastaStr="$2";;          # if is a reference genome
        -full-ref) fullRefFastaStr="$2";; # file with full length references
        -cluster) clustStr="$2";;         # users tells if to cluster or not
        -h-detail) print "%sn" "$helpStr"; exit;; # print more detalid help message
        -h) printf "%s\n" "$keyHelpStr"; exit;; # print help message

        # scoring
        -len-weight) lenWeightDbl="$2";;  # Weight of read length in scoring
        -mapq-weight) mapqWeightDbl="$2";;# Weight of mapping quality in scoring
        -ref-mapq) refFiltMapqDbl="$2";;  # min mapq in reference filtering step

        # racheting
        -min-con-dif) minDiffDbl="$2";;   # consensus must be x% different
        -min-con-perc-mis) minMisDbl="$2";; # consensuses must differ by x% mis
        -min-rachet-mapq) minRachetMapqDbl="$2";; # min read mapq for racheting
        -reads-per-rachet) numReadsForCluster="$2";; # number reads per batch
        -medaka-clust-ratio) maxMedakaClustRatio="$2";; # readsPerRachet * clustRatio = clusterSize
        -discard) discardFailedClustPercInt="$2";; # number reads to discard

        # polishing
        -mask-depth) minMaskDepthInt="$2";;
        -rounds-racon) roundsRaconInt="$2";;
        -medaka-model) medakaModelStr="$2";;

        # clustering input
        -min-reads) minClustReadsInt="$2";; # min number reads to keep a cluster
        -min-map-reads) minMapReadsInt="$2";; # min number reads to keep ref
        -min-perc-edges) minPercEdgesDbl="$2";; # min edges from % of reads
        -max-clust-read-dif) maxClustErrorRateDbl="$2";; # max diference beteween reads in clustering

        # Scoring input
        -max-dif) maxErrorRateDbl="$2";;  # max difference between query & ref
        -q-base) baseQInt="$2";;          # min Q-score to keep a bases
        -min-length) minReadLenInt="$2";; # min read length to keep aligment
        -max-length) maxReadLenInt="$2";; # max read length to keep aligment
        -mapq) minMapqDbl="$2";;     # min mapping quality to keep aligment
        -mean-q) minMeanQDbl="$2";;       # min mean Q-score to keep alignemnt
        -med-q) minMedainQDbl="$2";;      # min median Q-score to keep aligment
        -ins-A) maxAHomoLenInt="$2";;    # max A homopolymer size to keep insert
        -ins-T) maxTHomoLenInt="$2";;    # max T homopolymer size to keep insert
        -ins-G) maxGHomoLenInt="$2";;    # max G homopolymer size to keep insert
        -ins-C) maxCHomoLenInt="$2";;    # max C homopolymer size to keep insert
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

if [[ ! -f "$refFastaStr" ]]; then
    printf
        "%s is not a file, please provide fasta file with reference (-ref)\n" \
        "$refFastaStr";
    exit;
fi

if [[ ! -f "$fullRefFastaStr" ]]; then
    fullRefFastaStr="$refFastaStr";
fi # if no full length references provided

if [[ $minMapqDbl -lt 0 ]]; then
    printf "Mapping quality must be > -1\n";
    exit;
fi # check if the mapping qaulity is a valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Use reference to map, trim, & score reads
#   sec-3 sub-1: Filter out low quality reads
#   sec-3 sub-2: trim & bin reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Filter out low quality reads
#*******************************************************************************

# output: prefixStr--scores.tsv & prefixStr--filt-reads.fastq
bash "$scriptDirStr/filterReads.sh" \
        -f "$fastqFileStr" \
        -ref "$refFastaStr" \
        -t "$threadsInt" \
        -p "$prefixStr" \
        -len-weight "$lenWeightDbl" \
        -mapq-weight "$mapqWeightDbl" \
        -max-dif "$maxErrorRateDbl" \
        -q-base "$baseQInt" \
        -min-length "$minReadLenInt" \
        -max-length "$maxReadLenInt" \
        -mapq "$refFiltMapqDbl" \
        -mean-q "$minMeanQDbl" \
        -med-q "$minMedainQDbl" \
        -align-mean-q "$minAlignMeanQDbl" \
        -align-med-q "$minAlignMedainQDbl" \
        -ins-A "$maxAHomoLenInt" \
        -ins-T "$maxTHomoLenInt" \
        -ins-G "$maxGHomoLenInt" \
        -ins-C "$maxCHomoLenInt";

#*******************************************************************************
# Sec-3 Sub-2: Trim & bin reads
#*******************************************************************************

bash "$scriptDirStr/trimReads.sh" \
  -f "$prefixStr--filt-reads.fastq" \
  -t "$threadsInt" \
  -p "$prefixStr" \
  -ref "$refFastaStr";

rm "$prefixStr--filt-reads.fastq"; # no longer need

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Cluster reads
#    sec-4 sub-1: Check if should cluster, if not rename bins for next step
#    sec-4 sub-2: Start loop to cluster all fastq files
#    sec-4 sub-3: Get scores for reads in a bin
#    sec-4 sub-4: Run script to cluster reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-4 Sub-1: Check if should cluster, if not rename bins for next step
#*******************************************************************************

if [[ "$clustStr" != "yes" ]]; then
# if not clustering reads, then need to rename

  for strFastq in ./"$prefixStr--"*".fastq"; do
  # for all fastq files, insert "--cluster-NA" at end

    if [[ ! -f "$strFastq" ]]; then
      continue; # null case
    fi # if null case move on (loops exist)

    tmpInt="$(wc -l "$strFastq" | awk '{print $1}')";

    if [[ "$tmpInt" -lt "$minClustReadsInt" ]]; then
      rm "$strFastq";
      continue;
    fi # if to few reads to keep bin

    mv \
      "$strFastq" \
      "$(
        printf "%s" "$strFastq" |
        sed '
          s/\.\///;
          s/\./--cluster-NA./;
        '
      )";
  done # for all fastq files, insert "--cluster-0" at end
# if not clustering reads, then need to rename

#******************************************************************************
# Sec-4 Sub-2: Start loop to cluster all fastq files
#******************************************************************************

else
# Else I am clustering reads
  
  for strFastq in ./"$prefixStr--"*".fastq"; do
  # For all bins in the working directory
  
      if [[ ! -f "$strFastq" ]]; then
          continue;
      fi # no file, likey null case 
  
      # Get name of the reference that made the bin
      refStr="$(printf "%s" "$strFastq" | sed 's/.*--//; s/\.fastq//')";
  
      #************************************************************************
      # Sec-4 Sub-3: Get scores for reads in a bin
      #************************************************************************
  
      # Get the reads in this bin
      sed -n 's/^@//p;n;n;n;' \
              "$strFastq" \
          > "$prefixStr--tmp.grep";
  
      # grab out the scores of the reads in this bin
      grep -f "$prefixStr--tmp.grep" \
              "$prefixStr--scores.tsv" \
          > "$prefixStr--tmp--scores.tsv";
  
      tmpInt="$(wc -l "$prefixStr--tmp--scores.tsv" | awk '{print $1}')";
  
      if [[ "$tmpInt" -lt "$minClustReadsInt" ]]; then
          rm "$strFastq";
          continue;       # move to the next bin
      fi # If have to few reads in bin
  
      #************************************************************************
      # Sec-4 Sub-4: Run script to cluster reads
      #************************************************************************
  
      bash "$scriptDirStr/rachetAndCluster.sh" \
          -f "$strFastq" \
          -read-list "$prefixStr--tmp--scores.tsv" \
          -t "$threadsInt" \
          -p "$prefixStr--$refStr" \
          -discard "$discardFailedClustPercInt" \
          -mask-depth "$minMaskDepthInt" \
          -rounds-racon "$roundsRaconInt" \
          -medaka-model "$medakaModelStr" \
          -min-con-dif "$minDiffDbl" \
          -min-con-perc-mis "$minMisDbl" \
          -min-rachet-mapq "$minRachetMapqDbl" \
          -reads-per-rachet "$numReadsForCluster" \
          -medaka-clust-ratio "$maxMedakaClustRatio" \
          -min-reads "$minClustReadsInt" \
          -min-map-reads "$minMapReadsInt" \
          -min-perc-edges "$minPercEdgesDbl" \
          -max-dif "$maxClustErrorRateDbl" \
          -q-base "$baseQInt" \
          -min-length "$minReadLenInt" \
          -max-length "$maxReadLenInt" \
          -mapq "$minMapqDbl" \
          -mean-q "$minMeanQDbl" \
          -med-q "$minMedainQDbl" \
          -ins-A "$maxAHomoLenInt" \
          -ins-T "$maxTHomoLenInt" \
          -ins-G "$maxGHomoLenInt" \
          -ins-C "$maxCHomoLenInt";
  
    rm "$strFastq"; # remove none clustered reads file
  done # For all bins in the working directory
# else clustering reads
fi # check if clusteing reads or not

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Build the final consensus genomes
#    sec-5 sub-1: Check if null case & get names
#    sec-5 sub-2: Score reads
#    sec-5 sub-3: Pull out top reads & extract longest read
#    sec-5 sub-4: Polish my reads
#    sec-5 sub-5: Determine if should keep consensus (to similar to others?)
#    sec-5 sub-6: Assgin reads for consensus to a final bin
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-5 Sub-1: Check if null case, get ref name, extract untrimmed reads
#*******************************************************************************

for strFastq in ./"$prefixStr--"*"cluster-"*".fastq"; do
# For all kept clusters

    if [[ ! -f "$strFastq" ]]; then
        continue;      # loop will quit next round
    fi # if on null case (no file)

    # Get name of the reference that made the cluster
    refStr="$(printf "%s" "$strFastq" | sed 's/--cluster.*//; s/.*--//;')";

    # Get cluster-number
    clusterStr="$(printf "%s" "$strFastq" | sed 's/.*--//; s/\.fastq//')";

    # Extract untrimmed reads from orignal fastq file
    sed \
        -n \
        'p;n;n;n;' \
        "$strFastq" |
      "$scriptDirStr/fastqGrep" \
        -stdin-filt \
        -fastq "$fastqFileStr" \
      > "$prefixStr--tmp.fastq";

    #***************************************************************************
    # Sec-5 Sub-2: score reads
    #***************************************************************************

    # output: prefixStr--scores.tsv & prefixStr--filt-reads.fastq
    bash "$scriptDirStr/filterReads.sh" \
            -f "$prefixStr--tmp.fastq" \
            -ref "$fullRefFastaStr" \
            -t "$threadsInt" \
            -p "$prefixStr--tmp" \
            -len-weight "$lenWeightDbl" \
            -mapq-weight "$mapqWeightDbl" \
            -max-dif "$maxErrorRateDbl" \
            -q-base "$baseQInt" \
            -min-length "$minReadLenInt" \
            -max-length "$maxReadLenInt" \
            -mapq "$refFiltMapqDbl" \
            -mean-q "$minMeanQDbl" \
            -med-q "$minMedainQDbl" \
            -align-mean-q "$minAlignMeanQDbl" \
            -align-med-q "$minAlignMedainQDbl" \
            -ins-A "$maxAHomoLenInt" \
            -ins-T "$maxTHomoLenInt" \
            -ins-G "$maxGHomoLenInt" \
            -ins-C "$maxCHomoLenInt";

    #***************************************************************************
    # Sec-5 Sub-3: Pull out reads & extract longest read
    #***************************************************************************

    tmpInt="$(wc -l "$prefixStr--tmp--scores.tsv" | sed 's/[ \t].*//')";
    if [[ "$tmpInt" -lt "$minClustReadsInt" ]]; then
      rm "$strFastq";
      continue;
    fi # if bin had to few reads, discard

    # output: prefixStr--scores.tsv & prefixStr--filt-reads.fastq
    # get the top scoring read
    head -n+1 \
            < "$prefixStr--tmp--scores.tsv" |
        sed 's/[ \t].*//' |
        "$scriptDirStr/fastqGrep" \
            -stdin-filt \
            -fastq "$prefixStr--tmp.fastq" |
        sed -n '1s/^@/>/; 1,2p' \
        > "$prefixStr--tmp--longest-read.fasta";
        
    # get remaining reads
    sed -n "2,$maxMedakaReadsInt{s/[ \t].*//p;}" \
            < "$prefixStr--tmp--scores.tsv" |
         "$scriptDirStr/fastqGrep" \
            -stdin-filt \
            -fastq "$prefixStr--tmp.fastq" \
        > "$prefixStr--tmp--other-reads.fastq";

    #***************************************************************************
    # Sec-5 Sub-4: Polish my reads
    #***************************************************************************

    # Polish my reads
    bash "$scriptDirStr/buildConsensus.sh" \
       -i "$prefixStr--tmp--other-reads.fastq" \
       -I "$prefixStr--tmp--longest-read.fasta" \
       -t "$threadsInt" \
       -p "$prefixStr--$refStr--$clusterStr" \
       -c "$minMaskDepthInt" \
       -r "$roundsRaconInt" \
       -m "$medakaModelStr";

    #***************************************************************************
    # Sec-5 Sub-5: Determine if should keep consensus (to similar to others?)
    #***************************************************************************

    fileStr="";

    # trim consensus to only the region of interest
    bash "$scriptDirStr/trimReads.sh" \
        -f "$prefixStr--$refStr--$clusterStr--consensus.fasta" \
        -t "$threadsInt" \
        -p "$prefixStr--$refStr--$clusterStr--consensus" \
        -ref "$refFastaStr";
   mv "$prefixStr--$refStr--$clusterStr--consensus--"*".fasta" \
      "$prefixStr--$refStr--$clusterStr--consensus--trim.fasta";

    if [[ -f "$prefixStr--final--consensuses--trim.fasta" ]]; then
    # Check if I need to compare consenuses (past consensuses exist)
        fileStr="$( \
          bash "$scriptDirStr/compareConsensuses.sh" \
            -ref "$prefixStr--$refStr--$clusterStr--consensus--trim.fasta" \
            -f "$prefixStr--final--consensuses--trim.fasta" \
            -t "$threadsInt" \
            -min-diff "$minDiffDbl" \
            -min-mismatches "$minMisDbl").fastq"; 
    fi # Check if I need to compare consenuses (past consensuses exist)

    #***************************************************************************
    # Sec-5 Sub-6: Assgin reads for consensus to a final bin
    #***************************************************************************

    if [[ -f "$fileStr" ]]; then
        rm "$prefixStr--$refStr--$clusterStr--consensus.fasta" \
           "$prefixStr--$refStr--$clusterStr--consensus--trim.fasta";
    else
        fileStr="$prefixStr--$refStr--$clusterStr.fastq";

        # Add new consensus to list of consensuses kept
        cat "$prefixStr--$refStr--$clusterStr--consensus--trim.fasta" \
          >> "$prefixStr--final--consensuses--trim.fasta";

    fi # check if need to remove the new consensus

    # Add reads to the final assigned bin
    rm "$strFastq"; # No longer need
    cat "$prefixStr--tmp.fastq" >> "$fileStr";
done # for all kept clusters

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: get read counts & clean up
#    sec-6 sub-1: get read counts
#    sec-5 sub-2: clean up files
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-6 Sub-1: get read counts
#*******************************************************************************

# Get read counts per bin
wc -l "$prefixStr--"*"--cluster"*".fastq" |
  awk \
   'BEGIN{OFS="\t"};
    {print $1 / 4, $2;}; #number lines / lines per fastq entry, file name
  ' > "$prefixStr--counts.tsv";

#*******************************************************************************
# Sec-6 Sub-2: clean up files
#*******************************************************************************

if [[ "$clustStr" == "yes" ]]; then
  rm \
      "$prefixStr--tmp.grep" \
      "$prefixStr--final--consensuses--trim.fasta" \
      "$prefixStr--"*"--consensuses.fasta";
fi # if did clustering, need to remove some extra files

rm \
    -r \
    "$prefixStr--tmp.fastq" \
    "$prefixStr--tmp--scores.tsv" \
    "$prefixStr--tmp--filt-reads.fastq" \
    "$prefixStr--tmp--other-reads.fastq" \
    "$prefixStr--tmp--longest-read.fasta" \
    "$prefixStr--scores.tsv" \
    "$prefixStr--"*"racon"* \
    "$prefixStr--"*"--medaka" \
    "$prefixStr--"*"read-depth.txt";

exit;
