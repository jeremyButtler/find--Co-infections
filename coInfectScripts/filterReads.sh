#!/usr/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Filter reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: filterReads.sh
# Use: Filters input fastq using input fastq & scoreReads
# Input:
#   -f: reads to align               (Required)
#   -ref: reference to map reads to  [Required]
#   -p: prefix to name everything
# Output:
#   File: fastq file with kept reads ("$prefixStr--filt-reads.fastq")
#   File: prints score for each fastq file to "$prefixStr--scores.tsv"
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
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# user input variables
fastqFileStr="";                              # reads to get mapq's for
minMapqDbl="20";                               # min mapq to keep an aligment
threadsInt=3;                                  # number threads to use
refFastaStr="";                             # reference to map & trim reads with
prefixStr="out";                                # prefix to name everything

# scorin variables
lenWeightDbl=1.25;                          # Weight of read length in scoring
mapqWeightDbl=1;                          # Weight of mapping quality in scoring

# filtering parameters
maxErrorRateDbl=0.08;                        # max difference between read & ref
baseQInt=8;                                   # min Q-score for each base
minReadLenInt=600;                            # min read length to keep aligment
maxReadLenInt=1000;                           # max read length to keep aligment
minMeanQDbl=13;                              # min mean Q-score to keep aligment
minMedainQDbl=13;                         # min medain Q-score to keep alignment
minAlignMeanQDbl=13;                         # min aligned mean Q-score
minAlignMedainQDbl=13;                       # min aligned median Q-score

# max homopolymer sizes to keep inserts from
maxAHomoLenInt=1;
maxTHomoLenInt=1;
maxGHomoLenInt=1;
maxCHomoLenInt=1;

# script variables
scriptDirStr="$(dirname "$0")";                # directory with this script

# help message
helpStr="$(basename "$0") -f reads.fastq [options ...]
        -h: print this help message & exit
        -f: file to do all versus all mapping on [Required]
        -t: number of threads to use with minimap2 [Default: 3]
        -ref: fasta with single reference to map & trim reads with [Required]

        Read scoring:
        -len-weight: Weight of aligned read length in scoring [Default: 1.25]
            - Score is square root of read length
        -mapq-weight: Weight of mapping quality in scoring [Default: 1]

        Read filtering input:
        -max-dif: max difference between query & reference to keep aigment
            - Default: 0.08%
        -q-base: min Q-score to keep a mismatch or insertion
            - Default: 8
        -min-length: min read length to keep aligment
            - Default: 600
        -max-length: max read length to keep aligment
            - Default: 1000
        -mapq: min mapping quality to keep aligment
            - Default: 20
        -mean-q: min mean Q-score to keep alignemnt
            - Default: 20
        -med-q: min median Q-score to keep aligment
            - Default: 20
        -align-mean-q: min aligned mean Q-score
            - Default: 20
        -align-med-q: min aligned median Q-score
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
        -t) threadsInt="$2";;             # number of threads to use
        -p) prefixStr="$2";;              # prefix to name output files
        -ref) refFastaStr="$2";;          # if is a reference genome

        # calcualating final read scores
        -len-weight) lenWeightDbl="$2";;  # Weight of read length in scoring
        -mapq-weight) mapqWeightDbl="$2";;# Weight of mapping quality in scoring

        # Scoring input
        -max-dif) maxErrorRateDbl="$2";;  # max difference between query & ref
        -q-base) baseQInt="$2";;          # min Q-score to keep a bases
        -min-length) minReadLenInt="$2";; # min read length to keep aligment
        -max-length) maxReadLenInt="$2";; # max read length to keep aligment
        -mapq) minMapqDbl="$2";;     # min mapping quality to keep aligment
        -mean-q) minMeanQDbl="$2";;       # min mean Q-score to keep alignemnt
        -med-q) minMedainQDbl="$2";;      # min median Q-score to keep aligment
        -align-mean-q) minAlignMeanQDbl="$2";; # min aligned mean Q-score
        -align-med-q) minAlignMeanQDbl="$2";;  # min aligned median Q-score
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

if [[ ! -f "$refFastaStr" ]]; then
    printf \
		"%s is not a file, please provide fasta file with reference (-ref)\n" \
		"$refFastaStr";
    exit;
fi

if [[ $minMapqDbl -lt 0 ]]; then
    printf "Mapping quality must be > -1\n";
    exit;
fi # check if the mapping qaulity is a valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: align, filter, & score reads
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   # --sam-hit-only: does not output unmapped reads
   # --secondary=no: Only output the best alignment

minimap2 \
        --sam-hit-only \
        --secondary=no \
        --eqx \
		-ax map-ont \
		-t "$threadsInt" \
		"$refFastaStr" \
		"$fastqFileStr" |
	"$scriptDirStr/scoreReads" \
		-stdin \
		-min-q "$baseQInt" \
		-min-map-q "$minMapqDbl" \
		-min-read-length "$minReadLenInt" \
		-max-read-length "$maxReadLenInt" \
		-min-mean-q "$minMeanQDbl" \
		-min-median-q "$minMedainQDbl" \
		-min-aligned-mean-q "$minAlignMeanQDbl" \
		-min-aligned-median-q "$minAlignMedainQDbl" \
		-ins-A "$maxAHomoLenInt" \
		-ins-T "$maxTHomoLenInt" \
		-ins-G "$maxGHomoLenInt" \
		-ins-C "$maxCHomoLenInt" |
	awk -f "$scriptDirStr/filterScoreReadsOut.awk" \
			-v maxDiffDbl="$maxErrorRateDbl" |
	awk -v lenWeightDbl="$lenWeightDbl" \
			-v mapqWeightDbl="$mapqWeightDbl" \
			'BEGIN{OFS="\t"};
             {print lenWeightDbl * sqrt($4) + mapqWeightDbl * $3, $2}
    ' |
	sort -k 2 -r |
	uniq -f 1 |
    awk 'BEGIN{OFS="\t"}; {print $2, $1}' |
    sort -n -r -k 2 \
	> "$prefixStr--scores.tsv";
    # fliping around because unique will not just look at first feild

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Extract filtered reads to a fastq file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

awk '{print $1}' \
		< "$prefixStr--scores.tsv" |
	"$scriptDirStr/fastqGrep" \
			-stdin-filt \
			-fastq "$fastqFileStr" \
	> "$prefixStr--filt-reads.fastq";

exit;
