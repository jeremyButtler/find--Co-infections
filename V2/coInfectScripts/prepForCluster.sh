#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: Variable declerations
#   sec-2: Get and check user input
#   sec-3: Do an all versus all read mapping
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
readsFastqStr="";                             # reads to get mapq's for
minMapqDbl="13";                              # min mapq to keep an aligment
threadsInt=3;
maxErrorRateDbl=0.07;                       # max diference to keep alignemnt

# scoreReads parameters
baseQInt=13;                                  # min Q-score for each base
minReadLenInt=600;                            # min read length to keep aligment
maxReadLenInt=1000;                           # max read length to keep aligment
minMeanQDbl=20;                              # min mean Q-score to keep aligment
minMedainQDbl=20;                         # min medain Q-score to keep alignment
minAlignMeanQDbl=20;                         # min aligned mean Q-score
minAlignMedainQDbl=20;                       # min aligned median Q-score

# max homopolymer sizes to keep inserts from
maxAHomoLenInt=1;
maxTHomoLenInt=1;
maxGHomoLenInt=1;
maxCHomoLenInt=1;

# script variables
prefInt="$(date +"%Y%m%d%H%M%S")$((RANDOM))"; #prefix: date, time & random num.
intCnt=0;                                      # counter for loop
scriptDirStr="$(dirname "$0")";                # directory with this script

# help message
helpStr="$(basename "$0") -f reads.fastq [options ...]
        -h: print this help message & exit
        -f: file to do all versus all mapping on
        -max-dif: max difference between query & reference to keep aigment
        -q-base: min Q-score to keep a mismatch or insertion
        -min-length: min read length to keep aligment
        -max-length: max read length to keep aligment
        -mapq: min mapping quality to keep aligment
        -mean-q: min mean Q-score to keep alignemnt
        -med-q: min median Q-score to keep aligment
        -align-mean-q: min aligned mean Q-score
        -align-med-q: min aligned median Q-score
        -ins-A: max size of an A homopolymer size to keep an insertion
        -ins-T: max size of an A homopolymer size to keep an insertion
        -ins-G: max size of an A homopolymer size to keep an insertion
        -ins-C: max size of an A homopolymer size to keep an insertion
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
        -f) readsFastqStr="$2";;          # file for all versus all mapping
        -t) threadsInt="$2";;             # number of threads to use
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

    shift;                               # move to arugment
    shift;                               # move to next parameter
done # whild have cmd arguments

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$readsFastqStr" ]]; then
    printf "%s is not a file, please povide file with -f\n" "$readsFastqStr";
    exit;
fi # Check if input file was valid

if [[ $minMapqDbl -lt 0 ]]; then
    printf "Mapping quality must be > -1\n";
    exit;
fi # check if the mapping qaulity is a valid

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Do an all versus all read mapping
#    sec-3 sub-1: Find the number of reads in the file & loop though all reads
#    sec-3 sub-2: Make 2 tempory fastq files: 1: with only ref & 2: no ref read
#    sec-3 sub-3: Map reads to ref read & keep read name, ref name, & mapq
#    sec-3 sub-4: Clean up temporary files and exit
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Find the number of reads in the file & loop thruogh all reads
#*******************************************************************************

numReadsInt="$(wc -l "$readsFastqStr" | awk '{print $1}')";
numReadsInt="$((numReadsInt / 4))";
    # wc -l: gets the number of reads
    # awk: prints only the number
    # $( / 4()) divides the line count by 4 (4 lines per fastq entry)

printf "Mapping reads to read number: 0\n" >&2;

while [[ $intCnt -lt $numReadsInt ]]; do   # intCnt is index 0
# while there are references to map

    #***************************************************************************
    # Sec-3 Sub-2: Make 2 tempory fastq files: 1: with only ref & 2: no ref read
    #***************************************************************************

    awk -f "$scriptDirStr/fastqExtractReadNumber.awk" \
            -v readNumInt="$intCnt" \
            -v refFileStr="$prefInt--tmp-ref.fastq" \
            -v readFileStr="$prefInt--tmp-reads.fastq" \
            < "$readsFastqStr";

    #***************************************************************************
    # Sec-3 Sub-3: Map reads to ref read & keep read name, ref name, & mapq
    #***************************************************************************

    tput cuu 2 cuf 30 el >&2;                   # go up one line
    #tput cuu 1 cub 40 ed >&2;                   # go up one line
        # cuu 2: go back two lines
        # cuf 32: move cursor right 32 spaces (left is cub)
        # el clear to end of line 
        # &>2: redirect to stderr

    intCnt=$((intCnt + 1));          # Count the next read
    printf "%s\n" \ "$intCnt" \ >&2;

    # Using minimap2 with ava-ont settings 
    minimap2 -a \
			-k15 \
			-w5 \
			-e0 \
			-m100 \
			-r2k \
			-t "$threadsInt" \
			--eqx \
			"$prefInt--tmp-ref.fastq" \
			"$prefInt--tmp-reads.fastq" \
			2> /dev/null |
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
		awk -v maxDiffDbl="$maxErrorRateDbl" \
				-f "$scriptDirStr/filterScoreReadsOut.awk" |
        awk 'BEGIN{OFS="\t"}; {print $1, $2, $3}';   # filter out aligned length
done # while there are references to map

#*******************************************************************************
# Sec-3 Sub-4: Clean up temporary files and exit
#*******************************************************************************

# remove the temporary files
rm "$prefInt--tmp-ref.fastq" \
	"$prefInt--tmp-reads.fastq";

exit;
