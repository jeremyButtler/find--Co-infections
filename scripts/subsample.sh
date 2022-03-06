#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    section-1: variable declerations
#    section-2: get and check user input
#    section-3: subsample reads with filtlong
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: subsample.sh
# Use: Subsamples reads to input amount
# Input:
#    -i: fastq with reads to subsample (string)
#        Required
#    -s: How many reads to subsample (string)
#        Default: 300
#    -n: max read length to keep (integer)
#        Default: 0 (no not filter by max read length)
#    -a: min read length to keep (integer)
#        Default: 600
#    -q: min mean q-score to keep a read (integer)
#        Default: 13
#    -t: number of threads to use (integer)
#        Default: 4
#    -p: prefix to call everything (string)
#        Default: out
# Output:
#    Files: fastq file with reads (prefix--subsample.fastq)
#    File: tsv with the stats output from filtlong (prefix--filtlong-stats.tsv)
# Requires:
#    filtlong version 0.2.1
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-1: Variable declerations
#    sec-1 sub-1: variables to hold user input
#    sec-1 sub-2: script variables
#    sec-1 sub-3: help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: variables to hold user input
#*******************************************************************************

readsFastqStr=""; # holds path to reads to filter
subsampleToInt=300; # how many reads to subsample
maxReadLenInt=0; # max read length
minReadLenInt=600; # min read length
minQInt=13; # min mean q-score to keep
threadsInt=4;
prefixStr="out"; # what to name everthing

#*******************************************************************************
# Sec-1 Sub-2: script variables
#*******************************************************************************

scriptDirStr="$(dirname "$0")"; # directory this script is in
filtStatsTsv="filtlong-stats.tsv"; # file to hold stats from filtlong
outFastqStr="subsample.fastq"; # fastq with subsampled reads

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -i reads.fastq [options ...]
    -i: fastq with reads to subsample [Required]
    -s: How many reads to subsample [Default: 300 reads]
    -n: max read length allowed [Default: 0 (any length)]
    -a: min read length to keep a read [Default: 600 bases]
    -q: min mean q-score to keep a read [Default: 13]
    -t: number of threads to use [Default: 4]
    -p: prefix to call everything [Default: out]
"

# Section-2: get and check user input
#    sec-2 sub-1: get user input
#    sec-2 sub-2: check user input

while getopts 'i:s:n:a:q:t:p:h:z' argsStr;
do # loop through all user input
    case $argsStr in
        i) readsFastqStr="$OPTARG";;
        s) subsampleToInt="$OPTRAG";;
        n) maxReadLenInt="$OPTARG";;
        a) minReadLenInt="$OPTARG";;
        q) minQInt="$OPTARG";;
        t) threadsInt="$OPTARG";;
        p) prefixStr="$OPTARG";;
        h) printf "%s\n" "$helpStr"; exit;;
        :) printf "%s\n-%s argument not valid\n" "$helpStr" "${OPTARG}; exit;;
        ?) printf "%s\n-%s is not a paramter\n" "$helpStr" "${OPTARG}; exit;;
        z) printf "%s\n-z is not a paramter\n" "$helpStr"; exit;;
    esac
done # loop through all user input

# Sec-2 Sub-2: check user input

if [[ ! -f "$readsFastqStr" ]]; then
    printf "%s does not exist. Input fastq with reads using -i\n" "$readsFastqStr";
    exit;
elif [[ "$(printf "%s" "$readsFastqStr" | sed 's/.*\.fastq/.fastq/')" != ".fastq" ]]; then
    print "%s is not a fastq file\n" "$readsFastqStr";
    exit;
fi # check if reads is a valid file

# Sec-2 Sub-3: Set output variable names and make directories

filtStatsTsv="$prefixStr--$filtStatsTsv";
outFastqStr="$prefixStr--$outFastqStr";

# Section-3: subsample reads
#    sec-3 sub-1: Get stats with filtlong
#    sec-3 sub-2: Prepare the filtlong stats for subsampling

# Sec-3 Sub-1: Get stats with filtlong

if [[ ! -f "$filtStatsTsv" ]]; then
# if need to get stats from filtong

    # filter out the read lengths not in range
    awk -v minLenInt="$minReadLenInt" \
		-v maxLenInt="$maxReadLenInt" \
		'{
                headStr = $0; # Header line of a read
                getline seqStr; # sequence line for a read
                getline plusStr; 
                getline qStr; # q-score line of a read
    
                lenSeqInt = length(seqStr);
                if(maxLenInt == 0 || maxLenInt >= lenSeqInt && minLenInt <= lenSeqInt)
                {printf "%s\n%s\n%s\n%s\n", headStr, seqStr, plusStr, qStr};
            }' < "$readsFastqStr" \
		> "$prefixStr--tmp.fastq";
    filtlong --min_mean_q "$minQInt" \
		--verbose \
		"$prefixStr--tmp.fastq" \
		2> "$prefixStr--tmp.tsv" \
		1> /dev/null; # will subsample fastq with filtlong output
    rm "$prefixStr--tmp.fastq"; # no longer need
fi # if need to get stats with filtlong

# Sec-3 Sub-2: Prepare the filtlong stats file for subsampling

# sort stats by the score filtong assigned (this is really quick)
if [[ ! -f "$filtStatsTsv" ]]; then
# if have not formated filtlong output yet
    headStr="$(sed -n '/Read/p' < "$prefixStr--tmp.tsv")"; # save the header
    sed '1,/Read/d; /^$/d; /Output/d' \
		< "$prefixStr--tmp.tsv" |
		awk -v minQInt="$minQInt" '
                BEGIN{FS=OFS="\t"};
                {if($3 >= minQInt){print $0}}' |
		sort -k 5 -nr > "$prefixStr--tmp2.tsv";
    # sed: Deletes the non table lines
    # awk 1: removes reads with low q-scores
    # sort: Sorts enteries by the total score (fith column)

    cat <(printf "%s\n" "$headStr") \
		"$prefixStr--tmp2.tsv" \
		> "$filtStatsTsv";
    rm "$prefixStr--tmp2.tsv"; # no longer need
fi # if have no formated the filtlong output yet

# Sec-3 Sub-3: Subsample reads

# grab the top scoring reads (making a grep filter)
sed -n '2,$p' < "$filtStatsTsv" |
	head -n "$subsampleToInt" |
	awk '{print $1}' > "$prefixStr--tmp.tsv";

# subsample reads
grep -f "$prefixStr--tmp.tsv" \
	-A 3 \
	< "$readsFastqStr" |
	sed '/^--$/d' \
	> "$outFastqStr";
rm "$prefixStr--tmp.tsv";

exit;
