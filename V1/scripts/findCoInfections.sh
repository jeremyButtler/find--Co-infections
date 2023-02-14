#!/bin/bash

################################################################################
# Name: binAndBuidConsensus.sh
# Use: runs severals smaller scripts to bin, subsample, and buid a
#      consensus genome from input reads
# Input:
#    -i: fastq file with reads to bin (string)
#        Required
#    -r: fasta file with references to map reads to [the bins] (string)
#        Required
#    -p: prefix to name the output files (string
#        Default: out
#  Binning variables
#    -Q: Min mapping quality needed to keep a read (integer)
#        Default: 20
#    -q: Min q-score need to keep a read (integer)
#        Default: 13
#    -a: Min aligned length needed to keep a read (integer)
#        Default: 600
#    -k: Min percentage of reads needed to keep a bin (Double)
#        Default: 0.4 (0.4% of reads)
#    -K: Min number of reads to keep a bin (integer)
#        Default 100
#  Sub-sampling variables
#    -s: How many reads to sub-sample (string)
#        Default: 300
#    -n: max read length to keep (integer)
#        Default: 0 (no not filter by max read length)
#  Building CONSENSUS genome variables
#    -c: Min read depth to not mask bases in the consensus genome (integer)
#        Default: 30
#    -r: Number of rounds to run racon (integer)
#        Default: 1
#    -m: Model to use with medaka polishing (string)
#        Default: r941_min_high_g351
#  Checking consensus variables
#    -g: Min % difference between consensus genomes (double)
#        Default: 1% different
#    -x: Min % mismatches different between consensus genomes (double)
#        Default: 0.3%
#    -C: make a circularized database (integer: 1 or 0)
#        Default: 0
# Output:
#    consensus genomes
# Requires:
#    minimap2
#    samtools
#    bamtools # binning
#    filtlong # stats for subsampling
#    racon
#    medaka
#    blastn # check the consensus genomes (need to change out)
# Script depends:
#    bin.sh # bins reads
#    subsample.sh # subsamples reads
#    buildConsensus.sh # builds the consenus genomes 
#    checkConsensus.sh # checks if the consensus genomes were from misbinned reads (to similar)
#    makeBlastDataBase.sh # checkConsensus uses (need to find better way to get percent id)
#    selectConsensus.awk # part of checkConsensus
###############################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-1: Variable declerations
#    sec-1 sub-1: variables to hold user input
#    sec-1 sub-2: script variables
#    sec-1 sub-3: help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: variables to hold user input
#*******************************************************************************

fastqFileStr="";
refFastaStr=""; # references to map reads against (fasta)
prefixStr="out"; # what to call everything
threadsInt=4;

# binning varialbes
minMapQInt=20; # min mapping quality of a read
minQInt=13; # min q-score to keep a read
minAlignedLenInt=600; # min length to keep a read
minBinPercDbl=0.4; # keep bins with at least 0.8% of reads
minBinNumDbl=100; # keep bins with at least 50 reads

# Sub-sampling variables
maxReadLenInt=0; # see what happens with whole genomes
subsampleToInt=300; # how many reads to subsample

# build consensus genome variables
minMaskDepthInt=30; # min depth to not mask bases in conesnsus
roundsRaconInt=1; # number of rounds to run racon
medakaModelStr="r941_min_high_g351"; # models used with guppy

# check consensus genomes (NEED TO ADD IN as options)
minConDiffDbl="99";
minMismatchId="0.3"; # min percentage of mismatches needed to keep a second consensus
circDataBaseBool=0; # do not make a circularized database

#*******************************************************************************
# Sec-1 Sub-2: script variables
#*******************************************************************************

scriptDirStr="$(dirname "$0")"; # directory with this script
longestReadFileStr=""; # CHECK IF NEED TO DELETE

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -i reads.fastq -r references.fasta
    -i: fastq with reads to bin [Required]
    -r: fasta with references to map reads to (the bins) [Required]
    -p: prefix to name all output [Default: out]
    -t: number of threads to use [Default: 4]

  Binning parameters
    -k: min % of mapped reads per bin [Default: 0.4 (0.4% of reads)]
    -K: min number of reads needed to keep a bin [Default: 100]
    -Q: Min mapping quality need to keep a read [Default: 20]
    -q: min quality score needed to keep a read [Default: 13]
    -a: min aligned length need to keep a read [Default: 600]

  Sub-sampling parameters
    -s: How many reads to sub-sample [Default: 300 reads]
    -n: max read length allowed [Default: 0 (any length)]

  Building consensus genome parameters
    -c: Min read depth to not mask bases in the consensus genome (Default: 30)
    -r: Number of rounds to run racon (Default: 1)
    -m: Model to use with medaka polishing (Default: r941_min_high_g351)

  Checking consensus genome parameters:
    -g: Min % difference between consensus genomes [Default: 1% different]
    -x: Min % mismatches different between consensus genomes [Default: 0.3%]
    -C: If not 0 make a circular database when checking the consensus genomes
        [Default: 0 (do not circularize)]
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-2: get and check user input
#    sec-2 sub-1: get user input
#    sec-2 sub-2: check user input
#    sec-2 sub-3: make directorys to hold files
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: get user input
#*******************************************************************************

while getopts 'i:r:p:t:k:K:Q:q:a:s:n:c:r:m:x:g:C:h:z' argsStr;
do # loop through all user input
    case $argsStr in
        i) fastqFileStr="$OPTARG";;
        r) refFastaStr="$OPTARG";;
        p) prefixStr="$OPTARG";;
        t) threadsInt="$OPTARG";;
        k) minBinPercDbl="$OPTARG";;
        K) minBinNumDbl="$OPTARG";;
        Q) minMapQInt="$OPTARG";;
        q) minQInt="$OPTARG";;
        a) minAlignedLenInt="$OPTARG";;
        s) subsampleToInt="$OPTARG";;
        n) maxReadLenInt="$OPTARG";;
        c) minMaskDepthInt="$OPTARG";;
        r) roundsRaconInt="$OPTARG";;
        m) medakaModelStr="$OPTARG";;
        g) minConDiffDbl="$OPTARG";;
        x) minMismatchId="$OPTARG";;
        C) circDataBaseBool="$OPTARG";;
        h) printf "%s\n" "$helpStr"; exit;;
        :) printf "%s\n-%s has no argument\n" "$helpStr" "${OPTARG}"; exit;;
        ?) printf "%s\n-%s is not valid\n" "$helpStr" "${OPTARG}"; exit;;
        z) printf "%s\n-z is not a valid argument\n" "$helpStr" "${OPTARG}"; exit;;
    esac
done # loop through all user input

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$fastqFileStr" ]]; then
    printf "%s is not a valid file\n. Provide a fastq with reads (-i)\n" "$fastqFileStr";
    exit;
elif [[ "$(printf "%s" "$fastqFileStr" | sed 's/.*\(\.fastq\)$/\1/')" != ".fastq" ]]; then
    printf "-i %s is not a fastq file\n" "$fastqFileStr";
    exit;
fi # check if the reads file exists and is a fastq

if [[ ! -f "$refFastaStr" ]]; then
    printf "%s is not a file\n. Provide a fasta with references (-r)\n" "$refFastaStr";
    exit;
elif [[ "$(printf "%s" "$refFastaStr" | sed 's/.*\(\.fasta\)$/\1/')" != ".fasta" ]]; then
    printf "-r %s is not a valide fasta file\n" "$refFastaStr";
    exit;
fi # check if the references to bin with are a valid fasta file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-3: Run scripts to build my consensus genomes
#     sec-3 sub-1: sort reads into bins
#     sec-3 sub-2: Subsample reads
#     sec-3 sub-3: build the consenus genomes
#     sec-3 sub-4: Check and make sure no consensus are identical
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: sort reads into bins
#*******************************************************************************

printf "\nBinning reads\n\n";

bash "$scriptDirStr/bin.sh" \
	-i "$fastqFileStr" \
	-r "$refFastaStr" \
	-t "$threadsInt" \
	-p "$prefixStr" \
	-a "$minAlignedLenInt" \
	-q "$minQInt" \
	-k "$minBinPercDbl" \
	-K "$minBinNumDbl" \
	-Q "$minMapQInt";

mv "$prefixStr--bins/"*tsv ./; # get read count files

#*******************************************************************************
# Sec-3 Sub-2: Subsample reads
#*******************************************************************************

for strFastq in "$prefixStr--bins/"*.fastq; 

do # loop and subsample all fastq files

    if [[ ! -f "$strFastq" ]]; then
        continue;
    fi # if on the null case

    printf "\nSubsampling reads for %s\n\n" "$strFastq";

    refStr="$(printf "%s" "$strFastq" |
				sed "s/.*$prefixStr--//; s/\.fastq//; s/--.*//")";
    bash "$scriptDirStr/subsample.sh" \
		-i "$strFastq" \
		-n "$maxReadLenInt" \
		-t "$threadsInt" \
		-s "$subsampleToInt" \
		-q "$minQInt" \
		-a "$minAlignedLenInt" \
		-p "$prefixStr--$refStr";
done # loop and subsample all fastq files

mkdir "$prefixStr--subsampling";
mv *filtlong-stats.tsv *.fastq "$prefixStr--subsampling";

#*******************************************************************************
# Sec-3 Sub-3: build the consenus genomes
#*******************************************************************************

for strFastq in "$prefixStr--subsampling/"*.fastq; 
do # loop and build consensus genomes for all subsampled fastq files

    if [[ ! -f "$strFastq" ]]; then
        continue;
    fi # if on the null case

    printf "\nBuilding consensus for %s\n\n" "$strFastq";

    # get new prefix
    refStr="$(printf "%s" "$strFastq" |
				sed 's/.*\///; s/\.fastq//; s/--subsample.*//')";
    sed -n '1s/@/>/;
        1s/$/--top-hit/p;
        2p;' \
		< "$strFastq" \
		> "$refStr--top-hit.fasta"; # read to polish with racon

    bash "$scriptDirStr/buildConsensus.sh" \
		-i "$strFastq" \
		-I "$refStr--top-hit.fasta" \
		-t "$threadsInt" \
		-c "$minMaskDepthInt" \
		-r "$roundsRaconInt" \
		-m "$medakaModelStr" \
		-p "$refStr";
done # loop and build consensus genomes for all subsampled fastq file

mkdir "$prefixStr--build-consensus";
mv "$prefixStr"*medaka \
	"$prefixStr"*racon* \
	"$prefixStr"*depth.txt \
	"$prefixStr"*top-hit.fasta \
	"$prefixStr--build-consensus";

#*******************************************************************************
# Sec-3 Sub-4: Check and make sure no consensus are identical
#*******************************************************************************

printf "\nChecking consensus similarity\n\n" "$strFastq";
bash "$scriptDirStr/checkConsensus.sh" \
	-i "$prefixStr--number-reads.tsv" \
	-p "$prefixStr" \
	-a "$minAlignedLenInt" \
	-g "$minConDiffDbl" \
	-x "$minMismatchId" \
	-C "$circDataBaseBool";
