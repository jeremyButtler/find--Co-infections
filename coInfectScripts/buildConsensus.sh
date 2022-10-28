#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC
#     section-1: varialbes decleration
#     section-2: get and check user input
#     section-3: build the consensus genome
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name:
# Use:
# Input:
#    #REQUIRED variables#
#    -i: Fastq file with reads to sample (string)
#        Required
#    -I: Fasta file with the longest read (to polish with racon)
#        Required
#
#    #GENERAL optinal variables#
#    -p: Prefix to call the files and directories (string)
#         Default: out
#    -t: Number of threads to use, when multi-threading is an option (integer)
#       Default: 4
#
#    #CONSENSUS genome optinal variables#
#    -c: Min read depth to not mask bases in the consensus genome (integer)
#        Default: 30
#    -r: Number of rounds to run racon (integer)
#        Default: 1
#    -m: Model to use with medaka polishing (string)
#        Default: r941_min_high_g351
# Output:
# Requirments:
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-1: varialbes decleration
#    sec-1 sub-1: user input
#    sec-1 sub-2: script variables
#    sec-1 sub-3: help message variable
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: user input
#*******************************************************************************

# REQUIRED variables
fastqFileStr=""; # input fastq file with reads
longestReadFileStr="";
prefixStr="out"; # prefix to call files and directories
threadsInt=4; # number of threads to use

# CONSENSUS genome optinal variables
minMaskDepthInt=30; # min depth to not mask bases in conesnsus
roundsRaconInt=1; # number of rounds to run racon
medakaModelStr="r941_min_high_g351"; # fill in later (default super accuracy)

#*******************************************************************************
# Sec-1 Sub-2: script variables 
#*******************************************************************************

raconPrefixStr="racon-round"; # for easy racon renamine in loop
medakaDirStr="medaka"; # holds output from medaka
raconFileStr=""; # holds the output from racon
depthFileStr="read-depth.txt"; # holds the read depth for each base
maskFileStr="mask.bed"; # tells which bases to mask
consensusFileStr="consensus.fasta"; # output
intRound=1; # counter for a loop

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -i reads.fastq [options ...]
    REQUIRED variables
    -i: Fastq file with reads to sample (Required)
    -I: fasta with the longest read (to polish) (Required)

    GENERAL optinal variables
    -p: Prefix to call the files and directories (Default: out)
    -t: Number of threads to use, when multi-threading is an option (Default: 4)

    Consensus genome settings
    -c: Min read depth to not mask bases in the consensus genome (Default: 300)
    -r: Number of rounds to run racon (Default: 1)
    -m: Model to use with medaka polishing (Default: r941_min_high_g351)"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-2:Get and check user input
#    sec-2 sub-1: get user input
#    sec-2 sub-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: get user input
#*******************************************************************************

while getopts ':i:I:h:p:t:q:c:r:m:z' argsStr;
do # loop though all user input
    case $argsStr in
        i) fastqFileStr="$OPTARG";;
        I) longestReadFileStr="$OPTARG";;
        h) printf "$helpStr\n"; exit;;
        p) prefixStr="$OPTARG";;
        t) threadsInt="$OPTARG";;
        c) minMaskDepthInt="$OPTARG";;
        r) roundsRaconInt="$OPTARG";;
        m) medakaModelStr="$OPTARG";;
        :) printf "%s\n -%s needs an argument\n" "$helpStr" "${OPTARG}"; exit 1;;
        ?) printf "%s\n -%s is not valid\n" "$helpStr" "${OPTARG}"; exit 1;;
    esac
done # loop through all user input

#*******************************************************************************
# Sec-2 Sub-2: check user input
#*******************************************************************************

if [[ ! -f "$fastqFileStr" ]];
then # if the fastq file with reads does not exist
    printf "%s is not a valid file\n" "$fastqFileStr";
    exit;
fi # if the fastq file with reads does not exist

if [[ ! -f "$longestReadFileStr" ]]; 
then # if the user did not provide a read to polish with
    printf "Racon needs a read to polish with, please povide one \(-I\)\n";
    exit;
fi # if the user did not provide a read to polish with

# set up file names
raconPrefixStr="$prefixStr--$raconPrefixStr";
medakaDirStr="$prefixStr--$medakaDirStr";
depthFileStr="$prefixStr--$depthFileStr";
maskFileStr="$prefixStr--$maskFileStr";
consensusFileStr="$prefixStr--$consensusFileStr";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-3:
#    sec-3 sub-1: Polish the longest read with racon
#    sec-3 sub-2: Polish the racon consesnus with medaka
#    sec-3 sub-3: Get the read depths with samtools depth
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Polish the longest read with racon
#*******************************************************************************

if [[ ! -f "$raconPrefixStr-$roundsRaconInt.fasta" ]];
then # if need to polish with racon

    while [[ intRound -le "$roundsRaconInt" ]];
    do # loop though all rounds of racon
        printf "Polishing read with racon. Round: %s\n" "$intRound";
        # add the round number to the racon name
        raconFileStr="$raconPrefixStr-$intRound.fasta";

        minimap2 -ax map-ont "$longestReadFileStr" "$fastqFileStr" -t "$threadsInt" \
				> tmp.sam;
        racon -m 8 -x -6 -g -8 -w 500 -t "$threadsInt" "$fastqFileStr" tmp.sam "$longestReadFileStr"  \
				> "$raconFileStr";
        longestReadFileStr="$raconFileStr";
        intRound=$((intRound + 1)); # prep for the next round
    done # loop though all rounds of racon
    
    rm tmp.sam; # no longer need
fi # if need to polish with racon

#*******************************************************************************
# Sec-3 Sub-2: Polish the racon consensus with medaka
#*******************************************************************************

if [[ ! -f "$medakaDirStr/consensus.fasta" ]];
then # if need to make a consensus genome with medaka
    # medaka polish the consensus
    # eval & conda activate are needed to activate medaka in conda
    #eval "$(conda shell.bash hook)"; # installed medaka by miniconda
    #conda activate medaka;
    source ~/medaka/venv/bin/activate; # activate python virtual enviroment

    medaka_consensus -t "$threadsInt" \
		-i "$fastqFileStr" \
		-d "$longestReadFileStr" \
		-m "$medakaModelStr" \
		-o "$medakaDirStr";

    deactivate; # exit python virtal enviorment
    #conda deactivate; # done with medaka (for deactivating conda virtaul env)
fi # if I need to make a consensus genome with medaka

#*******************************************************************************
# Sec-3 Sub-3: Get read depths with sam tools
#*******************************************************************************

if [[ ! -f "$depthFileStr" ]]; 
then # if need to make the depth file
    minimap2 -ax map-ont \
		"$medakaDirStr/consensus.fasta" \
		"$fastqFileStr" \
		-t "$threadsInt" |
		samtools sort |
		samtools view -b > "tmp.bam";
    # get read depth and mask
    samtools depth -a -d 0 tmp.bam > "$depthFileStr"; 
    rm tmp.bam;
fi # if need to make the depth file

#*******************************************************************************
# Sec-3 Sub-4: Mask low read depth bases
#*******************************************************************************

if [[ ! -f "$consensusFileStr" ]];
then # if need to mask the consensus genome
    awk -v minDepthInt="$minMaskDepthInt" \
		-v nameStr="$prefixStr" '
        BEGIN{
             printf "1s/.*/>%s/;\n", nameStr; # print the header
             printf "2s/[ \\t]*//g;\n"; # remove with space (just in case)
         } # BEGIN block
         {if($3 < minDepthInt){printf "2s/./N/%i;\n", $2}} # mask base
        ' < "$depthFileStr" \
		> mask.sed;
     sed -f mask.sed \
		"$medakaDirStr/consensus.fasta" \
		> "$consensusFileStr";
    rm mask.sed;
fi # if need to mask the consensus genome

exit; # done
