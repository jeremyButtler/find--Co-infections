#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    section-1: variable declerations
#    section-2: get and check user input
#    section-3: map and bin reads
#    section-4: Filter bins by percentage and number. Also, add metadata
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: bin.sh
# Use: bins reads from a fastq using references in a fasta file
# Input:
#    -i: fastq file with reads to bin (string)
#        Required
#    -r: fasta file with references to map reads to [the bins] (string)
#        Required
#    -p: prefix to name the output files (string
#        Default: out
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
# Output:
#    A fastq file for each reference (bin) in the prefix--bins directory
# Requires:
#    minimap2 (to align reads)
#    samtools (convert sam to bam and filter reads)
#    bamtools (to bin reads)
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

fastqFileStr="";
refFastaStr=""; # references to map reads against (fasta)
prefixStr="out"; # what to call everything
threadsInt=4;
minMapQInt=20; # min mapping quality of a read
minQInt=13; # min q-score to keep a read
minAlignedLenInt=600; # min length to keep a read
minBinPercDbl=0.4; # keep bins with at least 0.8% of reads
minBinNumDbl=100; # keep bins with at least 50 reads

#*******************************************************************************
# Sec-1 Sub-2: script variables
#*******************************************************************************

binStr="";
binDirStr="bins"; # name of directory to store bins (prefix add later)
readCountFileStr="number-reads.tsv"; # holds the number of reads per bin
readPercFileStr="percent-reads.tsv"; # holds the percentage of reads per bin
misbinDirStr="miss-bins"; # holds bins with to few reads
keepBinBool=0; # Stores if a bin ment the min % threshold to keep

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -i reads.fastq -r references.fasta
    -i: fastq with reads to bin [Required]
    -r: fasta with references to map reads to (the bins) [Required]
    -p: prefix to name all output [Default: out]
    -t: number of threads to use [Default: 4]
    -k: min % of mapped reads per bin [Default: 0.4 (0.4% of reads)]
    -K: min number of reads needed to keep a bin [Default: 100]
    -Q: Min mapping quality need to keep a read [Default: 20]
    -q: min quality score needed to keep a read [Default: 13]
    -a: min aligned length needt to keep a read [Default: 600]
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

while getopts 'i:r:p:t:k:K:Q:q:a:h:z' argsStr;
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

#*******************************************************************************
# Sec-2 Sub-3: Set directory names and make directories
#*******************************************************************************

binDirStr="$prefixStr--$binDirStr";
readCountFileStr="$prefixStr--$readCountFileStr";
readPercFileStr="$prefixStr--$readPercFileStr";
misbinDirStr="$prefixStr--$misbinDirStr";

if [[ ! -d "$binDirStr" ]]; then
    mkdir "$binDirStr";
fi # if need to make the directory to store the bins

if [[ ! -d "$binDirStr/$misbinDirStr" ]]; then
    mkdir "$binDirStr/$misbinDirStr";
fi # if directory for holding bins with too few reads does not exist

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-3: bin references and get counts
#    sec-3 sub-1: bin references with minimap2 and bamtools
#    sec-3 sub-2: create fastqs
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: bin reference with minimap2 and bamtools
#*******************************************************************************

minimap2 -ax map-ont \
	-t "$threadsInt" \
	"$refFastaStr" \
	"$fastqFileStr" |
	samtools sort -T reads.tmp |
	samtools view \
		-F 0x04 \
		-b \
		-q "$minMapQInt" \
		--min-qlen "$minAlignedLenInt" \
		-e "avg(qual)>=$minQInt" \
		--threads "$threadsInt" \
	> "$prefixStr--bins/tmp.bam";

cd "$prefixStr--bins" || exit;


printf "binning reads with bamtools\n";
bamtools split -in "tmp.bam" -reference;
rm "tmp.bam"; # no longer need # make the bins

#*******************************************************************************
# Sec-3 Sub-2: Create fastqs
#*******************************************************************************

for strBam in ./*REF_*.bam;
do # loop, convert all bam bins to fastq bins

    if [[ ! -f "$strBam" ]]; then
        continue;
    fi # if the null case, just move to the next case

    # get the name of the bin
    binStr="$(printf "%s" "$strBam" | sed 's/.*REF_//; s/\.bam//;')";

    # convert the bam to a fastq
    samtools fastq "$strBam" > "$prefixStr--$binStr.fastq";

    rm "$strBam"; # no longer need
done # loop, convert all bam bins to fastq bins

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-4: Get bin counts and check if bin has enough reads to keep
#    sec-4 sub-1: Get the number of reads per bin
#    sec-4 sub-2: Get the percentage of reads per bin
#    sec-4 sub-3: check if bin has enough reads to keep (by percentage)
#    sec-4 sub-4: check if bins has enough reads to keep (by number)
#    sec-4 sub-5: add metadata
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-4 Sub-1: Get stats
#*******************************************************************************

wc -l *.fastq |
	awk '
        BEGIN{
            refRowStr = "";
            countRowStr = "";
            totalInt = 0; # keep track of all reads
        } # BEGIN block

        {if($2 != "total")
        { # if not on the total row, count the number of reads
            gsub(/\.fastq/, "", $2); # remove the fastq ending
            gsub(/.*--/, "", $2); # remove the prefix
            $1 = $1 / 4; # fastqs have 4 lines per read
            totalInt = totalInt + $1;
            refRowStr = refRowStr "\t" $2;
            countRowStr = countRowStr "\t" $1;
       }} # if not on the total row, count reads

        END{
            print "Total-reads" refRowStr;
            print totalInt countRowStr;
        } # END block, print out counts
    ' > "$readCountFileStr";

#*******************************************************************************
# Sec-4 Sub-2: Print out the percent reads per bin
#*******************************************************************************

awk 'BEGIN{
        totalInt = 0;
        percRowStr = ""; # holds read percentages
    } # BEGIN block

    {if(NR == 1){print $0;} # header row

     else
     { # else on the second row
         for(intBin = 1; intBin <= NF; intBin++)
         { # loop through all bins

             {if(intBin == 1){totalInt = $1}
              else{percRowStr = percRowStr "\t" ($intBin / totalInt) * 100;}
             } # check if on the total column or not

         } # loop through all bins
      } # else on the second (data) row
    } # check which row on

    END{print totalInt percRowStr} # print out the percenatges
' < "$readCountFileStr" > "$readPercFileStr";

#*******************************************************************************
# Sec-4 Sub-3: determine if bin has enough reads to keep
#*******************************************************************************

awk 'BEGIN{FS=OFS="\t"}; # makin temp file to loop through
    {if(NR == 1)
     { # if on the first row
         for(intBin = 2; intBin <= NF; intBin++)
         {binAryStr[intBin - 1] = $intBin}; # get the bin name
     } # if on the frist row
     else
     { # else on the second row
         for(intBin = 2; intBin <= NF; intBin++)
         {percAryStr[intBin - 1] = $intBin}; # get % reads in the bin
     } # else if on the second row
    } # check if on the first or secon row
    END{
        for(intBin = 1; intBin <= length(binAryStr); intBin++)
        {print binAryStr[intBin], percAryStr[intBin];} 
    }' < "$readPercFileStr" > tmp.tsv;

while IFS= read -r binStr;
do # loop though all bins to see if have enough reads
    keepBinBool="$(printf "%s" "$binStr" |
		awk -v minPercDbl="$minBinPercDbl"  '
            {if($2 < minPercDbl){print 0} else {print 1}}')";

    if [[ "$keepBinBool" -le 0 ]]; then
        mv "$prefixStr--$(printf "%s" "$binStr" | awk '{print $1}').fastq" "$misbinDirStr";
    fi # if bin has to few reads to keep
done < tmp.tsv;

#*******************************************************************************
# Sec-4 Sub-4: Check if the bin has the minimum number of reads
#*******************************************************************************

awk 'BEGIN{FS=OFS="\t"}; # makin temp file to loop through
    {if(NR == 1)
     { # if on the first row
         for(intBin = 2; intBin <= NF; intBin++)
         {binAryStr[intBin - 1] = $intBin}; # get the bin name
     } # if on the frist row
     else
     { # else on the second row
         for(intBin = 2; intBin <= NF; intBin++)
         {percAryStr[intBin - 1] = $intBin}; # get # reads in the bin
     } # else if on the second row
    } # check if on the first or secon row
    END{
        for(intBin = 1; intBin <= length(binAryStr); intBin++)
        {print binAryStr[intBin], percAryStr[intBin];} 
    }' < "$readCountFileStr" > tmp.tsv;

while IFS= read -r binStr;
do # loop though all bins to see if have enough reads
    keepBinBool="$(printf "%s" "$binStr" |
		awk -v minNumDbl="$minBinNumDbl"  '
            {if($2 < minNumDbl){print 0} else {print 1}}')";

    if [[ "$keepBinBool" -le 0 ]]; then
        mv "$prefixStr--$(printf "%s" "$binStr" | awk '{print $1}').fastq" "$misbinDirStr";
    fi # if bin has to few reads to keep
done < tmp.tsv;

#*******************************************************************************
# Sec-4 Sub-5: add the meta data to the files with bin counts
#*******************************************************************************

paste <(printf "Run-name\tFastq\tReferences\t\n%s\t%s\t%s\t" \
		"$prefixStr" \
		"$(basename "$fastqFileStr")" \
		"$(basename "$refFastaStr")") \
	"$readCountFileStr" \
	> tmp.tsv;
mv tmp.tsv "$readCountFileStr";

paste <(printf "Run-name\tFastq\tReferences\t\n%s\t%s\t%s\t" \
		"$prefixStr" \
		"$(basename "$fastqFileStr")" \
		"$(basename "$refFastaStr")") \
	"$readPercFileStr" \
	> tmp.tsv;
mv tmp.tsv "$readPercFileStr";

exit;
