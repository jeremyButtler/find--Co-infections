################################################################################
# Name: checkConsensus.sh
# Use: checks in the consensus genomes are to similar (keeps consensus genome
#      with the most reads).
# Input:
#    -i: File with read counts to add to the blast results (string)
#        Required (file made by bin.sh)
#    -p: Prefix to name everything (string)
#        Default: out
#    -a: min aligned length need to count a consensus genome hit 
#        Default: 600
#    -g: Min % difference between consensus genomes (double)
#        Default: 1% different
#    -x: Min % mismatches different between consensus genomes (double)
#        Default: 0.3%
#    -C: make a circularized database (integer: 1 or 0)
#        Default: 0
# Output:
#    Moves to similar consensus genome to a directory
#        - consensus-filter
#    Also outputs a blast report on the kept consensus genomes:
#        - kept-consensus.tsv
# Requires:
#    blastn # check the consensus genomes (need to change out)
# Script depends:
#    makeBlastDataBase.sh # checkConsensus uses (need to find better way to get percent id)
#    selectConsensus.awk # does the actual selection
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

# HAVE MIN ALIGNED LENGTH (USE MIN READ LENGTH)
readCountsFileStr="";
prefixStr="";
minConDiffDbl="99";
minMismatchId="0.3"; # min percentage of mismatches needed to keep a second consensus
minAlignedLenInt=600; # min length to keep a read
circDataBaseBool=0; # do not make a circularized database

#*******************************************************************************
# Sec-1 Sub-2: script variables
#*******************************************************************************

scriptDirStr="$(dirname "$0")"; # location of scripts (including this one)

#*******************************************************************************
# Sec-1 Sub-3: help message
#*******************************************************************************

helpStr="$(basename "$0") -i reads.fastq -r references.fasta
    -i: File with read counts to add to the blast results [Required]
    -p: prefix to name all output [Default: out]
    -g : Min % difference between consensus genomes [Default: 1% different]
    -x : Min % mismatches different between consensus genomes [Default: 0.3%]
    -a: min aligned length need to keep a read [Default: 600]
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

while getopts 'i:p:a:x:g:C:h:z' argsStr;
do # loop through all user input
    case $argsStr in
        i) readCountsFileStr="$OPTARG";;
        p) prefixStr="$OPTARG";;
        a) minAlignedLenInt="$OPTARG";;
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

if [[ ! -f "$readCountsFileStr" ]]; then
    printf "%s is not a valid file\n. Provide a fastq with reads (-i)\n" \
		"$readCountsFileStr";
    exit;
fi # check if the file with the read counts exists

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-3: Determine consensus quality
#     sec-3 sub-1: blast all consensus against each other to find % id
#     sec-3 sub-2: make sed script to add number of reads to the blast file
#     sec-3 sub-3: Add the number of reads for each conseus genome to blast results
#     sec-3 sub-4: add the percent mismatches to each hit and complete header
#     sec-3 sub-5: Remove miss-binned consensus genomes (high similarity to others)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: find consensus % ids to eacho ther with blast
#*******************************************************************************

# mainly here so I have something to check for miss-assigned bins
cat *.fasta | sed 's/> />/' > tmp.fasta # file with all consensus genomes

if [[ "$(wc -l tmp.fasta | sed 's/[ \t].*//')" -le 2 ]]; then
    rm tmp.fasta; 
    exit; # no need to sort, on one consensus genome
fi # if just one consensus genome

bash "$scriptDirStr/makeBlastDataBase.sh" tmp.fasta tmpData "$circDataBaseBool";
blastn -query tmp.fasta \
	-db tmpData/tmpData--database \
	-outfmt "6 qseqid sacc length pident bitscore mismatch gaps gapopen" \
	-max_target_seqs 5 |
	awk 'BEGIN{
            FS=OFS="\t";
            print "consensus\tmatch\tlength\tpercent-id\tbitscore\tmismatches\tgaps\tgapopen"
        } # BEGIN block
        {if($4 < 100){print $0}}
    ' > "$prefixStr--consensus-report.tsv";
rm -r tmp.fasta tmpData; # files no longer needed

#*******************************************************************************
# Sec-3 Sub-2: make sed script to add reads per bin to the blast results
#*******************************************************************************

# add in read depths not the best way, but is a quick way
awk 'BEGIN{
        FS=OFS="\t";
     } # BEGIN block
     {if(NR == 1)
      { # put each row into an array
          {for(intCol = 5; intCol <= NF; intCol++) # skip meta data (first four columns)
           {headerStr[intCol] = $intCol;}
          } # loop and put all columns in an array
      } # put each column of header in an array
      
      else
      { # else print out number
          {for(intCol = 5; intCol <= NF; intCol++) # skip meta data (first four columns)
           {{if($intCol > 0){print headerStr[intCol], $intCol;}}}
          } # loop and print out all references with more then 1 read
     }} # check if on the first row or not
' < "$readCountsFileStr" |
	sed -n '1,$s/\([^\t]*\).\(.*\)/\/^[^\\t]*\1\/{s\/\$\/\\t\2\/};/; p' \
	> tmp.sed;
    # awk puts all reads in a line format
    # sed makes a sed script the adds the reads to the end
    # sed output: /^[\t]*sequenceName/{s/$/\tnumberReads/};
    # first line output by awk is usless and messes up sed

#*******************************************************************************
# Sec-3 Sub-3: Add the read percentage to there consensus genomes
#*******************************************************************************

# add the read counts to the blast report
sed -f tmp.sed "$prefixStr--consensus-report.tsv" > tmp.tsv;

# add the number of reads for the subject sequence
sed 's/*/*\\t[^\\t]*/' < tmp.sed > tmp2.sed;
sed -f tmp2.sed < tmp.tsv > tmp2.tsv;

rm tmp.sed tmp2.sed tmp.tsv; # no longer needed

#*******************************************************************************
# Sec-3 Sub-4: Add the percentage of mismatches (between consensus) and complete header
#*******************************************************************************

awk 'BEGIN{FS=OFS="\t"}
     {if(NR > 1){print $0, 100 * $6 / $3} # print out percent of mismatches
      else{print $0, "query-read-number", "subject-read-number", "Mismatch-percent"}
     } # check if on a header or hit
    ' < tmp2.tsv > "$prefixStr--consensus-blast.tsv";
rm tmp2.tsv; # no longer need the temp file

#*******************************************************************************
# Sec-3 Sub-5: Filter out miss-bined consensus genomes (very similar to others)
#*******************************************************************************

# get the list of consensus genomes to keep
awk -f "$scriptDirStr/selectConsensus.awk" \
	-v minDiffDbl="$minConDiffDbl" \
	-v minMisDbl="$minMismatchId" \
	-v minLenInt="$minAlignedLenInt" \
	< "$prefixStr--consensus-blast.tsv" \
	> "$prefixStr--kept-consensus.tsv";

if [[ ! -d "$prefixStr--consensus-filter" ]]; then
    mkdir "$prefixStr--consensus-filter";
fi # if need to make the directory to store less usefull files

mv "$prefixStr--consensus-blast.tsv" "$prefixStr--consensus-filter";
mv *.fasta "$prefixStr--consensus-filter";

while IFS= read -r goodConStr;
do # loop and copy all good consensus 
    if [[ ! -f "$prefixStr--consensus-filter/$goodConStr--consensus.fasta" ]]; then
        continue;
    fi # if on the null case

    cp "$prefixStr--consensus-filter/$goodConStr--consensus.fasta" ./;
done < <(tail -n+2 "$prefixStr--kept-consensus.tsv" | sed 's/\t.*//');

exit;
