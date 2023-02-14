#!/bin/awk

################################################################################
# Name: scoreSameFileReads.awk
# Use: Scores the reads based on aligned length & mapping quality
# Input:
#    < file.sam: sam file with reads to score                    [Required]
#    -v lenWeightDbl: Weight to apply to the aligned length      [Default: 1]
#    -v mapqWeightDbl: Weight to apply to mapping quality        [Default: 1]
# Ouput:
#     stdout: tsv: query-name\treference-name\tmapping-quality\taligned-length\n
# Note: Sam file cigar needs to be in eqx from (minimap2 --eqx -a)
################################################################################

BEGIN{
    FS=OFS="\t";
    totalBasesInt = 0;             # How may matches, mismatches, & insertions
    lenEntryInt = 0;               # Number of kept cigar entries

    if(lenWeightDbl == 0){lenWeightDbl = 1};
    if(mapqWeightDbl == 0){mapqWeightDbl = 1};
}; # BEGIN block

{ # main brace
    if($1 ~ /^@/){next;};          # on header line
    if($6 == "*"){next;};          # the read did not map to reference

    sub(/^[0-9]*S/, "", $6);       # remove soft masks at start
    sub(/[0-9]*S$/, "", $6);       # remove soft masks at end
    gsub(/[0-9]*[DH]/, "", $6);    # remove deletions from count
    lenEntryInt = split($6, basesInt, /[IX=]/);  # get array of numbers

    for(intCnt = 1; intCnt <= lenEntryInt; intCnt++)
        totalBasesInt = totalBasesInt + basesInt[intCnt];

    # score read & output query, reference, score = sqrt(aligned length) + mapq
    print $1, $3, lenWeightDbl * sqrt(totalBasesInt) + mapqWeightDbl * $5;
        
    totalBasesInt = 0;               # For next entry
}; # main brace
