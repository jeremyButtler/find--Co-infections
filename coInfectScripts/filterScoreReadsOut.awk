#!/bin/awk

################################################################################
# Name: filterScoreReadsOut.awk
# Use: Filters stats output from scoreReads & decides wich alignments to keep
# Input:
#    < file.tsv:
#        tsv output from scoreReads [Required]
#    -v maxDiffDbl:
#        Max difference to keep (mismatches + indels)/aligned length
#            [Default: 0.07 (7%)]
# Output: 
#    stdout: tsv: reference-name\tquery-name\tscore\taligned-length\n
# Note: mapq handled by scoreReads
################################################################################

BEGIN{
    FS=OFS="\t";
    if(maxDiffDbl == 0){maxDiffDbl = 0.07;};
}; # BEGIN block

{ # main block
    # $1 read, $2 ref, $3 mapq, $4 len, $5 aligned len, $6 matches,
    # $7 kept matches, $8 mis, $9 ins, $10 del, $11 med Q,
    # $12 mean Q

    if($1 ~ /^[Rr]ead/){next;};     # if on a header row
    if($1 ~ /^$/){next;};           # if on a header row
    if($5 < 1){next;};              # no matches kept

    percErrorDbl = ($8 + $9 + $10) / $5;
        # (mis + ins + del) / aligned length
    if(percErrorDbl > maxDiffDbl){next;}; # if error rate to high

    print $2, $1, $3, $5;     # print ref name, read name, mapq, aligned length
} # main block
