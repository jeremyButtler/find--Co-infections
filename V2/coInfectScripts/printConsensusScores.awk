################################################################################
# Name: printConsensusScores.awk
# Use: prints out % mismatches & %difference for consensuses
# Input: 
#    < file.sam: with mapped consensus
#    -v headerBool: 1: prints header for file
# Output:
#    stdout: Query\treference\t%-mismatches\t%-difference\tmapq\tmismatches\t
#            indels-&-mismathes\tsequence-length
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: begin block, print out header if requested
#    sec-2: Find number of mismatches & print out scores
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: begin block, print out header if requested
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
    FS=OFS="\t";
    if(headerBool == 1)
    { # if printing out the header
        print "Query",
              "Reference",
              "perMis",
              "percDiff",
              "mapq",
              "totalError",
              "length";          
    } # if printing out the header
} # BEGIN block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Find number of mismatches & print out scores
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # main Get the percent mismatches & mapq of consensuses
    if($1 ~ /^@/){next;};      # if a header line
    if($6 == "*"){next;};      # if not mapped

    totalMisInt = 0;           # make sure 0 for this read
    lenSeqInt = length($10);   # get length of the sequence

    gsub(/[0-9]*[SDI=]/,"",$6); # remove all but mismatches
    gsub(/[0-9]*$/,"",$6);      # remove matches at end

    numMisInt = split($6, misAryInt, "X");

    for(intMis = 1; intMis <= numMisInt; intMis++)  # find number of mismatches
        totalMisInt = totalMisInt + misAryInt[intMis];

    gsub(/NM:i:/, "", $12);                  # remove tag for total score

    # print query, reference, %mismatches, mapq, mismatches, total errors
    print $1,
          $3,
          totalMisInt/lenSeqInt,
          $12/lenSeqInt,
          $5,
          totalMisInt,
          $12,
          lenSeqInt;
} # main block
