#!/bin/awk -f

################################################################################
# TOC
#    section-1: variable declerations and set defaults
#    section-2: print header and selected consensus to keep (and top hit)
#    section-3: END block, print the last kept hit (if there is one)
################################################################################

################################################################################
# Name: selectConsensus.awk
# Use: filters consensus genomes based on identiy to other consensus genomes
#      Also takes read count into consdieration
# Input:
#     File: tsv file output by blast (-outfmt "6 qseqid sacc length pident mismtach gap gapopen")
#           with the read counts for all references add
#        - Row one should be a header or empty (is ignored)
#           - Column one: query
#           - Column two: subject (blast database hit)
#           - Column three: the aligned length
#           - Column four: percent idenity of hit
#           - Column nine: number of reads assigned to the consensus
#           - Column ten: number of reads assigned to the subject
#           - Column eleven: Percent of mismtaches (100 * mismatches / alinged length)
#     -v minDiffDbl: Min difference needed to keep a consensus (double)
#         Default: 98.5
#     -v minMisDbl: Min number of mismatches to keep a consensus (double)
#         Default: 0.3
#     -v minLenInt: Min aligned length needed to keep a blast hit (integer)
#         Default: 600
# Output:
#     File: tsv with the top hits for each reference kept
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-1: Begin block
#     sec-1 sub-1: variable declerations and checks
#     sec-1 sub-2: set the default values if nothing provided
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: Begin block (variable declerations and checks
#*******************************************************************************

BEGIN{
    FS=OFS="\t";
    seqOnStr = "";
    moveOnBool = 0; # if sequence is dupicate or alreayd printed ignore
    bestPrintBool = 0; # just get a good print line

    #***************************************************************************
    # Sec-1 Sub-2: Set defaults (if not provided)
    #***************************************************************************

    {if(minDiffDbl == ""){minDiffDbl = 98.5}}; # default value
    {if(minMisDbl == ""){minMisDbl = 0.3}}; # default value
    {if(minLenInt == ""){minLenInt = 600}}; # set the min length to keep a hit
} # BEGIN block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Section-2: filter the reads
#    sec-2 sub-1: Moving to the next sequence or print the header
#    sec-2 sub-2: check if sequence meets criteria to print
#    sec-2 sub-3: print out a last match if there was one
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: print header or deal with changing sequences
#*******************************************************************************

{if(NR == 1){print $0} # on the header

 else
 { # else on the data rows

     {if(seqOnStr != $1)
      { # if on the next sequence
          seqOnStr = $1;
          moveOnBool = 0; # restarting
          bestPrintBool = 0;

         {if(keepSeqStr != "")
          { # if need to print out a best hit
              print keepSeqStr;
              keepSeqStr = "";
         }} # if need to print out a best hit

     }} # if on the next sequence

    #***************************************************************************
    # Sec-2 Sub-2: Check if sequence is good or bad
    #***************************************************************************

     {if($3 < minLenInt){next;} # if the hit is not worth keeping

      else if($4 <= minDiffDbl && $11 >= minMisDbl)
      { # if the sequence is unique and have not already printed
         {if(moveOnBool < 1){keepSeqStr = $0;}};
         moveOnBool = 2; # printed hit out just move on
      } # if the sequence is unique and have not already printed

     else if(moveOnBool <= 1 && $4 > minDiffDbl || $11 < minMisDbl && seqOnStr != $2)
     { # else if sequence is not unique and have not discarded or printed

         {if($9 > $10)
          { # if may be the best hit of the range
              {if(bestPrintBool == 0){keepSeqStr = $0;}} # if was the best hit print
              moveOnBool = 1; # only check to make sure is the best (most reads)
          } # if this may be the best hit (most reads)
    
          else
          { # else not a good match, move one
              keepSeqStr = "";
              moveOnBool = 2;
          }  # else not a good match, move on
         } # check if a good match or if should ignore
     } # else if the sequence is not unquie, decide if worth printing

     else if(moveOnBool == 1 && $4 <= minDiffDbl && $11 >= minMisDbl && bestPrintBool == 0)
     { # else if found a better line to print 
         keepSeqStr = $0;
         bestPrintBool = 1; # found the best print line
     } # else if found a better print line
    } # See if shoud keep the read
 } # else if on a data row
} # check if on a header or a sequence

#*******************************************************************************
# Sec-2 Sub-3: print out the last match if have one
#*******************************************************************************

END{{if(keepSeqStr != ""){print keepSeqStr}}}
