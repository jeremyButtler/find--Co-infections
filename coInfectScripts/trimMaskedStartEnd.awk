#!/bin/awk

#*******************************************************************************
# TOC:
#    fun-1: trimAndPrintReads: trim read down using cigar & print to stdout
#    begin: Set variables
#    main: find best alignment for read & trim using trimAndPrintRead (fun-1)
#    end: trim & print the last read with trimAndPrintRead (fun-1)
#*******************************************************************************

#*******************************************************************************
# Name: trimMaskedStartEnd.awk
# Use: Trims the masked staring bases and ending bases from a sequence. The 
#      input file should be a sam file in eqx format (minimap2 --eqx -a).
# Input:
#    < file.sam: sam file with sequences to trim and assing to ref
#    -v prefStr: Prefix to add to output files
# Output:
#    files: Each read is assigned to reference (prefStr--refence.fastq)
#*******************************************************************************

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fun-1: trimAndPrintReads: trim read down using cigar & print to stdout
# Use:
#   a cigar to trim soft masks off an end of an read
# Output:
#   Prints trimmed read to stdout
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function trimAndPrintRead(readIdStr, seqStr, qStr, cigStr, refIdStr, prefStr)
{ # trimAndPrintRead
  # make out file name

  if(qStr != "*")
      fileStr = prefStr "--" refIdStr ".fastq";
  else
      fileStr = prefStr "--" refIdStr ".fasta";

  #Split cigar and get the sequence lengths
  lenCigarInt = split(cigStr, cigarAryStr, "="); # get length of cigar
  lenSeqInt = length(seqStr); # get the length of the sequence

  if(cigarAryStr[lenCigarInt] ~ /S$/) # S=soft (H=hard, minimap2 did trim)
  { # If need to trim the end bases off
      # remove the notation for the masking, leaving the number
      sub(/S/, "", cigarAryStr[lenCigarInt]);

      # trim my sequence
      trimInt = lenSeqInt - cigarAryStr[lenCigarInt];
      seqStr = substr(seqStr, 1, trimInt);     # fastq line

      if(qStr != "*")
        qStr = substr(qStr, 1, trimInt);     # Q-score line
  } # If need to trim the end bases off

  if(cigarAryStr[1] ~ /S/) # S=soft trim (H=hard, minimap2 already did trim)
  { # If need to trim the starting bases off
      # remove the notation for the masking, leaving the number
      sub(/S.*/, "", cigarAryStr[1]); # .* removes correct bases

      # trim off the starting bases that are masked
      seqStr = substr(seqStr, cigarAryStr[1] + 1);    # sequence line

      if(qStr != "*")
        qStr = substr(qStr, cigarAryStr[1] + 1);    # Q-score line
  } # If need to trim the starting bases off

  if(seqStr == ""){return;} # nothing to print out

  if(qStr != "*")                   # is a fastq file
    printf "@%s\n%s\n+\n%s\n", readIdStr, seqStr, qStr >> fileStr;
  else                              # likely from fasta file
    printf ">%s\n%s\n", readIdStr, seqStr >> fileStr;
  return;
} # trimAndPrintRead

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BEGIN: set variables
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{ # BEGIN block
   FS="\t";
   if(prefStr == ""){prefStr = "out"};
}; # BEGIN block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main: find best alignment for read & trim using trimAndPrintRead (fun-1)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # awk main block
  if($1 ~ /^@/){next;}                    # Ignore headers
  if($3 == "*"){next;}                    # read did not map

  if($1 == readIdStr)
  { # if still on the old read, determine if better aligment to use

      if($5 > mapqInt)     # unsigned mapqs are 0 (best is always top)
      { # if this line has a better mapq
        mapqInt = $5;  # mapping quality
        cigStr = $6;   # cigar line
        refIdStr = $3;   # get name of reference mapped to
      }; # if this line has a better mapq
 
      if($10 != "*")
      { # if have a sequence & q-score line, but on an old read
        seqStr = $10;
        qStr = $11;
      }; # if have a sequence & q-score line, but on an old read

      next; # still have aligments to check and see if better
  } # if still on the old read, determine if better aligment to use

  if(readIdStr != 0)     # if not the first round
      trimAndPrintRead(readIdStr, seqStr, qStr, cigStr, refIdStr, prefStr); # using best alignment

  # Get the new reads information
  readIdStr = $1;
  refIdStr = $3;
  mapqInt = $5;  # mapping quality
  cigStr = $6;   # cigar line
  seqStr = $10;
  qStr = $11;
}; # awk main block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# END: trim & print the last read with trimAndPrintRead (fun-1)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

END{trimAndPrintRead(readIdStr, seqStr, qStr, cigStr, refIdStr, prefStr);};
