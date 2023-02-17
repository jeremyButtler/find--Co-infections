#!/bin/awk

################################################################################
# Name: getNumRefReads.awk
# Use:
#   Gets the number of reads simulated by badread for a single reference
# Input:
#   < file.fastq: Fastq file to get number of reference simulated reads from
#   -v refStr: Name of reference to get read counts for
# Output:
#   sdtout: prints number of reads for reference
################################################################################
{ # MAIN block
  if(NR % 4 != 1){next};   # not a header line
  #sub(/,.*/, "", $2); # Make sure only ref name in header
  if($2 ~ refStr){numRefInt++;}; # count number of reads from ref
} # MAIN block

END{
  if(numRefInt == "")
    print 0; # prevent blank prints
  else
    print numRefInt;
} # END block
