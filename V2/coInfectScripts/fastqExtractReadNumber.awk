#!/bin/awk

################################################################################
# Name: fastqExtractReadNumber.awk
# Use: Makes two fastq files, one with only the target read & one with every 
#      read except the target read
# Input:
#   -v readNumInt:
#     - The number off the read to extract (1st read if file, 2nd, ..., last?)
#     - Default: 0 (first read)
#     - Note: 0 index (so 1st read is 0, 2nd read is 1, 3rd read is 2, ...)
#   -v refFileStr
#     - What to name the file fram a single read
#   -v readFileStr
#     - What to name file with every read except the target read
# Output:
#   File: file with only the target read (refFileStr)
#   File: file with every read except the target read (readFileStr)
################################################################################

BEGIN{
    readNumInt *= 4; 
    readNumInt++
}; # BEGIN block

{ # main block
    if(NR == readNumInt)
        fileStr = refFileStr;      # file to print read to
    else
        fileStr = readFileStr;      # file to print read to

    # print the read to the target file
    print $0 > fileStr;
    getline;
    print $0 > fileStr;
    getline;
    print $0 > fileStr;
    getline;
    print $0 > fileStr;
} # main block
