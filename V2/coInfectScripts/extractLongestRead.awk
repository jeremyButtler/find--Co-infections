################################################################################
# Name: extractLongestRead.awk
# Use: Extracts the longest read from a fastq file
# Input:
#    < file.fastq: fasatq file with reads to separate out
#    -v filePrefixStr: prefix to name the output files
# Output:
#    File: fasta with longest read (filePrefix--longest-read.fasta)
#    File: fastq with the remaing reads (filePrefix--other-reads.fastq)
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: BEGIN block
#    sec-2: Extrac the full fastq entry
#    sec-3: If longer than current longest read, make longst & print old read
#    sec-4: Else is not the longest read, print out to other read file
#    sec-5: Print out longest read to longest read file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: BEGIN block
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
    if(filePrefixStr == ""){filePrefixStr = "out";};

    otherFile = filePrefixStr "--other-reads.fastq";
    longFile = filePrefixStr "--longest-read.fasta";
} # BEGIN block

{ # main block

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-2: Extrac the full fastq entry
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    headerStr = $0;
    getline seqStr;
    getline spacerStr;
    getline qStr;

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-3: If longest read, print out old longest read to other file & recored
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    lenSeqInt = length(seqStr);

    if(lenSeqInt > longestReadInt)
    { # if is the current longest read
        longestReadInt = lenSeqInt;

        if(longestHeaderStr != "")
        { # If not the first round (no longest read yet)
            printf "%s\n%s\n%s\n%s\n",
                longestHeaderStr,
                longestSeqStr,
                longestSpacerStr,
                longestQStr > otherFile;
        } # If not the first round (no longest read yet)

        # record the longest entry
        longestHeaderStr = headerStr;
        longestSeqStr = seqStr;
        longestSpacerStr = spacerStr;
        longestQStr = qStr;
    } # if is the current longest read

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-4: Else is not the longest read, print out to other read file
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    else
    { # else is not the longest read
        printf "%s\n%s\n%s\n%s\n",
            headerStr,
            seqStr,
            spacerStr,
            qStr > otherFile;
    } # else is not the longest read
} # main block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Print out longest read to longest read file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

END{
    sub(/^@/, ">", longestHeaderStr);
    printf "%s\n%s",
        longestHeaderStr,
        longestSeqStr > longFile;
} # END block
