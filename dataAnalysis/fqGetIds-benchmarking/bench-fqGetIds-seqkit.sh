#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: Variable declerationS, check input, check output
#    sec-2: Get read counts
#    sec-3: Loop though each percentage (hard coded) & build filter file
#    sec-4: Run blank cases to ensure to create more stability
#    sec-5: Run time trials
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerationS, check input, check output
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

fastqStr="$1";
numRepInt="$2";    # number of replicates to do
readsInFastqInt=0; # holds number of reads in fastqFile
intRep=1;
numReadsInt=0;   # number reads extracted
maxThreads=3;    # max number of threads to test
statsFileStr="fqGetIds-seqkit-benchmark.tsv";
#sedCmdStr="p;n;n;n;n;n;"; # for Illumina (fqGetIds needs)
sedCmdStr="p;n;n;n;"; # for Nanopore

if [[ ! -f "$fastqStr" ]]; then
    printf \
        "No fastq provided (first argument)\n";
fi # if no ids provided

if [[ "$numRepInt" == "" ]]; then
    numRepInt=10;
fi # if user did not provide a number of replicates to do

if [[ ! -f "$statsFileStr" ]]; then
    { # merge printf commands
        printf "Program\tfastqFile\treadsInFastq\textractionSize";
        printf "\treadsExtracted\tReplicate\tThreads\tElapsedTime";
        printf "\tCPUKernTime\tCPUserTime\tMaxResidentMemoryInKb";
        printf "\tPercentCPU\n";
    } > "$statsFileStr";
fi # if need to create the stats file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get read counts
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

readsInFastqInt="$(
    awk '
        { # MAIN
            intReads++;         # header
            numSeqLinesInt = 0; # number of lines in sequence entry

            # Nanopore entries are one line, but Illumina likes to use two line
            # So I need a more flexible system
            while($0 !~ /^\+/)
            { # while are are sequence entries
                numSeqLinesInt++;
                getline;
            } # while are are sequence entries

            getline;            # move to fastq entry

            for(intQ = 0; intQ < numSeqLinesInt; intQ++)
                getline; # move to header entry
        }; # MAIN block

        END{print intReads;};
    ' < "$fastqStr" \
 )"; # get the number of reads in the fastq file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Loop though each percentage (hard coded) & build filter file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ -f "tmp-test-benchmark.filt" ]]; then
    rm "tmp-test-benchmark.filt"; # make sure awk has new temporary filke
fi # make sure not test file


for intPer in 1 10 100 250 350 500 600 750 850 1000; do
# For the tested percentages

    numReadsInt="$(
        awk \
            -v numIdsToKeepInt="$(((intPer * readsInFastqInt) /1000))" \
            -v totalReadsI="$readsInFastqInt" \
            '
	    BEGIN{
	        totalIdsI = 0; # Will hold number of read ids I read in
	        modValI = int(totalReadsI / numIdsToKeepInt);
	    };

            { # MAIN
                if(keptIdsInt > numIdsToKeepInt)
                    exit;

                totalIdsI++;  # count number of reads I have read in

                if(totalIdsI % modValI == 0)
                { # if is a read I want to keep
                    # clean up read id for seqkit (fastqGrep can handle this)
                    sub(/^@/, "", $0); # remove @ marking fastq header
                    sub(/ .*/, "", $0); # remove extra information in header

                    keptIdsInt++;  # keep track of how many ideas I kept
                    print $0 >> "tmp-test-benchmark.filt";  # print out the header
                } # if is a read I want to keep

                numSeqLinesInt = 0;
                getline;            # get off the header

                # Nanopore entries are one line, but Illumina has used
                # two lines for entries. So I need a more flexible system
                while($0 !~ /^\+/)
                { # while are are sequence entries
                    numSeqLinesInt++;
                    getline;
                } # while are are sequence entries

                for(intQ = 0; intQ < numSeqLinesInt; intQ++)
                    getline; # move to header entry
            }; # MAIN block

            END{print keptIdsInt}; # make sure user gets counts
          ' < "$fastqStr" \
    )"; # make filter file & get counts

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-4: Run blank cases to ensure to create more stability
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # Run a few through away cases, so measuring best time usage
    ../fqGetIds \
        -f "tmp-test-benchmark.filt" \
        -fastq "$fastqStr" \
      > "tmp--benchmark--test.fastq";

     seqkit \
       grep \
       -j 1 \
       -f "tmp-test-benchmark.filt" \
       "$fastqStr" \
     > "tmp--benchmark--test.fastq";

     seqkit \
       grep \
       -j 2 \
       -f "tmp-test-benchmark.filt" \
       "$fastqStr" \
     > "tmp--benchmark--test.fastq";

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-5: Run time for se kit trials
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    intRep=1;

    # Set up meta data for the stats file
    metaStr="$fastqStr	$readsInFastqInt	$numReadsInt";

    while [[ "$intRep" -le "$numRepInt" ]]; do
    # For all desiered replicates

        iThread=1;

        while [[ "$iThread" -le "$maxThreads" ]]; do
        # while have threads to test
           /usr/bin/time \
               -f "%e\t%S\t%U\t%M\t%P" \
               -o "tmp-time.tsv" \
             seqkit \
                grep \
                -f "tmp-test-benchmark.filt" \
                -j "$iThread" \
                "$fastqStr" \
             > "tmp--benchmark--test.fastq"

           timeStr="$(cat "tmp-time.tsv")";

           numReadsI="$(
               sed -n 'p;n;n;n;' "tmp--benchmark--test.fastq" | wc -l
           )"; # find number of reads extracted
               # seqkit converts the Illumina fastq 6 read entries into
               # 4 read entries, so do not need a sed command

           printf "seqkit\t%s\t%s\t%s\t%s\t%s\n" \
               "$metaStr" \
               "$numReadsI" \
               "$intRep" \
               "$iThread" \
               "$timeStr" \
             >> "$statsFileStr";
           
           iThread=$((iThread + 1));
        done # while have threads to test

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-6: Get time and stats for fqGetIds run with speed complie
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        /usr/bin/time \
            -f "%e\t%S\t%U\t%M\t%P" \
            -o "tmp-time.tsv" \
          ../fqGetIds \
            -f "tmp-test-benchmark.filt" \
            -fastq "$fastqStr" \
          > "tmp--benchmark--test.fastq";

        timeStr="$(cat "tmp-time.tsv")";
        numReadsI="$(
            sed -n "$sedCmdStr" "tmp--benchmark--test.fastq" | wc -l
        )"; # find number of reads extracted

        printf "fqGetIds\t%s\t%s\t%s\t%s\t%s\n" \
            "$metaStr" \
            "$numReadsI" \
            "$intRep" \
            "1" \
            "$timeStr" \
          >> "$statsFileStr";
       
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-7: Get time and stats for fqGetIds with memory compile
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        /usr/bin/time \
            -f "%e\t%S\t%U\t%M\t%P" \
            -o "tmp-time.tsv" \
          ../fqGetIdsMem \
              -f "tmp-test-benchmark.filt" \
              -fastq "$fastqStr" \
          > "tmp--benchmark--test.fastq";

        timeStr="$(cat "tmp-time.tsv")";
        numReadsI="$(
            sed -n "$sedCmdStr" "tmp--benchmark--test.fastq" | wc -l
        )"; # find number of reads extracted

        printf "fqGetIdsMem\t%s\t%s\t%s\t%s\t%s\n" \
            "$metaStr" \
            "$numReadsI" \
            "$intRep" \
            "1" \
            "$timeStr" \
          >> "$statsFileStr";

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-?: Clean up
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        intRep=$((intRep + 1));
    done

    rm "tmp-test-benchmark.filt"; # make sure awk has new temporary file
done # For the tested percentages

rm "tmp--benchmark--test.fastq";
rm "tmp-time.tsv";

exit;
