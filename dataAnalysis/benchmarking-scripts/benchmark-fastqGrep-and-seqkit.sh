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

if [[ ! -f "$fastqStr" ]]; then
    printf \
        "No fastq provided (first argument)\n";
fi # if no ids provided

if [[ "$numRepInt" == "" ]]; then
    numRepInt=10;
fi # if user did not provide a number of replicates to do

if [[ ! -f "fastqExtractNoGmpBenchmark.tsv" ]]; then
    { # merge printf commands
        printf \
          "Program\tfastqFile\treadsInFastq\treadsExtracted\tReplicate"
        printf \
          "\tThreads\tElapsedTime\tCPUKernTime\tCPUserTime";
         printf \
          "\tMaxResidentMemoryInKb\tPercentCPU\n";
    } > "fastqExtractBigNumBenchmark.tsv";
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


for intPer in 1 10 25 35 50 60 75 85 100; do
# For the tested percentages

    numReadsInt="$(
        awk \
            -v numIdsToKeepInt="$(((intPer * readsInFastqInt) / 100))" \
            '
            { # MAIN
                if(keptIdsInt > numIdsToKeepInt)
                    exit;

                keptIdsInt++;  # keep track of how many ideas I kept

                # clean up read id for seqkit (fastqGrep can handle this)
                sub(/^@/, "", $0); # remove @ marking fastq header
                sub(/ .*/, "", $0); # remove extra information in header
                print $0 >> "tmp-test-benchmark.filt";  # print out the header

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
    ./fastqGrep \
        -f "tmp-test-benchmark.filt" \
        -fastq "$fastqStr" \
      > "tmp--benchmark--test.fastq";

     seqkit \
       grep \
       -f "tmp-test-benchmark.filt" \
       -j 1 \
       "$fastqStr" \
     > "tmp--benchmark--test.fastq";

     seqkit \
       grep \
       -f "tmp-test-benchmark.filt" \
       -j 2 \
       "$fastqStr" \
     > "tmp--benchmark--test.fastq";

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-5: Run time trials
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    while [[ "$intRep" -le "$numRepInt" ]]; do
    # For all desiered replicates
        /usr/bin/time \
            -f "seqKit\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t1\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
          seqkit \
            grep \
            -f "tmp-test-benchmark.filt" \
            -j 1 \
            "$fastqStr" \
          > "tmp--benchmark--test.fastq";
        
        /usr/bin/time \
            -f "seqKit\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t2\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
          seqkit \
            grep \
            -f "tmp-test-benchmark.filt" \
            -j 2 \
            "$fastqStr" \
          > "tmp--benchmark--test.fastq";
        
        /usr/bin/time \
            -f "seqKit\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t4\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
          seqkit \
            grep \
            -f "tmp-test-benchmark.filt" \
            -j 4 \
            "$fastqStr" \
          > "tmp--benchmark--test.fastq";
        
        /usr/bin/time \
            -f "seqKit\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t7\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
          seqkit \
            grep \
            -f "tmp-test-benchmark.filt" \
            -j 7 \
            "$fastqStr" \
          > "tmp--benchmark--test.fastq";
        
        /usr/bin/time \
            -f "fastqGrep-hash\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t1\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
        ./fastqGrep \
            -f "tmp-test-benchmark.filt" \
            -fastq "$fastqStr" \
          > "tmp--benchmark--test.fastq";
        
        /usr/bin/time \
            -f "fastqGrep-tree\t$fastqStr\t$readsInFastqInt\t$numReadsInt\t$intRep\t1\t%e\t%S\t%U\t%M\t%P" \
            -o "fastqExtractBigNumBenchmark.tsv" \
            -a \
        ./fastqGrep \
            -f "tmp-test-benchmark.filt" \
            -fastq "$fastqStr" \
            -no-hash \
          > "tmp--benchmark--test.fastq";
    
        intRep=$((intRep + 1));
    done # loop through all replicates

    rm "tmp-test-benchmark.filt"; # make sure awk has new temporary file
done # For the tested percentages

rm "tmp--benchmark--test.fastq";
exit;
