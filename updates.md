# Use

Holds the change log of all updates I have have applied to find
  co-infections since March 2023.

The version number for find co-infections V3 now includes the date it
  was last updated. The pattern in the version number is 3.yearMonthDay.
  So, the version number 3.20230309 would mean that find co-infections
  was last updated on March 29th, 2023.

# Standing problems:

1. Racon (-enable-racon) errors out on occasion (planning on looking at
   after I finsh the alignment step for the majority consensus step).
2. fqGetIds and Illumina reads (see bottom of page, less likely to be
   resolved).

# Changes

## 2023-03-26

1. Debugged alignSeq, which currently is a Needleman-Wunsch alignment
   only.
    - See supplemental programs for information on and how to install.

## 2023-03-22:

1. This is mainly an updated for adding new figures and the results of
   my benchmarking with the dataset used to benchmark the ASHURE
   pipeline.
2. Changed the required 0.3% of reads to keep a cluster to 1%.
2. Also added in alignSeq.c and its supporting code in V3. This code
   is not complete, but is here as a back up.
   - Its goal is to add a Needleman Wunsch aligner that I can use in the
     majority consensus step. This step will use the best read to start
     a consensus and then align the next read to the consensus and add
     in all insertions. This re-alignment step will repeat until a
     consensus is built. I am hoping this will reduce the number of
     deletion errors the majority consensus step has.
   - Currently this aligner is overrunning the memory bounds on large
     alignments and needs some
     more work.
   - It also is having a problem flagging mismatches.
   - It only keeps two rows of scores (the previous scores and the 
     current row) and a matrix of two bit elements having on direction
     for each score (including discarded scores). This reduces the
     memory usage by a good amount but selects results in only a single
     path being selected (like the Hirschbergs (a divide and conquer
     Needleman Wunsch algorithm with even less memory usage). This could
     store all four possible directions (match, mismatch, indel,
     deletion) by using a matrix of four bit elements.

## 2023-03-15:

1. Fixed an issue with find co-infections V3 crashing on ASHURE dataset
   A when -skip-bin was combined with -pimers ASHURE-primers.fasta.
   - This was an issue with findXReads balance not being reset and
     freeing bigNum structs.
   - Fixed by changing the stack declarations of readInfo structures
     in findXBestReads (function used with extractTopReads) with heap
     declarations.

## 2023-03-14:

1. Changed version number to reflect the version number and the date
   of when it was last pushed to git hub.
   - Pattern is 3.yearMonthDay.
   - For example: 3.20230309 would mean that this version of find
     co-infections was last updated on March 29th, 2023.
2. Added in a primer trimming step to find co-infections
   - Trims any primer mappings in reads.
   - It will print out multiple entries when removing the primer splits
     the read. Each entry will be taged with "--D1_19-number", were
     the number is an integer value that represents the number of times
     this read was split.
3. Allowed my filtering of reads by length and mean/median quality
   score, primer trimming step, and selecting the top reads steps to 
   be compiled as separate programs (see supplemental-programs.md for
   more details).
   - These may not be very useful, but are here if you want them.
4. Fixed an issue were the best read could be added back in to the 
   fastq file twice (created an infinite loop).
5. Fixed an issue were the read to read mapping settings were being
   given to the read to consensus subsampling step. This did result in
   an infinite loop when the consensus could not bin all reads.
6. Fixed an issue with cluster counts over 9 being backwards (01 instead
   of 10).
7. -keep-unmapped-reads setting has been added to trimSamFile to allow
   the user to keep unmapped reads.
   - This is not very useful in general, but I had some issues with
     trimSam removing unmapped reads when I was trimming the reads in 
     the rolling circle replication dataset used to benchmark the
     ASHURE pipeline.
9. Found that the majority consensuses step truncated some consensuses
   by 200 to 400 base pairs in the ASHURE data set.
   - It looks like it was due to supplemental alignments being counted
     as separate reads. I was able to mostly get around this by trimming
     reads by primers and ignoring supplemental alignments.
10. Added in additional settings to set a minimum consensus length
   (-min-con-length), the minimum aligned length to keep reads during
   subsampling (-min-read-read-map-length & -min-read-con-map-length).
11. Binning parameters are no longer printed to the log when the binning
    step is skipped.

## 2023-03-05:

1. Fixed several minor errors in fqGetIds that would cause long
   overflows and also updated how it ran.
   - Their will be graphs comparing fqGetIds to seqkit in the
     supplemental programs document.
   - Fixed an issue were 4000 looked the same as 400 or 40 or 4 by
     change 0 to 16 in my big number. This makes each hex character
     (0-9 & a-f) take 5 bits instead of 4 bits, but also prevents
     pulling out non-targeted reads.
   - Changed the limb size from an long to an integer, so that the 
     sum of the big number could be stored in a long and used for
     hashing and comparisions.
2. Unresolved issues with fqGetIds:
   - fqGetIds does not know the difference between 186:235 and 18:6235.
     This could be fixed by setting ":" to 17, but will increase the 
     memory usage for Illumina data and may increase time.
   - fqGetIds may have a hard time handling Illumina fastq files that
     are merged together (one file containing multiple fastq files from
     different sequencers).
