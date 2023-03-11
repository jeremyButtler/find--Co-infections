# Use

Holds the change log of all updates I have have applied to find
  co-infections since March 9th, 2023.

The version number for find co-infections V3 now includes the date it
  was last updated. The pattern in the version number is 3.yearMonthDay.
  So, the version number 3.20230309 would mean that find co-infections
  was last updated on March 29th, 2023.

# Changes

## 2023-03-10:

1. Changed version number to reflect the version number and the date
   of when it was last pushed to git hub.
   - Pattern is 3.yearMonthDay.
   - For example: 3.20230309 would mean that this version of find
     co-infections was last updated on March 29th, 2023.
2. Fixed an issue were the best read could be added back in to the 
   fastq file twice (created an infinite loop).
3. Fixed an issue were the read to read mapping settings were being
   given to the read to consensus subsampling step. This did result in
   an infinite loop when the consensus could not bin all reads.
4. Fixed an issue with cluster counts over 9 being backwards (01 instead
   of 10).
5. -keep-unmapped-reads setting has been added to trimSamFile to allow
   the user to keep unmapped reads.
   - This is not very useful in general, but I had some issues with
     trimSam removing unmapped reads when I was trimming the reads in 
     the rolling circle replication dataset used to benchmark the
     ASHURE pipeline (this data set should give me a good idea of how
     accuratly I can build a consensus).
6. Found that the majority consensuses steps did truncate some
   consensuses by 200 to 400 base pairs in the ASHURE data set.
   - It looks like in this case it was due to supplemental alignments
     being counted as separate reads. I was able to mostly get around
     for trimmed datasets this by ignoring all supplemental alignments.
     Their still some truncation, but it is know under 80 base pairs.
   - This problem has not been fixed for running the untrimmed reads 
     from the ASHURE dataset through find co-infections. I am working to
     add in a primer trimming step to reduce this.
7. Added in additional settings to set a minimum consensus length
   (-min-con-length), the minimum aligned length to keep reads during
   subsampling (-min-read-read-map-length & -min-read-con-map-length).
   - These settings are currently not being output to the log. I will
     try to get around to this after I have made some more progress.

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

