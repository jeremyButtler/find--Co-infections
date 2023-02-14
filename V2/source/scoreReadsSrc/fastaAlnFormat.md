# What is the fastqAln format: #

The fastqAln format is an binary fastq file that has been aligned to at least
  one reference. Each sequence in the fastqAln format stores the nucleotide of
  the reference sequence used, the nucleotide of the original sequence, if the
  base should be masked, if the base is not aligned to the reference
  (soft masked), and if the base was a match, mismatch, or indel. In addition
  this format also keeps the quality entry of the fastq file. So, the size of a
  fastq aligned file is the same as a fastq file, but the amount of information
  stored is greater.

This file format allows me to speed up my clustering step, by allowing me to 
  call minimap2 only once to align reads to a reference genome. I can then mask
  out low quality bases, but still keep the original nucleotide. The not soft 
  clipped regions of these reads can the be compared in an all versus all
  alignment without having to call minimap2 multiple times. I can even trim the
  reads before hand to reduce the amount of time ignoring soft masked regions.
  Finally, I can convert a fastqAln file back to a fastq file.

# fastqAln format: #

## Header entry: ##

The header starts with a number with the number of characters in the header (
  this includes deliminators), followed by a tab and the read id. Everything
  after the read id is space deliminated. The other information in the header 
  is the mapped reference, flag marking if reference bases used (1 yes, 0 no),
  the sam file flag, and mapping quality. Other entries can be added after and
  should be separated with spaces. The header ends with an tab.

Header: 15	07e KY800362 1 0 30   
Header datatypes: unsigned-long char[] char[] char unsinged-int unsigned-int    

## Sequence entry: ##

The sequence entry starts with the number of characters in the sequence entry
  (this includes the tab at the start and end) and a tab. Next is a character
  for each base in the sequence followed by a tab to end the sequence entry.
  Each character in the sequence entry has eight possible flags.

WARNING: This format is unable to record anonymous bases. So, if anonymous
  bases are present in the reference, there will be errors. To reduce this
  bases that map to anonymous bases will be masked.

| Flag |              0                 |                1                 |
|:----:|:------------------------------:|:--------------------------------:|
|  1   | Reference base is a transition | Reference base is a transversion |
|  2   |        Is not an mismatch^1^   |          Is an mismatch^1^       |
|  4   |        Is not an deletion      |          Is an deletion          |
|  8   |        Is not an insertion     |          Is an insertion         |
|  16  |          Is not masked         |            Is masked             |
|  32  |          Not soft masked       |            Is soft masked        |
|  64  |          Purine                |            Pyrimidine            |
|  128 |          A or G                |            C or T                |

Table: fastqAln format. ^1^ will also be flagged if for deletions to tell if
  reference base was an Purine or Pyrimidine. For deletions. 0 is Purine, 1 is
  Pyrimidine. Also, flag 1 changes from transition  and transervsion to 
  A and G or T and C for deletions. Flag 8 is set to insertions to ensure that
  newline (10), carriage return (13), and tab (9) are never used. For insertions
  (8 set to 1) the deletion, mismatch, and reference base flags are set 0.

Sequence entry: 20	18 characters here	
Sequence datatypes: unsigned-long	char[]  

Note: Reverse complement reads are reverse complement before printing (so
  printed in the reference forward direction). This is so all sequences can
  be compared easily (no worrys about if sequences are reverse complements). Pay
  attention to the sam file flag if you care about reverse complement bases.

## Quality score entry: ##

The quality score (Q-score) entry starts with a number telling the number of
  characters in the q-score line (includes deliminiators) followed by a tab. The
  q-score entry comes next as a series of characters and ends in a tab. The only
  difference from the fastq q-socre entry is that deletions are marked with a
  space (fastq file does not show deletions).

Q-score entry: 20	18 q-scores here12	
Q-score datatypes: unsigned-long	char[]  

# Programs I used to process fastqAln formats: #

You likely have an interest in this format if you have made it this far. So here
  is a list of programs I am working on or have finished for version three
  of findCoInfections that work with fastqAln files.

  1. samToFastqAln:
    - Converts a sam file to fastqAln file.
    - Requires an eqx cigar entry (23X10=10I1=) (minimap2 --eqx).
    - Only keeps mapped alignments with sequences. So, only one entry will be 
      kept for a multiple reference alignment.
    - Masks bases that map to anonymous bases in the reference. 
    - Recommended prep: minimap2 --eqx --MD -ax setting ref.fasta reads.fastq.
    - Using a fasta file will result in no Q-score entry "2\t\t".
  2. maskFastqAln:
    - Mask low quality bases in a fastqAln file
    - Low quality bases are determined by q-score and/or reference homopolymer
      size for indels (deletions do not have quality scores).
  3. scoreFastqAln:
    - Puts hamming distance after mapping quality in each sequence header and
      a score (generally square root, but may make flexible)
    - Ignores masked bases.
    - Also puts: number of kept mismatches, indels, matches, mean Q-score, 
      median Q-score, mean aligned Q-score, and medain aligned Q-score in the
      header.
  4. fastqAlnFilter:
    - Removes reads beneath input min thresholds from scores found using
      scoreFastqAln.
  5. trimFastqAln:
    - Removes softmasking at start and end of reads.
  6. binFastqAlnByRef:
    - Splits fastqAln file into several smaller fastqAln files, named after
      the reference a read mapped to.
  7. fastqAlnToFastq:
    - Converts fastqAln file to fastq file.
    - Allows masking to be applied as N's (default is no).
  8. allVesusAllHam:
    - Gets percent hamming distance by comparing all reads mapping to the same
      reference in a fastqAln file.
    - NOTE TO SELF: MAKE THIS HANDLE MULTI REF FILES
  9. fastqAlnGrep:
    - Extract reads by read id in a fastqAln file.
  10. fastqAlnToConsensus:
      - Builds a simple majority consensus from a fastqAln file with only one
        reference. Any maksed reads are ignored in the consensus building.
      - The hope is that a majority consensus will provide a higher quality 
        "read" to racon + medaka. This hopefully will result in more accurate
         results.
  11. fastqAlnDetectFalseRec:
    - Compares a consensus to the reads it was built from to see if the
      consensus is a merging of two different variants.
    - The idea is most false positives in this pipeline are do to misbinned
      reads being merged together.
    - The consensus and reads must be in fastqAln format and be aligned to the
      same reference.
