# Use:

Finds co-infections in nanopore sequenced reads.

## Requirements:

  - Required
    1. GCC (to compile findCoInft).
    2. minimap2 (https://github.com/lh3/minimap2)
  - Optional:
    1. racon (https://github.com/isovic/racon)
    2. medaka (https://github.com/nanoporetech/medaka)
      - Can be installed by conda or python virtual enviroments
      - For python virtual env install; the folder with medaka needs to
        be located in the users home directory (~/medaka).

## Install:

cd V3;
make;
mv findCoInft /path/to/install;
chmod a+x /path/to/install/findCoInft;

## Run:

findCoInft -fastq reads.fastq -ref refferences.fasta [other options...]

## Some quick options:
  - -h:
    - Show the main help message (their are multiple)
  - -prefix:                                                     [Out]
    - Prefix to add to file names.
  - -threads:                                                    [3]
    - Number of threads to use with Minimap2/Racon/Medaka.
  - -min-reads-per-bin:                                          [100]
    - Minimum number of reads to keep a cluster or bin.
  - -max-reads-per-con:                                          [300]
    - Max number of of reads to use in building a consensus.
  - -enable-racon:                                               [No]
    - Use Racon in building consensuses.
  - -enable-medaka:                                              [No]
    - Use Medaka in building consensuses.
  - -model:                                         [r941_min_high_g351]
    - Model to use with Medaka_consensus

## Testing:

Currently benchmarking findCoInfections version three, and will update
  this section as data comes in. For version two testing results, see
  README.md in V2.

## How find co-infections version three works:

The first step in find co-infections involves binning the reads with the
  supplied references. Reads that have mapping qualities under 20 for
  the best reference, have multiple mappings (probably chimeras), did
  not map to any reference, have low mean or median Q-score (under 13),
  or are larger or smaller than the specified read length (600 to 1000)
  are removed. Any Read that passes the filtering step are then assigned
  to the bin representing the best reference they mapped to.

After binning comes the clustering step, were the best read for a bin
  is selected by mapping quality, integer median Q-score, and length.
  All reads in the bin are then mapped to this read and any reads that
  have a mapping quality under 10, more than a 2% difference in SNPs, or
  a 3% total difference are discarded. The top three hundred reads are
  then selected from the kept reads by median Q-score and length. The
  selected read and its top three hundred mapped reads are then used
  to build a consensus.

This consensus is then rebuilt with a new set of top three hundred
  reads. The top three hundred reads are found by mapping all reads in
  the bin to the consensus and then removing any read that has a mapping
  quality under 13, differs from the consensus by 1.5% SNPs, or is 2%
  different than the consensus. The best three hundred reads are then 
  selected by mapping quality, median Q-score, and read length. The
  consensus and top three hundred reads are then used to build yet
  another consensus. This rebuilding of the consensus step with a new
  set of reads is repeated twice.

Next the reads in the bin are mapped to the current consensus to see
  which reads belong in the consensus's cluster. False positive clusters
  are removed from each bin by comparing the consensuses in each bin to
  each other. Clusters that have similar consensuses (under 1%
  difference in SNPs or 3% in total difference) are considered the same
  and the cluster with the fewest reads is merged into the cluster with
  the most reads (the consensus is deleted for the merged cluster).

After the clustering step, consensuses between bins are compared with
  the same method as at the end of the clustering step. Clusters that
  are found to have similar consensuses are again merged with the same
  method.

![Diagram of how version three works, for asci flow chart see 
  diagrams/v3FlowDiagram.ditaa](diagrams/v3FlowIdagram.svg)

The default consensus building step builds a majority consensus using
  the top read or consensus and the top three hundred reads. In the
  first part a linked list (sequence) with the references bases is made,
  with bases from reads having quality scores under ten being set to
  nothing. The top read or consensus is then mapped to the top three
  hundred reads. 

Each mapped read is used to get a count of the number of matches, SNPs,
  and insertions for each position. Each match or SNP needs a Q-score of
  at lest ten to be counted, while insertions need a Q-score of five.
  Any new SNPs are inserted into the linked list as an alternative base,
  while new insertions are inserted as a new base. However, an
  alternative base is used if their is an insertion already at the
  position.

```
                 |a 10 reads|   <- Alternative bases
                     | 
   |t 1 read |   |c  2 reads|                  |insertion t 2 reads|
       |             |                                  |
   |a 2 reads|---|t  1 reads|---|a 10 reads|---|insertion a 1 read |
Position: 1          2              3                 ?
Total: 3 reads    13 reads         10 reads         4 reads

 An example of what the sequence linked list might look like. The first
   base would be discarded, because 3/15 = 20% support, which is less
   than the required 40% support for non-insertions. The insertion would
   also be discarded, because (2 + 1)/15 = 20% support, which is less
   than the required 30% support for insertions. Their is no alignment
   step for insertions.
```

Insertion errors in the sequences are found by finding the total number
  of reads supporting a position (including alternative bases).
  Non-insertion positions that have less than 40% support are removed 
  and insertion positions that have less than 30% are removed. The base
  with the most supporting reads is kept, while the other alternative
  bases are removed. Leaving the final consensus.

## OS's/platforms tested on

  - OpenBSD (no Medaka)
  - Debian
  - Raspberry pi 1B with a 32 bit Debian install (no Medaka)
    - The independent install of fqGetIds did crash on this, but
      findCoInft ran fine.
    - I should fix this issue at some point.

## Possible plans:

1. Add in full consensus building step to version three and maybe some
   masking for poorly supported bases. The majority consensus step will
   remove these bases, so I am not to worried about this.
2. Put consensus building steps in separate file, so can make it into a
   separate C program.
3. Maybe multi-thread fqGetIds (in V2 was fastqGrep). Likely will not
   happen.
