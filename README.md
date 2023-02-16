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
    - Show the main help message.
    - Their are enough parameters in find co-infections to 
      make a single help message daunting. To get around this
      I split the parameters into multiple help messages by
      catagory. The main help message shows what I think are
      the most useful parameters and all commands to show
      the other help messages.
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

## Advice on fine tuning this find co-infections:

Find co-infections allows you to change most of the parameters at every
  step, which means that their are a lot of parameters that can be
  changed. My advice is to get several good data sets that you know the
  answers to and pick a couple key parameters to fiddle around with.

Another tip of advice, make sure that some of your tests do not have
  references that match the sequences in the sample. That way you can
  test a realistic case.

## Testing:

Currently benchmarking findCoInfections version three, and will update
  this section as data comes in. For version two testing results, see
  README.md in V2. One note I will make about version two results is
  that I caught 18 cases of the same pair of references being run twice
  in all reads, version one, version two, and version two with
  clustering.  I think these came in because my computer shut down
  during testing at one point and though I caught most of the files at
  the time, but missed the 18 files that got duplicates from.

### Brief Methods:

To benchmark find co-infections we simulated reads for pairs of
  references using badread on default settings (for a list of reads
  see tbl:interGenotype in tables.md). The similarities between the two
  references in each reference pair ranged between, 97% similar
  (97.6% to 96.9%), 98% similar (97.7% to 98.2%), and 95%
  (98.6% to 98.9%). The percent similarity between consensuses was found
  by mapping the references to each other with minimap2 -a
  (output sam file) and then dividing the total error (NM:i tag) by the
  length of the sequence. 

Some of the reads in our reference pairs existed in our database, while
  others did not. To reduce the impact of exact matches in our database
  influencing our result, we ensure that only ten of thirty reference
  pairs had both references in the database for reference pairs that
  were at or under 98% similar (98.5% always has one reference out)
  (see tbl:interGenotype in tables.md). We also had ten reference
  pairs that were at least 1% different than any references in the
  database (see tbl:interGenotype in tables.md). For the final ten
  reference pairs we ensured that one reference was in the database,
  while one reference was at least 1% different from any reference in
  the database (see tbl:interGenotype in tables.md).

For all reference pairs we simulated 50% reads for reference pairs or
  80% of reads for the major strain and 20% of reads for the minor
  strain. We also simulated an additional set of reads that reversed the
  major and minor strain for each reference pair that had only one
  reference in the database. This resulted in 30 sets of simulated reads
  when the reference pairs were 98.5% similar and 40 sets of simulated
  reads when the reference pairs were 98% or 97% similar.

### Results:

![Figure 1: Number of references correctly made when 50% of reads came
  from the major strain. The horizontal line shows the maximum number of
  consensuses per test](figures/Similarity--20000--v3--50-50--SNPs.svg)

![Figure 2: Number of references correctly made when 80% of reads came
  from the major strain. The horizontal line shows the maximum number of
  consensuses per test](figures/Similarity--20000--v3--80-20--SNPs.svg)

Figure 1 and 2 are showing the number of correct identified consensuses
  and extra (false positives) consensuses when the major strain
  contributed 50% of reads (Figure 1) and 80% of reads (figure 2). In
  all cases find co-infections V3 is detecting less false positives and
  close to if not more correct consensuses than Longshot, Clair 3, or
  the previous versions of find co-infections. Find co-infections V3
  also has less false positive SNPs than any other tested pipeline.

For indels we found that find co-infections V3 had more indel free
  consensuses than any other tested pipeline
  (figures/Similarity--20000--v3--80-20--indels.svg). However, V3 also
  has a few consensuses with higher indel counts than Clair 3, Longshot,
  and the previous versions of find co-infections, so their are times V3
  will have worse indel accuracy then other programs
  (figures/Similarity--20000--v3--80-20--indels.svg).

### Take away:

We have shown that find co-infections V3 is better than find
  co-infections V2 at detecting co-infections in simulated datasets.
  It is hard to say if this improvement is due to a better ability to
  handle reads simulated by Badread or is an actual improvement.

## Weaknesses:

Find co-infections has only been tested using data simulated from
  badread (Wick 2019) on default settings, which does produce low
  quality reads. However, these simulated reads can be used to build
  more accurate consensuses than consensuses built with real reads of
  the same quality
  (https://github.com/rrwick/Badread/tree/main/comparison).
  This could pose a problem for find co-infections, since find
  co-infections relies on accurate SNP calls in consensuses to identify
  co-infections and to hone in on a good set of reads for building the
  consensuses used in the consensus comparison steps.

Find co-infections has only been tested on data sets that contain two
  references per test. This means that I am unsure how find
  co-infections will act when their are more than two viral strains in a
  sample or when their are less than two viral strains in a sample.
  However, find co-infections version three did not detect any false 
  positives in the similarity data set, so it will probably not detect
  false positives consensuses when their is only one viral strain in a
  sample.

Find co-infections is designed to run on data sets that have very deep
  read depths, such as the read depths produced by PCR.

I need to test a somatic variant caller to round out the variant caller
  selection. Not likely to be done, but is something I missed.

```
Wick RR. Badread: simulation of error-prone long reads.
  Journal of Open Source Software. 2019;4(36):1316.
  doi:10.21105/joss.01316.
```

## How find co-infections version three works:

### An overview of the general steps

The first step in find co-infections involves binning the reads with the
  supplied references. Reads that have mapping qualities under 20 for
  the best reference,  did not map to a reference, have low mean or
  median Q-score (under 13), or are larger or smaller than the specified
  read length (600 to 1000) are removed. Also supplementary and
  secondary alignments for any read are ignored. Any Read that passes
  the filtering step are then assigned to the bin representing the best
  reference they mapped to.

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
  diagrams/v3FlowDiagram.ditaa](diagrams/v3FlowDiagram.png)

### The consensus building steps

Find co-infections can build consensuses using an inbuilt consensus
  builder, with Racon, or Medaka. It can also polish the consensus from
  one consensus building method with another method. The first
  consensus, is built using the inbuilt consensus builder (this can be
  turned off). This consensus can be polished with four rounds of
  Racon or Medaka_consensus, when Racon is not used. The consensus built
  or polished by Racon is then polished by Medaka_consensus if Medaka is
  being used. I do not think Medaka can use a fastq file as a reference,
  so you should plan on at least using Racon or the inbuilt consensus
  builder with Medaka.

The default and inbuilt consensus building step builds a majority
  consensus using the top read or consensus and the top three hundred
  reads. It using the Q-scores assigned to each base to determine if a
  match, SNP, or insertion should be kept. Deletions are not recorded,
  but can be found at the end based on the percentage of read support
  each position has.

The first step is to convert the reference sequence into a linked list
  (sequence) that can support alternative bases (SNPs) and insertions.
  Each position is recorded, but only bases that have quality scores
  over ten are counted for the read support. For consensuses, all bases
  are considered to have a quality score of at least ten.

After building the linked list the reads are mapped to the reference
  read or consensus for scoring. Each SNP or match that has a quality 
  score of ten or more are counted, while each insertion that has a 
  quality score of five or more is counted. For new SNPs at a position,
  an new alternative base is made, while for new insertions a new base
  is added to the sequences linked list. Alternative bases are add to
  insertions when the position supporting position already has an
  insertion. Multiple insertions are added to the list when their is
  more than one insertion in a read at a single position.

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

Insertion errors in the sequences are found in the linked sequence list
  by finding the total number of reads supporting a position. This
  includes alternative bases. Non-insertion positions that occurred in
  less than 40% of reads are removed and insertion positions that
  occurred in less than 30% of reads are removed. The base with the most
  supporting reads is kept, while the other alternative bases are
  removed. Once completed, we are left with the final consensus.

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
4. Probably should document how to interface with the small programs I
   made (trimSamFile, scoreReads, and fqGetIds). May never be done.

## Thanks:

Their are many people who made this project possible.

  - To my family for their support as a debugged. In particular my dad,
    Jeff Buttler, who sometimes made suggestions for how to improve my
    code. He introduced me to using look up tables in my code.
  - Danielle Wrenn: Who suggested using sam flags to remove
    supplementary alignments, which allowed multi threading of minimap2
    for the binning step.
  - Devin Drown: Who was one of my mentors through graduate school and
    is still a mentor. He helped in the development my graphs, helped in
    the benchmarking for my pipeline, and gave me advice.
  - Matthew Redlinger: Who is a skilled bioinformation and was always
    willing to give advise.
  - Eric Bortz: Who was one of my mentors in graduate school and still
    is one of my mentors. He and Matt gave me the encouragement to start
    making find co-infections.
  - The entire Drown lab at the University of Alaska Fairbanks, who
    viewed and gave feedback on early reports of find co-infections.
  - The entire Bortz lab at the University of Alaska Anchorage, who also
    viewed and gave feedback on my presentations.
  - The EMME group at the University of Alaska Fairbanks and the
    University of Alaska Anchorage. For being a group I could present
    my results to.
  - To Molly Murphy: Who was one my graduate committee and provided
    valuable feedback on my work during my defense and thesis writing 
    stages.
  - To Naoki Takebayasi: Who was one my graduate committee and provided
    valuable feedback on my work during my defense and thesis writing
    stages.
  - To Nataliia Rudova and everyone else who set up the experimental 
    design and did the sequencing for the first chapter of my thesis
    (See doi.org/10.3390/v14050924 for a full list). In this project I
    was introduced to co-infections and the problem of detecting
    co-infections between different vial strains in noisy reads. This
    was the reason I made version one of find co-infections.
  - To Leigh et al., who I never have met or interacted with, but who
    published a paper (http://dx.doi.org/10.3390/v12080801) that set up
    the ground work and gave me the idea of how to make my first version
    of find co-infections.
  - To Ryan Wick, who I again never met or interacted with, but made
    the read simulation software I used to test find co-infections with. 
  - Finally the Universities of Fairbanks and Anchorage Alaska, who 
    provided me with TA ships while I was in graduate school. During
    this time I wrote and tested version one and two of find
    co-infections.
