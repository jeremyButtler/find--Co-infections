# Use:

Finds co-infections in nanopore sequenced reads.

I updated the license to CCO because it has the same terms
  and conditions as the unlicense (none), but is more
  acknowledged than the unlicense. 

## Notes:

The install for fqGetIds was removed, due seqkit version
  2.x being faster and using less memory than fqGetIds.
  fqGetIds was benchmarked against seqkit version 0.15,
  which did not have the speed or memory upgrades of
  version 2.

## Requirements:

  - Required
    1. GCC (to compile findCoInft).
    2. minimap2 (https://github.com/lh3/minimap2)
  - Optional:
    1. racon (https://github.com/lbcb-sci/racon)
    2. medaka (https://github.com/nanoporetech/medaka)
      - Can be installed by conda or python virtual enviroments
      - For python virtual env install; the folder with medaka needs to
        be located in the users home directory (~/medaka).

## Install:

```
cd V3;
make;
mv findCoInft /path/to/install;
chmod a+x /path/to/install/findCoInft;
```

## Run:

findCoInft -fastq reads.fastq -ref refferences.fasta [other options...]

## Some quick options:
  - -h:
    - Show the main help message.
    - Also shows commands for additional help messages.
  - -prefix:                                                     [Out]
    - Prefix to add to file names.
  - -threads:                                                    [3]
    - Number of threads to use with Minimap2/Racon/Medaka.
  - -min-reads-per-bin:                                          [100]
    - Minimum number of reads to keep a cluster or bin.
  - -max-reads-per-con:                                          [300]
    - Max number of of reads to use in building a consensus.
  - -extra-consensus-steps:                                      [2]
      - Number of times to rebuild the consensus using a
        new set of best reads.
  - -min-read-length:                                           [600]
    - Discard reads with read lengths under this setting
      during the binning step.
    - Minimum read length to keep a read when binning.
    - Length is compared after the trimming step.
  - -max-read-length:                                          [1000]
    - Discard reads with read lengths over this setting
      during the binning step.
    - Length is compared after the trimming step.
  - -min-median-q:                                             [10]
    - Minimum read median quality score needed to keep a
      read in inital filtering steps.
  - -min-mean-q:                                               [10]
    - Minimum read mean quality score needed to keep a read
      in inital filtering steps.
  - -min-per-reads:                                        [0.01=1%]
    - Minimum percentage of total clustered reads needed
      to keep a cluster.
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

For more detail on my methods (still brief) see dataAnaylsis/README.md.
  For the similarity test I simulated reads with Badread using --seed 
  1026 and --quantity 20000x for references. The reference pairs used
  to simulate reads can be found in the
  dataAnalysis/PCV2--similarity--reference-pairs/ directory.

The number of false positive SNPs (SNP errors) and false positive indels
  (indel errors) between the built consensuses and their reference was
  found by mapping the consensus to the reference pair use to simulate
  the reads used to build them with minimap2 --eqx -a (output sam file).
  The number of SNP errors and indel errors was then found with an awk
  script (dataAnalysis/benchmarking-scirpts/getMisIndelCnt.awk). This
  awk script counted the number of mismatches in the cigar entry and
  found the number of indel errors by subtracted the total errors
  (NM:i tag) from the number of SNP errors.

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
  consensuses for simulated reads than any other tested pipeline
  (figures/Similarity--20000--v3--80-20--indels.svg). However, when I
  run find co-infections V3 on real porcine circovirus type 2 reads
  (From doi.org/10.3390/v14050924), I noticed that most consensuses had
  at least one to three indels to their manually curated sequences,
  is higher than the number of V1 sequences with indels
  (data not shown). This shows that simulated data does not completely
  capture the error rates.

  | UP9 Num. sequences | Num. insertions | Num. deletions | Num. indels |
  |:------------------:|:---------------:|:--------------:|:-----------:|
  |         7          |      0          |       0        |      0      |
  |        10          |      0          |       1        |      1      |
  |         2          |      1          |       0        |      1      |
  |         3          |      0          |       2        |      2      |
  |         1          |      1          |       1        |      2      |
  |         3          |      0          |       3        |      3      |

  Table showing the number insertions and deletions consensuses built
    using UP9 PCV2 sequences from Ukraine. Consensus were built using V3
    and were comapred to their manually curated V1 counterparts. Most of
    the un-manually curated V1 sequences had no indels to their manually
    curated sequences. (Sequences were in doi.org/10.3390/v14050924).

Also in the UP9 dataset, I found that find co-infections V3 detected 
  the same number of samples without having co-infections as V1.
  However, for samples with co-infections, V3 would often detect
  additional co-infections that had support from less than 1% of reads.
  This is worrying, because these could be false positives only seen
  at deeper read depths. I would advise removing any co-infections that
  have less than 1% of read support
  (See dataAnalaysis/all-kept-read-counts.tsv for my counts).

### Real data testing methods

find co-infections version three was run on datasets A, B, and A and B
  combined that Baloglu et al. (2021) used to benchmark the ASHURE
  pipeline (see dataAnalysis/README.md for more details). To avoid
  building a reference database I trimmed the reads down that mapped to
  one of the references for the 50 OTUs in the datasets. Supplemental
  mappings were discarded, but unmapped reads were kept
  (trimSamFile -keep-unmapped-reads). I also ran a set of separate tests
  were the reads were trimmed based on the input primers and not by the
  references. Scores were found by mapping the reads to their refrences
  and scoreReads. Any scores discared by scoreReads were marked NA.

Find co-infections was used with various combinations of Racon
  (-enable-racon), Medaka (-enable-medaka) -model r941_min_high_g303,
  and the the majority consensus step
  (disabled by -disable-majority-consensus). Primers were added with
  -primers primers.fasta when reference trimming was not done. In all
  cases the binning step was skipped with -skip-bin. The maximum read
  length was set to 750 base pairs (-max-read-len 750) and a cluster
  needed at least 0.3% of all clustered reads to be kept
  (-min-per-reads 0.003).

Baloglu, B., Chen, Z., Elbrecht, V., Braukmann, T., MacDonald, S. and
  Steinke, D. (2021). A workflow for accurate metabarcoding using
  nanopore MinION sequencing. Methods Ecol Evol. 2021;12:794-804.
  https://doi.org/10.1111/2041-210X.13561

### Real data testing results

![Stacked bar chart showing that Racon detected the most false positive
  consensus, had the most innacurate consensus for SNPs. Medaka
  improved the results from both Racon and the majority consensus step.
  Finally it also shows that the reference trimming reduced the number
  of false positives](figures/ASHURE-bench--snps.svg)

The chart above shows the number of false positive and correctly
  detected consensus find co-infections detected. The solid line
  indicates the total number of OTUs, while the dashed line indicates
  the number of OTUs that have over 100 reads without supplemental
  mappings. In each name maj standrs for majority consensus, rac stands
  for Racon, med sntads for Medaka, Ref stands for reference trimming 
  and Primer stands for primer trimming.

We can see that Racon detected the most false positives and had a higher
  number of SNP errors than either medaka or the majority consensus.
  Medaka did well alone, but had even better results when it was
  combined with the majority consensus step.

We also see a reduced chance of detecting false positives for Medaka and
  the majority consensus steps when the reads were trimmed by the
  reference. In all cases we found at least one consensus that did not
  map to a reference. However at least one of these NAs does have some
  support for being a real consensus (see dataAnaylsis/README.md for
  my evidence).

For indels we found that Medaka performed the best and that Racon
  performed the worst, with the majority consensus being a close second
  to Racon (see figures/ASHURE-bench--indels).

![Dot plot showing that few false positives detected by Medaka and the
  majority consensus step could be removed by requiring 1% of reads
  be present](figures/ASHURE-bench--percReads.svg)

We found that by requiring 1% of reads instead of 0.3% of reads would
  reduce the number of false positives detected for Medaka and majority
  consensus down to one unmapped consensus, which is likely the one
  consensus that might be a genuine OTU. However, this reduction does
  come with a large decrease in the number of detected OUTs.

### Take away:

We have shown that find co-infections V3 is better than find
  co-infections V2 at detecting co-infections in simulated datasets.
  However, we are unsure how well find co-infections will work when real
  reads are used. It is possible that improvements in version three may 
  only apply to simulated reads.

We also found that find co-infections did have similar results to 
  find co-infections version one for the porcine circovirus reads 
  version one was used on (doi.org/10.3390/v14050924), however, we
  also found that version three does produce more consensuses with
  indels and may have a problem with false positives at deep read
  depths. One problem is that this is not a ground truth dataset and
  thus, we have no idea what the real answer is for any of these
  samples.

We have tested find co-infections version three on a real
  dataset that was basecalled with guppy version 3.2.2 and found that 
  with the right settings it could produce reliable results. However,
  the sensitivity is reduced compared to other pipelines like ASHURE. 
  One weakness is that this dataset as OTUs that are 15% different and
  so it did not test find co-infections ability to detect similar
  co-infections.

We can also say from the real dataset that Medaka is the best option 
  and that the majority consensus step could provide some improvement
  in consensus accuracy. However, it does reduce the consensus length
  down a bit (see figures/ASHURE-bench--lenDiff.svg). Also not mentioned
  previously is that the primer trimming step does sometimes add some
  bases at the end (see figures/ASHURE-bench--extraConLen.svg).
 
## Weaknesses:

Find co-infections is set up for reads of at least 700 base pairs. Since
  I am using percentages, these settings will scale up. However, these
  settings will likely not scale down for smaller read lengths, such as
  100 base pairs.

Find co-infections has only been tested using data simulated from
  badread (Wick 2019) on default settings, which does produce low
  quality reads. However, these simulated reads can be used to build
  more accurate consensuses than consensuses built with real reads of
  the same quality
  (https://github.com/rrwick/Badread/tree/main/comparison).
  This could pose a problem for find co-infections, since find
  co-infections relies on accurate SNP calls in consensuses to identify
  co-infections and to hone in on a good set of reads for building the
  consensuses used in the consensus comparison steps. At the very least
  is means our accuracy results are optimistic.

Find co-infections has only been tested on data sets that contain two
  references per test. This means that I am unsure how find
  co-infections will act when their are more than two viral strains in a
  sample or when their are less than two viral strains in a sample.
  However, find co-infections version three did not detect any false 
  positives in the similarity data set, so it will probably not detect
  false positives consensuses when their is only one viral strain in a
  sample.

Find co-infections is designed to run on data sets that have very deep
  read depths, such as the read depths produced by PCR. So it will 
  probably not work on datasets not amplified by methods like PCR.

I need to test a somatic variant caller to round out the variant caller
  selection. Not likely to be done, but is something I missed.

The outy flow cells, with their rereading of reads will likely make
  find co-infections obsolete when they come out for public use
  (hopefully these will be released soon).

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
  - Raspberry pi 2B with a 32 bit Debian install (no Medaka)

## Possible plans:

1. Add in full consensus building step to version three and maybe some
   masking for poorly supported bases. The majority consensus step will
   remove these bases, so I am not to worried about this.
2. Maybe multi-thread fqGetIds (in V2 was fastqGrep). Likely will not
   happen.
3. Probably should document how to interface with the small programs I
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
    the benchmarking for my pipeline, and gave me advice. One big thing
    that helped was introducing me to the ASHURE data set.
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
  - To Baloglu et al. (2021), who I again never meant, but were willing
    to post their mock community on the SRA, which allowed me to test
    find co-infections on real data.
  - To the people (never met) who coded baba
    (https://baba.sourceforge.net/). It gave me a great visual on how
    the Needleman Wunsch algorithim worked. Their were many other
    sources, but this was the one that was the most useful to me.
  - I should also should thank the people how coded Minimap2, Racon,
    and Medaka, since these are integrated into find co-infections.
  - Finally the Universities of Fairbanks and Anchorage Alaska, who 
    provided me with TA ships while I was in graduate school. During
    this time I wrote and tested version one and two of find
    co-infections.
