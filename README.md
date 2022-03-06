## Co-infection detection scripts ##

These scripts can detect co-infections in Nanopore data, when a database of
 references is provided.

The amount of difference between references in your database will determine how
 sensitive the pipeline is to co-infections. If you provide a database were the
 references are at least 4% different, then you will not be able to detect
 co-infections that are 3% different.

However, with increased sensitivity comes the increased chance of detecting
 noise. More similar references, means more chance to track errors in the reads
 **WORK ON THE DESCRIPTION**. Allowing more reads to be miss-mapped to similar
 reference or keeping of low quality reads, due to have a closely matching
 reference in the database.

I found that a database with 2% difference between references worked well, so
 long as there are good filters in place. These settings are my default settings
 in my scripts.

### How to run ###

bash findCoInfections.sh -i reads.fastq -r references.fasta [options ...]

Some options you may be likely to change:

1. -p: prefix to call the output
2. -a: min read length to keep (also min aligned length for blastn)
3. -n: max read length to keep
4. -m: model to used to basecall the reads (currently r941_min_high_g351)
5. -t: number of threads to use
6. -h: print out all parameters you can change (help message)

### How the pipeline works ###

![Pipeline diagram](figures/co-infection--pipeline-figure.svg]

### Performance ###

We tested the pipeline by simulating reads for a co-infection. Our references
used for the simulation came from the PCV2 capsid gene. The capsid gene is ~ 700
bp long and has a similarity between 90% to 100%. Our references were split into
pairs, with each pair having a similarity between 98 to 98.5%.

Reads were simulated with badread, using a mean accuracy of 90% and 20000x read
depth. We simulated reads for co-infections having the minor strain at 50%, 5%,
and 1%. Co-infections were detected by inputting  the simulated reads (-i) and 
a reference database with 2% difference (-r) through my pipeline. Non-default
optional settings were: MAPQ of 20 or 30 (-Q 20 or -Q 30), do not check
differences mismatches between consensus genomes (-x 0), a max percent identity
between references of 99% (-g 1), and 0 reads to keep a bin (-K 0).

command line: 

```
bash binAndBuildSubtypes \
	-f reads.fastq \
	-r database.fasta \
	-p prefix \
	-Q 20 or 30 \
	-x 0 \
	-g 1 \
	-K 0;
```
			
We found that my input settings did detect co-infections at least half of the
time (Figure 1). However, there were times were we also missed co-infections or
detected an extra co-infection, which was from noisy reads (Figure 1).

[!Figure 1: Number of consensus genomes detected in my pipeline.
](../two-percent-results/Num-con-graph.tiff)

We simulated what would happen if we had removed bins will less then 100 reads
and consensus genomes with less then 0.3% of mismatches in R. To simulate the
percent of mismatches between consensus genomes mapping to the same reference we
multiplied the number of mismatches for the most error prone consensus genome
for each reference by 2.

We found that removing bins with 100 reads combined with removing consensus
genomes with less then 0.3% mismatches removed all noise, while having only a 
slight increase in missed consensus genomes (Figure 2). We, also noticed that a
MAPQ of 20 was able to detect more co-infections then a MAPQ of 30, when the
minor strain was at 1% (Figure 2).

[!Figure 2: Number of consensus genomes detected in my pipeline when bins had
 at least 100 reads and there was at least 0.3% mismatches different between 
 consensus genomes.
](../two-percent-results/Num-con-graph.tiff)

### Some future questions and directions ###

Some things that need further exploration:

1. Effect of read depth on co-infection detection
    - Deeper read depths may allow higher numbers of noisy reads to bin
2. Testing if requiring 0.4% of mapped reads per bin is needed.
3. Finding an ideal mapping quality
4. Having higher database similarities

Some things to add to my pipeline:

1. Only do the % similarity and % mismatch checks between consensus genomes on
   the region of interest. Right know it does these checks on the entire
   consensus.
2. Does building a consensus of consensus improve the read accuracy.
    - Idea here is sub-sample 1000 read depth bins to build multiple consensus.
    - Then make a majority rules consensus genome with samtools
    - Currently accuracy improves little past 300 reads (link in papers)
3. Can we use some of the tricks in NanoQ to detect variants?
    - Replace low quality and low Q-score bases with N (unless is another
      anonymous base [like d]) [per kept bin]
    - Cluster reads by mismatches, ignoring indels. Also, will ignore indels
      that are unique to only a few reads.
    - Build a consensus genome from the clusters
