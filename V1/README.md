## Co-infection detection scripts ##
          
I am currently working on testing the second version of this pipeline.
One problem I found is that by not trimming reads to the references
can result in some false postives being kept. I think it is due to
the pipeline including untrimmed regions in the comparision. The next
version will trim reads and consensuses.

This pipeline detects co-infections in Nanopore reads using a database
 of references. References should be a region of the genome with high
 variation, instead of the full genome. The pipeline maps the reads to the best 
 reference and then splits the reads into bins by reference. Each bin is
 polished to make a consensus for each sequence in the co-infection. For a
 flow chart showing all steps in this pipeline see Figure 1.

The consensuses can be longer than the region provided, allowing for a region 
 of 600 to 700 bp to bin reads from larger, 2kb amplicons. However, the pipeline
 does not guarantee that a read covering the entire 2kb region will be used in
 the polishing step. Manual steps may be needed if the user wants consensuses
 longer than the region of interest.

### How to run ###

bash findCoInfections.sh -i reads.fastq -r references.fasta [options ...]

Some arguments you may want to change:

1. -p: prefix to call the output
2. -a: min read length to keep (also min aligned length for blastn)
3. -n: max read length to keep
4. -m: model to used to basecall the reads (currently r941_min_high_g351)
5. -t: number of threads to use
6. -h: print out all parameters you can change (help message)

### Requirements: ###

1. minimap2
2. samtools version 1.14: Needs to be able to run this command
    - samtools view --min-qlen number -e "avg(qual)>=number" -q number
3. bamtools
4. filtlong
5. racon
6. medaka
    - Should be installed by miniconda.
    - If not change lines 179, 180, and 188 in buildConsensus.sh to your command
      to activate and deactivate medaka.
7. blastn and makeblastdb (in same blast install package from ncbi)

### How the pipeline works ###

<figure>
    <img
      src="figures/co-infection--pipeline-figure.svg"
      width=60%
      height=60%
    > <!--show my pipeline flow digram image-->
    <figcaption>
      <em>
        <b>Figure 1:</b> Flow chart showing each step in our pipeline.
      </em>
    </figcaption>
<figure>
    
&nbsp;
    
### Performance ###

We tested our pipeline by simulating reads from 237 pairs of porcine circovirus
 type 2 capsid genes, which are roughly 700bp long. Each pair of capsid genes
 were from a single subtype and were 98% to 98.5% similar to each other. We had
 at least one capsid pair for each PCV2 subtype. Reads were simulated with
 badread, using a mean accuracy of 90% and 20000x read depth. For each capsid 
 pair, we varied the percentage of reads from the minor variant. Our levels 
 were 50%, 5%, and 1% of reads from the minor variant.

Command line: 

```
bash findCoInfections.sh \
	-f reads.fastq \
	-r capsid.fasta \
	-p prefix \
	-Q 20 or 30 \
	-x 0 \
	-g 1 \
	-K 0;
```

<figure>                             <!--make a Figure-->
    <image src="figures/Num-con-graph.png" 
       width=50%
       height=50%
    >
    <figcaption>             <!--Add a caption to the Figure-->
      <em>                   <!--<em> for italics, <b> for bold-->
        <b>Figure 2:</b>
        Number of co-infections detected for each reference pair.
      </em>
    </figcaption>
</figure>
    
&nbsp; <!--Add an empty line-->
    
We found that our initial pipeline settings detected co-infections at least half
 of the time (Figure 2). However, some of time we missed co-infections or
 detected
 extra co-infections, which were from noisy reads (Figure 2). The number of
 missed co-infections or extra co-infections increased as the percentage of
 reads from the minor variant decreased (Figure 2).

&nbsp; <!--Add an empty line-->
    
<figure>
    <img src="figures/Num-con-depth100-misPerc0_3-filter-graph.png" 
       width=50%
       height=50%
     >
    <figcaption>
        <em>
          <b>Figure 3:</b> Number of co-infections detected for each reference
          pair when we require at least 100 reads per consensus and a minimum
          difference of 0.3% mismatches between all consensuses.
         </em>
    <figcaption>
</figure>
    
&nbsp;
        
We simulated what would happen if we had removed bins will fewer than 100 reads
 and consensus genomes with fewer than 0.3% of mismatches in R. To simulate the
 percent of mismatches between consensus genomes mapping to the same reference
 we multiplied the number of mismatches for the most error prone consensus
 genome in a reference pair by 2.

We found that removing bins with 100 reads combined with removing consensus
 genomes with less than 0.3% mismatches removed all noise (extra co-infections),
 while having only a slight increase in missed co-infections (Figure 3). Also, a
 MAPQ of 20 detected more co-infections than a MAPQ of 30 (Figure 3).

### Accuracy ###

<figure>
  <img src="figures/Num-con-depth100-misPerc0_3-mismatch-filter-graph.png"
    width=50%
    height=50%
  > <!--image settings-->
  <figcaption>
    <em>
      <b>Figure 4:</b> Number of consensuses with one or more mismatches. The
      horizontal bar shows the number of consensuses from minor variants when
      95% or 99% of reads were from the major variant. The number of consensus
      detected in each consensuses pair is shown on the x-axis as 1 or 2. The
      high error rate in the consensus from 99-1 1 was due to aggressive
      filtering on the major variant.
    </em>
  </figcaption>
<figure>

&nbsp;

We looked at the error rate in our consensuses by looking at the number of
 indels (Figure 5) and the number of mismatches (Figure 4) in each consensus.
 We found mismatches were rare (13 out of 237 consensuses). Half of 
 the consensuses with mismatches were from reference pairs that detected only
 one variant and had 50% of reads from the minor variant (50-50 1, Figure 4).

Mismatches in consensuses built from reference pairs that detected only one
 variant were mostly from the reference pairs that had 50% of reads from the
 minor variant (Figure 4). The only exception is one case that had 99% of reads
 from the major consensus (99-1 1, Figure 4). However, in the 99-1 1 case, the
 major variant lost most of its reads to aggressive filtering. Which, resulted
 in the minor variant having 272 of the 949 reads in the bin used to make the
 consensus. In both cases, mismatches are likely from the polisher merging both
 variants into one consensus.

A quarter of consensuses with mismatches were from consensus pairs that
 detected two variants and had 99% of reads from the major variant 
 (99-1 2, Figure 4). All 99-1 2 consensuses with mismatches were from the minor
 variant (Figure 4.). Sometimes, the minor variant had fewer than 200 reads
 in a bin and had some misbinned, noisy reads from the major variant.
 Allowing mismatches from lower read depths or merging of the misbinned major
 variant reads with the correctly binned minor variant reads.
 

<figure>
    <img
      src="figures/Num-con-depth100-misPerc0_3-indel-filter-graph.png"
      width=50%
      height=50%
    > <!--img: insert the indel graph image-->
    <figcaption>
      <em> <!--start italics-->
        <b>Figure 5:</b>  Number of consensuses with at least one indel. The
        horizontal bar shows the number of consensuses from minor variants when
        95% or 99% of reads were from the major variant. The number of
        consensus detected in each consensuses pair is shown on the x-axis as 
        1 or 2.
      </em> <!--end italics--> 
    </figcaption>
</figure>

&nbsp;

We found that most consensuses (300 of 474) had at least one indel when the 
 minor variant contributed 50% of reads (Figure 5). Most of the difference
 between the ratios of major and minor variants were due to an increase in
 missed co-infections (Figure 5). The missed co-infections explain the reduced
 the number of minor variant consensuses seen as the percent of reads from the
 major variant increases.

### Some future questions and directions ###

1. Right now, the similarity between consensuses is found by comparing the
   entire consensus, instead of just the region of interest. This is ok when
   the amplicons only amplify the region of interest, but may cause missed
   co-infections if amplicons are larger than the region of interest.
2. Use a full genome reference to detect and use longer reads for polishing.
   Allowing the user to not have to do manual steps to ensure consensuses are
   longer than the region of interest.
