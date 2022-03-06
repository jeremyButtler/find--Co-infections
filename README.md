## Co-infection detection scripts ##
          
This pipeline detectects co-infections in Nanopore reads using a database
 of references. The database should not conatain full genome reads, but
 instead have regions that distiguish variants for your virus.

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

<figure>
    <img src="figures/co-infection--pipeline-figure.svg" width=50% height=50%>
    <figcaption><em><b>Figure 1:</b>
        Flow digram showing each step in our pipeline
        </em>
    </figcaption>
<figure>
    
&nbsp;
    
### Performance ###

We tested the pipeline by simulating reads for pairs of porcine circovirus type 2
 capsid genes (around 700bp long). Each pair of capsid genes were 98% to 98.5%
 similiar.Reads were simulated with badread, using a mean accuracy of 90% and 20000x
 read depth. We simulated reads for co-infections having the minor strain at 50%, 5%,
 and 1%.

command line: 

```
bash findCoInfections.sh \
	-f reads.fastq \
	-r database.fasta \
	-p prefix \
	-Q 20 or 30 \
	-x 0 \
	-g 1 \
	-K 0;
```

<figure>                             <!--make a figure-->
    <image src="figures/Num-con-graph.png" 
       width=50%
       height=50%
    >
    <figcaption>             <!--Add a caption to the figure-->
        <em><b>Figure 2:</b> <!--<b> is to bold text, <em> for italics-->
            Number of co-infections detected for each reference pair.
    </em></figcaption>
</figure>
    
&nbsp; <!--Add an empty line-->
    
We found that my input settings did detect co-infections at least half of the
time (Figure 2). However, there were times were we also missed co-infections or
detected an extra co-infection, which was from noisy reads (Figure 2).

&nbsp; <!--Add an empty line-->
    
<figure>
    <img src="figures/Num-con-depth100-misPerc0_3-filter-graph.png" 
       width=50%
       height=50%
     >
    <figcaption>
        <em><b>Figure 3:</b>
            Number of co-infections detected for each reference pair when
            we require at least 100 reads per consensus and a differnce of 
            0.3% mismatches between all consensus.
         </em>
    <figcaption>
</figure>
    
&nbsp;
        
We simulated what would happen if we had removed bins will less then 100 reads
 and consensus genomes with less then 0.3% of mismatches in R. To simulate the
 percent of mismatches between consensus genomes mapping to the same reference we
 multiplied the number of mismatches for the most error prone consensus genome
 in a reference pair by 2.

We found that removing bins with 100 reads combined with removing consensus
 genomes with less then 0.3% mismatches removed all noise, while having only a 
 slight increase in missed co-infections (Figure 3). Also, a MAPQ of 20 was able to 
 detect more co-infections then a MAPQ of 30 (figure 3).
        
### Some future questions and directions ###

Some things that need further exploration:

1. Effect of read depth on co-infection detection
    - Deeper read depths may allow higher numbers of noisy reads to bin
2. Finding an ideal mapping quality

Some things to add to my pipeline:

1. Only do the % similarity and % mismatch checks between consensus genomes on
   the region of interest. Right know it does these checks on the entire
   consensus.
2. Right know this is limited to your data base. It would be nice to do
   some clustering after binning to detect variants that were missed in
   the binning step.
    - Replace low quality and low Q-score bases with N (unless is another
      anonymous base [like d]) [per kept bin]
    - Cluster reads by mismatches, ignoring indels. Also, will ignore indels
      that are unique to only a few reads.
    - Build a consensus genome from the clusters
