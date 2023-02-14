# Use:

Used to find co-infections in nanopore sequenced reads.

## Requirements:

1. minimap2 (https://github.com/lh3/minimap2)
2. samtools (https://github.com/samtools)
  - You can typically install this from a repository.
3. racon (https://github.com/isovic/racon)
4. medaka (https://github.com/nanoporetech/medaka)
  - Install as a python virtual environment in the users home directory (~)
  - Change to conda by uncommenting (remove # symbol) lines 176, 177, and 188 in
    coInfectScripts/buildConsensus.sh. You will also need to comment out
    (add a # at start) lines 179 and 187.
    

## Install:

make;

### Install in specific location

mv findCoInfections.sh coInfectScripts /path/to/install;

chmod -R a+x /path/to/install/coInfectScripts /path/to/install/findCoInfections.sh;

## Run:

bash findCoInfections.sh -f reads.fastq -ref reference-database.fastq -p prefix

## Some quick options:

  - -f: Fastq file to search for co-infections in                    [Required]

  - -ref: Reference to filter, trim, & bin reads with                [Required]
    - Use to bin reads and trim final consensuses for
      filtering. Both the trimmed and full length
      consensuses are kept.

  - -full-ref: Full length reference to align reads to               [None]
    - Used to pick the best read when building final consensuses

  - -t: Number of threads to use                                     [3]

  - -p: prefix for output bins & consensuses                         [out]

  - -h: Print help message

  - -h-detail: Print the extended help message

## Testing:

This may not be the best analysis, but is just a quick run through of what I 
  found and how findCoInfections compares to other pipelines. I will be updating
  this with a comparison to Nano-Q in a bit. It also is not the best writing,
  but it will do for now.

We compared findCoInfections to Longshot and Clair3 using simulated
  co-infections. Reads for the major strain and minor strain in each simulated 
  co-infection were simulated using references for porcine circovirus type 2
  ORF2 and badread simulate v 0.2.0 with a --seed of 1026 and a --quantity of
  4000x, 15000x, 40000x for strains that were from different genotypes and
  20000x reads for strains that were 1.6% to 3% different. Our pipeline was 
  tested with a reference database containing references that were all at least
  2% different. We also tested a sparse database that only had four references.
  For Clair3 and Longshot, we used the major strain, minor strain,
  consensus genome, and a reference that was at least 3% different as the 
  reference for variant calling.

![
  Figure-1:
  Number of false positive SNPs in the consensus genomes built using Clair3.
  Major, minor, consensus, and distant indicates the reference used for variant
  calling.
](../figures/Genotype--ClairRefTest.png)

We found that the accuracy of Clair3 was reliant on the input reference genome.
  With more distant reference genomes resulting in consensuses with more 
  false positive SNPs then using the consensus genome or the reference one of
  the strains was simulated from (Figure-1). This reliance on a close reference
  resulted in at least one consensuses having many false positive SNPs when
  the major and minor strain were from different genotypes (Figure-1). However,
  this effect was reduced as the major and minor strains became more similar
  (See figures/Similarity--20000--ClairRefTest.png).

We found similar results for Longshot (see figures/Genotype--LongshotRefTest.png
  and figures/Similarity--20000--LongshotRefTest.png)

![
  Figure-2:
  Number of correct consensus identified by Longshot, Clair3, and
  findCoInfections. New, and old indicates if the newer (version 2) or
  older version of findCoInfections was used. None indicates that all reads
  were used to build a consensus.  While cluster indicates the new pipeline was
  used with clustering. The horizontal line indicates the actual number of
  consensuses.
](../figures/Similarity--20000--correctConsensusCount.png)

We compared how well Longshot, Clair3, and findCoInfections could detect
  co-infections. We found that findCoInfections often outperformed Longshot and
  Clair in detecting co-infections (Figure-2). Also, we found that the
  consensuses built by findCoInfections often had fewer false positive SNPs
  then the consensuses built by Longshot and Clair3 (Figure-2).

![
  Figure-3:
  Number of incorrectly identified consensuses. New and old indicate if the 
  new or older pipeline was used. New-1.4-50 indicates the percent
  difference between variants (1.4) and the percent of reads from the major
  variant (50%).
](../figures/Similarity--20000--wrongConsensusCount.png)

We also looked at the number of false positive consensuses findCoInfections
  built. We found that The newer version of findCoInfections built fewer
  false positive consensuses (only 2) than the version of findCoInfections
  (Figure-3). However, many of the false positives from the older version of
  findCoInfections had no false positive SNPs (Figure-3). Suggesting that the
  read trimming step of the new findCoInfections improved its ability to detect
  false positive consensuses.

![
  Figure-4:
  Number of correct consensus identified by findCoInfections when a sparse 
  database is used. New, and old indicates if the newer (version 2) or
  older version of findCoInfections was used. None indicates that all reads
  were used to build a consensus.  While cluster indicates the new pipeline was
  used with clustering. The horizontal line indicates the actual number of
  consensuses.
](../figures/Sparse--Similarity--20000--correctConsensusCount.png)

We tested the effect of the database by running findCoInfections with a database
  containing only four references. We found that the number of detected
  co-infections decreased for all versions of findCoInfections (Figure-4).
  However, at the highest major and minor strain ratios, we found that
  clustering did improve the detection of co-infections (Figure-4). However, 
  this improvement greatly decreased with as the major strain contributed more
  reads and as the major and minor strain became more similar (Figure-4; Also
  see figures/Sparse--Genotype--correctConsensusCount.png).

![
  Figure-5:
  Number of incorrectly identified consensuses when a sparse database was used.
  New and old indicate if the new or older pipeline was used. Clustering
  indicates if clustering was used with version 2 of findCoInfections.
  New-1.4-50 indicates the percent difference between variants (1.4) and the
  percent of reads from the major variant (50%).
](../figures/Sparse--Similarity--20000--wrongConsensusCount.png)

Finally, we looked at how well the new and old version of findCoInfections
  could avoid false positive consensuses when a sparse database was used. We
  found that the new findCoInfections without clustering detected the fewest 
  false positives (1), followed by clustering, and then the old
  findCoInfections (Figure-5). With almost all of clusterings false positives
  being when the major strain contributed only 50% of reads and was 3% different
  then the minor strain (Figure-5). However, the old pipeline did not follow
  any pattern (Figure-5).

## Take away:

We have shown that findCoInfections can detect co-infections more often then 
  either Longshot or Clair3. We have also shown that those consensuses are 
  not effected by the reference genome. This results in fewer false positive 
  SNPs in at least one consensus when the major and minor strains are more 
  distant. Finally, we have shown that the current version of findCoInfections
  without clustering will not detect many false positive consensuses, even
  when the database is very limited. However, when the database is sparse, 
  the ability of findCoInfections to detect co-infections is limited.

## Future directions:

Working on paper and building the next version of findCoInfections. 
