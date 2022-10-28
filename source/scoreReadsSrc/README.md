# scoreReads:

## Use:

Find the number of mismatches and insertions above a certain quality threshold
  in a sam file. Also finds the mean q-score, mean aligned q-score,
  median q-score, median aligned q-score, length, and aligned length. Insertions
  are also filtered by their homopolymer length. While, deletions are only
  filtered by the length of their homopolymer (no q-score).

The sam file needs to have an eqx cigar entry (minimap2 --eqx). The stats,
  query id, and reference ids of alignments not ignored (beneath min threshold)
  are printed to stdout in a tab delimited format.

## Install:

make
mv scoreReads /path/to/install
chmod a+x /path/to/install/scoreReads

## How to use:

scoreReads -f file.sam > stats.tsv

minimap2 --eqx -ax map-ont ref.fasta reads.fastq | scoreReads -stdin > stats.tsv

Score reads does not output the sam file, so the sam file is lost when -stdin
  is used. Only use -stdin when you do not want to keep the sam file.

## Additional parameters:

1. -min-q: Min Q-score to keep alignment [Default 13]
2. -min-map-q: min mapping quality to keep alignment [Default 20]
3. -min-read-length: min read length to keep read [Default 600]
4. -max-read-length: max read length to keep read (0 is any) [Default 1000]
5. -ref-del: Reference to use in comparing deletions [Default: None]
  - This turns deletions into indels
6. -h or any non-valid parameter to print error message.

## Notes:

I am not sure how far I trust the mean and median Q-scores, since they often
  seem a high for data simulated by badread. I have double checked my
  calculations and can see nothing wrong.

Sadly this code is a mess and needs to be cleaned up at some point. However, I
  will not take the time to do this, since I am splitting the functions of 
  score reads into several smaller programs in my next version of
  findCoInfections. I will also be making a unique format to work with these
  modules. See the fastqAlnFormat.md for an idea of my future plans for
  scoreReads.
