# Purpose

Some of the steps for find co-infections version three have been set up
  so that they can be compiled separately. In some cases these are 
  specialized tasks that may not exists or are here to eliminate a
  dependency (they are still specialized). I am not sure how useful
  others will find these programs, but I have included them in case
  someone does find one of these programs useful. The goal of this
  readme is to show you how to build and use these programs.

### BuildCon

#### Overview of buildCon

BuildCon runs the consensus building method used in find co-infections.
  It can run the inbuilt consensus builder, Racon, or Medaka. It finds
  a good quality read to build a consensuses with by looking for the 
  read with the  highest median integer Q-score. The consensus is then
  built using a multi round polishing approach. You can also provide a
  tsv (-stats) made using scoreReads (mentioned later) to select the
  best read by mapping quality to a reference. 

You can also build consensus by polishing a reference (-ref). However,
  this only means that buildCon will attempt to build a consensus with
  the reference first. It will switch to hunting for the best read if
  a consensus could not be built with the reference.

BuildCon uses Minimap2 and can also use Racon and Medaka.

The help message is called with -h.

#### How to build and run buildCon

```
cd V3
make buildCon
mv buildCon /path/to/install
chmod a+x /path/to/install/buildCon

builCon -fastq reads.fastq
```

### ScoreReads

#### Overview of scoreReads

Score reads is designed to count the number of matches, SNPs,
  insertions, and deletions for every read in a sam file. It has two
  counts, one for total matches, SNPs, insertions, and deletions and
  the other for the number of filtered matches, SNPs, insertions, and
  deletions. The filtered output is selected by the bases q-score for
  everything but deletions. Insertions and deletions are also filtered
  by the size of the homopolymer they are in. For deletions this is a
  guess, using the largest neighboring homopolymer.

Additional stats extracted or are found include the mapping quality,
  median Q-score, mean Q-score, length, aligned (non soft masked bases)
  length, aligned median Q-score, and aligned mean Q-score.

Score reads will also filter out reads with low mapping qualities
  (under 20), low Q-scores (under 13), and to large or to small by
  lengths.

Score reads requires an eqx cigar in the sam file. This can be made
  with minimap2 using --eqx. I would also advise removal of secondary
  alignments (--secondary=no), since these will result in inaccurate
  results when the primary alignment is in the opposite direction of the
  secondary alignment. This is due to score reads reading the sequence
  in the forward direction when it should be reading it backwards.

The help message is called with -h.

#### How to build and run scoreReads

```
cd V3
make scoreReads
mv scoreReads /path/to/install
chmod a+x /path/to/install/scoreReads

scoreReads -file reads.sam > read-stats.tsv
minimap2 --secondary=no --eqx -a ref.fasta reads.fastq | scoreReads -stdin > stats.tsv
```

### fqGetIds

FqGetIds is a program that extracts fastq entries using a provided set
  of read ids. It runs on less memory than seqkit and at times can be
  a bit faster than seqkit for extracting reads by id from a fastq
  files.

fqGetIds combines a hashing table for fast read id access and then an
  AVL tree to deal with collisions. For hashing only hex characters 
  (0 to 9 and a to f) are used in read id comparison. This works well
  for Nanopore data, were the read ids are hex based, but may have some
  rare issues when dealing with Illumina data. Each read id is converted
  to a big number (multiple integers) and the sum of all integer values
  is found. The hash is found by using Knuths multiplicative hash on the
  sum. For comparisons the sum is first compared and then each element
  in the big number is compared.

For benchmarking I used the Illumina data set originally used to
  benchmark seqkit (dataset_C.fastq) and a data set with 6 million 
  Nanopore PCR reads (around 600 to 700 base pairs). The testes were 
  run on a raspberry pi 1B with a 32 bit Debian linux instal (slow IO)
  and a computer with a solid state harddrive (fast IO) running a 64 bit
  Debian OS. For the raspberry pi I had to reduce read depths to around
  800000 reads to avoid memory issues (max of 1Gb).

![
  Figure showing fqGetIds is a as fast or faster than seqkit for
  Nanopore data and Illumina data
](figures/fqGetIds-seqkit-full-bench--elapsed-time.svg)

The elapsed time figure (first figure) is showing that fqGetIds is 
  taking as much or less time to extract reads from a fastq file with
  Nanopore or Illumina read ids. The difference becomes more dramatic
  as the percentage of extracted reads increases.

![
  Figure showing fqGetIds uses less memory for Nanopore and Illumina
  fastq files
](figures/fqGetIds-seqkit-full-bench--memory.svg)

The memory figure (second figure) shows that fqGetIds is using less
  memory than seqkit at all read depths to extract reads.

One thing I should note is that fqGetIds currently only compares the hex
  based characters (0 to 9 and a to f) in the read ids. This is ok for
  Nanpore data, which only uses hex based characters in their read ids,
  but might be a problem for Illumina data. One problem is that it can
  not tell the difference between 86:123 and 861:23. I did not have this
  problem in my benchmarking tests, but it still is possible that this
  problem might occur. Fixing this would likely result in slower times
  and more memory usage for Illumina data, but not Nanopore data.

Also I am not sure how well fqGetIds would work if you merged Illumina
  fastq files, since the Illumina sequencer id and flow cell ids include
  non-hex characters.

FqGetIds was first seen as fastqGrep in find co-infections. It is mainly
  here to remove additional dependencies.

The help message is called with -h.

#### How to build fqGetIds without multithreading

```
cd V3
make fqGetIds
mv fqGetIds /path/to/install
chmod a+x /path/to/install/fqGetIds
```

#### how to run fqGetIds

```
fqGetIds -fastq reads.fastq -f ids.txt > extractedReads.fastq
cat reads.fastq | fqGetIds -stdin-fastq -f ids.txt > extractReads.fastq
cat ids.txt | fqGetIds -fastq reads.fastq -stdin-filt > extReads.fastq
```

### BinReads

BinReads is the binning method used in find co-infections version three.
  It bins a set of reads to a reference and provides quality filters to
  control which reads should be kept or discarded. It will also remove
  any bin with a low number of reads and provide a list of how many 
  reads were binned to each reference.

BinReads uses Minimap2 to map reads to their references.

The help message is called with -h.

#### How to build and run binReads

```
cd V3
make binReads
mv binReads /path/to/install
chmod a+x /path/to/install/binReads

binReads -fastq reads.fastq -ref references.fasta
```

### TrimSamFile

TrimSamFile removes the softmasking at the start and end of every
  alignment in a sam file. It does this by trimming the cigar, sequence,
  and q-score entries. It is not much of a program, but it is a needed
  step for find co-infections to work. I am unaware of any other program
  that does this particular task.

The help message is called with -h.

#### How to build and run trimSamFile

```
cd V3
make trimSam
mv trimSamFile /path/to/install
chmod a+x /path/to/install/trimSamFile

trimSamFile -file reads.sam > trimmed.sam
minimap2 -a ref.fasta reads.fastq | trimSamFile -stdin > trimmed.sam
```

### trimPrimers

TrimPrimers uses a set of primers to trim primer sequences from reads.
  TrimPrimers will print out multiple fastq entries for a single read
  when the trimming of primers would split a read into multiple
  sequences. In this case each entry is marked with a --D1_19-x, were x
  is an integer that represents the number of splits the read had.

TrimPrimers uses minimap2 to map primers to a set of reads. This can
  be skipped by providing a paf file with mappings.

The help message can be called with -h.

#### How to build and run trimPrimers

```
cd V3
make trimPrimers
mv trimPrimers /path/to/install
chmod a+x /path/to/install/trimPrimers

trimPrimers -fastq reads.fastq -primers primers.fasta > reads-trim.fastq
trimPrimers -fastq reads.fastq -paf mappings.paf > reads-trim.fastq
minimap2 -k5 -w1 -s 20 -P primers.fasta reads.fastq | trimPrimers -fastq reads.fastq -stdin-paf > reads-trim.fastq
```

### filterReads

FilterReads removes reads that are under a minimum length, over a
  maximum length, under a minimum median quality score, and under a 
  minimum quality score.

The help message can be called with -h.

#### Building and using filterReads

```
cd V3
make filterReads
mv filterReads /path/to/install
chmod a+x /path/to/install/filterReads

filterReads -fastq reads.fastq > filtered-reads.fastq
trimPrimers -fastq reads.fq -paf map.paf | filterReads -stdin > file.fastq
```

### extractTopReads

ExtractTopReads extracts the highest scoring top X (default 300) reads
  from a fastq file. It finds top reads by using length and median
  quality scores. When a reference is used it will also extract the top
  reads by using mapping qualities.

One warning is that the structures used with fqGetIds are predeclared
  and stored on the stack instead of the heap, which means that their is
  less memory to use. So this program may error out (not enough memory)
  if you are extracting to many reads (probably over 10000).

The help message can be called with -h.

#### How to build and use extractTopReads

```
cd V3
make extractReads
mv extractTopReads /path/to/install
chmod a+x /path/to/install/extractTopReads

extractTopReads -fastq reads.fastq > top-reads.fastq
```
