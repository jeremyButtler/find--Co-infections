## Purpose

Find co-infections version three is composed of five smaller programs
  that if desired, can be complied separately. The purpose of this 
  readme is to give you an idea of what these programs are and how to
  compile them

### BuildCon

#### Overview of buildCon

BuildCon runs the consensus building method used in find co-infections.
  It can run the inbuilt consensus builder, Racon, or Medaka. It finds
  a good quality read to build a consensuses with by looking for the 
  read with the  highest median integer Q-score. The consensus is then
  built using a multi round polishing approach. You can also provide a
  tsv (-stats) made using scoreReads (mentioned later) to select the
  best read by mapping quality to a reference. 

  You can also build a consensus by polishing a reference (-ref).

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
  median Q-score, mean Q-score, length, alignened
  (non soft masked bases) length, aligned median Q-score, and
  aligned mean Q-score.

Score reads will also filter out reads with low mapping qualities
  (under 20), low Q-scores (under 13), and to large or to small by
  lengths.

Score reads requires an eqx cigar in the sam file. This can be made
  with minimap2 using --eqx. I would also advise removal of secondary
  alignments (--secondary=no), since these will result in inaccurate
  results when the primary alignment is in the opposite direction of the
  secondary alignment. This is due to score reads reading the sequence
  in the forward direction when it likely should be reading it
  backwards.

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
  of read ids. It runs on less memory than seqkit, but also runs slower
  than seqkit when you are extracting less than 33% of nanopore
  sequenced reads. When extracting over 33% of the file it is a toss up,
  with fqGetIds doing better as read depth increases. It does a bit
  better on Illumina reads than seqkit, but this is due to processing
  the read ids backwards instead of the normal forwards direction.

FqGetIds was first seen as fastqGrep in find co-infections. It is mainly
  here to remove additional dependencies.

The help message is called with -h.

#### How to build and run fqGetIds

```
cd V3
make fqGetIds
mv fqGetIds /path/to/install
chmod a+x /path/to/install/fqGetIds

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
  step for find co-infections to work and I am unaware of any other
  programs that do this particular task.

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
