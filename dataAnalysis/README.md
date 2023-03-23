## Purpose

This directory has the data, scripts, reference pairs, and references I
  used to benchmark find co-infections version three. For all tests I
  simulated reads using Badread v0.2.0. It also includes the
  benchmarking scripts and data I used to benchmark fqGetIds.

## ASHURE dataset benchmarking

The ASHURE dataset I used to benchmark find co-infections was an mock
  community I download from NCBI's SRA. It was originally used to
  benchmark how well the ASHURE pipeline could detect OTUs in Nanopore
  sequenced, PCR amplified reads. It was created and submitted to the
  SRA by Baloglu et al. (2021). The dataset has two separate sequencing
  runs (A and B) that contain the same 50 OUTs, but were prepared
  differently. Both datasets were amplified with rolling circle
  replication and were basecalled with guppy version 3.2.2
  (See Baloglu et al. (2021) for further details). The rolling circle 
  replication did make analysis more tricky for my pipeline.

My ASHURE dataset benchmarking scrips and results can be found in
  ASHURE-dataset-benchmarking.

Baloglu, B., Chen, Z., Elbrecht, V., Braukmann, T., MacDonald, S. and
  Steinke, D. (2021). A workflow for accurate metabarcoding using
  nanopore MinION sequencing. Methods Ecol Evol. 2021;12:794-804.
  https://doi.org/10.1111/2041-210X.13561

### Files in ASHURE-dataset-benchmarking

1. ASHURE-SRA-accensions-for-dataset-A.txt
    - Accession numbers to download dataset A with NCBI prefetch
    - ```while IFS= read -r lineStr; do prefetch "$lineStr"; fastq-dump "$lineStr"; done < ASHURE-SRA-accensions-for-dataset-A.txt```
2. ASHURE-SRA-accensions-for-dataset-B.txt
    - Accession numbers to download dataset A with NCBI prefetch
    - ```while IFS= read -r lineStr; do prefetch "$lineStr"; fastq-dump "$lineStr"; done < ASHURE-SRA-accensions-for-dataset-B.txt```
3. ASHURE-dataset-primers.fasta
    - Primers to use with find co-infections for benchmarking the ASHURE
      pipeline (from Baloglu et al. 2021).
4. ASHURE-dataset-references.fasta
    - References for the 50 OTUs in ASURE (from Baloglu et al. 2021).
    - Each reference is marked by the haplotype ASHURE Detected, not
      the accession numbers.
5. ASHURE-graphs.r
    - Makes graphs (in ../figures) used in this repository.
    - Is hardcoded for ASHURE-full-20230314-trim--duplicate-marked.tsv
    - Run: Rscript ASHURE-graphs.r
6. ASHURE-full-20230314-trim--duplicate-marked.tsv
    - File with benchmarking results with failed runs removed and
      duplicates marked
7. ASHURE-full-20230314-trim--stats.tsv
    - File output by benchASHUREDataset.sh
8. benchASHUREDataset.sh
    - Runs find co-infections for a single benchmark
    - Run: ```benchASHUREDataset.sh -fastq reads.fastq -ref refs.fasta -prefix output-name```
        - Needs trimSamFile (cd V3; make trimSam; mv trimSamFile ../)
        - Needs trimPrimers (cd V3; make trimPrimers; mv trimPrimers ../)
        - Needs findCoInft (cd V3; make findCoInft; mv findCoInft ../)
        - Needs scoreReads (cd V3; make scoreReads; mv scoreReads ../)
        - Needs GNU time (/usr/bin/time)
    - Commands:
        - -primers file.fasta: primers to trim with
        - -fastq file.fastq: Reads to find co-infections in
        - -refs ref.fasta: References to trim (no primers) and score
          clusters with
        - -max-read-len: Maximum read length to keep a read
        - -min-read-len: Minimum read length to keep a read
        - -min-perc-reads: Minimum percentage of clustered reads needed
          to keep a cluster
        - -enable-racon: Use racon
        - -enable-medaka: Use medaka
        - -disable-maj: Disable the majority consensus step
        - -rnds-polish: Number of times to rebuild the consensus
        - -keep-files: Keep fastq and fasta files when done
            - (discard)
9. extra-OTU--majority-consensus-only.fasta
    - An extra OTU my pipeline detected in the ASHURE data set.
    - It was 98% similarity to *Leptolegina* sp. (BOLD:AAX5717, seq id.
      HQ708212, file: ASHURE-Extra-OTU.fasta).
    - Supplemental table 1. in Baloglu et al. 2021 had 18 Illumina reads
      map to this reference with 98% similarity. The Illumina reads 
      only coverd half of this reference.
    - This support does not prove this extra OTU is reall, but does
      raise the possiblity.


## fqGetIds benchmarking

All fqGetIds benchmarkings scripts and data are in the
  fqGetIds-benchmarking directory.

1. The bench-fastqGrep-seqkit.sh bash script contains the commands I
   used to benchmark fqGetIds against seqkit. The fastq file is provided
   as the first argument and the number of replicates as the second.
   For Illumina data you will want to uncomment line 24
   (`sedCmdStr="p;n;n;n;n;"`) and comment out line 25
   (`sedCmdStr="p;n;n;n;"`). Otherwise it will look like fqGetIds is
   getting the wrong read counts (Illumina has six lines per fastq,
   entry, while Nanopore has four lines per fastq entry).
   - To run script: `bash bench-fastqGrep-seqkit.sh reads.fastq 10`
2. fqGetIds-bench-good-and-pi.tsv contains that data I used to benchmark
   fqGetIds and seqkit.
   - Note some of read extractions sizes for seqkit are off, due to me
     not realizing that seqkit converted the six line Illumina
     entries into four line entries. This was fixed for the Pi run, but
     not the good/fast IO run.
3. fqGetIdsSeqkitBench.r contains the commands I used to make the
   figures I posted for benchmarking. The file and commands are
   hard coded, so you will have to modify this script for your own uses.
   This will mainly be in section one (Sec-1) and section two (Sec-2).

## Read simulation and reference pairs used in read simulation:

### How reads were simulated:

All reads were simulated using the twoRefRunBadread.sh script in the
  benchmarking-scripts directory (all my scripts that simulate reads
  call this script). To ensure consistency across all tests I used the
  same seed (1026) with badread. This seed is the default setting for
  twoRefFunBadread.sh.

### Labeling system:

All reference pair files are labeled with the percent of reads from the
  major strain, percent of reads from the minor strain, and the
  references used to make them. The format is
  %majorReads-%minorReads-majorReference-minorReference.fasta.

### Hepatitis C data set:

The reference pairs I used to simulate co-infections for the HCV
  (hepatitis C) test can be found in HCV-reference-pairs. The
  individualReferences folder in HCV-reference-pairs has a fasta file
  for each reference used to make the reference pairs.

This data set has four different major/minor strain ratios. The first
  is 50/50 (50% of reads from the major strain), then 75/25
  (75% of reads from the major strain), then 95/5 (95% of reads from the
  major strain), and finally 99/1 (99% of reads from the major strain).

Badread commands for this test:
   - badread simulate --quantity 40000x --seed 1026 --reference file.fasta
   - badread simulate --quantity 15000x --seed 1026 --reference file.fasta
   - badread simulate --quantity 4000x --seed 1026 --reference file.fasta

|  Major   | Major % Diff. from Db. |  Minor   | Minor % Diff. from Db. | Major/Minor Difference |
|:--------:|:----------------------:|:--------:|:----------------------:|:----------------------:|
| GQ370235 |         1.59%          | GQ370237 |         1.59%          |          1.98%         |
| GQ370234 |         1.03%          | GQ370236 |         1.51%          |          1.90%         |
| GQ370229 |         0.63%          | GQ370233 |         1.43%          |          2.06%         |
| GQ370226 |         1.75%          | GQ370232 |         1.83%          |          1.90%         |
| GQ370220 |         1.75%          | GQ370231 |         1.67%          |          2.06%         |
| GQ370228 |         1.83%          | GQ370230 |         1.83%          |          2.14%         |
| GQ370224 |         1.19%          | GQ370227 |         1.83%          |          2.14%         |
| GQ370223 |         1.51%          | GQ370225 |         1.59%          |          1.83%         |
| GQ370127 |         1.35%          | GQ370131 |         1.19%          |          1.90%         |
| GQ370123 |         1.90%          | GQ370128 |         1.03%          |          2.14%         |
| GQ370237 |         1.59%          | GQ370309 |         1.11%          |         10.87%         |
| GQ370233 |         1.43%          | GQ370323 |         1.90%          |         10.56%         |
| GQ370309 |         1.11%          | GQ370235 |         1.59%          |         10.87%         |
| GQ370323 |         1.90%          | GQ370234 |         1.03%          |         10.71%         |
| GQ370309 |         1.11%          | GQ370319 |         1.35%          |         10.79%         |
| GQ370236 |         1.51%          | GQ370309 |         1.11%          |         10.79%         |
| GQ370309 |         1.11%          | GQ370230 |         1.83%          |         10.63%         |
| GQ370232 |         1.83%          | GQ370309 |         1.11%          |         10.71%         |
| GQ370223 |         1.51%          | GQ370323 |         1.90%          |         10.71%         |
| GQ370323 |         1.90%          | GQ370225 |         1.59%          |         10.71%         |

Table:
  Major and Minor is the reference used to simulate reads for the major
    or minor strain
  % Diff. from Db. is the percent difference from the closest reference
    in the HCV database.

### Porcine cirovirus type 2 similarity test:

The goal of this test was to test how well find co-infections could
  respond to very similar co-infections. The similarities tested were
  97%, 98%, and 98.5%. Each percentage had 30 reference pairs, with ten
  of those reference pairs having both references the database (except
  for the 98.5%, wich was to similar to have this). The remaining 20
  references pairs either had only one reference in the data base or no
  references in the database. I simluated reads twice for the ten
  references pairs with one reference in the database. To ensure
  different results I swapped the major and minor strain. This double
  simulation is not shown in the table, but is why each similarity level
  has an additional ten files.

This data set, with the one reference in database duplicates, has a
  30 reference pairs for the 98.5% similarity test and 40 reference
  pairs for the 98% and 97% datasets.

This data set has two different major/minor strain ratios. The first
  is 50/50 (50% of reads from the major strain) and the second is
  80/20 (80% of reads from the major strain).

This data set has two different major/minor strain ratios. The first
  is 50/50 and the second is 80/20.

Badread command for this test:
   - badread simulate --quantity 20000x --seed 1026 --reference file.fasta

| Reference-1 | Ref-1 DB %id | Reference-2 | Ref-2 DB %id | Ref-1 Ref-2 %id |
|:-----------:|:------------:|:-----------:|:------------:|:---------------:|
|  FJ233908   |     98.6%    |  HQ202949   |      100%    |       98.6%     |
|  FJ644562   |     98.6%    |  KY806003   |      100%    |       98.6%     |
|  MH509733   |     98.6%    |  LC008137   |      100%    |       98.6%     |
|  LC278346   |     98.6%    |  KX828216   |      100%    |       98.6%     |
|  HQ202950   |     98.6%    |  JF683403   |      100%    |       98.6%     |
|  KX641126   |     98.6%    |  KY656098   |      100%    |       98.6%     |
|  LC310740   |     98.7%    |  KX828216   |      100%    |       98.7%     |
|  KT867860   |     98.7%    |  KX510064   |      100%    |       98.7%     |
|  JN382189   |     98.7%    |  JN382187   |      100%    |       98.7%     |
|  AB512130   |     98.7%    |  KT868187   |      100%    |       98.7%     |
|  AB361574   |     99.1%    |  KX098762   |      100%    |       98.0%     |
|  AB512130   |     98.7%    |  JF683403   |      100%    |       98.0%     |
|  MT423827   |     98.7%    |  KX169329   |      100%    |       98.0%     |
|  LC310740   |     98.7%    |  JF683403   |      100%    |       97.9%     |
|  EU450613   |     98.3%    |  KX828216   |      100%    |       98.0%     |
|  KP081538   |     98.7%    |  KY806003   |      100%    |       98.0%     |
|  KR868575   |     98.9%    |  KX098762   |      100%    |       97.9%     |
|  KT868437   |     98.3%    |  KX828216   |      100%    |       98.0%     |
|  KX169308   |     98.3%    |  KY656098   |      100%    |       97.9%     |
|  KX641126   |     98.6%    |  MF314288   |      100%    |       97.9%     |
|  JF690912   |     98.7%    |  KX169329   |      100%    |       96.6%     |
|  MH509735   |     98.2%    |  LC004750   |      100%    |       96.5%     |
|  JN382189   |     98.7%    |  KP231118   |      100%    |       97.0%     |
|  KU697031   |     98.6%    |  JN382175   |      100%    |       96.4%     |
|  JF690912   |     98.7%    |  MF314288   |      100%    |       96.6%     |
|  HQ202950   |     98.6%    |  HQ202949   |      100%    |       96.9%     |
|  KP081538   |     98.7%    |  FJ804417   |      100%    |       96.7%     |
|  KP081541   |     98.4%    |  KP081542   |      100%    |       96.5%     |
|  KU697049   |     98.9%    |  MF142266   |      100%    |       96.9%     |
|  LC278346   |     98.6%    |  HQ202948   |      100%    |       96.6%     |
|  KT867952   |     98.9%    |  AB462385   |      98.9%   |       98.6%     |
|  KU697227   |     98.7%    |  DQ629127   |      98.7%   |       98.6%     |
|  KX298474   |     98.5%    |  KC514991   |      98.6%   |       98.6%     |
|  KY656102   |     98.7%    |  KP231128   |      98.4%   |       98.6%     |
|  LC278333   |     98.9%    |  KR868575   |      97.9%   |       98.6%     |
|  MF737379   |     98.9%    |  GQ358997   |      98.9%   |       98.6%     |
|  MT769305   |     98.9%    |  GQ358997   |      98.9%   |       98.6%     |
|  LC278333   |     98.9%    |  EU450613   |      98.0%   |       98.6%     |
|  KY656019   |     98.2%    |  KC514991   |      98.6%   |       98.9%     |
|  LC310740   |     98.7%    |  KR868575   |      97.9%   |       98.6%     |
|  MH465420   |     98.9%    |  EU755377   |      98.9%   |       98.0%     |
|  MG813260   |     98.4%    |  JQ866919   |      98.9%   |       97.9%     |
|  MK005838   |     98.3%    |  KM624030   |      98.9%   |       97.9%     |
|  KU697031   |     98.6%    |  JQ866919   |      98.9%   |       98.0%     |
|  KT868365   |     98.3%    |  KR868575   |      97.9%   |       98.0%     |
|  MK504383   |     98.4%    |  KR868575   |      97.9%   |       98.0%     |
|  KT867958   |     98.8%    |  KR868575   |      97.9%   |       98.0%     |
|  KX169298   |     98.8%    |  KP231128   |      98.4%   |       98.0%     |
|  LC278333   |     98.9%    |  HQ202951   |      98.4%   |       97.9%     |
|  LC278348   |     98.9%    |  AB512130   |      98.0%   |       97.9%     |
|  LC310732   |     98.3%    |  EU057188   |      98.6%   |       97.0%     |
|  LC278346   |     98.6%    |  AF465211   |      98.4%   |       96.7%     |
|  LC278346   |     98.6%    |  AB512130   |      98.0%   |       97.0%     |
|  LC278346   |     98.6%    |  KM924366   |      98.9%   |       97.2%     |
|  KX641126   |     98.6%    |  GQ358997   |      98.9%   |       97.0%     |
|  KX169298   |     98.7%    |  KP231145   |      98.4%   |       97.0%     |
|  KX169298   |     98.7%    |  AF201311   |      98.7%   |       97.4%     |
|  KT868374   |     98.3%    |  DQ629114   |      98.7%   |       97.6%     |
|  KY810321   |     98.4%    |  HQ202951   |      98.4%   |       97.2%     |
|  MK504383   |     98.4%    |  HQ202951   |      98.4%   |       97.4%     |
|  KY806003   |     100%     |  KC514972   |      100%    |       97.7%     |
|  KY806003   |     100%     |  EU755373   |      100%    |       97.7%     |
|  KY806003   |     100%     |  KC620537   |      100%    |       97.9%     |
|  KX828216   |     100%     |  JF683403   |      100%    |       97.9%     |
|  KT868187   |     100%     |  JF683403   |      100%    |       97.9%     |
|  KY940535   |     100%     |  HM038034   |      100%    |       98.0%     |
|  KX828216   |     100%     |  HQ202949   |      100%    |       98.0%     |
|  KT868187   |     100%     |  HQ202949   |      100%    |       98.0%     |
|  MH341484   |     100%     |  JX512855   |      100%    |       98.0%     |
|  MF142266   |     100%     |  JF683403   |      100%    |       98.2%     |
|  KY806032   |     100%     |  KC620537   |      100%    |       96.9%     |
|  KX510064   |     100%     |  HQ202949   |      100%    |       96.9%     |
|  KX169329   |     100%     |  KP081542   |      100%    |       96.9%     |
|  KX098762   |     100%     |  HQ202948   |      100%    |       96.9%     |
|  MN735211   |     100%     |  KJ094606   |      100%    |       96.9%     |
|  KX169329   |     100%     |  KJ139962   |      100%    |       97.0%     |
|  MH341484   |     100%     |  AY556477   |      100%    |       97.0%     |
|  MF589543   |     100%     |  EF524523   |      100%    |       97.0%     |
|  MF314288   |     100%     |  KC515014   |      100%    |       97.0%     |
|  MF314288   |     100%     |  JQ181595   |      100%    |       97.0%     |

Table:
  References used in percent similarity testing.
  Reference-1 and Reference-2 are the two references used in reference
    pair.
  Ref-1 DB %id: is the % similarity between reference-1 and its most
    similar reference in the database.
  Ref-2 DB %id: is the % similarity between reference-2 and its most
    similar reference in the database.
  Ref-1 Ref-2 %id: Is percent similarity between reference-1 and
    reference-2

### Porcine circovirus type 2 genotype data set:

This dataset has reads from different genotypes mixed together. Their
  were two references used to simulate each genotype. In all cases only
  one of the two references was used with each genotype. The reference
  chosen in a genotype pair was chosen on a whim to ensure some
  level of randomness. I also make a reference pair with both references
  from the same genotype. The labeling for these reference pairs is a 
  little different, in that I used the genotype instead of reference in
  the file name.

This data set has four different major/minor strain ratios. The first
  is 50/50 (50% of reads from the major strain), then 75/25
  (75% of reads from the major strain), then 95/5 (95% of reads from the
  major strain), and finally 99/1 (99% of reads from the major strain).

Badread commands for this test:
   - badread simulate --quantity 40000x --seed 1026 --reference file.fasta
   - badread simulate --quantity 15000x --seed 1026 --reference file.fasta
   - badread simulate --quantity 4000x --seed 1026 --reference file.fasta

| Major Variant | Major Genotype | Minor Variant | Minor Genotype |
|:-------------:|:--------------:|:-------------:|:--------------:|
|   HM038034    |      a         |   KX828216    |       a        |
|   HM038034    |      a         |   KU697268    |       b        |
|   KX828216    |      a         |   KC515014    |       d        |
|   HM038034    |      a         |   LC008137    |       f        |
|   KY806003    |      b         |   KU697268    |       b        |
|   KY806003    |      b         |   KC515014    |       d        |
|   KY806003    |      b         |   LC004750    |       f        |
|   KY806003    |      b         |   KP420197    |       g        |
|   KC515014    |      d         |   MN735211    |       d        |
|   KC515014    |      d         |   LC008137    |       f        |
|   KC515014    |      d         |   JX512856    |       g        |
|   LC004750    |      f         |   LC008137    |       f        |
|   LC004750    |      f         |   JX512856    |       g        |
|   KP420197    |      g         |   JX512856    |       g        |

Table:
  List of reference pairs used in genotype co-infections test.

## Databases

The references I used to detect co-infections can be found in dataBases.
  Each database is a fasta file with multiple references used to bin
  the simulated reads. For the porcine circovirus type 2 (PCV2)
  databases each reference came from open reading frame 2 (encodes th
  capsid protein) and is around 700 bp long. For the Hepatitis C
  database each reference targets the first 1260 bases of the E2 gene.

### Porcine circovirus type 2 similarity database:

PCV2--2-percent-different-database.fasta has my complete porcine
  circovirus type 2 database, which has 66 references. Each reference
  is at least 2% different from all other references in the database.
  This database is a bit out of date, but will probably work ok for
  real studies.

| Accession | Accession | Accession | Accession | Accession | Accession |
|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
| AY556477  | HQ202944  | JX512856  | KP081547  | KX169329  | MF142266  |
| AY864814  | HQ202948  | KC514972  | KP081548  | KX510057  | MF314288  |
| DQ856569  | HQ202949  | KC515014  | KP081549  | KX510064  | MF589543  |
| EF524518  | HQ591379  | KC620532  | KP081550  | KX828216  | MH341484  |
| EF524523  | JF690911  | KC620537  | KP081551  | KY656098  | MH465430  |
| EF619037  | JN133304  | KJ094603  | KP231118  | KY806003  | MH465431  |
| EU057186  | JF927978  | KJ094606  | KP231129  | KY806032  | MK005850  |
| EU136711  | JN382175  | KJ139962  | KT819159  | KY947570  | MN735211  |
| EU755373  | JN382187  | KP081542  | KT868187  | KY940535  |           |
| FJ804417  | JQ181595  | KP081543  | KX098762  | LC004750  |           |
| HM038034  | JX512855  | KP081546  | KX169322  | LC008137  |           |

Table:
  List of references used to make the 2% reference database.

### Porcine circovirus type 2 genotype database:

PCV2--sparse--database.fasta has my sparse data base, which consists of
  four references. I tried to match each reference to a different
  genotype, however, I am not sure how well I did. I can at least say
  genotype a and b are represented.

| Reference | Minimum Difference | Maximum Difference | Genotype |
|:---------:|:------------------:|:------------------:|:--------:|
| HQ591379  |        2.85%       |          9.40%     |   b/f    |
| HQ202944  |        2.28%       |          9.22%     |    a     |
| KX169298  |        1.28%       |          9.26%     |    d     |
| KJ094600  |        3.69%       |          9.36%     |   f/b    |

Table:
  Reference set used to detect co-infections in sparse reference set
    test.
  Minimum Difference is the closest reference used in the similarity
    test and genotype test.
  Maximum Difference is most different reference used in the similarity
    test and genotype test.
  Genotype is the closest genotype to the reference.

### Hepatitis C database:

HCV--2-percent-different--database.fasta has my database for
  benchmarking my Hepatitis C tests. Each reference is at least two
  percent different than any other reference in the data base. This is
  not a great database and is missing much of the Hepatitis C diversity.
  I would not advice using this database for Hepatitis C, studies.

| Reference | Reference | Reference |
|:---------:|:---------:|:---------:|
| GQ370316  | GQ370176  | GQ370322  |
| GQ370239  | GQ370167  | GQ370321  |
| GQ370238  | GQ370326  | GQ370142  |
| GQ370132  | GQ370325  | GQ370136  |
| GQ370161  | GQ370324  | GQ370134  |

Table:
  Genbank accession numbers for each reference in our reference database
    used in HCV testing.

## Scripts

R-scripts-and-data has the scripts I used to make the graphs for version
  three of find co-infections and the data sheet I used with the script.

benchmarking-scripts has the scripts I used for benchmarking. It is a
  a mess.
