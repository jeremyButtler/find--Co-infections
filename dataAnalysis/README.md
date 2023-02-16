## Purpose

This directory has the data, scripts, reference pairs, and references I
  used to benchmark find co-infections version three.

The references I used to simulate co-infections involving similar viral
  strains can be found in similarityTestReferencePairs. I will be adding
  a folder for each test I ran, so expect to see more as more data comes
  out.

The references I used to detect co-infections can be found in dataBases.
  These are the references I ran with find co-infections.
  PCV2SparseCp.fasta has one reference from four different genotypes of
  PCV2 and was the v3-sparse database (same for V2 figures).
  PCV2TwoPercentDifferenceCp.fasta has my complete database, which has
  66 references. Each reference is at least 2% different from all other
  references in the database.

RScriptsAndData has the scripts I used to make the graphs for version
  three of find co-infections and the data sheet I used with the script.

benchmarkingScripts has the scripts I used for benchmarking. It is a bit
  of a mess.
