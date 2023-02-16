## Overview

Here is a quick guide to how to run some of the scripts in these
  directories. These scripts are here so that you can see how I tested
  my data and replicate my results.

## Find Co-infection V3 Testing:

To test version three of my pipeline you will want to run either
  runV3Test.sh, runSingleV3Test.sh or runFqV3Test.sh. All of these
  pipelines will require fqGetIds to run. runFqV3Test also requires
  scoreReads to run.

  ```
  cd ../V3/;
  make fqGetIds;
  make scoreReads;
  mv fqGetIds scoreReads ../dataAnalysis;
  chmod a+x ../dataAnalysis/fqGetIds ../dataAnalysis/scoreReads;
  ```

  - runV3Test.sh:
    - Run:
      - bash runV3Test fasta-directory database.fasta prefix read-depth
    - takes in a directory of fasta files and simulates
      reads reads for each fasta file in the directory. It then runs
      find co-infections version three on the simulated file and gets
      the stats. After running find co-infections it deletes the
      simulated reads.
  - runSingleV3Test.sh:
    - Single run of runV3Test. Type runV3Test.sh for a help message.
  - runFqV3Test:
    - Used to run Version three with fastq files simulated using
      twoRefRunBadread.sh
      (bash twoRefRunBadread.sh -h for more information)
    - bash runFqV3Test -h for more information
    - Gives a bit more data than previous tests.
    - Note the reference directory must have the name of both references
      the fasta file. See similarityTestReferencePairs for some
      examples.

## Other pipelines:

Other options include runClairTest.sh for running a docker install of
  Clair three, runLongShotTest.sh to benchmark Longshot, and
  runNanoQTest.sh to benchmark NanoQ.

  - For Clair 3:
    - bash runClairTest.sh reference.sh reference-database.fasta read-depth
    - Clair three must be installed as hkubal/clair3:latest or
      you will have to change line 79 of runClair3.sh.
      - hkubal was a misspelling on my part when installing Clair 3.
  - For Longshot:
    - bash runLongShotTest.sh reference.sh reference-database.fasta read-depth
  - For NanoQ:
    - bash runNanoQ.sh reference.sh reference-database.fasta read-depth
    - Make sure Nano-Q is installed at ~/Nano-Q or /usr/local/bin/Nano-Q
    - You will need to also need the updates I made for NanoQ. See
      jeremybuttler/Nano-Q--wrapper/Nano-Q--upgrades
    - You will also need Racon for this script. This was my way of
      getting around Nano-Q sometimes building consensuses having only
      N's.

## V2 and V1

Their are scripts in this directory to run version one and two of find
  co-infections (runTest.sh and runSingleTest.sh). You will need Medaka
  and Racon installed on your machine. For Medaka you may need to
  install it by both miniconda and python virtual environments. I am not
  sure if V1 used conda or python env.

```
Set up for benchmarking version one and two of find co-infections:
   cd ../V2;
   make;
   cd ../;
   cp -r V2 dataAnalysis/findCoInfections;
   chmod -R a+x dataAnalysis/findCoInfections;

   cp -r V1 dataAnalysis/oldFindCoInfections;
```
    

  - runTest.sh:
    - bash runTest.sh reference-directory reference-database.fasta prefix
    - For read depth: set readDepthInt (line 29) to the target depth.
    - This will run V1, V2, and build a consensus from all reads.
  - runSingleTest.sh:
    - Run bash runSingleTest.sh -h for more information.
    - This will run V1, V2, and build a consensus from all reads.

## For Detecting Duplicates:

You can run detectDupRef.awk to detect duplicates in the data output by
  all scripts, except for the V3 scripts. Their are some issues with
  column alignments their.

  - awk -f detectDupRef.awk data.tsv
