# Use:

Finds co-infections in nanopore sequenced reads.

## Requirements:

  - Required
    1. GCC (to compile findCoInft).
    2. minimap2 (https://github.com/lh3/minimap2)
  - Optional:
    1. racon (https://github.com/isovic/racon)
    2. medaka (https://github.com/nanoporetech/medaka)
      - Can be installed by conda or python virtual enviroments
      - For python virtual env install; the folder with medaka needs to
        be located in the users home directory (~/medaka).

## Install:

cd V3;
make;
mv findCoInft /path/to/install;
chmod a+x /path/to/install/findCoInft;

## Run:

findCoInft -fastq reads.fastq -ref refferences.fasta [other options...]

## Some quick options:
  - -h:
    - Show the main help message (their are multiple)
  - -prefix:                                                     [Out]
    - Prefix to add to file names.
  - -threads:                                                    [3]
    - Number of threads to use with Minimap2/Racon/Medaka.
  - -min-reads-per-bin:                                          [100]
    - Minimum number of reads to keep a cluster or bin.
  - -max-reads-per-con:                                          [300]
    - Max number of of reads to use in building a consensus.
  - -enable-racon:                                               [No]
    - Use Racon in building consensuses.
  - -enable-medaka:                                              [No]
    - Use Medaka in building consensuses.
  - -model:                                         [r941_min_high_g351]
    - Model to use with Medaka_consensus

## Testing:

Currently benchmarking findCoInfections version three, and will update
  this section as data comes in. For version two testing results, see
  README.md in V2.

## OS's/platforms tested on

  - OpenBSD (no Medaka)
  - Debian
  - Raspberry pi 1B with a 32 bit Debian install (no Medaka)
    - The independent install of fqGetIds did crash on this, but
      findCoInft ran fine.
    - I should fix this issue at some point.

## Future directions:

1. Add in full consensus building step to version three and maybe some
   masking for poorly supported bases. The majority consensus step will
   remove these bases, so I am not to worried about this.
2. Put consensus building steps in separate file, so can make it into a
   separate C program.
3. Maybe multi-thread fqGetIds (in V2 was fastqGrep). Likely will not
   happen.
