# Use:

Extracts reads from fastq file using read id.

# Reason here:

fastqGrep reduces the number of dependencies needed for findCoInfections.
  However, fastqGrep is also slightly slower than seqkit. This is not a problem
  because, the clustering step and initial read mapping steps will always take
  much larger amounts of time then the read extraction steps
  (under one minute. See comparison of gmp fastqGrep and seqkit in
  libGMPVersioN folder). So, taking slightly more time then seqkit will have
  little effect on noticed pipeline performance.

The libGMPVersioN of fastqGrep is faster then the default version, though still
  slower then seqkit. However, it also requires the dev files for the gmp
  library and thus, is not the default install. I am only using the gmp library
  for string to number conversions and large number comparisons, which could
  easily be done without gmp. I am planning to remove this dependency later in
  the next version of findCoInfections, which I just started working on.

If you wish to use the libgmp version install the libgmp dev files for your OS.
  For debain this can be done with apt-get install libgmp-dev. Then move the 
  contents of libGMPVersion to this directory.

# How to use:

1. Providing input by file
  - fastqGrep -f idList.txt -fastq reads.fastq > out.fastq
2. Providing read id list by stdin
  - fastqGrep -stdin-filt -fastq reads.fastq > out.fastq
3. Providing fastq file by stdin
  - fastqGrep -f idList.txt -stdin-fastq > out.fastq
4. Other paramaters:
  - -v to extract reads not in idList.txt (same parameter for seqkit and grep)
  - -V print version number

# How to install:

make;
mv fastqGrep /path/to/install;
chmod a+x /path/to/install;

# How fastqGrep works:

fastqGrep uses a hash table combined with a balanced tree (using AVL method) to
   identify matching read ids. The hash table reduces the best case look up time
   to O(1), while the balanced tree ensures that the worst case look up time is
   around, but not quite O(nlog(n)). This reduces the impact of collisions in 
   hashing and is likely why the gmp version of the pipeline is faster then
   seqkit when extracting all reads (see if figures in libGMPVersion folder).

The string version (default) finds the hash by multiplying all characters in a
  c-string and then taking the modulo of the array size. While the gmp version
  flips the c-string, discards any non-hex digits (not 0-9 or a-f) and converts
  the c-string to a big number. The hash is found using Knuths multiplicate hash
  on the 64 most significant bits of the big number. The non-hashed big number
  is also used to detect if the read id matches any input read ids.
