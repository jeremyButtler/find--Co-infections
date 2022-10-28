# Use:

Extracts reads from fastq file using read id.

# Reason here:

fastqGrep reduces the number of dependencies needed for findCoInfections.
  However, fastqGrep is also slightly slower than seqkit. This is not a problem
  because, the clustering step and initial read mapping steps will always take
  much larger amounts of time then the read extraction steps

The regular version of fastqGrep converts the hex characters in read ids to
  large numbers. However, in the process it ignores and letter that is not
  0 to 9 or a to f. This works well for nanopore sequence reads, which have
  read id's in hex and somewhat well for Illuminia, which only has a few non-hex
  characters for machine and flow cell ids.

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
   around, but not quite O(nlog(n)).

The string version (default) finds the hash by multiplying all characters in a
  c-string and then taking the modulo of the array size. While the default
  version reverses the c-string, discards any non-hex digits
  (not 0 to 9 or a to f) and converts the c-string to a big number. The hash is
  then found using Knuths multiplicate hash on the 64 most significant bits of
  the big number. The non-hashed big number is also used to detect if the read
  id matches any input read ids.
