# Use

This is an fqGetIds that include ":" and the letters "A to V" in its big
  number conversion step. This should reduce some of the potential
  issues that could in theory happen with Illumina, but will also
  increase the memory usage (it is storing more characters from the
  name). To include U to Z would require a 6 bit system, which would
  require even more memory usage.

For some odd reason this version is a bit slower for extracting a few
  reads from a deep fastq file with Nanopore reads. It should not be and
  I am not sure what is happening.

```
make
cp fqGetIdsIll /path/to/install/fqGetIdsIll
chmod a+x /path/to/install/fqGetIdsIll

fqGetIdsIll -f extract.ids -fastq reads.fasta > out.fastq
```


