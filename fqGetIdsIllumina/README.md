# Use

This is an fqGetIds that include ": & G to V" in its big number
  conversion step. This should reduce some of the potential issues that
  could in theory happen with Illumina. For some odd reason this version
  is a bit slower for extracting a few reads from a deep fastq file with
  Nanopore reads. It should not be, so I am not sure what is happening.

```
make
cp fqGetIdsIll /path/to/install/fqGetIdsIll
chmod a+x /path/to/install/fqGetIdsIll

fqGetIdsIll -f extract.ids -fastq reads.fasta > out.fastq
```


