# Reason for this seemingly duplicate code

This version of fqGetIds is here because the version complied with
  find co-infections V3 is slow when extracting a small percentage of
  reads Nanopore reads from a very deep fastq file (6 million reads) for
  some odd reason. I would reduce down the extra files by adding in the
  one function it uses in fqAndFaFun.c, but found that it actually
  increased the time usage when extracting a low percentage of Nanopore
  reads from a large (6 million read) fast file somehow.

At this point I am not sure if it is just my imagination or something 
  deeper with the version of GCC I used, were the extra files causes
  some other optimization to kick in. Either way I do not know GCC flags
  enough to figure this out, so you have multiple different versions.

This is the version of fqGetIds I used for benchmarking.

```
make
mv fqGetIds /place/to/install/fqGetIdsFast
chmod a+x /place/to/install/fqGetIdsFast

fqGetIds -f filter.ids -fastq reads.fastq
```

Their are is some more information and graphs in the supplemental
  programs.
