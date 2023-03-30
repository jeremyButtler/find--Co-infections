AlignSeq has a Needleman Wunsch and Waterman Smith alignment algorithm
  that I was originally hoping to use to improve my majority consensus
  step. This did not work out very well, so I scrapped the idea, but
  kept the alignment programs here in case any one finds them useful. I
  was hoping to code in a Hirschberg, but at this point I have no reason
  to try.

Both alignments uses a two bit direction matrix, which results in at
  least 4x less memory than a traditional Needleman Wunsch. During
  scoring it decides the direction to travel back and ignores any other
  potentially equal paths. However, this is only compared to other
  Needleman Wunsch and Waterman Smith alignments. When this is compared
  to a Hirschberg, which uses a linear instead of N^2 memory usage, it
  is very memory heavy and slow.

AlignSeq can also do a Waterman Smith alignment to get a single answer.
  However this Waterman Smith only keeps track of a single best entry
  (the one that is nearest to the bottom right of the matrix), so it
  will not give multiple alignments. A Hirschberg is a better option.
  One way to make this deal with multiple sub alignments would be to
  record the position of the best score per base, start of with the best
  score and then explore the paths of any unmapped bases. It would take
  (n + m) * constant amount of memory to do, but would allow for more
  sub alignments. This is not currently in the algorithm an I am
  mentioning this here to record the thought, in case I or any one else
  will ever need it.

AlignSeq is only memory efficient when compared to other Needleman
  Wunsch or Waterman Smith alignments ((N^2)/4 instead of N^2). However,
  it will never be as efficient has a Hirschberg alignment, which
  operates in linear space (less than 5mb to align two 200kb sequences).

AlignSeq does not use decimals, so if you want decimals for the gap
  extension penatly you will have to multiply all scores by 10.

Note: This program is more set up for noisy reads. This means that the
  gap opening has a lighter penalty than the gap extenstion, which means
  that this program favors single gaps instead of longer gaps. For a
  reference or consensus alignment you will want to set the gap open
  penalty higher and lower the gap extension penalty.

The help message can be called with -h. For an exmaple of the scoring
  matrix see V3/scoring-matrix.txt.

#### How to build alignSeq

```
cd V3
make alignSeq
mv alignSeq /path/to/install
chmod a+x /path/to/install/alignSeq

# For a global alignment (Needleman Wunsch)
alignSeq -query query.fasta -ref ref.fasta > alignment.aln

# For a single local alignment (Waterman Smith)
alignSeq -use-water -query query.fasta -ref ref.fasta > alignment.aln
```
