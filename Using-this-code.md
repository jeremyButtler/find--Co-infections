# Use

This document is here to give you an idea of how you could use various
  parts of this code for you own project. In it I go over my
  documentation method and data structures and functions that might be
  of interest. I will not cover everything.

## The defualt settings

Most of the default settings, include the commands that run Minimap2,
  Racon, and Medaka can be found in defaultSettings.h.
  
## The samEntry structer and its functions

The samEntry structure is a structure you will often seen throughout my
  code. It has a buffer named sameEntryCStr, which is a c-string stored
  on the heap for reading a line from a sam file. Additional pointers
  (reference id, query id, cigar entry, sequence entry, q-score entry
  all point to the buffer. It also includes variables to hold stats or
  are used to find stats, such as Q-scores (values, histograms, and
  totals), lengths, mapping qualities, flags, position, and error
  counts. See samEntryStruct.h for a list of variables to manipulate.

Some useful functions include:

- blankSamEntry
    - This will set all non-buffer pointers to 0 and set all stats to
      0 (this includes the histograms).
- blankReadStats
   - This is like blankSamEntry, but will not blank any pointers.
- initSamEntry
   - Initializes a samEntry structure. It should only be called once.
       - Sets buffer to 0
       - Sets buffer length to 0
       - Calls blankSamEntry
    - Input: pointer to samEntry structure
- freeHeapSamEntry
    - Frees a samEntry structure that was allocated on the heap
    - Input: pointer to samEntry structure
- freeStackSamEntry
    - Frees a samEntry structure that was on the stack
    - Input: pointer to samEntry structure
- readSamLine
    - Reads in a single line from a sam file
    - This will resize the buffer (samEntryCStr) when the buffer is to
      small to hold the entire line.
    - Input: pointer to samEntry structure
    - Input: pointer to sam file
- printSamStats
    - Prints the stats in a sam file to tsv file
    - Does not print flag or position
    - Input: pointer to samEntry structure
    - Input: pointer to tsv file to print to
- samToFq
    - Prints sam entry to a fastq file as a fastq entry
    - Input: pointer to samEntry structure
    - Input: pointer to fastq file to print to

## The readStat structer

The readStat structure holds the stats printed out a samEntry structure,
  but it does not include the Q-score histograms or buffer. It is here
  to read in files that were printed out from the samEntry structure.

Some useful functions include:

- blankReadStat
   - Set all variables in readStat to 0
   - Input: pointer to readStat structure
- readStatsFileLine
   - Read in a single line from the tsv file made using a samEntry
     structure.
   - Input: pointer to readStat structure
   - Input: pointer to tsv file to read from
- printReadStat
   - Print out the stats from a readStat structure.
   - Input: pointer to readStat structure
   - Input: pointer to tsv file to print to
- samEntryToReadStat
   - Copy stats from a samEntry structure to a readStat structure.
   - Input: pointer to readStat structure to copy to
   - Input: pointer to samEntry structure to copy from

## Some additional structures

Some additional structures you might encounter across multiple functions
  are the minAlnStats and readBin structures. The minAlnStats structure
  holds the minimum thresholds to keep a mapping and can be found in
  minAlnStatsStruct.c/h. It is initiated with a blankMinStats function
  for read to reference mappings, blankMinStatsReadRead for read to 
  read mappings, blankMinStatsReadCon for read to consensus mappings,
  and blankMinStatsConCon for consensus to consensus comparison
  settings.

The readBin structure (findCoInftBinTree) is for managing the files the
  go with a single bin or cluster. It holds the reference id (I often do
  not update this for clusters), fastq file name, name of the stats file
  from score reads, the name of the best read, the name of the fastq
  file with the top reads, and the consensus name. It also stores a
  count of the number of reads in the bin or cluster. It is made with
  makeReadBin and freed with freeReadBin.

You can free a tree or list of readBin's with freeBinTree, convert a
  readBin tree to a list call cnvtBinTreeToList, merge readBins with
  mergeBins, delete the files in a bin with binDeleteFiles, or remove a
  bin from a linked list (not a tree) with rmBinFromList. These
  functions are at the bottom of findCoInftBinTree.c/h.

## Some useful functions that go with the samEntry structure

- trimSamEntry (in trimSam.c)
    - Trims of the soft masked regions of a sam file.
    - Note: this will not trim secondary and unmapped alignments.
    - Returns: 0 if trimmed the entry
    - Input: pointer to samEntry structure with sequence to trim.
- FindQScores (in FCIStatsFun.c)
    - Finds median and mean Q-scores for the sam file entry in a
      samEntry structure.
    - Input: pointer to a samEntry structure
- readRefFqSeq: (fqAndFaFun.c/h)
    - Reads in a sequence from a fastq file into a samEntry structure.
      This is used for getting a reference for scoreAln (next function).
    - Input: fastq file to get read from
    - Input: samEntry structure to hold the reference sequence and
      Q-score
    - Input: 1 to not remove spaces and new lines from the sequence and
      q-score entries, 0 to remove extra spaces and new lines from the
      sequence and q-score entries. For non-space removal, the length
      is found by the number of read in characters, so extra newlines
      will through this off.
- scoreAln:
    - Gets the stats for an sam alignment. This only gets the aligned
      median and aligned mean Q-scores.
    - Required that your sam entry have an --eqx cigar entry.
    - Input: minAlnStats structure, with min base Q-scores and maximum
      hompolymer lengths. (see minAlnStatsStruct.h for the structure).
    - Input: samEntry structure with mapping to score.
    - Input: samEntry structure with reference sequence (0 to ignore the
      reference).
    - Input: 1 if the reference has a q-score entry, otherwise 0.
    - Input: 1 if you want to use the reference only for deletions,
      otherwise 0.
- filterReads (fqAndFaFun.c/h)
    - Filters out reads with low Q-scores and to small or to long
      lengths.
    - Makes a new fastq file with filtered reads.
    - Input: c-string with fastq file name to filter
    - Input: c-string with new fastq file name to hold filtered reads
    - Input: samEntry structure to work with
    - Input: minAlnStats with minimum length, maximum length, minimum
      mean Q-score, and minimum median Q-score.
- binReads (binReadsFun.c/h)
    - Bins reads by reference (calls minimap2). Each read is binned to
      a file named from the input prefix and the assigned reference
      name.
    - You can expect a fastq file and a stats file from scoreReads for
      each bin.
    - Returns: An AVL tree of readBin structures (has reference name and
      binned fastq file name).
    - Input: Name of fastq file with reads to bin
    - Input: Name of file with references to bin reads to
    - Input: Prefix to add to the name of each files name
    - Input: Number of threads to have minimap2 use (c-string)
    - Input: 1: remove any read with supplementary alignments; 0 do not
        - Note: The actual supplemental alignments are discared either
          way. This just determines if the primary alignment is kept.
    - Input: 1: trim the reads to the reference; 0 do not
    - Input: samEntry structure to work with
    - Input: samEntry structure to work with when removing any read with
      supplementary alignments.
    - Input: minAlnStats structure with minimum stats to keep a read.
      Almost all stats are used (Think aligned stats are the only ones
      not used).
    - Input: Variable to hold the error report (2 = open FILE error,
      4 & 8 = write file error, 64 for memory error).

## Consensus building

The consensus building functions are in buildConFun.h. The overall
  method is called by buildCon. Single rounds of consensus building are
  called by buildSingleCon. The individual consensus building functions
  are simpleMajCon, buildConWithRacon, and medakaPolish. Each of the 
  individual functions have their own structure for settings and storing
  the output consensuses length (majConStruct, raconStruct, and
  medakaStruct). These structures are stored in the conBuildStruct
  structure, which also has the more global settings for buildCon.

You can set up default settings for the conBuildStruct by calling
  initConBuildStruct. You can turn a consensus building method on or
  of by setting the useXBl in majConSet (majconStruct),
  raconSet (raconStruct), or medakaSet (medakaStruct) structures of the
  conBuildStruct.

- buildCon:
    - Builds a multi-round consensus using the specified consensus
      building methods.
    - Warning this function will change the input fastq file, so make
      sure it is a copy (copyFile in fqandFaFun.c will do this).
    - Input: readBin structure with the reads to build a consensus for.
    - Input: a reference sequence to start out with or 0 for finding
      the best read.
    - Input: Number of threads to run minimap2 with (as c-string)
    - Input: conBuildStruct with the consensus building settings.
    - Input: samEntry structure for reading in a sequence
    - Input: samEntry structure for extracting reads
    - Input: minAlnStats, minimum requirements to keep a read when
      hunting for the best reads (read to read mapping).
    - Input: minAlnStats, minimum requirements to keep a read when
      hunting for the best reads (read to consensus mapping).

## fqGetIds

This is not a detailed list and more can be found by browsing the files,
  however, this should give you a vague if of how to interface with
  this.

- fqGetIdsSearch.c has the driver functions for fqGetIds
- fqGetIdsFqFun.c and fqAndFaFun.c have the fastq functions for fqGetIds
- fqGetIdsAVLTree.c has the functions for managing the tree
- fqGetIdsStructs.c has the functions for interacting with the
  structures.

fqGetIds is a bit more difficult to work with and involves keeping track
  of several smaller variables. These include the number used in hashing
  (majicNumULng), the array size as a power of two (digPerKeyUChar),
  a hash table, the size of the hash table, and a stack of readInfo
  structures (used with readNodeStack). It uses a structure called
  readInfo, which has a bigNum structure with the read id and the
  pointers to turn it into an AVL tree.

The first step is to build a linked list of readInfo nodes that are only
  connected by the readInfo->rightChild pointer. Each readInfo structure
  can be made by calling makeReadInfoStruct(read Id, length). This will
  make a readInfo structure with the input read id as a big number. This
  structure can be freed by calling freeReadInfoStruct(&readInfo).

Your readInfo stack can be created by the commands bellow. Once created
  you only need to pass the stack around, but do not need to deal with
  it (including freeing). I am using it in this fashion to avoid speed
  issues of creating it each time.

```
struct readNodeStack stack[256];
stack[0].readNode = 0;
stack[255].readNode = 0;
```

You can also write new read ids to an existing bigNum structure by
  calling  strToBackwardsBigNum(bigNum pointer, new ID, &lenght of ID).
  This structure can then be freed with
  freeBigNumStruct(&bigNum pointer).

After you list is built you can then build the hash table using
  readListToHash (See fun-6 in fqGetIdsHash.c). To search the hash table
  call
  findReadInHashTbl(bigNum pointer, &majic number, digPerKeyUChar, hash table).
  To free the hash table call
  freeHashTbl(& hash table, size of table, stack)

  For examples of using fqGetIds on a fastq file see the extractReads
  function (fun-3 fqGetIdsSearch.c). The fastq file functions are
  located in fqGetIdsFqFun.c and fun-7 in (moveToNextFastqEntyr) in
  fqAndFqFun.c

## Misc

These functions probably have little use, but might of use in some
  special cases

- getNumReadsInFq (fqAndFaFun.c/h)
    - Find the number of reads in a fastq file.
    - Input: Pointer to fastq file (FILE).

## AlignSeq

This is the documentation for alignSeq, which uses a two bit
  Needleman Wunsch alignment algorithm. Most of the functions can be
  found in alignments.c/h. It also requires cStrToNumberFun.c/h and
  defaultSettings.h (sec-9 holds the gap penalties and scoring matrix
  settings).
   
   - Uses the alnSet structure for storing alignment settings.
       - gapStartPenaltyI: penatly for opening new gaps.
       - gapExtendPenaltyI: penalty for extending a gap
       - It also has the scoring matrix and variables that set the
         priorities for alignments
   - initAlnSet:
       - Initalizes the alnSet structure to default values.
       - Input pointer to an alnSet structure to initalize.
   - readInScoreFile:
       - Takes in a score matrix file (V3/scoring-matrix.txt) and
         updates the scores.
       - Input: Pointer to alnSet structures to add scores to.
       - Input: scoreFILE, FILE pointer to file to get scores from.
   - setBasePairScore
       - Sets the alignment matrix score for a single pair of bases
       - Input: Pointer to the query base
       - Input: Pointer to the reference base
       - Input: Score to add to the scoring matrix
       - Input: Pointer to alnSet structer with the scoring matrix
   - WatermanSmithAln
       - Performs a Waterman Smith alignment on a pair of sequences.
       - This version only finds a single answer and so is not really
         that great (a Hirschberg Waterman Smith is much better).
       - The only difference for input and output from the
         NeedleManWunschAln (next entry) is that the alignment array
         also has soft mask flags (defSoftRefFlag, defSoftQueryFlag, or
         defSoftRefFlag + defSoftQueryFlag).
   - NeedleManWunschAln
       - Performs a Needleman Wunsch alignment on a pair of sequences
       - Returns: an alignment array with flags marking if the base is
                  an match, snp, insertion or deltion
           - This array needs to be freeded (on heap)
           - Use: defBaseFlag (snp/match), defInsFlag, defDelFlag to
             identify the error type.
       - Input: c-string with the query sequence to align
       - Input: Starting point on the query sequence (index 1)
       - Input: Ending point on the query sequence (index 1)
       - Input: c-string with the reference sequence to align
       - Input: Starting point on the reference sequence (index 1)
       - Input: Ending point on the query reference (index 1)
       - Input: Pointer to alnSet structure with settings
       - Input: lenErrAryUI, will hold the length of the alignment array
       - Input: Score of the alignment
   - cnvtAlnErrToSeq:
       - Coverts an alignment error array to an aligned sequence
       - The aligned sequence is stored on heap and needs to be freeded
       - Input: c-string with sequence
       - Input: starting position on sequence
       - Input: 1: is a query sequence, 0 is a reference sequence
       - Input: error array (from NeedleManWunschAln) with alignment
       - lenErrAryUI: length of the error array
   - cnvtAlnErrAryToLetter
       - Converts an error array to have a human readable format (I is
         insertion, = is match, X is SNP, D is deletion, s is a query
         soft mask, P is a reference soft mask, and S is both reference
         and query soft mask). This will only count matches if a
         reference and query sequence is provided, otherise it will
         print X for SNPs and matches.
       - Input: Reference sequence used in the alignment.
         Used to detect matches (use 0 to ignore matches).
       - Input: Query sequence used in the alignment.
         Used to detect matches (use 0 to ignore matches).
       - Input: error Array (from NeedleManWunschAln) to convert to
         human readable format.
   - readFaSeq (from sequence.c)
       - Reads in a fasta sequence and modifies input arrays to hold the
         input.
       - Arrays will be reallocated as needed (so plan on putting them
         on the heap and freeing them as needed).
       - Input: FILE pointer to fasta file to read sequneces from
       - Input: buffer to hold the header entry in the fasta file
       - Input: length of header entry (0 will allocate more memory)
       - Input: buffer to hold the sequence entry in the fasta file
       - Input: length of sequence entry (0 will allocate more memory)
       - Input: Variable to hold the length of the inptu sequence
   - twoBit functions:
       - Found in twoBirArrays.c/h
       - These functions are here for working with two bit arrays.

## Primer trimming

The primer trimming functions are from the trimPrimers\*.c/h files and
  also take a lot of code from fqGetIds. The trimPrimers function, which
  does the work can be found in trimPrimersSearch.c/h. TrimPrimers calls
  makeReadPrimList to get the primer cooridinates (readPrimHash.c/h),
  readPrimListToHash (trimPrimerHash.c/h) to make a hash table or
  readPrimListToTree (trimPrimersAvlTree.c/h) to make an AVL tree, 
  extractAndTrimReads (trimPrimersSearch.c/h) to extract the reads, 
  and freeReadPrimHashST (trimPrimersHash.c/h) to free the tree/hash
  table.

I am only going to go over trimPrimers, but mentioned the other
  functions if you want to build your own driver functions.

-trimPrimers:
    - Trims the primer regions off reads in a fastq file and splits each
      read into multiple reads if needed.
    - Input: Fasta file name with primers to trim from reads (only used
      if their was no input paf file).
    - Input: paf file name with primer coordinates (0 to skip)
    - Input: 1; paf file is from stdin; 0 is from a file name
    - Input: Fastq file name with reads to trim
    - Input: File name of output file ('\0' for stdout)
    - Input: Number of threads to use (c-string)
    - Input: 1;Use hashing in the read extraction; 0 do not use hashing.

## read extraction

the readExtract.c/h holds the functions to extract a targeted number of
  reads from a fastq file. Most of the functions are for single reads,
  however, findBestXReads is designed to extract the top X best reads
  from a file.

- findBestXReads:
    - Extracts the top X best reads from a file using the mapping
      quality, then integer median Q-score, and then length to select
      reads.
    - Input: Number of reads to extract
    - Input: Variable to hold the number of reads extracted
    - Input: Number of threads for minimap2 (as c-string)
    - Input: 1: select best reads using mapping quality; 0 ignore
      mapping quality.
    - Input: minAlnStats structure holding the stats needed to even
      keep a read. All minimum percentage variables are used as well as
      minimum length and mapping quality are used for mapping qualities.
      Median and mean Q-score replace mapping qualities when mapping
      quality is ignored.
    - Input: samEntry structure for reading in alignments or sequences
    - Input: samEntry structure to hold a fastq reference (input 0 to
      ignore this)
    - Input: readBin structure with fastq file to extract reads from
    - Input: 1: do not use a reference (defaults to Q-scores); 0 use the
      reference to map reads.
    - Input: 1; build the name for the output fastq from the input fastq 
      file; 0 Use the top reads file for the name. (Use 1).
      
# Documentation method

My documentation style may not be the best style, but if you want to 
  look though my code you will need to have an idea of how I document
  my code. It has evolved since I first stared find co-infections and
  even as I was writing version three of find co-infections, so you
  may notice some differences between here and some of my code.

Most files with start out wit ha header that gives an idea of what the
  code will be used for. It may also include a list of C dependencies
  that you will have to use. It does not always include system calls,
  but does include standard c libraries and my own code. Unfortunately
  I have not kept everything as up to date as it should be, so their
  may be some missing dependencies. This block is distinguished by
  /\*#####

```
/*######################################################################
# Name: someCode
# Use:
#   - This code does something
# Includes:
#   - Included in the .h file
#   o Included by an include in the .h file
# Note:
#   - Some of these blocks may have a dependency tree, so the bullet
#     format may not hold
######################################################################*/
```

The next block is the SOF (start of functions), SOP (start of program),
  or in some older code the TOC (table of contents block). This block
  is distinguished by header with /\*>>>> (older) or a /\*~~~~ (newer)
  header. In this block are the function numbers (fun-x), function name,
  and a brief description of what the function does. You can navigate
  to the function by searching for "Fun-x TOC", were x is the function
  number. This block is not in the header files.

```
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF:
'   fun-1 doSomething:
'     o This does something, but I am not sure what
'   fun-2 doNothing:
'     o This is an infinite loop
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
```

Above every function is my output bock, which tells what my functions
  modifies or returns. It will also mention any warnings or notes I have
  put down. In newer code it is distinguished by /\*-----, while in
  older code it looks like the first block (/\*#####).

```
/*---------------------------------------------------------------------\
| Output:
|   - Returns:
|     o 1 if succeeded  
|     o 0 if failed
| Notes:
|   - Something
\---------------------------------------------------------------------*/
```

For my functions I enter every parameter on a separate line. After
  each parameter is a brief description of what it might do (this may
  be a bit cryptic at times). The ) ending the function definition also
  has a brief description of what the function does.

```
void doSomething(
    char inC,   /*Charter to do something to*/
    char *outC, /*Character to modify*/
) /*This does something to a inC and puts change outcome in outC*/
```

The next block is my function TOC block, which lists what each section
  of the function does. For functions that are to small for sections I
  use Fun-x TOC: Sec-1 Sub-1: functionName. In older versions this block
  used ">>>".

```
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' Fun-x TOC: function name
'   o sec-1: variable declerations
'   o sec-2: do something
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

Or if their are no sections

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' Fun-x TOC: Sec-1 Sub-1: function name
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
```

The next block is my section block, which is used to mark sections of
  code and tells what that section of code does. It will also have a
  list of subsections in the section when there are subsections.

```
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Fun-1 Sec-1: This code does something
^   o sec-1 sub-1: the frist part of this section
^   o sec-1 sub-2: the second part of this section
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

Or if no subsections

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Fun-1 Sec-1: This code does something
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

```

The final block is my subsection block. It just list what the goal of
  the subsection is.

```
/**********************************************************************\
* Fun-1 Sec-1 Sub-1: This subsection does something
\**********************************************************************/
```

You will notice that I am using lower case when I am listing things in a
  table of contents, but am using upper case for the sections. This 
  allows a case sensitive search to quickly navigate to a particular
  function, section of a function, or subsection of a function.

Finally all statements that use braces will always have a duplicate
  description or name (if a function) at the start and end of the braced
  statement. This is here so I can easily tell were a brace goes.

```
if(x == y)
{ /*If x is equal to y then subtract and return z*/
    z = x - y; /*Not a good example, I would normally just return this*/
    return z;
} /*If x is equal to y then subtract and return z*/
```
