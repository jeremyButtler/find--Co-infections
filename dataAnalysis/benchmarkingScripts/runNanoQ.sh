#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: Variable declerations
#    sec-2: Get and check user input
#    sec-3: Run Nano-Q
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: runNanoQ.sh
# Use:
#    Is a wrapper to run Nano-Q (https://github.com/PresetonLeung/Nano-Q)
# Input:
#    -fastq:
#        - fastq file with reads to detect varaints in
#    -ref:
#         - fasta file with references to align reads to
#    -p:
#         - prefix to name everything
#    -racon:
#        - Use racon to build a consensus instead of Nano-Q
#    -h:
#        - Print help message (has all possible commands)
# Output:
#    file: prefix.bam used with Nano-Q
#    Diretory: prefix, holding output from Nano-Q
# Run:
#    bash runNanoQ.sh -fastq reads.fastq -ref refs.fasta
# Requires
#    Nano-Q, see help message for possible install locations
#    racon, Only if using -racon option
#    samtools
#    minimap2
#    awk
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#    sec-1 sub-1: User input variables
#    sec-1 sub-2: Variables specific to script
#    sec-1 sub-3: Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-1 Sub-1: User input variables
#*******************************************************************************

readsFastqStr=""; # Fastq file for NanoQ to detect variants in
refFastaStr="";   # fasta file with references
minMapqDbl=20;    # minimum mapping quality to keep an aligned read
minLenInt=600;    # minimum aligned length to keep a read
prefixStr="Nano-q-out"; # prefix for file names
numThreadsInt=3;     # number threads to use (and subprocesses Nano-Q luanches)
cutOffOffsetInt=50;  # how much to offset cut off value for Nano-Q
pathToNanoQStr="";   # path to nanoq.py
useRaconBool=0;     # use racon to build consensuses
jumpIntervalInt=0;          # number reads to process per cluster

# Variables for Nano-Q

codonStartInt=1;       # frist codon in coding region
readLenCutOffInt=0;    # max read legnth??
nanoQBaseQInt=5;       # base Q-score to keep base
maxHamDistInt=234;     # maximum haming distance to group reads in cluster
minReadsPerClustInt=30; # minimum number of reads to keep cluster
singleArgsStr="";      # holds the single argument options for Nano-Q

#*******************************************************************************
# Sec-1 Sub-2: Variables specific to script
#*******************************************************************************

numRefsInt=0; # number of references in fasta file
jumpIntervalIntervalInt=0;     # jump interval are the number of reads Nano-Q caclulates
                       # the hamming distance at a time. Each jump interval
                       # launches a subrocess. If poorly choosen this will
                       # launch so many process that it will lock up the syste
                       # For example a jump distance of 10 for 20000 reads will
                       #   launch 2000 subprocess. However, you system will
                       #   kill Nano-Q before this, due to too much RAM usage
scriptDirStr="$(dirname "$0")"; # location of this script

#*******************************************************************************
# Sec-1 Sub-3: Help message
#*******************************************************************************

helpStr="$(basename "$0") -fastq reads.fastq -ref refs.fasta
  Input:
    -fastq: Fastq file of reads to find variants in         [Required]
    -ref: Fasta file with references to align reads against [Required]
    -p: Prefix for file names                               [Nano-q-out]
    -racon: Use racon to build a cosensus for each cluster  [0 = no]
            of reads.
        - Sets -kc (keep clusters) to ON.
    -min-q: Minimum mapping quality to keep read            [20]
    -min-len: Minimum aligned length to keep read           [600]
    -path-to-nanoq: Path to nanoq.py. If not provided checks
        - $(pwd)/Nano-Q/nano-q.py
        - $(dirname "$0")/nano-q.py
        - $(dirname "$0")/Nano-Q/nano-q.py
        - $HOME/Nano-Q/nano-q.py
        - /usr/local/bin/Nano-Q/nano-q.py
    -offest-cut-off: Number to offset Nano-Q read cut off by [50]
        - Nano-Q read cut off (-l) is found by
          length of shortest reference - offset-cut-off
        - If user provides -l, then this paramter is ignored
    -t: Number of threads to use                            [3]
        - Sets Nano-Q's jumpt distance to
           total mapped reads / number threads.
          This prevents Nano-Q from launching over a
          1000 subprocess for PCR data. However, this will
          not always launch three subprocess.
    -h, --h, -help, or --help: print this help message

    Variables for Nano-Q that take input
    -c: Starting codon position in the reference genome     [1]
    -l: Length cut off for read size                   [shortest reference - 50]
        - Nano-Q will discard sequences that are under
           this length, due to not trimming these sequences.
           My best guess is that Nano-Q assumes there
           will be primers at the end.
    -j: Number of reads to process per thread               [1000]
        - If not supplied calculated based on number mapped
	  reads.
	- Warning: for old Nano-Q this can launch more 
	  process then number threads specified.
	- For my scribbles.pxy file, this luanches -t
          (number threads) subprocess per cluster + main
          python process
    -q: Quality threshold cut off for each base             [5]
    -ht: Maximum hamming distant to group reads in cluster  [234]
    -mc: Minimum number of reads to keep a cluster          [30]

    Variables telling Nano-Q to turn on a feature (Do not provide arguments) 
    -d: Draw dendrogram to keep help determine HD-Threshod  [OFF]
        (-ht) cutoff
    -hd: Keep the hamming distance files (can be large)     [OFF]
    -kc: Keep the clustered reads (can be large)            [OFF]
  Output:
    Nano-Q makes a file called Results with output
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get and check user input
#    sec-2 sub-1: Get user input
#    sec-2 sub-2: Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: Get user input
#*******************************************************************************

while [ $# -gt 0 ]; do
# while their is user input to check
    # Check if their is an agrument (parameters start with -)
    if [[ "$(printf "%s" "$2" | sed 's/^-.*/-/')" == "-"  || "$2" == "" ]]; then
    # if parameter has no argument
        if [[ "$singleArgsStr" != "" ]]; then
            singleArgsStr="$singleArgsStr ";
        fi # if need to add space between single arguments

        case $1 in             # checking single argument input for Nano-Q
            -racon) useRaconBool=1;  # Use racon to build a consensuses (con)
                    singleArgsStr="$singleArgsStr-kc";; # need reads for con
            -d) singleArgsStr="$singleArgsStr-d";;
            -hd) singleArgsStr="$singleArgsStr-hd";;
            -kc) singleArgsStr="$singleArgsStr-kc";;
            ?) printf \
                   "%s\n%s has no arguments\n" \
                   "$helpStr" \
                   "$1"; exit;;
       esac
       shift;    # get of parameter to next set of input
       continue; # move on in loop (avoid double shift later)
   fi # if parameter has no argument

   case $1 in
       -fastq) readsFastqStr="$2";;                 # fastq file with reads
       -ref) refFastaStr="$2";;                     # fasta file with references
       -min-q) minMapqDbl="$2";;                    # min mapq to keep read
       -min-len) minLenInt="$2";;                   # min length to keep read
       -offest-cut-off) cutOffOffsetInt="$2";;      # offset for Nano-Q cut off
       -t) numThreadsInt="$2";;                        # number threads to use
       -p) prefixStr="$2";;                         # prefix for file names
       -path-to-nanoq) pathToNanoQStr="$2";;        # path to nanoq program

       # Nano-Q variables
       -c) codonStartInt="$2";;                    # frist codon position
       -j) jumpIntervalInt="$2";;                  # num reads to process/thread
       -l) readLenCutOffInt="$2";;                 # max read legnth??
       -q) nanoQBaseQInt="$2";;                    # base Q-score to keep base
       -ht) maxHamDistInt="$2";;                  # max ham dist to keep cluster
       -mc) minReadsPerClustInt="$2";;             # min reads to keep cluster

       # various help commands to print the help message
       -h) printf "%s\n" "$helpStr"; exit;;
       --h) printf "%s\n" "$helpStr"; exit;;
       -help) printf "%s\n" "$helpStr"; exit;;
       --help) printf "%s\n" "$helpStr"; exit;;
       ?) printf \                                  # invalid input
              "%s\n-%s is not a vaid paremeter\n" \
              "$helpStr" \
              "$1"; exit;;
   esac

   shift;                                           # move to argument
   shift;                                           # move to next parameter
done # whild there is user input to check

#*******************************************************************************
# Sec-2 Sub-2: Check user input
#*******************************************************************************

if [[ ! -f "$readsFastqStr" ]]; then
# if user did not input a valid fastq file
    printf \
        "%s is not a valid fastq file. Input fastq file with reads (-fastq)\n" \
        "$readsFastqStr";
    exit;
fi # if user did not input a valid fastq file

if [[ ! -f "$refFastaStr" ]]; then
# if user did not input a valid fastq file
    printf \
        "%s is not a valid file. Input fasta file with reference(s) (-ref)\n" \
        "$refFastaStr";
    exit;
fi # if user did not input a valid fastq file

if [[ ! -f "$pathToNanoQStr" ]]; then
# If the user did not provide a path to nano-q.py
    if [[ -f "Nano-Q/nano-q.py" ]]; then
    # if no nano-q.py path provided check some default locatoins
        pathToNanoQStr="Nano-Q";

    elif [[ -f "$scriptDirStr/nano-q.py" ]]; then
        pathToNanoQStr="$scriptDirStr";

    elif [[ -f "$scriptDirStr/Nano-Q/nano-q.py" ]]; then
        pathToNanoQStr="$scriptDirStr/Nano-Q";

    elif [[ -f "$HOME/Nano-Q/nano-q.py" ]]; then
        pathToNanoQStr="$HOME/Nano-Q";

    elif [[ -f "/usr/local/bin/Nano-Q/nano-q.py" ]]; then
        pathToNanoQStr="/usr/local/bin/Nano-Q";

    else
        printf \
            "Can not find nano-q.py. Please provide path with -path-to-nanoq\n";
    fi # if no nano-q.py path provided check some default locatoins
fi # check if can find path to Nano-Q

# Remove duplicate -kc commands. Otherwise Nano-Q sets off
singleArgsStr="$( \
    printf \
      "%s" \
      "$singleArgsStr" | 
    awk '
      { # MAIN
         printf "%s", $1; # print out the first argument (want to keep)

         if($1 == "-kc")
             intKc = 1;  # is a -kc, make so can remove all other -kc
         else
             intKc = 0; # Mark that I still need to find the first -kc

          for(intCol = 2; intCol <= NF; intCol++)
          { # loop through all feilds
              if($intCol != "-kc")
              { # If was a parameter I want to keep
                  printf " %s", $intCol; # not -kc, so ok to print out
                  continue;
              } # If was a paramter I want to  keep

              if(intKc < 1)
              { # if the first -kc
                      printf " %s", $intCol; # print the -kc, since first time
                      intKc = 1;            # mark so no other -kc are printed
              } # if the first -kc
          } # loop through all feilds
      }; # MAIN block
    ' \
)"; # Get the number of times -kc was supplied

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Run Nano-Q
#    sec-3 sub-1: Make indexed bamfile for Nano-Q
#    sec-3 sub-2: Find the cut off length of Nano-Q if user did not supply one
#    sec-3 sub-3: Find jump interval (number threads) & number of refs supplied
#    sec-3 sub-4: Run Nano-Q
#    sec-3 sub-5: Run Racon (if user wanted)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: Make indexed bamfile for Nano-Q
#*******************************************************************************

if [[ ! -f "$prefixStr.bam" ]]; then
# If need to make a bam file
    minimap2 \
        --MD \
        --secondary=no \
        -ax map-ont \
        -t "$numThreadsInt" \
        "$refFastaStr" \
        "$readsFastqStr" |
      awk \
        -v numThreadsInt="$numThreadsInt" \
        -v prefStr="$prefixStr" \
        ' # Get the maximum number of reads mapped to a single reference
          # Uses the Mapq value to figure this out.
        BEGIN{
           intRef = 0;
           if(prefStr == "")
               prefStr = "out";
           if(numThreadsInt == "")
               numThreadsInt = 1;
           outFileStr = prefStr "--count.txt";
        }; # BEGIN block

        { # MAIN block
            print $0; # print out samfile for samtools

            if($1 ~ /^@/)
                next;                           # Is a header, move on
        
            if($5 != "*" && $5 > 0)
            { # if this is the best matching reference for the read
                if(intRef == 0)
                { # If this is the frist kept alignment
                    intRef = 1;                 # Number of references with mapped reads
                    refAryStr[intRef] = $3;     # Name of reference one
                    refCntAryInt[intRef] = 1;   # Number of reads mapped to reference 1
                    next;                       # move onto the next alignment
                }; # If this is the first kept alignment 
        
                for(intCnt = 1; intCnt <= intRef; intCnt++)
                { # for all mapped references, check if matches current reference
                    if(refAryStr[intCnt] == $3)
                    { # If was a match, incurment count
                        refCntAryInt[intCnt]++;     # Incurment count for mathching ref
                        break;                      # found match, no need to look further
                    }; # If was a match, incument count
        
                    if(intCnt == intRef)
                    { # If there was not matching reference (is new)
                        intRef++;                   # adding in new reference
                        refAryStr[intRef] = $3;     # get new reference name
                        refCntAryInt[intRef] = 1;   # set counter for new ref to 1
                        break;
                    }; # If there was no matching reference (is new)
                }; # for all mapped references, check if matches current reference
            } # if this is the best matching reference for the read
        }; # MAIN block

        END{ # END block
            for(intCnt = 1; intCnt <= intRef; intCnt++)
            { # for all references with reads, check which has most reads
                if(refCntAryInt[intCnt] > mostMappedReadsInt)
                    mostMappedReadsInt = refCntAryInt[intCnt];
            }; # for all references with reads, check which has most reads
        
            printf "%i", mostMappedReadsInt / numThreadsInt > outFileStr; # print out the max read count/ref
        }; # END block
        ' |
      samtools sort \
        -@ "$numThreadsInt" \
        -T tmp.txt |
      samtools view \
        -@ "$numThreadsInt" \
        -F 0x04 \
        -b \
        -q "$minMapqDbl" \
        --min-qlen "$minLenInt" \
      > "$prefixStr.bam";

      # the awk script is massive and adds in a bit more processing time,
      # but gets the maximum number of reads mapped to a single reference.
      # This will allow setting of a better jump setting when multiple
      # references are used.
      # I put hte awk script here instead of in a file to avoid having to
      # have to move two files around.
fi # If need to make a bam file

if [[ ! -f "$prefixStr.bam.bai" ]]; then
# if need to index the bam file
    samtools index \
        -@ "$numThreadsInt" \
        "$prefixStr.bam";   # Index bamfile for Nano-Q
fi # if need to index the bam file

#*******************************************************************************
# Sec-3 Sub-2: Find the cut off length of Nano-Q if user did not supply one
#*******************************************************************************

if [[ "$readLenCutOffInt" -lt 1 ]]; then
    readLenCutOffInt="$(
      awk \
          -v offSetInt="$cutOffOffsetInt" \
          '
          BEGIN{shortestLenInt = 1000000000;}

          { # MAIN BLOCK
              if($1 ~ /^>/)
              { # if on a header line
                  if(seqStr == "")
                      next;
                  tmpLenInt = length(seqStr);  # find previous sequences length

                  if(tmpLenInt < shortestLenInt)
                      shortestLenInt = tmpLenInt;   # find shortest length

                  seqStr = "";                      # reset for next sequence
                  next;                             # move to the sequence
              } # if on a header line

              seqStr = seqStr $0;
          } # MAIN BLOCK

          END{
              tmpLenInt = length(seqStr);  # find previous sequences length

              if(tmpLenInt < shortestLenInt)
                  shortestLenInt = tmpLenInt;   # find shortest length

              print shortestLenInt - offSetInt;  # print out cut off value
              # Subtracting 50 so Nano-Q will be force to trim every sequence
              # Otherwise does not print out sequences
          } # END BLOCK
          ' \
        < "$refFastaStr" \
    )";
fi # if user did not provide a read cut off value

#*******************************************************************************
# Sec-3 Sub-3: Find jump interval (number threads) and number of refs supplied
#*******************************************************************************

# Get the size of jump interval for a given thread
#  (otherwise Nano-Q launches to many subrocess)

if [[ "$jumpIntervalInt" -lt 1 ]]; then
    jumpIntervalIntervalInt="$(cat "$prefixStr--count.txt")";
fi # if have to find the number of jumps to do

if [[ "$jumpIntervalIntervalInt" == "" ]]; then
    printf \
        "No reads in %s were mapped, please provide a new set of references\n" \
        "$readsFastqStr";
    exit;
fi # if no reads were mapped

# -@: is number of threads for samtools view
# -c: count number of reads
# -F 4: Only count mapped reads (Really is ignore unmapped reads)

numRefsInt="$( \
  sed \
    -n \
    '/^>/p;' \
    < "$refFastaStr" |
  wc -l |
  awk '{print $1}'
)"; # finder number of references the user input

#*******************************************************************************
# Sec-3 Sub-4: Run Nano-Q
#*******************************************************************************

if [[ "$singleArgsStr" == "" ]]; then
# if runing Nano-Q without any of the single argument parameters
    python3 \
        "$pathToNanoQStr/nano-q.py" \
        -b "$prefixStr.bam" \
        -c "$codonStartInt" \
        -l "$readLenCutOffInt" \
        -nr "$numRefsInt" \
        -q "$nanoQBaseQInt" \
        -j "$jumpIntervalIntervalInt" \
        -ht "$maxHamDistInt" \
        -ct "$minReadsPerClustInt";
# if runing Nano-Q without any of the single argument parameters
else
# else running Nano-Q with at least one single argument command
    python3 \
        "$pathToNanoQStr/nano-q.py" \
        -b "$prefixStr.bam" \
        -c "$codonStartInt" \
        -l "$readLenCutOffInt" \
        -nr "$numRefsInt" \
        -q "$nanoQBaseQInt" \
        -j "$jumpIntervalIntervalInt" \
        -ht "$maxHamDistInt" \
        -ct "$minReadsPerClustInt" \
        "$singleArgsStr";
# else running Nano-Q with at least one single argument command
fi # check if runing Nano-Q without any of the single argument parameters

# Change the results file name to the prefix name
mv \
    "Results" \
    "$prefixStr";

#*******************************************************************************
# Sec-3 Sub-5: Run Racon (if user wanted)
#   - In my testing I found that Nano-Q seems to like to put out consensuses
#     with only N's. My best guess is that their is enough varaition in my
#     bases to not meet the minimum criteria for keeping a base. This may be
#     do to read quality of badread simulated reads.
#   - Nano-Q only builds a dumb consensus. So, I should get similar or better 
#     results using racon. This will create a consensuses without solid N's
#*******************************************************************************

# Remove unneeded files
rm \
    "$prefixStr.bam" \
    "$prefixStr.bam.bai" \
    "$prefixStr--count.txt";

if [[ "$useRaconBool" -lt 1 ]]; then
    exit;
fi # if not running racon, then finshed

cd "$prefixStr"; # move into directory with clusters

# move the consensus Nano-Q built to a separate directory
mkdir "$prefixStr--Nano-Q-build-consensuses";
mv \
    *ConsensusFinal.fa \
    "$prefixStr--Nano-Q-build-consensuses";

# Run racon
for strClust in ./Clusters/*.fa; do
# For all clusters Nano-Q found, make a consensus

    if [[ ! -f "$strClust" ]]; then
        continue; # is null case, just let loop finish
    fi # if is the null case

    # Grab the top 300 reads to build a consensues
    head \
        -n 2 \
        < "$strClust" \
      > "$prefixStr--first-read.fasta"; # the first read is polished

    sed \
        -n \
        '3,600p; # print out the 299 reads after the first read' \
        < "$strClust" \
      > "$prefixStr--other-reads.fasta"; # Used to polish first read

    minimap2 \
        -ax map-ont \
        -t "$numThreadsInt" \
        "$prefixStr--first-read.fasta" \
        "$prefixStr--other-reads.fasta" \
      > "$prefixStr--tmp.sam";

    racon \
        -m 8 \
        -x -6 \
        -g -8 \
        -w 500 \
        -t "$numThreadsInt" \
        "$prefixStr--other-reads.fasta" \
        "$prefixStr--tmp.sam" \
        "$prefixStr--first-read.fasta" \
      > "$( \
            printf \
                "%s" \
                "$strClust" |
              sed '
                s/^\.\/Clusters\///; # remove the ./Clusters/ from the file name
                s/\.fa/_Consensus.fa/; # add _Consensus to the file name
              ' \
      )"; # build consensus & output to file (named ref_cluster#_Consensus.fa)

    # Remove the unneeded files
    rm \
        "$prefixStr--tmp.sam" \
        "$prefixStr--first-read.fasta" \
        "$prefixStr--other-reads.fasta";
done # For all clusters Nano-Q found, make a consensus

exit;
