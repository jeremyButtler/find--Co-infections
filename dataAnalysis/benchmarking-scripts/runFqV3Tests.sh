#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#   sec-1: variable declerations
#   sec-2: Get and (check user input; Not done)
#   sec-3: Fun findCoInft & get stats for consensuses
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

########################################################################
# Name:
# Use:
#   simlulates read for input pair of references with badread & runs
#   findCoInfections scripts on simulated data
# Input:
#   -fastq-dir:
#     o directory of fastq files simulated with               [Required]
#       twoReadRunBadread.sh
#   -ref-dir:                                                 [Required]
#     o directory of references used to make the
#       fastq files (needs ref names in file name)
#   -map-ref:                                                 [Required]
#     o Fasta file to use in binning reads
#   -p:                                                       [out]
#     o prefix to name output files
# Output:
#   File: with stats for each pipeline         ["$prefixStr--stats.tsv"]
########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#    sec-1 sub-1: varaibles holding user input
#    sec-1 sub-2: script variables
#    sec-1 sub-3: Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#***********************************************************************
# Sec-1 Sub-1: Variables holding user input
#***********************************************************************

fastqDirStr="";              # fasta file with references
refDbFaStr="";               # references to map reads to
refDirStr="";
prefixStr="out";

#***********************************************************************
# Sec-1 Sub-2: script varaibles
#***********************************************************************

dashLnStr=""; # dash line for separating tests
majVarStr=""; # major variant
minVarStr=""; # minor variant
numMajVarReadsInt=0; # number of major variant reads
numMinVarReadsInt=0; # number of reads from minor variant
numJunkInt=0;        # number of junk sequences
numRandInt=0;        # number of random sequences
pLineStr="";         # Line to print out
refStatStr="";       # holds number of errors in reference
mapRefStr="";        # reference consensus mapped to in reference library
refFaStr="";         # reference file used to simulate reads
scriptDirStr="$(dirname "$0")"; # this script's directory

#***********************************************************************
# Sec-1 Sub-3: Help message
#***********************************************************************

helpStr="$(basename "$0") -fasta file.fasta -map-ref references.fasta -ref-dir 
  Use:
    simlulates read for input pair of references with badread & runs
    findCoInfections scripts on simulated data
  Input:
    -fastq-dir:
      o directory of fastq files simulated with               [Required]
        twoReadRunBadread.sh
    -ref-dir:                                                 [Required]
      o directory of references used to make the
        fastq files (needs ref names in file name)
    -map-ref:                                                 [Required]
      o Fasta file to use in binning reads
    -p:                                                       [out]
      o prefix to name output files
  Output:
    File: with stats for each pipeline               [prefix--stats.tsv]
"; # help message

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get user input
#   sec-2 sub-1: get user input
#   sec-2 sub-2: check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#***********************************************************************
# Sec-2 Sub-1: get user input
#***********************************************************************

while [ $# -gt 0 ]; do
# While there is input to read

  if [[ "$2" == "" ]]; then
    printf "%s\n%s has no arguments\n" \
        "$helpStr" \
        "$1";
    exit;
  fi # if argument is blank

  case $1 in
    -fastq-dir) fastqDirStr="$2";;
    -ref-dir) refDirStr="$2";;
    -map-ref) refDbFaStr="$2";;
    -p) prefixStr="$2";;
    -h) printf "%s\n" \
            "$helpStr" &&
          exit;;
   ?) printf "%s\n%s is not valid\n" \
          "$helpStr" \
          "$1" &&
        exit;;
  esac

  shift;  # move to parameter
  shift;  # move to next argument
done # while their are user arguemnts to check

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Run badread & find co-infections
#   sec-3 sub-1: Print out header & loop through fastq files
#   sec-3 sub-2: get metadata & referenence names from the fastq file
#   sec-3 sub-3: Build consensuses using v3 pipeline
#   sec-3 sub-4: Extract full length reads, with ref meta data for bin
#   sec-3 sub-5: Get number of major & minor strain reads in bin
#   sec-3 sub-6: Get the stats & mappped reference for the bin
#   sec-3 sub-7: Print out all extracted data & clean up
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#***********************************************************************
# Sec-3 Sub-1: Print out header & loop through fastq files in directory
#***********************************************************************

{ # stats header print block
  printf "pipeline\tpercId\tmajorRef\tminorReft\tnoMajorRef\t"
  printf "\tnoMinorRef\tnoBinMajorRef\tnoBinMinorRef\tnoBinJunk";
  printf "\tnoBinRand\ttime\tmemory\tmappedRef\tacutalRef";
  printf "\tmismatches\tindels\tscoreRef\tlength\talnLength\tmatches"
  printf "\tsnps\tkeptIns\tkeptDel\tins\tdel\n"
} > "$prefixStr--stats.tsv" # stats header print block

dashLnStr="~~	~~	~~	~~	~~";
dashLnStr="$dashLnStr	~~	~~	~~	~~";
dashLnStr="$dashLnStr	~~	~~	~~	~~	~~";
dashLnStr="$dashLnStr	~~	~~	~~	~~	~~	~~";
dashLnStr="$dashLnStr	~~	~~	~~	~~	~~";

for strFq in ./"$fastqDirStr"*".fastq"; do
# for all fastq files

    if [[ ! -f "$strFq" ]]; then
        continue;                # null case
    fi

    # so can tell when miss entries
    printf "%s\n" "$dashLnStr" >> "$prefixStr--stats.tsv";

    #*******************************************************************
    # Sec-3 Sub-2: Get meta data & reference names from fastq file
    #*******************************************************************

    pLineStr="$(
      awk \
        'BEGIN{FS="_"; OFS="\t"};
         {print $2, $3, $4, $5, $6}
      ' <(printf "%s" "$strFq") \
    )"; # get meta data (id, major ref, minor ref, read maj, reads min)

    # get the major and minor strain/variants references name
    # have to remove -ORF2, because can mess up next step
    majVarStr="$(\
        printf "%s" "$pLineStr" |
          awk '{sub(/-ORF2/, "", $2); print $2}'\
    )";
    minVarStr="$(\
        printf "%s" "$pLineStr" |
          awk '{sub(/-ORF2/, "", $3); print $3}'\
    )";

    # get the fasta file with the references
    refFaStr="$(\
        find "$refDirStr"*"$majVarStr"*"$minVarStr"*".fasta" |
        sed -n '/.fasta/{p;q;};' \
    )";

    #*****************************************************************
    # Sec-3 Sub-3: Build consensuses using v3 pipeline
    #*****************************************************************

    /usr/bin/time \
        -f "%e\t%M" \
        -o "$prefixStr--time.tsv" \
      "$scriptDirStr/../findCoInft" \
          -fastq "$strFq" \
          -ref "$refDbFaStr" \
          -prefix "$prefixStr--v3";

    #*******************************************************************
    # Sec-3 Sub-4: Extract full length reads, with ref meta data for bin
    #*******************************************************************

    for strCon in ./"$prefixStr--v3--"*"consensus.fasta"; do
    # for all kept consensuses, get stats
    
      if [[ ! -f "$strCon" ]]; then
        continue; # null case
      fi # if the null case
    
      conFastqStr="$( \
        printf "%s" "$strCon" |
        sed 's/\.\///; s/--consensus\.fasta/.fastq/;' \
      )"; # find the fastq file used to build ref

      # Get fastq file with metadata
      sed -n 'p;n;n;n;' < "$conFastqStr" |
        "$scriptDirStr/../fqGetIds" \
          -stdin-filt \
          -fastq "$strFq" \
        > "$prefixStr--v3--bin.fastq";
    
      #*****************************************************************
      # Sec-3 Sub-5: Get number of major & minor strain reads in bin
      #*****************************************************************

      numMajVarReadsInt="$( \
        awk \
           -f "$scriptDirStr/getNumRefReads.awk" \
           -v refStr="$majVarStr" \
           < "$prefixStr--v3--bin.fastq" \
      )"; # find number of major variant reads
    
      numMinVarReadsInt="$( \
        awk \
          -f "$scriptDirStr/getNumRefReads.awk" \
          -v refStr="$minVarStr" \
          < "$prefixStr--v3--bin.fastq" \
      )"; # Find the number of reads in bin from minor variant
    
      numJunkInt="$( \
        awk \
          -f "$scriptDirStr/getNumRefReads.awk" \
          -v refStr="junk_seq" \
          < "$prefixStr--v3--bin.fastq" \
      )"; # Find the number of reads in bin from minor variant
      
      numRandInt="$( \
        awk \
          -f "$scriptDirStr/getNumRefReads.awk" \
          -v refStr="random_seq" \
          < "$prefixStr--v3--bin.fastq" \
      )"; # Find the number of reads in bin from minor variant

      #*****************************************************************
      # Sec-3 Sub-6: Get the stats & mappped reference for the bin
      #*****************************************************************
    
      refStatStr="$( \
        bash "$scriptDirStr/getBestAlignment.sh" \
            "$strCon" \
            "$refFaStr" \
      )"; # get the errors for the consensus
    
      mapRefStr="$( \
        bash "$scriptDirStr/getBestAlignment.sh" \
            "$strCon" \
            "$refDbFaStr" |
        awk '{print $1}' \
      )"; # Get the closest reference in reference database

      # get extra stats
      refScoreStr="$( \
          minimap2 --eqx -a "$refFaStr" "$strCon" |
          awk '{if($10 != "*"){print $0}}' |
          "$scriptDirStr/../scoreReads" -stdin |
          awk 'BEGIN{OFS="\t"};
	       {
                getline; # move past the header

                sub(/.*/, "", $1); # remove query id
                sub(/.*/, "", $3); # remove mapq
                sub(/.*/, "", $7); # remove kept matchs

                # remove the Q-score entries
                sub(/.*/, "", $11);
                sub(/.*/, "", $12);
                sub(/.*/, "", $13);
                sub(/.*/, "", $14);

                # remove totalSNP (no-qscore so all snps kept)
                sub(/.*/, "", $15);

                oldLine = $0;
                numSnpsI = $8;

                getline;
                sub(/.*/, "", $1); # remove query id
                sub(/.*/, "", $3); # remove mapq
                sub(/.*/, "", $7); # remove kept matchs

                # remove the Q-score entries
                sub(/.*/, "", $11);
                sub(/.*/, "", $12);
                sub(/.*/, "", $13);
                sub(/.*/, "", $14);

                # remove totalSNP (no-qscore so all snps kept)
                sub(/.*/, "", $15);

                # figure out wich match was better by number of snps
                if(numSnpsI < $8)
                    print oldLine;
                else
                    print $0;
                exit;
              }' \
      )"; # get the scoreReads output (for dels/ins)
    
      #*****************************************************************
      # Sec-3 Sub-7: Print out all extracted data & clean up
      #*****************************************************************

      # print out data for this consensus
      printf "V3\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$pLineStr" \
          "$numMajVarReadsInt" \
          "$numMinVarReadsInt" \
          "$numJunkInt" \
          "$numRandInt" \
          "$(cat "$prefixStr--time.tsv")" \
          "$mapRefStr" \
          "$refStatStr" \
          "$refScoreStr" \
        >> "$prefixStr--stats.tsv";
    done # for all kept consensuses, get stats

    rm "$prefixStr--v3"*; # remove all files pipeline made
    rm "$prefixStr--time.tsv"; # Just to be safe
done # for all fastq files

exit;
