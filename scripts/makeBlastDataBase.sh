#!/bin/bash

# Name: makeBlastDataBase.sh
# Use: makes a blast data base from an input fasta file
# Input:
#     $1: fasta file with references (nuclotide only)
#     $2: database name
#     $3: If 0 do not create a circulaized database
# Output:
#     Directory with the blast data base
# Requires:
#    makeblastdb from ncbi


if [[ ! -f "$1" ]];
then # if not a valid file
    printf "%s does not exist\n" "$1";
    printf "usage: bash makeBlastDataBase.sh reference.fasta data-base-name\n"
    exit;
fi # if not a valid file

mkdir "$2";

if [[ "$3" -gt 0 ]]; 
then # if allowing data base to handle circular sequences (duplicate sequence)
    printf "Cirucularing %s (doubling each genome)\n" "$1";
    awk '{if($1 ~ /^>/)
          { # if found a header line

              {if(seqStr == ""){print $1;} # just the header
               else
               { # else have a sequence to print out
                  gsub("\r", "",  seqStr); # had issue with ^M in file
                  printf "\n%s%s\n%s", seqStr, seqStr, $1; # print doubled sequence line
               } # else have a sequence to print out
              } # check if printing just header or previous sequence

              seqStr = ""; # clean for the next sequence
           } # if found a header line

           else{seqStr = seqStr $1}
          } # check if sequence line or header

          END{
            gsub("\r", "",  seqStr); # had issue with ^M in file
            printf "\n%s%s", seqStr, seqStr; # print the last sequence out
          } # END block
    ' < "$1" |
      sed '/^$/d' |
	sed -n '/>/{N; # Move past the header to the sequence
        s/[ \t]*//g; # reomve all white space
        p;}' \
		> tmp.fasta
        # sed removes any white spaces in non header lines
    makeblastdb -in tmp.fasta -out "$2/$2--database" -dbtype nucl -parse_seqids;
    rm tmp.fasta
    exit; 
fi # if wanted to allow data base to handle circular genomes (double length)

makeblastdb -in "$1" -out "$2/$2--database" -dbtype nucl -parse_seqids;
exit;
