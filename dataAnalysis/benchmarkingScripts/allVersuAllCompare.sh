intReads=0; # read on
seqStr="$1"; # references to use with minimap2
minSimDbl="$2"; # min similarty to keep alignment
maxSimDbl="$3"; # min similarty to keep alignment

if [[ "$maxSimDbl" -lt "$minSimDlb" ]]; then
    maxSimDbl=100;  # 100%
fi # if user did not specify a max simulatity range

numReadsInt="$( \
 sed \
    -n \
    '/^>/p;' \
    "$seqStr" |
  wc -l |
  awk '{print $1}' \
)"; # Get the number of reads in the query file

while [[ "$intReads" -lt "$numReadsInt" ]]; do
# While their are reads to compare
    sed \
        -n \
	"$((intReads * 2 + 1)), $((intReads * 2 + 2)) p;" \
        "$seqStr" \
      > "tmp-ref.fasta";

    head \
        -n $((intReads * 2)) \
        < "$seqStr" \
      > "tmp-other.fasta";

    tail \
        -n+$((intReads * 2 + 3)) \
        < "$seqStr" \
      >> "tmp-other.fasta";

    minimap2 \
        --secondary=yes \
        -a \
        "tmp-ref.fasta" \
        "tmp-other.fasta" |
      awk \
        -v minSimDbl="$minSimDbl" \
        -v maxSimDbl="$maxSimDbl" \
	'
        { # MAIN block 
            if($1 ~ /^@/)
                next;
            sub(/NM:i:/, "", $12);

            percDifDbl = 100 * $12 / 1260;

            if(percDifDbl >= minSimDbl && percDifDbl <= maxSimDbl)
                print $1, $3, 100*$12/1260; # if in read keeping range
        } # MAIN block
      '

    intReads=$((intReads + 1));
done # while their are reads to compare

rm \
    "tmp-ref.fasta" \
    "tmp-ref.fasta";
