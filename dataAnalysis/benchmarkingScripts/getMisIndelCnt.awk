# Name: getMisIndelCnt.awk
# Use:
#   Gets the number of mismatches & indels for sequence
# Input:
#   < file.sam: Sam file to get scores from
# Output:
#   reference\tmismatches\tindels

{ # MAIN
  if($1 ~ /^@/){next;};                # if header of sam file
  if($5 == 0){next;};             # if mapping quality under 0 (not best)

  misCigStr = $6;                      # cigar line
  sub(/[0-9]*$/, "", misCigStr);       # remove matches at end
  gsub(/[0-9]*[=DISH]/, "", misCigStr); # rem indels, matches, softmasks

  # Find number of mismatches
  lenMisInt = split(misCigStr, misAryInt, "X");
  totalMisInt = 0;

  for(intMis = 1; intMis <= lenMisInt; intMis++)
    totalMisInt = totalMisInt + misAryInt[intMis];

  sub(/NM:i:/, "", $12); # remove tag from total errors
  printf "%s\t%i\t%d\n", $3, totalMisInt, $12 - totalMisInt;
    # $12 is indels + mismatches. $12 - totalMisInt = total indels
} # MAIN
