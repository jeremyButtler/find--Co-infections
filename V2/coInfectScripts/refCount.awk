BEGIN{
  intClust = 0;

  if(fileStr == "")
    fileStr = "out--ref-counts.tsv";
};
{ # MAIN
 if(readStr != $1)
 { # if moving to next read
    if(readStr != "")
      printf "%s\t%s\n", readStr, intClust >> fileStr;

    intClust = 1;
    readStr = $1;
    print $0;
    next;
 }; # if moving to next read

  intClust++
  print $0;
}; # MAIN

END{printf "%s\t%s\n", readStr, intClust >> fileStr;};
