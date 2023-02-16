#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: Columns in expected dataset
#    sec-2: Begin block
#    sec-3: Check if test is duplicate
#    sec-4: END block, print out file test
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

################################################################################
# Name: detectDupRef.awk
# Use:
#    - Detects duplicate references in a series of my co-infection pipeline
#      tests and prints TRUE for duplicate with the most errors. NA for no
#      reference
# Input:
#    file.tsv: 
#        - Test results with duplicates, should be organized by pipeline & test
# Output:
#    Prints: file.tsv to screen with; FALSE, TRUE, or NA at end
################################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Sec-1 Sub-1: Columns in expected dataset
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# $1 = pipeline used
# $2 = seed used with badread
# $3 = percecent id
# $4 = reference for major variant
# $5 = reference for minor variant
# $6 = number of reads made for major variant
# $7 = number of reads made for minor variant
# $8 = number of major variant reads binned
# $9 = number of minor variant reads binned
# $10 = number of junk reads binned
# $11 = number of random reads binned
# $12 = time pipeline took to run in seconds
# $13 = max resident memory pipeline used
# $14 = mapped reference from database or the variant number
# $15 = minor or major variant reference (consensus mapped closest to)
# $16 = mismatches
# $17 = indels
# $18 = will be duplicated (not in file)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Sec-2 Sub-1: BEGIN block
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
    OFS = "\t";
    numConInt = 0;     # keeps track of number of consensus output by test
} # number of matches for each test

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Check if test is duplicate
#    sec-1 sub-1: Check if on the header
#    Sec-3 sub-2: Check if moving onto another set of tests / pipeline
#    sec-3 sub-3: Check if have duplicate references
#    sec-3 sub-4: Record values for the new consensus
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


{ # main block

    #***************************************************************************
    # Sec-3 Sub-1: Check if on the header
    #***************************************************************************

    if(NR == 1)
    { # if was the header row
        print $0, "Duplicate";  # print the header
        next;      # move onto the data
    } # if was the header row

    #***************************************************************************
    # Sec-3 Sub-2: Check if moving onto another set of tests / pipeline
    #***************************************************************************

    if(pipeStr!=$1 || majorRefStr!=$4 || minorRefStr!=$5 || totalReads!=$6 + $7)
    { # if is a new test / pair or references
        for(intCon = 1; intCon <= numConInt; intCon++)
            print testResultAry[intCon], dupAryStr[intCon];

        pipeStr = $1;           # pipeline used
        majorRefStr = $4;       # reference for major variant
        minorRefStr = $5;       # reference for minor varaint
        totalReads = $6 + $7;   # total major and minor variant reads made
        testResultAry[1] = $0;  # record the line for latter printing out
        misAryInt[1] = $16;
        indelAryInt[1] = $17;

        if($15 != "")
            dupAryStr[1] = "FALSE"; # assume is not a duplicate, until proven
        else
            dupAryStr[1] = "NA";    # Their is no reference

        mappedRefAryStr[1] = $15
        numConInt = 1;         # reset for new set of test.
        next;                 # move onto the next consensus
    } # if is a new test / pair or references

    #***************************************************************************
    # Sec-3 Sub-3: Check if have duplicate references
    #***************************************************************************

    dupAryStr[numConInt + 1]= "FALSE"; # Assume this consensus is not a duplicate

    if($15 != "")
    { # if the consensus mapped to the major or minor variant reference
        for(intCon = 1; intCon <= numConInt; intCon++)
        { # for all mapped references
            if(mappedRefAryStr[intCon] == $15)
            { # if there is a duplicate reference
                # check which consensus to flag as duplicate
                if(dupAryStr[intCon] == "TRUE")
                    continue; # wait to flag duplicate till compare mismatches

                if(misAryInt[intCon] > $16)
                    dupAryStr[intCon] = "TRUE"; # other cosensus more mismatches

                else if(misAryInt[intCon] < $16)
                { # else if this consenuss has more mismatches
                    dupAryStr[numConInt+1]= "TRUE";
                    break; # no need to continue
                } # else if this consenuss has more mismatches

                else if(indelAryInt[intCon] > $17)
                    dupAryStr[intCon] = "TRUE"; # other cosensus > indels

                else if(indelAryInt[intCon] < $17)
                { # else if, this consensus has more indels
                    dupAryStr[numConInt + 1]= "TRUE";
                    break;  # no need to continue
                } # else if, this consensus has more indels
 
                else
                { # else both cosensus are equal, flag this one as duplicate
                    dupAryStr[numConInt + 1]= "TRUE"; # both equal, flag this
                    break;  # no need to continue
                } # else both cosensus are equal, flag this one as duplicate
            } # if there is a duplicate reference
        } # for all mapped references
    } # if the consensus mapped to the major or minor variant reference

    #***************************************************************************
    # Sec-3 Sub-4: Record values for the new consensus
    #***************************************************************************

    numConInt++;
    testResultAry[numConInt] = $0;  # record the line for latter printing out
    mappedRefAryStr[numConInt] = $15;
    misAryInt[numConInt] = $16;
    indelAryInt[numConInt] = $17;

} # main block

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Sec-4 Sub-1: END block, print out file test
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

END{
    for(intCon = 1; intCon <= numConInt; intCon++)
        print testResultAry[intCon], dupAryStr[intCon];
} # END block

