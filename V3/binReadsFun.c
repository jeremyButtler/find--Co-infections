/*######################################################################
# Use:
#   o Holds functions related to read binning.
######################################################################*/

#include "binReadsFun.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' binReadsFun SOF:
'   fun-1 binReads:
'     o Bins read to a set of references.
'   fun-2 binReadToCon:
'     o Bins reads to a consensus to from a cluster.
'     o This differes from binReads in that it is extracting reads from
'       the former bin and it does not produce a stats file.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*---------------------------------------------------------------------\
| Output:
|   - Modifies:
|     o errUC to hold error output
|       - 1 for success
|       - 2 for file error
|       - 4 for unable to create a fastq file for a bin
|       - 8 for unable to create a stats file for a bin
|       - 64 for memory allocation error
|   - Creates:
|     o A fastq file with reads for each bin
|     o A stats file with the stats from scoreReads for each bin
|   - Returns:
|     o Tree of bins that the reads mapped into
\---------------------------------------------------------------------*/
struct readBin * binReads(
    char *fqPathCStr,        /*Fastq file to bin*/
    char *refsPathCStr,      /*References to bin with*/
    char *prefixCStr,
    char *threadsCStr,       /*Numbe of threads to use with minimap2*/
    char rmSupAlnBl,         /*Remove supplementary alignments*/
    char trimBl,             /*1: trim reads, 0: do not*/
    struct samEntry *newSam, /*Holds minimap2 output*/
    struct samEntry *oldSam, /*Holds previous line of minimap2 output*/
    struct minAlnStats *minStats,
    unsigned char *errUC     /*Reports any errors*/
) /*Bin reads with a set of references*/
{ /*binReads*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: binReads
    '   fun-1 sec-1: Variable declerations
    '   fun-1 sec-2: Run minimap2 & read in first line of output
    '   fun-1 sec-3: Check if first line is valid
    '   fun-1 sec-4: Get past the header & read in/score the first entry
    '   fun-1 sec-5: Bin reads
    '   fun-1 sec-6: Print out the last read
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-1: Varaible declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char minimap2CMDCStr[2048];
    char binFileCStr[256];
    char statFileCStr[256];
    char refIdCStr[256];
    char *tmpCStr = 0;
    char dupBL = 0;        /*Was the last read a duplicate*/
    uint8_t zeroUChar = 0;    /*For when I need to pass a 0 as a pointer*/
    char funErrUC = 0;     /*Holding err output from called functions*/
    unsigned char printStatsHeadUC = 0;/*tells to print stat file head*/

    struct samEntry *tmpSam = 0;       /*For swaping newSam and oldSam*/
    struct samEntry *samZeroStruct = 0; /*Points were avoiding ref*/

    struct readBin *tmpBin = 0;
    struct readBin *binTree = 0;

    struct readBinStack binStack[200]; /*Stack for read bin AVL tree*/

    FILE *stdinFILE = 0;
    FILE *fqBinFILE = 0;
    FILE *statFILE = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Run minimap2 & read in first line of output
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Running minimap2 with on thread so that I can detect duplicate
      entries. Otherwise minimap2 will have no order for ouput
      mappings*/
    tmpCStr = cStrCpInvsDelm(minimap2CMDCStr, minimap2CMD);

    /*I can only use one thread when remove all reads with supplemental
      alignments, but can use multiple when ignoring supplemental
      alignments*/
    if(rmSupAlnBl & 1)
        tmpCStr = cpParmAndArg(tmpCStr, "-t", "1");
    else
        tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);

    tmpCStr = cpParmAndArg(tmpCStr, refsPathCStr, fqPathCStr);

    stdinFILE = popen(minimap2CMDCStr, "r"); /*run minimap2*/

    blankSamEntry(oldSam); /*Remove old stats in sam file*/
    funErrUC = readSamLine(oldSam, stdinFILE);
        /*get the first line from minimap2*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Check if first line is valid
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(!(funErrUC & 1))
    { /*If an error occured*/
        pclose(stdinFILE);    /*No longer need open (due to error*/

        if(!(funErrUC & 64))
        { /*If errored out*/
            *errUC = 2;      /*If was not a memory allocation error*/
            return 0;
        } /*If errored out*/

        else
        { /*If had memory allocation error*/
            *errUC = 64;
            return 0;
        } /*If had memory allocation error*/
    } /*If an error occured*/

    if(*(oldSam->samEntryCStr) != '@')
    { /*If their is no header line, minimap2 likely errored out*/
        pclose(stdinFILE);
        *errUC = 2;
        return 0;
    } /*If minimap2 did not produce a header, it likely error out*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-4: Get past the header & read in/score the first entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(funErrUC & 1)
    { /*While not past the first header*/
        if(*oldSam->samEntryCStr == '@')
        { /*If was a header*/
            blankSamEntry(oldSam); /*Remove old stats in sam file*/
            funErrUC = readSamLine(oldSam, stdinFILE); /*read new line*/
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        if(trimBl & 1)
            funErrUC = trimSamEntry(oldSam); /*Trim the read*/
        else if(oldSam->flagUSht & (2048 | 256 | 4))
            funErrUC = 1 << 2; /*Let next step know to ignore*/

        if(funErrUC >> 2)
        { /*If entry did not have a sequence, discard*/
            /*make the new alignment (not dupicate) the old alignment*/
            tmpSam = oldSam;
            oldSam = oldSam;
            oldSam = tmpSam;

            /*Read the next sam entry*/
            blankSamEntry(oldSam); /*Remove old stats in sam file*/
            funErrUC = readSamLine(oldSam, stdinFILE);
            continue;
        } /*If entry did not have a sequence, discard*/

        else
            funErrUC = 1; /*So next loop fires*/

        findQScores(oldSam); /*Find the Q-scores*/

        scoreAln(
            minStats,      /*thesholds for read to reference map*/
            oldSam,
            samZeroStruct, /*Not using reference for scoring*/ 
            &zeroUChar,    /*Not using reference, so no Q-score*/
            &zeroUChar     /*Not using reference, so no deletions*/
        );

        if(rmSupAlnBl & 1)
        { /*If removing all reads with supplemental alignments*/
            blankSamEntry(newSam); /*Remove old stats in sam file*/
            funErrUC = readSamLine(newSam, stdinFILE);
        } /*If removing all reads with supplemental alignments*/

        else 
        { /*Else ignoring supplemental alignments*/
            newSam = oldSam;
            blankReadStats(newSam); /*Remove old stats in sam file*/
        } /*Else ignoring supplemental alignments*/

        break; /*Found the start of the second entry*/
    } /*While not past the first header*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-5: Bin reads
    ^   fun-1 sec-5 sub-1: Detect if supplmental & check if removing
    ^      reads with supplemental alignments
    ^   fun-1 sec-5 sub-2: Trim and score sam file alignment
    ^   fun-1 sec-5 sub-3: Set up bin file names
    ^   fun-1 sec-5 sub-4: Check if can open the bin & stats file
    ^   fun-1 sec-5 sub-5: Append reads to file & add to tree
    ^   fun-1 sec-5 sub-6: Print out the stats and fastq entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-1 Sec-5 Sub-1: Detect if supplmental & if removing all reads
    *    with supplementals, remove
    \******************************************************************/

    while(funErrUC & 1)
    { /*While their is a samfile entry to read in*/
        dupBL = 0; /*So that I know I can print out the last read*/

        if(*newSam->samEntryCStr == '@')
        { /*If was a header*/
            blankSamEntry(newSam); /*Remove old stats in sam file*/
            funErrUC = readSamLine(newSam, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        if(newSam->flagUSht & (2048 | 256 | 4))
        { /*If is supplemental (2048), secondary (256), or no map (4)*/
            blankSamEntry(newSam); /*Remove old stats in sam file*/
            funErrUC = readSamLine(newSam, stdinFILE);

            if((rmSupAlnBl & newSam->flagUSht) & 2048)
            { /*If is a supplemental alignment & remove supplementals*/
                dupBL = 1; /*Mark that a duplicate was detected*/

                while(newSam->flagUSht & (2048 | 256 | 4))
                { /*While the entries are duplicates*/
                    /*Remove previous reads stats & get next alignment*/
                    blankSamEntry(newSam);
                    funErrUC = readSamLine(newSam, stdinFILE);

                    if(!(funErrUC & 1))
                        break;      /*end of file or other error*/
                } /*While the entries are duplicates*/

                blankSamEntry(oldSam);

                if(funErrUC & 1)
                { /*If need to grab the next line still*/
                    funErrUC = readSamLine(oldSam, stdinFILE);
                    dupBL = 0;
                } /*If need to grab the next line still*/

                tmpSam = newSam;
                newSam = oldSam;
                oldSam = tmpSam;
            } /*If is a supplemental alignment & remove supplementals*/

            continue; /*Move check the next alignment*/
        } /*If is supplemental (2048), secondary (256), or no map (4)*/
    
        /**************************************************************\
        * Fun-1 Sec-5 Sub-2: Trim and score sam file alignment
        \**************************************************************/

        /*Convert & print out sam file entry*/
        if(trimBl & 1)
            funErrUC = trimSamEntry(newSam); /*If trimming reads*/

        findQScores(newSam); /*Find the Q-scores*/

        scoreAln(
            minStats, /*thesholds for read to reference map*/
            newSam,
            samZeroStruct, /*Not using reference for scoring*/ 
            &zeroUChar,    /*Not using reference, so no Q-score*/
            &zeroUChar     /*Not using reference, so no deletions*/
        );

        if(rmSupAlnBl & 1)
            tmpSam = oldSam;
        else
            tmpSam = newSam;

        if(
          checkRead(minStats, tmpSam) == 0 ||
          !(checkIfKeepRead(minStats, tmpSam) & 1)
        ) { /*If the read is under the min quality, discard*/
            if(rmSupAlnBl & 1)
            { /*If doing a chimera removal*/
                /*make the new alignment (not dupicate) old alignment*/
                tmpSam = newSam;
                newSam = oldSam;
                oldSam = tmpSam;
            } /*If doing a chimera removal*/

            /*Remove old stats, not keeping and get next alinment*/
            blankSamEntry(newSam);
            funErrUC = readSamLine(newSam, stdinFILE);
            continue;
        } /*If the read is under the min quality, discard*/

        /**************************************************************\
        * Fun-1 Sec-5 Sub-3: Set up bin file names
        \**************************************************************/

        /*Grab the reference id*/
        cStrCpInvsDelm(refIdCStr, tmpSam->refCStr);

        /*Build the bin file name*/
        tmpCStr = cStrCpInvsDelm(binFileCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--");
        tmpCStr = cStrCpInvsDelm(tmpCStr, refIdCStr);
        cStrCpInvsDelm(tmpCStr, ".fastq"); /*Add in fastq ending*/

        tmpCStr = cStrCpInvsDelm(statFileCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--");
        tmpCStr = cStrCpInvsDelm(tmpCStr, refIdCStr);
        cStrCpInvsDelm(tmpCStr, "--stats.tsv");/*Add stats file ending*/

        /**************************************************************\
        * Fun-1 Sec-5 Sub-4: Check if can open the bin & stats file
        \**************************************************************/

        statFILE = fopen(statFileCStr, "a");

        if(statFILE == 0)
        { /*If can not open the stats file*/
            freeBinTree(&binTree);
            pclose(stdinFILE);
            *errUC = 4;
            return 0;
        } /*If can not open the stats file*/

        fqBinFILE = fopen(binFileCStr, "a");

        if(fqBinFILE == 0)
        { /*If can not open the fastq binning file*/
            freeBinTree(&binTree);
            pclose(stdinFILE);
            fclose(statFILE);
            *errUC = 8;
            return 0;
        } /*If can not open the stats file*/

        /**************************************************************\
        * Fun-1 Sec-5 Sub-5: Append reads to file & add to tree
        \**************************************************************/

        tmpBin =
            insBinIntoTree(
                refIdCStr,
                binFileCStr,       /*Fastq file for the bin*/
                statFileCStr,      /*Stats file for the bin*/
                &binTree, /*Root of bin tree*/
                binStack  /*Stack to use in rebalencing tree*/
        ); /*Find or add bin to tree*/

        if(tmpBin == 0)
        { /*If a memory error occured*/
            freeBinTree(&binTree);
            fclose(statFILE);
            fclose(fqBinFILE);
            pclose(stdinFILE);
            *errUC = 64;
            return 0;
        } /*If a memory error occured*/

        /**************************************************************\
        * Fun-1 Sec-5 Sub-6: Print out the stats and fastq entry
        \**************************************************************/

        if(tmpBin->numReadsULng == 1) /*This is a new bin*/
            printStatsHeadUC = 1;  /*Add header to stats file*/
        else
            printStatsHeadUC = 0; /*Else do not print the header*/

        /*Print out the old sam entry (is not a duplicate)*/
        /*Add sequence and stats to their files*/
        samToFq(tmpSam, fqBinFILE); /*Print sequence to fastq file*/

        /*Print the stats to its bin file*/
        printSamStats(tmpSam, &printStatsHeadUC, statFILE);

        fclose(fqBinFILE);
        fclose(statFILE);

        fqBinFILE = 0;  /*So program knows that no file is open*/
        statFILE = 0;   /*So program knows that no file is open*/
            
        blankSamEntry(tmpSam); /*Remove old stats in sam file*/
 
        if(rmSupAlnBl & 1)
        { /*If doing a chimera removal*/
            /*Swap pointers around (the new alignent becomes the old)*/
            tmpSam = newSam;
            newSam = oldSam;
            oldSam = tmpSam;
        } /*If doing a chimera removal*/

        /*Read in the next line*/
        funErrUC = readSamLine(newSam, stdinFILE);
    } /*While their is a samfile entry to read in*/

    pclose(stdinFILE);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-6: Print out the last read
    ^   fun-1 sec-6 sub-1: Check if should keep last read
    ^   fun-1 sec-6 sub-2: Prepare the fastq & stats file names
    ^   fun-1 sec-6 sub-3: Check if can open fastq & stats file
    ^   fun-1 sec-6 sub-4: Add read to the tree of bins
    ^   fun-1 sec-6 sub-5: Print out the read & its stats to the bin
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-1 Sec-6 Sub-1: Check if should keep last read
    \******************************************************************/

    if(rmSupAlnBl & 1 && dupBL == 0)
    { /*If I have a final read to print out*/
        if(checkRead(minStats, oldSam) == 0 ||
          !(checkIfKeepRead(minStats, oldSam) & 1)
        ) return binTree; /*If is a low quality read*/

        /**************************************************************\
        * Fun-1 Sec-6 Sub-2: Prepare the fastq & stats file names
        \**************************************************************/

        /*Grab the reference id*/
        cStrCpInvsDelm(refIdCStr, oldSam->refCStr);

        /*Build the bin file name*/
        tmpCStr = cStrCpInvsDelm(binFileCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--");
        tmpCStr = cStrCpInvsDelm(tmpCStr, refIdCStr);
        cStrCpInvsDelm(tmpCStr, ".fastq"); /*Add in fastq ending*/

        tmpCStr = cStrCpInvsDelm(statFileCStr, prefixCStr);
        tmpCStr = cStrCpInvsDelm(tmpCStr, "--");
        tmpCStr = cStrCpInvsDelm(tmpCStr, refIdCStr);
        cStrCpInvsDelm(tmpCStr, "--stats.tsv");/*Add stats file ending*/

        /**************************************************************\
        * Fun-1 Sec-6 Sub-3: Check if can open fastq & stats file
        \**************************************************************/

        statFILE = fopen(statFileCStr, "a");

        if(statFILE == 0)
        { /*If can not open the stats file*/
            freeBinTree(&binTree);
            pclose(stdinFILE);
            *errUC = 4;
            return 0;
        } /*If can not open the stats file*/

        fqBinFILE = fopen(binFileCStr, "a");

        if(fqBinFILE == 0)
        { /*If can not open the fastq binning file*/
            freeBinTree(&binTree);
            pclose(stdinFILE);
            fclose(statFILE);
            *errUC = 8;
            return 0;
        } /*If can not open the stats file*/

        /**************************************************************\
        * Fun-1 Sec-6 Sub-4: Add read to the tree of bins
        \**************************************************************/

        tmpBin =
            insBinIntoTree(
                refIdCStr,
                prefixCStr,       /*Fastq file for the bin*/
                statFileCStr,      /*Stats file for the bin*/
                &binTree, /*Root of bin tree*/
                binStack  /*Stack to use in rebalencing tree*/
        ); /*Find or add bin to tree*/

        if(tmpBin == 0)
        { /*If a memory error occured*/
            freeBinTree(&binTree);
            fclose(statFILE);
            fclose(fqBinFILE);
            pclose(stdinFILE);
            *errUC = 64;
            return 0;
        } /*If a memory error occured*/

        /**************************************************************\
        * Fun-1 Sec-6 Sub-5: Print out the read & its stats to the bin
        \**************************************************************/

        /*Print out the old sam entry (is not a duplicate)*/
        /*Add sequence and stats to their files*/
        samToFq(oldSam, fqBinFILE); /*Print sequence to fastq file*/

        /*Print the stats to its bin file*/
        printSamStats(oldSam, &printStatsHeadUC, statFILE);

        fclose(fqBinFILE);
        fclose(statFILE);
    } /*If I have a final read to print out*/

    return binTree;
} /*binReads*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|        - 1: if succeded
|    Modifies:
|        - fastq in binClust->fqPathCStr to be the fastq for the cluster
|        - fastq in binTree->fqPathCStr to not have clustered reads
\---------------------------------------------------------------------*/
uint8_t binReadToCon(
    const uint8_t *clustUChar,      /*Cluster on*/
    struct readBin *binTree,        /*Bin working on*/
    struct readBin *binClust,       /*Bin to assign reads & consensus*/
    struct samEntry *samStruct,     /*To hold temporary input*/
    struct minAlnStats *minStats,   /*Min stats needed to keep a read*/
    char *threadsCStr            /*Number threads to use with Minimap2*/
) /*Maps reads to consensus and keeps reads that meet user criteria*/
{ /*binReadToCon*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC: binReadToCon
    '    fun-2 sec-1: Variable declerations
    '    fun-2 sec-2: Set defaults & run minimap2
    '    fun-2 sec-3: Make temporary files & open the old stats file
    '    fun-2 sec-4: Check each minimap2 alignment to see if keeping
    '    fun-2 sec-5: Clean up and rename files
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t errUChar = 0;
    uint8_t oneUChar = 1;
    uint8_t zeroUChar = 0;
    uint8_t headBool = 0; /*Tells if frist round in stats file*/

    char minimap2CmdCStr[2048];  /*Holds minimap2 command to run*/
    char *tmpCStr = 0;
    char *tmpStatsCStr = "2023-01-17-1239-stats-tmp-01928327456123.tsv";
    char *tmpFqCStr = "2023-01-17-1239-fastq-tmp-019283274561234.fastq";

    struct samEntry *zeroSam = 0; /*Just to tell no reference struct*/

    FILE *tmpStatsFILE = 0; /*Stats keeping*/
    FILE *clustFILE = 0;/*Holds reads that mapped to the consensuses*/
    FILE *otherBinFILE = 0;/*Holds reads that did not map*/
    FILE *stdinFILE = 0;   /*Holds minimap2 output*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Set defaults & run minimap2
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    binTree->numReadsULng = 0;  /*Reseting size after binning*/
    binClust->numReadsULng = 0; /*For counting number reads in bin*/

    tmpCStr = cStrCpInvsDelm(minimap2CmdCStr, minimap2CMD);
    tmpCStr = cpParmAndArg(tmpCStr, "-t", threadsCStr);
    cpParmAndArg(
        tmpCStr,
        binClust->consensusCStr,
        binTree->fqPathCStr
    ); /*Add the file names to the minimap2 command*/

    stdinFILE = popen(minimap2CmdCStr, "r");

    /*Remove the old stats data in the structures*/
    blankSamEntry(samStruct);

    /*Read First line so can check if errored out*/
    errUChar = readSamLine(samStruct, stdinFILE);

    if(*samStruct->samEntryCStr != '@')
    { /*If their is no header*/
        pclose(stdinFILE);
        return 2; /*Minimap2 failed*/        
    } /*If their is no header*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-3: Make temporary files & open the old stats file
    ^     fun-2 sec-3 sub-1: open the temporary files & bin stat file
    ^     fun-2 sec-3 sub-2: Build the clusters fastq name
    ^     fun-2 sec-3 sub-3: read the first line of the stats file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-2 Sec-3 Sub-1: open the temporary files & bin stat file
    \******************************************************************/

    tmpStatsFILE = fopen(tmpStatsCStr, "w"); /*Open the temp file*/

    otherBinFILE = fopen(tmpFqCStr, "w"); /*file for discarded reads*/

    /******************************************************************\
    * Fun-2 Sec-3 Sub-2: Build the clusters fastq name
    \******************************************************************/

    /*Copy bin fastq name to cluster fastq name*/
    strcpy(binClust->fqPathCStr, binTree->fqPathCStr);

    tmpCStr = binClust->fqPathCStr;

    while(*tmpCStr != '\0')
        ++tmpCStr;

    tmpCStr -= 6; /*get to '.' in ".fastq"*/

    /*Add --*/
    *tmpCStr = '-';
    ++tmpCStr;
    *tmpCStr = '-';
    ++tmpCStr;

    strcpy(tmpCStr, "cluster-");
    tmpCStr += 8;

    /*Add the cluster number to the fastq file name*/
    tmpCStr = uCharToCStr(tmpCStr, *clustUChar);
    strcpy(tmpCStr, ".fastq"); /*Add fastq ending to cluster fastq*/

    clustFILE = fopen(binClust->fqPathCStr, "w"); /*make cluster fastq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-4: Check each minimap2 alignment to see if keeping
    ^    fun-2 sec-4 sub-1: Check if is a header
    ^    fun-2 sec-4 sub-2: Score read
    ^    fun-2 sec-4 sub-3: Check if Score meets requirements
    ^    fun-2 sec-4 sub-4: Read next line from minimap2 & stats file
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /******************************************************************\
    * Fun-2 Sec-4 Sub-1: Check if is a header
    \******************************************************************/

    while(errUChar & 1)
    { /*While their is a samfile entry to read in*/

        if(*samStruct->samEntryCStr == '@')
        { /*If was a header*/
            /*Read in next entry*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*If was a header*/

        if(samStruct->flagUSht & 256 || samStruct->flagUSht & 2048)
        { /*If was a secondary or supplementary alignement, ignore*/
            /*Read in next entry*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*If was a secondary or supplementary alignement, ignore*/

        findQScores(samStruct); /*Find the Q-scores*/

        if(samStruct->flagUSht & 4)
        { /*Make sure the read mapped to something*/
            samToFq(samStruct, otherBinFILE);
            printSamStats(samStruct, &headBool, tmpStatsFILE);
            ++binTree->numReadsULng; /*Update total scores in bin*/

            /*Read in next entry*/
            blankSamEntry(samStruct);
            errUChar = readSamLine(samStruct, stdinFILE);
            continue; /*Is a header line, move to next line in file*/
        } /*Make sure the read mapped to something*/

        /**************************************************************\
        * Fun-2 Sec-4 Sub-2: Score read
        \**************************************************************/

        scoreAln(
            minStats,
            samStruct,
            zeroSam,    /*Do not use reference for dels (no q-score)*/ 
            &oneUChar,   /*Mapped read has Q-score*/
            &zeroUChar   /*Reference has no q-score entry*/
        ); /*Score the alignment*/

        /**************************************************************\
        * Fun-2 Sec-4 Sub-3: Check if Score meets requirements
        \**************************************************************/

        if(
            samStruct->mapqUChar < minStats->minMapqUInt ||
            !(checkIfKeepRead(minStats, samStruct) & 1)
        ) { /*If read does not belong in this cluster*/
            /*Print out fastq & stats to the temporary files
              These will be made into the bins fastq files later*/
            samToFq(samStruct, otherBinFILE);
            printSamStats(samStruct, &headBool, tmpStatsFILE);
                /*Need the Q-scores for future clustering steps 
                  Other stats not big deal
                  I would like to save the orginal stats, but the order
                    is different from the fastq file
                */

            ++binTree->numReadsULng; /*Update total scores in bin*/
        } /*If read does not belong in this cluster*/

        else
        { /*else teh read belongs to the cluster*/
            samToFq(samStruct, clustFILE); /*Print read to cluster fq*/
            ++binClust->numReadsULng; /*Adding another read to the bin*/
        } /*else teh read belongs to the cluster*/

        /**************************************************************\
        * Fun-2 Sec-4 Sub-4: Read next line from minimap2 & stats file
        \**************************************************************/

        /*Read in the next line*/
        blankSamEntry(samStruct);
        errUChar = readSamLine(samStruct, stdinFILE);
    } /*While their is a samfile entry to read in*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-5: Clean up and rename files
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    pclose(stdinFILE);
    fclose(clustFILE);
    fclose(otherBinFILE);
    fclose(tmpStatsFILE);

    /*Remove the bins old fastq & stats files*/
    remove(binTree->fqPathCStr);
    remove(binTree->statPathCStr);

    /*Assign the temporary fastq & stats files to the bin*/
    rename(tmpFqCStr, binTree->fqPathCStr);
    rename(tmpStatsCStr, binTree->statPathCStr);

    return 1; /*No errors*/
} /*binReadToCon*/
