/*##############################################################################
# Name: ScoreReadsScoringWrappers
# Use: holds wrappers that check read scores and call score read functions
##############################################################################*/

#include "scoreReadsScoringWrappers.h"

/*##############################################################################
# ScoreReadsScoringWrappers TOC:
#   fun-3: scoreReadsFile: get & print read scores with input from a file
#   fun-4: scoreReadsStdin: get & print read scores with input from stdin
#   fun-5: refScoreReadsStdin: use reference when scoring reads with stdin input
#   fun-6: refScoreReadsFile: use reference when scoring reads with file input
#   fun-7: initForScoring: initate struct and allocate memory for scoreReads Var
#   fun-8: scores alignment of a single aligment in a sam file
#   fun-9: Uses a reference & scores a single aligment in a sam file
#   fun-10: checkRefEntry: Checks if reference meets min requirments
##############################################################################*/

/*##############################################################################
# Name: scoreReadsFile
# Use: Scores reads using input from a sam file
# Input:
#    samFile: Sam file the reads are in (FILE)
#    maxLineLenULng: Max line length in sam file (unsigned long)
#    minReadStats: Struct with Min stats to keep a read (minStats pointer)
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Modifies: Closes samFile
##############################################################################*/
void scoreReadsFile(FILE *samFile,
                    unsigned long maxLineLenULng,
                    minStats *minReadStats)
{ /*scoreReadsFile*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3: scoreReads
    #    fun-3 sec-1: Declare variables & initalize structers
    #    fun-3 sec-2: Score each read in the file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Declare variables & initalize structers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *samLineCStr = 0,            /*Line of samfile working on*/
         *oldSamLineCStr = 0,         /*Hold the older samfile line with seq*/
         printHeadChar = 1;           /*Tells to print out header*/

    /*Structers are from scoreReadsStructers.h*/
    samEntry samStruct,              /*samfile stats and points to seq*/
             oldSamStruct;           /*samfile stats and points to seq*/

    initForScoring(&samStruct, &samLineCStr, maxLineLenULng);
    initForScoring(&oldSamStruct, &oldSamLineCStr, maxLineLenULng);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Score each read in the file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    while(fgets(samLineCStr, maxLineLenULng, samFile)) /*IO lib inFileIO.h*/
         scoreAligment(samLineCStr,     /*line with sam entry (alignment)*/
                       &samStruct,      /*will have index of samLineCStr*/
                       oldSamLineCStr,  /*line with last sam entry*/
                       &oldSamStruct,   /*holds index of last sam entry*/
                       minReadStats,    /*min stats to keep an alignment*/
                       &printHeadChar   /*1: print header, 0 do not*/
        ); /*Get and print the score for a single aligment*/
 
    fclose(samFile);
    free(samLineCStr);
    free(oldSamLineCStr);
    return;
} /*scoreReadsFile*/

/*##############################################################################
# Name: scoreReadsStdin
# Use: Scores reads using input from stdin
# Input:
#    maxLineLenULng: Max line length in sam file (unsigned long)
#    minReadStats: Struct with Min stats to keep a read (minStats pointer)
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Modifies: Closes samFile
##############################################################################*/
void scoreReadsStdin(unsigned long maxLineLenULng, minStats *minReadStats)
{ /*scoreReadsStdin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4: scoreReads
    #    fun-4 sec-1: Declare variables & intalize structers
    #    fun-4 sec-2: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Declare variables & initalize structers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *samLineCStr = 0,            /*Line of samfile working on*/
         *oldSamLineCStr = 0,         /*Hold the older samfile line with seq*/
         printHeadChar = 1;           /*Tells to print out header*/

    /*Structers are from scoreReadsStructers.h*/
    samEntry samStruct,              /*samfile stats and points to seq*/
             oldSamStruct;           /*samfile stats and points to seq*/

    initForScoring(&samStruct, &samLineCStr, maxLineLenULng);
    initForScoring(&oldSamStruct, &oldSamLineCStr, maxLineLenULng);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(fgets(samLineCStr, maxLineLenULng, stdin)) /*IO lib inFileIO.h*/
        scoreAligment(samLineCStr,     /*line with sam entry (alignment)*/
                      &samStruct,      /*will have index of samLineCStr*/
                      oldSamLineCStr,  /*line with last sam entry*/
                      &oldSamStruct,   /*holds index of last sam entry*/
                      minReadStats,    /*min stats to keep an alignment*/
                      &printHeadChar   /*1: print header, 0 do not*/
        ); /*Get and print the score for a single aligment*/

    free(samLineCStr);
    free(oldSamLineCStr);
    return;
} /*scoreReadsStdin*/

/*##############################################################################
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Frees: seqCStr & qCStr in refEntry struct
##############################################################################*/
char refScoreReadsStdin(
    unsigned long maxLineLenULng,   /*Max line length to read in sam file*/
    struct minStats *minReadStats,  /*Min thresholds to keep an aligment*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    struct samEntry *refEntry       /*Refence to compare read to*/
) /*Score each alignment using the provided reference & cigar*/
{ /*refScoreReadsStdin*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 TOC:
    #    fun-5 sec-1: Declare variables & intalize structers
    #    Fun-6 Sec-2: initate variables & allocate memory
    #    fun-5 sec-3: Find reference Q-scores & length & decide if should keep
    #    fun-5 sec-5: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1: Declare variables & initalize structers
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *samLineCStr = 0,            /*Line of samfile working on*/
        *oldSamLineCStr = 0,         /*Hold the older samfile line with seq*/
        printHeadChar = 1;           /*Tells to print out header*/

    /*Structers are from scoreReadsStructers.h*/
    struct samEntry
        samStruct,                   /*samfile stats and points to seq*/
        oldSamStruct;                /*samfile stats and points to seq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: initate variables & allocate memory
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initForScoring(
        &samStruct,
        &samLineCStr,
        maxLineLenULng
    ); /*Initate samStruct to have all variables set to 0*/

    if(samLineCStr == 0)
    { /*If malloc failed to grab memory*/
        fprintf(
            stderr,
            "Not enough memory: fun-6: refScoreReadsFile: scoreReads.c 658\n"
        );
        return 0;
    } /*If malloc failed to grab memory*/

    initForScoring(
        &oldSamStruct,
        &oldSamLineCStr,
        maxLineLenULng
    ); /*Set all values in oldSamStruct to 0*/

    if(oldSamLineCStr == 0)
    { /*If malloc failed to grab memory*/
        free(samLineCStr);
        fprintf(
            stderr,
            "Not enough memory: fun-6: refScoreReadsFile: scoreReads.c 673\n"
        );
        return 0;
    } /*If malloc failed to grab memory*/


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-3: Find reference Q-scores & length & decide if should keep
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(checkRefEntry(refEntry, minReadStats) == 4)
    { /*If refernce is beneath min qaulity*/
        free(samLineCStr);
        free(oldSamLineCStr);
        return 4;
    } /*If refernce is beneath min qaulity*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-4: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(fgets(samLineCStr, maxLineLenULng, stdin))
    { /*While there is a line in stdin to score*/
        refScoreAligment(
            refEntry,       /*Struct with reference sequence & Q-score line*/
            refForDelChar,   /*Only using reference for deletions?*/
            samLineCStr,     /*line with sam entry (alignment)*/
            &samStruct,      /*will have index of samLineCStr*/
            oldSamLineCStr,  /*line with last sam entry*/
            &oldSamStruct,   /*holds index of last sam entry*/
            minReadStats,    /*min stats to keep an alignment*/
            &printHeadChar   /*1: print header, 0 do not*/
        ); /*Get and print the score for a single aligment*/
    } /*While there is a line in stdin to score*/

    free(samLineCStr);
    free(oldSamLineCStr);
    return 2;
} /*refScoreReadsStdin*/

/*##############################################################################
# Output:
#    stdout: Prints the scores of the reads to stdout
#    Frees: seqCStr & qCStr in refEntry struct
#    Closes: samFile
#    Returns: 4: If reference beneath min quality
#             2: if sucess
#             0: If malloc failed
##############################################################################*/
char refScoreReadsFile(
    FILE *samFile,                  /*Sam file with aligments to score*/
    unsigned long maxLineLenULng,   /*Max line length to read in sam file*/
    struct minStats *minReadStats,  /*Min thresholds to keep an aligment*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    struct samEntry *refEntry       /*Refence to compare read to*/
) /*Score each alignment using the provided reference & cigar*/
{ /*refScoreReadsFile*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 TOC:
    #    fun-6 sec-1: Declare variables
    #    fun-6 sec-2: initaite variables & allocate memory
    #    fun-6 sec-3: Find reference Q-scores & length & decide if should keep
    #    fun-6 sec-4: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: Declare variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char
        *samLineCStr = 0,            /*Line of samfile working on*/
        *oldSamLineCStr = 0,         /*Hold the older samfile line with seq*/
        printHeadChar = 1;           /*Tells to print out header*/

    /*Structers are from scoreReadsStructers.h*/
    struct samEntry
        samStruct,                   /*samfile stats and points to seq*/
        oldSamStruct;                /*samfile stats and points to seq*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: initaite variables & allocate memory
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    initForScoring(
        &samStruct,
        &samLineCStr,
        maxLineLenULng
    ); /*Initate samStruct to have all variables set to 0*/

    if(samLineCStr == 0)
    { /*If malloc failed to grab memory*/
        fprintf(
            stderr,
            "Not enough memory: fun-6: refScoreReadsFile: scoreReads.c 774\n"
        );
        fclose(samFile);
        return 0;
    } /*If malloc failed to grab memory*/

    initForScoring(
        &oldSamStruct,
        &oldSamLineCStr,
        maxLineLenULng
    ); /*Set all values in oldSamStruct to 0*/

    if(oldSamLineCStr == 0)
    { /*If malloc failed to grab memory*/
        free(samLineCStr);
        fclose(samFile);
        fprintf(
            stderr,
            "Not enough memory: fun-6: refScoreReadsFile: scoreReads.c 799\n"
        );
        return 0;
    } /*If malloc failed to grab memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Find reference Q-scores & length & decide if should keep
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(checkRefEntry(
           refEntry,
           minReadStats
       ) == 4
    ) { /*If refernce is beneath min qaulity*/
        free(samLineCStr);
        free(oldSamLineCStr);
        fclose(samFile);
        return 4;
    } /*If refernce is beneath min qaulity*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-4: Score each read provided by stdin
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(
        fgets(
            samLineCStr,
            maxLineLenULng,
            samFile
        )
    ) { /*While there is a line in stdin to score*/
        refScoreAligment(
            refEntry,       /*Struct with reference sequence & Q-score line*/
            refForDelChar,   /*Only using reference for deletions?*/
            samLineCStr,     /*line with sam entry (alignment)*/
            &samStruct,      /*will have index of samLineCStr*/
            oldSamLineCStr,  /*line with last sam entry*/
            &oldSamStruct,   /*holds index of last sam entry*/
            minReadStats,    /*min stats to keep an alignment*/
            &printHeadChar   /*1: print header, 0 do not*/
        ); /*Get and print the score for a single aligment*/
    } /*While there is a line in stdin to score*/

    free(samLineCStr);
    free(oldSamLineCStr);
    fclose(samFile);
    return 2;
} /*refScoreReadsFile*/

/*##############################################################################
# Output:
#    Modifies: samStruct to have variables set to 0
#    Allocates: memory to samLIneCStr (if malloc failes, returns 0)
##############################################################################*/
void initForScoring(samEntry *samStruct, /*struct to set variables to 0 in*/
                    char **samLineCStr,   /*char pointer to allocate memory to*/
                    unsigned long maxLineLenULng /*size of memory to allocate*/
) /*initate structs and allocate memory for scoreReads variables*/
{ /*intalizeForScoring*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1 TOC: initate struct and allocate memory for scoreReads Var
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blankSamFileEntry(samStruct); /*blank my stats structure*/

    *samLineCStr = malloc(sizeof(char) * (maxLineLenULng + 1));

    if(samLineCStr == 0)
        return;

    *(*samLineCStr + maxLineLenULng) = '\0'; /*make into null string*/

    return;
} /*intalizeForScoring*/

/*##############################################################################
# Use: scores and determins if should keep a signle aligment in a sam file
# Output: 
#    stdout: if alignment in sam file is kept, prints aligment stats to stdout
#    Modifies: When move to new set of alignments, modifies oldSamEntry and 
#        oldSamLineCStr to have the new entry
##############################################################################*/
void scoreAligment(char *samLineCStr,     /*line with sam entry (alignment)*/
                  samEntry *samStruct,    /*will have index of samLineCStr*/
                  char *oldSamLineCStr,   /*line with last sam entry*/
                  samEntry *oldSamStruct, /*holds index of last sam entry*/
                  minStats *minReadStats, /*min stats to keep an alignment*/
                  char *printHeadChar     /*1: print header, 0 do not*/
)
{ /*scoreAligment*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 TOC: scores alignment of a single aligment in a sam file
    #    fun-8 sec-1: Check if is a valid aligment (not header/no reference )
    #    fun-8 sec-2: Check if new read, set ptrs to sequence & Q-score entries
    #    fun-8 sec-3: check if alignment meets min thresholds before scoring
    #    fun-8 sec-4: score the read
    #    fun-4 sec-5: check if alignemnt scores meet min thresholds
    #    fun-4 sec-6: Print out the alignemnt score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1: Check if is a valid aligment (not header/no reference )
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*samLineCStr == '@')
        return; /*This is a header line*/

    /*Remove old stats in sam file*/
    blankSamFileEntry(samStruct); /*from scoreReadsStructers.h*/

    /*Find entry lines (query, flag, reference, mapq, cigar, seq, q-score*/
    if(indexSamEntry(samLineCStr, samStruct) < 1)
        return; /*invalid sam file entry*/
        /*only resets Q-score and sequence ptr when entry != '*'*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-2: Check if new read, set ptrs to sequence & Q-score entries
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*(*samStruct).seqCStr) != '*')
    { /*If starting a new sequence*/
        (*samStruct).readLenULng = findQScores((*samStruct).qCStr,
                                               &(*samStruct).meanQDbl,
                                               &(*samStruct).medianQDbl
        ); /*Find the Q-scores for the entire entry*/
      
        deepCpSamFileEntry(samStruct,
                           oldSamStruct,
                           samLineCStr,
                           oldSamLineCStr
        ); /*Copy the entry into a new structer and C-string*/
    } /*If starting a new sequence*/

    else /*Not a new seqence, copy sequence data*/
        cpSamFileEntry(oldSamStruct, samStruct);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-3: check if alignment meets min thresholds before scoring
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((*samStruct).flagUInt & 4)
         return; /*move onto the next entry (not mapped)*/

    if((*oldSamStruct).mapqUInt < (*minReadStats).minMapqUInt)
        return; /*move onto the next entry (has to low mapq)*/
            /*Minimap2 only gives MAPQ for best match*/

    if((*samStruct).medianQDbl < (*minReadStats).minMedianQDbl &&
       *(*samStruct).qCStr !='*')
        return; /*move onto the next entry (median Q-score to low)*/

    if((*samStruct).meanQDbl < (*minReadStats).minMeanQDbl &&
       *(*samStruct).qCStr !='*')
        return; /*move onto the next entry (mean Q-score to low)*/

    if((*samStruct).readLenULng > (*minReadStats).maxReadLenULng &&
       (*minReadStats).maxReadLenULng > 0)
        return; /*move onto the next entry (read is to long)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-4: score the read
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Score reads and get the aligned stats*/
    scoreRead(minReadStats, samStruct);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-5: check if alignemnt scores meet min thresholds
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((*samStruct).medianAligQDbl < (*minReadStats).minAlignedMedianQDbl &&
       *(*samStruct).qCStr !='*')
        return; /*move onto the next entry (aligned median Q-score low)*/

    if((*samStruct).meanAligQDbl < (*minReadStats).minAlignedMeanQDbl &&
       *(*samStruct).qCStr != '*')
        return; /*move onto the next entry (aligned mean Q-score to low)*/

    if((*samStruct).readAligLenULng < (*minReadStats).minReadLenULng)
        return; /*Read to to short*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-6: Print out the alignemnt score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    stdoutReadScore(samStruct, printHeadChar); /*print entry stats*/

    return;
} /*scoreAligment*/

/*##############################################################################
# Output: 
#    stdout: if alignment in sam file is kept, prints aligment stats to stdout
#    Modifies: When move to new set of alignments, modifies oldSamEntry and 
#        oldSamLineCStr to have the new entry
##############################################################################*/
char refScoreAligment(
    struct samEntry *refEntry,     /*Holds refence sequence & Q-score lines*/
    char refForDelChar,             /*Marks if only using ref for deltions*/
    char *samLineCStr,             /*line with sam entry (alignment)*/
    struct samEntry *samStruct,    /*will have index of samLineCStr*/
    char *oldSamLineCStr,          /*line with last sam entry*/
    struct samEntry *oldSamStruct, /*holds index of last sam entry*/
    struct minStats *minReadStats, /*min stats to keep an alignment*/
    char *printHeadChar     /*1: print header, 0 do not*/
) /*Scores alignments & after scoring determins if should keep alignment*/
{ /*scoreAligment*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 TOC: scores alignment of a single aligment in a sam file
    #    fun-9 sec-1: Check if is a valid aligment (not header/no reference )
    #    fun-9 sec-2: Check if new read, set ptrs to sequence & Q-score entries
    #    fun-9 sec-3: check if alignment meets min thresholds before scoring
    #    fun-9 sec-4: score the read
    #    fun-4 sec-5: check if alignemnt scores meet min thresholds
    #    fun-4 sec-6: Print out the alignemnt score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: Check if is a valid aligment (not header/no reference )
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*samLineCStr == '@')
        return 1; /*This is a header line*/

    /*Remove old stats in sam file*/
    blankSamFileEntry(samStruct); /*from scoreReadsStructers.h*/

    /*Find entry lines (query, flag, reference, mapq, cigar, seq, q-score*/
    if(indexSamEntry(samLineCStr, samStruct) < 1)
        return 1; /*invalid sam file entry*/
        /*only resets Q-score and sequence ptr when entry != '*'*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Check if new read, set ptrs to sequence & Q-score entries
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*(*samStruct).seqCStr) != '*')
    { /*If starting a new sequence*/
        (*samStruct).readLenULng = findQScores((*samStruct).qCStr,
                                               &(*samStruct).meanQDbl,
                                               &(*samStruct).medianQDbl
        ); /*Find the Q-scores for the entire entry*/
      
        deepCpSamFileEntry(samStruct,
                           oldSamStruct,
                           samLineCStr,
                           oldSamLineCStr
        ); /*Copy the entry into a new structer and C-string*/
    } /*If starting a new sequence*/

    else /*Not a new seqence, copy sequence data*/
        cpSamFileEntry(oldSamStruct, samStruct);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: check if alignment meets min thresholds before scoring
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((*samStruct).flagUInt & 4)
         return 1; /*move onto the next entry (not mapped)*/

    if((*oldSamStruct).mapqUInt < (*minReadStats).minMapqUInt)
        return 1; /*move onto the next entry (has to low mapq)*/
            /*Minimap2 only gives MAPQ for best match*/

    if(
        samStruct->medianQDbl < minReadStats->minMedianQDbl &&
        *samStruct->qCStr !='*'
    ) /*If their is a Q-score line & Q-score is under min requrirments*/
        return 1; /*move onto the next entry (median Q-score to low)*/

    if(samStruct->meanQDbl < minReadStats->minMeanQDbl &&
       *samStruct->qCStr !='*')
        return 1; /*move onto the next entry (mean Q-score to low)*/

    if((*samStruct).readLenULng > (*minReadStats).maxReadLenULng &&
       (*minReadStats).maxReadLenULng > 0)
        return 1; /*move onto the next entry (read is to long)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-4: score the read
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Score reads and get the aligned stats*/
    refScoreRead(minReadStats, samStruct, refEntry, refForDelChar);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-5: check if alignemnt scores meet min thresholds
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((*samStruct).medianAligQDbl < (*minReadStats).minAlignedMedianQDbl &&
       *(*samStruct).qCStr !='*')
        return 1; /*move onto the next entry (aligned median Q-score low)*/

    if((*samStruct).meanAligQDbl < (*minReadStats).minAlignedMeanQDbl &&
       *(*samStruct).qCStr != '*')
        return 1; /*move onto the next entry (aligned mean Q-score to low)*/

    if((*samStruct).readAligLenULng < (*minReadStats).minReadLenULng)
        return 1; /*Read to to short*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-6: Print out the alignemnt score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    stdoutReadScore(samStruct, printHeadChar); /*print entry stats*/

    return 1;
} /*scoreAligment*/

/*##############################################################################
# Output:
#    Modifes: refEntry: to have read length and Q-scores of reference
#    Returns: 4: If reference is beneath min stats
#             2: If reference is ok
#             0: If reference does not have Q-score line to check
##############################################################################*/
char checkRefEntry(
    struct samEntry *refEntry,
    struct minStats *minReadStats
) /*Checks if reference meets the min requirments provided by the user*/
{ /*checkRefEntry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 TOC:
    #    fun-10 sec-1: Get Q-scores & check if Q-score entry
    #    fun-10 sec-2: Check if references meets min median Q-score
    #    fun-10 sec-3: Check if references meets min mean Q-score
    #    fun-10 sec-4: Check if references meets the min length requirment
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-10 Sec-1: Get Q-scores & check if Q-score entry
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refEntry->qCStr == 0)
        return 0;             /*Nothing to check*/

    refEntry->readLenULng = findQScores(
                               refEntry->qCStr,
                               &refEntry->meanQDbl,
                               &refEntry->medianQDbl
    ); /*Find the Q-scores for the entire entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Check if references meets min median Q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refEntry->meanQDbl < minReadStats->minMedianQDbl)
    { /*If reference is beneath min quality*/
        fprintf(
            stderr,
            "Reference has lower medain Q-score (%f) than min median Q (%f)\n",
            refEntry->meanQDbl,
            minReadStats->minMeanQDbl
        ); /*Tell user reference is to short*/
        return 4;
    } /*If reference is beneath min quality*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: Check if references meets min mean Q-score
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refEntry->meanQDbl < minReadStats->minMeanQDbl)
    { /*If reference is beneath min quality*/
        fprintf(
            stderr,
            "Reference has lower mean Q-score (%f) than min mean Q (%f)\n",
            refEntry->meanQDbl,
            minReadStats->minMeanQDbl
        ); /*Tell user reference is to short*/
        return 4;
    } /*If reference is beneath min quality*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-4: Check if references meets the min length requirment
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refEntry->readLenULng < minReadStats->minReadLenULng)
    { /*If reference is to short*/
        fprintf(
            stderr,
            "Reference is smaller (%lu) than the min read length of %lu\n",
            refEntry->readLenULng,
            minReadStats->minReadLenULng
        ); /*Tell user reference is to short*/
        return 4;
    } /*If reference is to short*/

    return 2;
} /*checkRefEntry*/
