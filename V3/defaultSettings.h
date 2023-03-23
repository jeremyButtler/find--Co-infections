/*---------------------------------------------------------------------\
| Use:
|    o Holds the default settings and global definitions for
|      find co-infections and smaller pipelines used by
|      find co-infections
| Note:
|    Note all settings here are used, but in some cases are here for
|    settings structure default values with functions.
\---------------------------------------------------------------------*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOD: Start of definitions
'    sec-1: Defaults dealing with variable sizes
'    sec-2: System commands
'    sec-3: General settings
'    sec-4: Specific settings for separate consensus building steps
'    Sec-5: Read to reference mapping and scoreRead default settings
'    sec-6: read to read mapping settings
'    sec-7: read to consensus mapping settings
'    sec-8: consensus to consensus mapping settings
'    sec-9: Majority consensus and alignment settings
'      - The base scoring matrix is long (due to anonymous base support)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef DEFAULTSETTINGS_H
#define DEFAULTSETTINGS_H

/**********************************************************************\
* Sec-1: variable sizes
\**********************************************************************/

#define bitsInFlt sizeof(float) << 3 /*Number bits in a float*/
#define bitsInSht sizeof(short) << 3 /*Number of bits in a short*/
#define Q_ADJUST 33 /*offest to get q-score of 0*/
#define MAX_Q_SCORE 94 /*highest possible Q-score*/

/**********************************************************************\
* Sec-2: System commands
\**********************************************************************/

#define minimap2CMD "minimap2 --eqx --secondary=no -a -x map-ont"
#define raconCMD "racon -m 8 -x 6 -g -8 -w 500"

#define medakaCMD "/bin/bash -c \"source ~/medaka/venv/bin/activate;"
#define medakaCMDEnd "; deactivate\"" /*Telling to quite enviroment*/
#define medCondCMD "eval \"$(conda shell.bash hook)\"; conda activate medaka;"
#define medCondCMDEnd "; conda deactivate"

/*Command for mapping reads to primers*/
#define defMinimap2PrimCMD "minimap2 -k5 -w1 -s 20 -P"

/**********************************************************************\
* Sec-3: General settings
\**********************************************************************/

/*Gereral find co-infections settings*/
#define defVersion 3.20230323  /*The Version number of this program*/
    /*Format is version.yearMonthDay*/
#define defPrefix "out"      /*Default prefix to use*/
#define defThreads "3"       /*Default number of threads to use*/
#define rmReadsWithSupAln 0
    /*1: remove reads with supplementary aligments*/
#define defReadsPerCon 300
    /*Maximum number of reads to use to build a consensus*/
#define minReadsPerBin 100 
    /*Minimum number reads to keep a bin or build a consensus*/
#define defMinPercReads 0.01 /*(1% of all clustered reads)*/
  /*Minimum percentage of all clustered reads needed to keep a cluster*/

#define defSkipBinBl 0      /*Do not skip binning step*/
#define defSkipClustBl 0    /*Do not skip clustering step*/

#define defNumPolish 2      /*Number of times to rebuild the consensus*/
#define defMinConLen 500     /*consusens must be at least 500bp*/

/**********************************************************************\
* Sec-4: Specific settings for separate consensus building steps
\**********************************************************************/

#define defUseMajCon 1    /*1: Use majority in consnesus building*/
/*#define majConMinBaseQ 10  Min Q-score to keep base in majority con*/
#define majConMinBaseQ 7 /*Min Q-score to keep base in majority con*/
/*#define percBasesPerPos 0.4 majority consensus (-maj-con-min-bases)*/
#define percBasesPerPos 0.35 /*majority consensus (-maj-con-min-bases)*/
#define majConMinInsQ 5   /*Min Q-score to keep insertion in majority*/
#define percInsPerPos 0.3 /*% of supporting reads to keep insertion*/

#define defUseRaconCon 0  /*1: Racon in consensus building*/
#define defRoundsRacon 4    /*Default number of racon rounds/consensus*/

#define defUseMedakaCon 0 /*1: Use medaka in consensus building*/
#define defMedakaModel "r941_min_high_g351" /*Model to use with medaka*/
#define defCondaBl 0      /*Default no, but my code will autofind*/

/**********************************************************************\
* Sec-5: Read to reference mapping settings/scoreRead default settings
\**********************************************************************/

#define readRefMinPercSNPs 0.07 /*Min % of SNPs to keep a read*/
#define readRefMinPercDiff 1    /*Min % difference to keep a read*/
#define readRefMinPercDels 1    /*Min % of deletions to keep a read*/
#define readRefMinPercInss 1    /*Min % of insertions to keep a read*/
#define readRefMinPercIndels 1  /*Min % of indels to keep a read*/

#define readRefMapq 20          /*Min mapq to keep read mapped ref*/
#define readRefMinBaseQ 10      /*Min Q-score needed to keep a base*/

/*read to reference mapping quality thresholds*/
#define readRefMinMedQ 10
#define readRefMinMeanQ 10
#define readRefMinAlnMedQ 10
#define readRefMinAlnMeanQ 10

#define readRefMinReadLen 500   /*Min length to keep a read*/
    /*Also the mininum length to keep a consensus length*/
#define readRefMaxReadLen 1000  /*Maximum length to keep a read*/

/*read to reference mapping settings for insertion homopolymers*/
/*Maximum homopolymer size to keep an insertion in*/
#define readRefMaxInsAHomo 1
#define readRefMaxInsTHomo 1
#define readRefMaxInsCHomo 1
#define readRefMaxInsGHomo 1

/*read to reference mapping settings for deletion homopolymers*/
/*Maximum homopolymer size to keep an deletion in*/
#define readRefMaxDelAHomo 0
#define readRefMaxDelTHomo 0
#define readRefMaxDelCHomo 0
#define readRefMaxDelGHomo 0

/**********************************************************************\
* Sec-6: read to read mapping settings
\**********************************************************************/

/*Read mapped to read filters*/
#define readReadMapq 10          /*Min mapq to keep read mapped ref*/
#define readReadMinBaseQ 10      /*Min Q-score needed to keep a base*/

#define readReadMinPercSNPs 0.021 /*Min % of SNPs to keep a read*/
#define readReadMinPercDiff 0.03 /*Min % difference to keep a read*/
#define readReadMinPercDels 1     /*Min % of deletions to keep a read*/
#define readReadMinPercInss 1     /*Min % of insertions to keep a read*/
#define readReadMinPercIndels 1   /*Min % of indels to keep a read*/

/*read to reference mapping settings for insertion homopolymers*/
/*Maximum homopolymer size to keep an insertion in*/
#define readReadMaxInsAHomo 1
#define readReadMaxInsTHomo 1
#define readReadMaxInsCHomo 1
#define readReadMaxInsGHomo 1

/*read to reference mapping settings for deletion homopolymers*/
/*Maximum homopolymer size to keep an deletion in*/
#define readReadMaxDelAHomo 0
#define readReadMaxDelTHomo 0
#define readReadMaxDelCHomo 0
#define readReadMaxDelGHomo 0
#define readReadMinReadLen 500   /*Min length to keep a read*/

/*THSE SETTINGS WILL HAVE LITTE AFFECT IN findCoInft, MOVE ON TO NEXT
  SECTION. They are here for setting default settings*/

#define readReadMinMedQ 13
#define readReadMinMeanQ 13
#define readReadMinAlnMedQ 13
#define readReadMinAlnMeanQ 13

/*Maximum length to keep a read*/
#define readReadMaxReadLen 0  /*0 is for any length*/

/**********************************************************************\
* Sec-7: read to consensus mapping settings
\**********************************************************************/

/*THESE SETTINGS WILL HAVE AN AFFECT*/
/*read to consensus quality thresholds*/
#define readConMapq 20          /*Min mapq to keep read mapped ref*/
#define readConMinBaseQ 10      /*Min Q-score needed to keep a base*/

/*Read mapped to consensus filters*/
#define readConMinPercSNPs 0.015 /*Min % of SNPs to keep a read*/
#define readConMinPercDiff 0.02  /*Min % difference to keep a read*/
#define readConMinPercDels 1     /*Min % of deletions to keep a read*/
#define readConMinPercInss 1     /*Min % of insertions to keep a read*/
#define readConMinPercIndels 1   /*Min % of indels to keep a read*/

/*read to reference mapping settings for insertion homopolymers*/
/*Maximum homopolymer size to keep an insertion in*/
#define readConMaxInsAHomo 1
#define readConMaxInsTHomo 1
#define readConMaxInsCHomo 1
#define readConMaxInsGHomo 1

/*read to reference mapping settings for deletion homopolymers*/
/*Maximum homopolymer size to keep an deletion in*/
#define readConMaxDelAHomo 0
#define readConMaxDelTHomo 0
#define readConMaxDelCHomo 0
#define readConMaxDelGHomo 0
#define readConMinReadLen 500  /*Min length to keep a read*/

/*THESE SETTINGS WILL HAVE LTTILE AFFECT, FOR findCoIfnt MOVE TO NEXT
  SECTION. They are here for setting default settings*/
#define readConMinMedQ 13
#define readConMinMeanQ 13
#define readConMinAlnMedQ 13
#define readConMinAlnMeanQ 13

/*Maximum length to keep a read*/
#define readConMaxReadLen 0 /*0 is any length*/

/**********************************************************************\
* Sec-8: consensus to consensus mapping settings
\**********************************************************************/

/*THESE WILL HAVE AN AFFECT Consensus mapped to consensus filters*/
#define conConMinPercSNPs 0.01  /*Min % of SNPs to keep a consensus*/
#define conConMinPercDiff 0.03  /*Min % difference to keep a consensus*/
#define conConMinPercDels 1     /*% of deletions to keep a consensus*/
#define conConMinPercInss 1     /*% of insertions to keep consensus*/
#define conConMinPercIndels 1   /*Min % of indels to keep a consensus*/

/*Maximum homopolymer size to keep an insertion in*/
#define conConMaxInsAHomo 1
#define conConMaxInsTHomo 1
#define conConMaxInsCHomo 1
#define conConMaxInsGHomo 1

/*Maximum homopolymer size to keep an deletion in*/
#define conConMaxDelAHomo 0
#define conConMaxDelTHomo 0
#define conConMaxDelCHomo 0
#define conConMaxDelGHomo 0

/*THESE SETTINGS WILL BE IGNORED and are here for setting default
  settings*/
#define conConMapq 0           /*Min mapq to keep con mapped ref*/
#define conConMinBaseQ 0       /*Min Q-score needed to keep a base*/

#define conConMinMedQ 0
#define conConMinMeanQ 0
#define conConMinAlnMedQ 0
#define conConMinAlnMeanQ 0

#define conConMinReadLen 0   /*Min length to keep a con*/
#define conConMaxReadLen 0  /*Maximum length to keep a con (is is all)*/

/**********************************************************************\
* Sec-9: Majority consensus and alignment settings
\**********************************************************************/

// alignment kmer settings
#define defPercOverlap 0.1 // 10% overlap between chunks in kmer mapping
#define defMaxChunks 1000  // Work with no more than 1000 chunks at once
#define defMinChunkLen 100 
     // Stop kmer mapping when chunks are at 100 basepairs
#define defMinPercKmer 0.80
    // Reference & query chunks must share 80% of kmers
#define defMinScore 0.90   // Alignments must be 90% similar to keep

// Aligment Needleman-Wunsch settings
#define defGapStartPenalty 0  // adding an indel has a slight penatly
#define defGapExtendPenalty -3 // adding an indel has a slight penatly

// Scoring matrix is EDNAFULL or something close to it. I think I have
// the anonymous bases with the correct scores.

// Needleman-Wunsch alignment scoring matrix (non-anonymous)
#define defAToA 5        // Penalty for an query A to reference A
#define defAToT -4       // Penalty for an query A to reference T
#define defAToG -4       // Penalty for an query A to reference G
#define defAToC -4       // Penalty for an query A to reference C

#define defTToA -4       // Penalty for an query T to reference A
#define defTToT 5        // Penalty for an query T to reference T
#define defTToG -4       // Penalty for an query T to reference G
#define defTToC -4       // Penalty for an query T to reference C

#define defGToA -4       // Penalty for an query G to reference A
#define defGToT -4       // Penalty for an query G to reference T
#define defGToG 5        // Penalty for an query G to reference G
#define defGToC -4       // Penalty for an query G to reference C

#define defCToA -4       // Penalty for an query C to reference A
#define defCToT -4       // Penalty for an query C to reference T
#define defCToG -4       // Penalty for an query C to reference G
#define defCToC 5        // Penalty for an query C to reference C

// Needleman-Wunsch alignment scoring matrix anonymous (A)
#define defAToW 1        // Penalty for an query A to reference W (AT)
#define defAToS -4       // Penalty for an query A to reference S (CG)
#define defAToM 1        // Penalty for an query A to reference M (AC)
#define defAToK -4       // Penalty for an query A to reference K (GT)
#define defAToR 1        // Penalty for an query A to reference R (AG)
#define defAToY -4       // Penalty for an query A to reference Y (CT)
#define defAToB -4       // Penalty for an query A to reference B (CGT)
#define defAToD -1       // Penalty for an query A to reference D (AGT)
#define defAToH -1       // Penalty for an query A to reference H (ACT)
#define defAToV -1       // Penalty for an query A to reference V (ACG)
#define defAToN -2       // Penalty for an query A to reference N (ACGT)
#define defAToX -2       // Penalty for an query A to reference X (ACGT)

#define defWToA 1        // Penalty for an query W to reference A (AT)
#define defSToA -4       // Penalty for an query S to reference A (CG)
#define defMToA 1        // Penalty for an query M to reference A (AC)
#define defKToA -4       // Penalty for an query K to reference A (GT)
#define defRToA 1        // Penalty for an query R to reference A (AG)
#define defYToA -4       // Penalty for an query Y to reference A (CT)
#define defBToA -4       // Penalty for an query B to reference A (CGT)
#define defDToA -1       // Penalty for an query D to reference A (AGT)
#define defHToA -1       // Penalty for an query H to reference A (ACT)
#define defVToA -1       // Penalty for an query V to reference A (ACG)
#define defNToA -2       // Penalty for an query N to reference A (ACGT)
#define defXToA -2       // Penalty for an query X to reference A (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (C)
#define defCToW -4       // Penalty for an query C to reference W (AT)
#define defCToS 1        // Penalty for an query C to reference S (CG)
#define defCToM 1        // Penalty for an query C to reference M (AC)
#define defCToK -4       // Penalty for an query C to reference K (GT)
#define defCToR -4       // Penalty for an query C to reference R (AG)
#define defCToY 1        // Penalty for an query C to reference Y (CT)
#define defCToB -1       // Penalty for an query C to reference B (CGT)
#define defCToD -4       // Penalty for an query C to reference D (AGT)
#define defCToH -1       // Penalty for an query C to reference H (ACT)
#define defCToV -1       // Penalty for an query C to reference V (ACG)
#define defCToN -2       // Penalty for an query C to reference N (ACGT)
#define defCToX -2       // Penalty for an query C to reference X (ACGT)

#define defWToC -4       // Penalty for an query W to reference C (AT)
#define defSToC 1        // Penalty for an query S to reference C (CG)
#define defMToC 1        // Penalty for an query M to reference C (AC)
#define defKToC -4       // Penalty for an query K to reference C (GT)
#define defRToC -4       // Penalty for an query R to reference C (AG)
#define defYToC 1        // Penalty for an query Y to reference C (CT)
#define defBToC -1       // Penalty for an query B to reference C (CGT)
#define defDToC -4       // Penalty for an query D to reference C (AGT)
#define defHToC -1       // Penalty for an query H to reference C (ACT)
#define defVToC -1       // Penalty for an query V to reference C (ACG)
#define defNToC -2       // Penalty for an query N to reference C (ACGT)
#define defXToC -2       // Penalty for an query X to reference C (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (G)
#define defGToW -4       // Penalty for an query G to reference W (AT)
#define defGToS 1        // Penalty for an query G to reference S (CG)
#define defGToM -4       // Penalty for an query G to reference M (AC)
#define defGToK 1        // Penalty for an query G to reference K (GT)
#define defGToR 1        // Penalty for an query G to reference R (AG)
#define defGToY -4       // Penalty for an query G to reference Y (CT)
#define defGToB -1       // Penalty for an query G to reference B (CGT)
#define defGToD -1       // Penalty for an query G to reference D (AGT)
#define defGToH -4       // Penalty for an query G to reference H (ACT)
#define defGToV -1       // Penalty for an query G to reference V (ACG)
#define defGToN -2       // Penalty for an query G to reference N (ACGT)
#define defGToX -2       // Penalty for an query G to reference X (ACGT)

#define defWToG -4       // Penalty for an query G to reference W (AT)
#define defSToG 1        // Penalty for an query G to reference S (CG)
#define defMToG -4       // Penalty for an query G to reference M (AC)
#define defKToG 1        // Penalty for an query G to reference K (GT)
#define defRToG 1        // Penalty for an query G to reference R (AG)
#define defYToG -4       // Penalty for an query G to reference Y (CT)
#define defBToG -1       // Penalty for an query G to reference B (CGT)
#define defDToG -1       // Penalty for an query G to reference D (AGT)
#define defHToG -4       // Penalty for an query G to reference H (ACT)
#define defVToG -1       // Penalty for an query G to reference V (ACG)
#define defNToG -2       // Penalty for an query G to reference N (ACGT)
#define defXToG -2       // Penalty for an query G to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (T)
#define defTToW 1        // Penalty for an query T to reference W (AT)
#define defTToS -4       // Penalty for an query T to reference S (CG)
#define defTToM -4       // Penalty for an query T to reference M (AC)
#define defTToK 1        // Penalty for an query T to reference K (GT)
#define defTToR -4       // Penalty for an query T to reference R (AG)
#define defTToY 1        // Penalty for an query T to reference Y (CT)
#define defTToB -1       // Penalty for an query T to reference B (CGT)
#define defTToD -1       // Penalty for an query T to reference D (AGT)
#define defTToH -1       // Penalty for an query T to reference H (ACT)
#define defTToV -4       // Penalty for an query T to reference V (ACG)
#define defTToN -2       // Penalty for an query T to reference N (ACGT)
#define defTToX -2       // Penalty for an query T to reference X (ACGT)

#define defWToT 1        // Penalty for an query T to reference W (AT)
#define defSToT -4       // Penalty for an query T to reference S (CG)
#define defMToT -4       // Penalty for an query T to reference M (AC)
#define defKToT 1        // Penalty for an query T to reference K (GT)
#define defRToT -4       // Penalty for an query T to reference R (AG)
#define defYToT -4       // Penalty for an query T to reference Y (CT)
#define defBToT -1       // Penalty for an query T to reference B (CGT)
#define defDToT -1       // Penalty for an query T to reference D (AGT)
#define defHToT -1       // Penalty for an query T to reference H (ACT)
#define defVToT -4       // Penalty for an query T to reference V (ACG)
#define defNToT -2       // Penalty for an query T to reference N (ACGT)
#define defXToT -2       // Penalty for an query T to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (W)
#define defWToW -1       // Penalty for an query W to reference W (AT)
#define defWToS -4       // Penalty for an query W to reference S (CG)
#define defWToM -2       // Penalty for an query W to reference M (AC)
#define defWToK -2       // Penalty for an query W to reference K (GT)
#define defWToR -2       // Penalty for an query W to reference R (AG)
#define defWToY -2       // Penalty for an query W to reference Y (CT)
#define defWToB -3       // Penalty for an query W to reference B (CGT)
#define defWToD -1       // Penalty for an query W to reference D (AGT)
#define defWToH -1       // Penalty for an query W to reference H (ACT)
#define defWToV -3       // Penalty for an query W to reference V (ACG)
#define defWToN -1       // Penalty for an query W to reference N (ACGT)
#define defWToX -1       // Penalty for an query W to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (S)
#define defSToW -4       // Penalty for an query S to reference W (AT)
#define defSToS -1       // Penalty for an query S to reference S (CG)
#define defSToM -2       // Penalty for an query S to reference M (AC)
#define defSToK -2       // Penalty for an query S to reference K (GT)
#define defSToR -2       // Penalty for an query S to reference R (AG)
#define defSToY -2       // Penalty for an query S to reference Y (CT)
#define defSToB -1       // Penalty for an query S to reference B (CGT)
#define defSToD -3       // Penalty for an query S to reference D (AGT)
#define defSToH -3       // Penalty for an query S to reference H (ACT)
#define defSToV -1       // Penalty for an query S to reference V (ACG)
#define defSToN -1       // Penalty for an query S to reference N (ACGT)
#define defSToX -1       // Penalty for an query S to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (M)
#define defMToW -2       // Penalty for an query M to reference W (AT)
#define defMToS -2       // Penalty for an query M to reference S (CG)
#define defMToM -1       // Penalty for an query M to reference M (AC)
#define defMToK -4       // Penalty for an query M to reference K (GT)
#define defMToR -2       // Penalty for an query M to reference R (AG)
#define defMToY -2       // Penalty for an query M to reference Y (CT)
#define defMToB -3       // Penalty for an query M to reference B (CGT)
#define defMToD -3       // Penalty for an query M to reference D (AGT)
#define defMToH -1       // Penalty for an query M to reference H (ACT)
#define defMToV -1       // Penalty for an query M to reference V (ACG)
#define defMToN -1       // Penalty for an query M to reference N (ACGT)
#define defMToX -1       // Penalty for an query M to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (K)
#define defKToW -2       // Penalty for an query K to reference W (AT)
#define defKToS -2       // Penalty for an query K to reference S (CG)
#define defKToM -4       // Penalty for an query K to reference M (AC)
#define defKToK -1       // Penalty for an query K to reference K (GT)
#define defKToR -2       // Penalty for an query K to reference R (AG)
#define defKToY -2       // Penalty for an query K to reference Y (CT)
#define defKToB -1       // Penalty for an query K to reference B (CGT)
#define defKToD -1       // Penalty for an query K to reference D (AGT)
#define defKToH -3       // Penalty for an query K to reference H (ACT)
#define defKToV -3       // Penalty for an query K to reference V (ACG)
#define defKToN -2       // Penalty for an query K to reference N (ACGT)
#define defKToX -2       // Penalty for an query K to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (R)
#define defRToW -1       // Penalty for an query R to reference W (AT)
#define defRToS -1       // Penalty for an query R to reference S (CG)
#define defRToM -1       // Penalty for an query R to reference M (AC)
#define defRToK -2       // Penalty for an query R to reference K (GT)
#define defRToR -1       // Penalty for an query R to reference R (AG)
#define defRToY -4       // Penalty for an query R to reference Y (CT)
#define defRToB -3       // Penalty for an query R to reference B (CGT)
#define defRToD -1       // Penalty for an query R to reference D (AGT)
#define defRToH -1       // Penalty for an query R to reference H (ACT)
#define defRToV -3       // Penalty for an query R to reference V (ACG)
#define defRToN -1       // Penalty for an query R to reference N (ACGT)
#define defRToX -1       // Penalty for an query R to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (Y)
#define defYToW -2       // Penalty for an query Y to reference W (AT)
#define defYToS -2       // Penalty for an query Y to reference S (CG)
#define defYToM -2       // Penalty for an query Y to reference M (AC)
#define defYToK -2       // Penalty for an query Y to reference K (GT)
#define defYToR -4       // Penalty for an query Y to reference R (AG)
#define defYToY -1       // Penalty for an query Y to reference Y (CT)
#define defYToB -2       // Penalty for an query Y to reference B (CGT)
#define defYToD -3       // Penalty for an query Y to reference D (AGT)
#define defYToH -1       // Penalty for an query Y to reference H (ACT)
#define defYToV -3       // Penalty for an query Y to reference V (ACG)
#define defYToN -1       // Penalty for an query Y to reference N (ACGT)
#define defYToX -1       // Penalty for an query Y to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (B)
#define defBToW -3       // Penalty for an query B to reference W (AT)
#define defBToS -1       // Penalty for an query B to reference S (CG)
#define defBToM -3       // Penalty for an query B to reference M (AC)
#define defBToK -1       // Penalty for an query B to reference K (GT)
#define defBToR -3       // Penalty for an query B to reference R (AG)
#define defBToY -1       // Penalty for an query B to reference Y (CT)
#define defBToB -1       // Penalty for an query B to reference B (CGT)
#define defBToD -2       // Penalty for an query B to reference D (AGT)
#define defBToH -2       // Penalty for an query B to reference H (ACT)
#define defBToV -2       // Penalty for an query B to reference V (ACG)
#define defBToN -1       // Penalty for an query B to reference N (ACGT)
#define defBToX -1       // Penalty for an query B to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (D)
#define defDToW -1       // Penalty for an query D to reference W (AT)
#define defDToS -3       // Penalty for an query D to reference S (CG)
#define defDToM -3       // Penalty for an query D to reference M (AC)
#define defDToK -1       // Penalty for an query D to reference K (GT)
#define defDToR -1       // Penalty for an query D to reference R (AG)
#define defDToY -3       // Penalty for an query D to reference Y (CT)
#define defDToB -2       // Penalty for an query D to reference B (CGT)
#define defDToD -1       // Penalty for an query D to reference D (AGT)
#define defDToH -2       // Penalty for an query D to reference H (ACT)
#define defDToV -2       // Penalty for an query D to reference V (ACG)
#define defDToN -1       // Penalty for an query D to reference N (ACGT)
#define defDToX -1       // Penalty for an query D to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (H)
#define defHToW -1       // Penalty for an query H to reference W (AT)
#define defHToS -3       // Penalty for an query H to reference S (CG)
#define defHToM -1       // Penalty for an query H to reference M (AC)
#define defHToK -3       // Penalty for an query H to reference K (GT)
#define defHToR -3       // Penalty for an query H to reference R (AG)
#define defHToY -1       // Penalty for an query H to reference Y (CT)
#define defHToB -2       // Penalty for an query H to reference B (CGT)
#define defHToD -2       // Penalty for an query H to reference D (AGT)
#define defHToH -1       // Penalty for an query H to reference H (ACT)
#define defHToV -2       // Penalty for an query H to reference V (ACG)
#define defHToN -1       // Penalty for an query H to reference N (ACGT)
#define defHToX -1       // Penalty for an query H to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (V)
#define defVToW -3       // Penalty for an query V to reference W (AT)
#define defVToS -1       // Penalty for an query V to reference S (CG)
#define defVToM -1       // Penalty for an query V to reference M (AC)
#define defVToK -3       // Penalty for an query V to reference K (GT)
#define defVToR -1       // Penalty for an query V to reference R (AG)
#define defVToY -3       // Penalty for an query V to reference Y (CT)
#define defVToB -2       // Penalty for an query V to reference B (CGT)
#define defVToD -2       // Penalty for an query V to reference D (AGT)
#define defVToH -2       // Penalty for an query V to reference H (ACT)
#define defVToV -1       // Penalty for an query V to reference V (ACG)
#define defVToN -2       // Penalty for an query V to reference N (ACGT)
#define defVToX -1       // Penalty for an query V to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (N)
#define defNToW -1       // Penalty for an query N to reference W (AT)
#define defNToS -1       // Penalty for an query N to reference S (CG)
#define defNToM -1       // Penalty for an query N to reference M (AC)
#define defNToK -1       // Penalty for an query N to reference K (GT)
#define defNToR -1       // Penalty for an query N to reference R (AG)
#define defNToY -1       // Penalty for an query N to reference Y (CT)
#define defNToB -1       // Penalty for an query N to reference B (CGT)
#define defNToD -1       // Penalty for an query N to reference D (AGT)
#define defNToH -1       // Penalty for an query N to reference H (ACT)
#define defNToV -1       // Penalty for an query N to reference V (ACG)
#define defNToN -1       // Penalty for an query N to reference N (ACGT)
#define defNToX -1       // Penalty for an query N to reference X (ACGT)

// Needleman-Wunsch alignment scoring matrix anonymous (X)
#define defXToW -1       // Penalty for an query X to reference W (AT)
#define defXToS -1       // Penalty for an query X to reference S (CG)
#define defXToM -1       // Penalty for an query X to reference M (AC)
#define defXToK -1       // Penalty for an query X to reference K (GT)
#define defXToR -1       // Penalty for an query X to reference R (AG)
#define defXToY -1       // Penalty for an query X to reference Y (CT)
#define defXToB -1       // Penalty for an query X to reference B (CGT)
#define defXToD -1       // Penalty for an query X to reference D (AGT)
#define defXToH -1       // Penalty for an query X to reference H (ACT)
#define defXToV -1       // Penalty for an query X to reference V (ACG)
#define defXToN -1       // Penalty for an query X to reference N (ACGT)
#define defXToX -1       // Penalty for an query X to reference X (ACGT)

#endif
