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
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef DEFAULTSETTINGS_H
#define DEFAULTSETTINGS_H

/**********************************************************************\
* Sec-1: variable sizes
\**********************************************************************/

#define bitsInFlt sizeof(float) << 3 /*Number bits in a float*/
#define bitsInSht sizeof(short) << 3 /*Number of bits in a short*/

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
#define defVersion 3.20230313  /*The Version number of this program*/
    /*Format is version.yearMonthDay*/
#define defPrefix "out"      /*Default prefix to use*/
#define defThreads "3"       /*Default number of threads to use*/
#define rmReadsWithSupAln 0
    /*1: remove reads with supplementary aligments*/
#define defReadsPerCon 300
    /*Maximum number of reads to use to build a consensus*/
#define minReadsPerBin 100 
    /*Minimum number reads to keep a bin or build a consensus*/
#define defMinPercReads 0.003 /*(0.3% of all clustered reads)*/
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

#endif
