library("ggplot2")
library("data.table")
library("ggpubr")
library("viridis")
# uses tidyr::replace_na

#*******************************************************************************
# Sec-0 Sub-1:
#  Output: till: Saves output graph as a tiff
#*******************************************************************************
saveGraph = function(
  nameStr           # name of graph to save
) # Use: saves a graph using ggsave function
{ # saveGraph
    ggsave(paste(nameStr, ".svg", sep = ""), 
           device = "svg", # save as tiff file
           dpi = 300,
    ); # ggsave (save the graph)
} # saveGraph


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Read in my data
statsDT =
    setDT(         # data frames are a bit easier to manipulate
        read.csv(
            "ASHURE-full-20230314-trim--duplicate-marked.tsv",
            sep = "\t",
            header = TRUE
));
graphDT = NULL; # for graphing

# variables to hold each graph
graphPerc = NULL;
graphSnp = NULL;
graphIndel = NULL;
graphDiff = NULL;
graphTime = NULL;
graphRam = NULL;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Set up names
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Replace unkown aligned/reference lengths with 0
statsDT$AlignedReadLength =
    tidyr::replace_na(data = statsDT$AlignedReadLength, 0);

statsDT$referenceLength =
    tidyr::replace_na(data = statsDT$referenceLength, 0);

statsDT$readLength = tidyr::replace_na(data = statsDT$readLength, -1);

# Reduce the name of medakas model (high & fast combined with version)
statsDT$model = gsub("r941_min_", "", statsDT$MedakaModel);
statsDT$model = gsub("_g", "-", statsDT$model);

# Get the file used (A, B, AB)
statsDT$dataset = gsub(".*--", "", statsDT$fqFile);
statsDT$dataset = gsub(".fastq", "", statsDT$dataset);

# Get if primer triming was done
statsDT$primers = gsub(".*primers.fasta", "Primer", statsDT$primers);
statsDT$primers = gsub(".*Ref", "Ref", statsDT$primers);

# Set up the consensus building names
statsDT$test =
    sapply(
        statsDT$MajorityConsensus,
        FUN = (function(x) ifelse(x == TRUE, "maj", "NA"))
); # Detect if the majority consensus was used

statsDT$test =
    paste(
       statsDT$test, 
       sapply(
           statsDT$RaconUsed,
           FUN = (function(x) ifelse(x == TRUE, "rac", "NA"))
       ), # Add the racon name
       sep = "-"
); # Detect if racon was used

# Recored the consensus buildig method used with medaka
statsDT$noMedCol = statsDT$test;

statsDT$test =
    paste(
       statsDT$test, 
       sapply(
           statsDT$MedakaUsed,
           FUN = (function(x) ifelse(x == TRUE, "med", "NA"))
       ), # Add the racon name
       sep = "-"
); # Detect if racon was used

statsDT$test = paste(statsDT$test, statsDT$primer, sep = "-");
# Set up the name for each test
#statsDT$test = paste(statsDT$test, statsDT$model, sep = "-");
statsDT$percID =
    100 * (1 - ((statsDT$SNPs + statsDT$Insertions + statsDT$Deletions)
                 / statsDT$AlignedReadLength
               )
);

statsDT$qScore = 10 * log10(statsDT$percID);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: Graph SNPs
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT = statsDT;

# make duplicate colum more easy to understand
graphDT$Duplicate = gsub("TRUE", "False +", graphDT$Duplicate);
graphDT$Duplicate = gsub("FALSE", "TRUE +", graphDT$Duplicate);

graphDT$test =
    factor(
        graphDT$test,
        levels = 
            c(
                "NA-rac-NA-Ref",
                "NA-rac-NA-Primer",
                "maj-rac-NA-Ref",
                "maj-rac-NA-Primer",
                "NA-rac-med-Ref",
                "NA-rac-med-Primer",
                "maj-rac-med-Ref",
                "maj-rac-med-Primer",
                "NA-NA-med-Ref",
                "NA-NA-med-Primer",
                "maj-NA-med-Ref",
                "maj-NA-med-Primer",
                "maj-NA-NA-Ref",
                "maj-NA-NA-Primer"
        ) # catagories
);

graphDT$dataset = factor(graphDT$dataset, levels = c("A", "B", "AB"));


graphDT$charSNP =
    sapply(
        graphDT$SNPs,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(is.na(x), "zz NA",
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x))))
); # keep mismatches beneath 4, that way it is easier to compare

graphSnp = 
    ggplot(
        data = graphDT,
        aes(
            x = test, # consensus building method
            fill = charSNP
        ) # make the graph
);

graphSnp = 
   graphSnp +
   stat_count(position = position_stack(reverse = TRUE)) +
   facet_grid(
       cols = vars(dataset),
       rows = vars(Duplicate)
   ) +
   scale_fill_viridis(
       discrete = TRUE,
       direction = -1,
       option = "D",
       name = "False SNPs"
   ) + 
   geom_hline(aes(yintercept = DectableReferences), linetype="dashed") +
   geom_hline(aes(yintercept = TotalReferences)) +
   xlab("Consensus bulding method") +
   ylab("Number of consensuses") +
   theme_pubr() +
   theme(axis.text.x = element_text(angle = 90));

saveGraph("ASHURE-bench--snps");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Graph indels
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT$charIndel =
    sapply(
        graphDT$Insertions + graphDT$Deletions,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(is.na(x), "zz NA",
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x))))
); # keep mismatches beneath 4, that way it is easier to compare

graphIndel = 
    ggplot(
        data = graphDT,
        aes(
            x = test, # consensus building method
            fill = charIndel
        ) # make the graph
);

graphIndel = 
   graphIndel +
   stat_count(position = position_stack(reverse = TRUE)) +
   facet_grid(
       cols = vars(dataset),
       rows = vars(Duplicate)
   ) +
   scale_fill_viridis(
       discrete = TRUE,
       direction = -1,
       option = "D",
       name = "False Indels"
   ) + 
   geom_hline(aes(yintercept = DectableReferences), linetype="dashed") +
   geom_hline(aes(yintercept = TotalReferences)) +
   xlab("Consensus bulding method") +
   ylab("Number of consensuses") +
   theme_pubr() +
   theme(axis.text.x = element_text(angle = 90));

saveGraph("ASHURE-bench--indels");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-6: Graph length diff
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT$lenDiffI = graphDT$AlignedReadLength - graphDT$referenceLength;

graphDiff =
    ggplot(
        data = graphDT[graphDT$Duplicate != "False +",],
        #data = graphDT,
        aes(
            x = test,
            y = lenDiffI,
            col = test
        )
);

graphDiff =
    graphDiff +
    geom_point(
        size = 3,
        alpha = 0.5,
        position = position_jitter(height = 0, width = 0.3)
    ) +
   facet_grid(
       cols = vars(dataset),
       #rows = vars(Duplicate)
   ) +
   scale_color_viridis(
        discrete = TRUE,
        option = "D",    # Default pallete
        name = "False positive"
    ) +
   xlab("Consensus bulding method (False +'s removed)") +
   ylab("Aligned consensus length - reference length") +
   theme_pubr() +
   theme(
       axis.text.x = element_text(angle = 90), # rotate x-axis labels
       legend.position = "none"                # remove the legend
   );

saveGraph("ASHURE-bench--lenDiff");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Graph time
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT$elapsedTime = graphDT$elapsedTime / 60; # convert seconds to min

graphTime =
    ggplot(
        data = graphDT,
        aes(
            x = test,
            y = elapsedTime,
            col = test
        )
);

graphTime =
    graphTime +
    geom_point(
        size = 3,
        alpha = 0.5,
        position = position_jitter(height = 0, width = 0.3)
    ) +
   facet_grid(
       cols = vars(dataset)#,
       #rows = vars(Duplicate)
   ) +
   scale_color_viridis(
        discrete = TRUE,
        option = "D",    # Default pallete
        name = "Pipeline"
    ) +
   xlab("Consensus bulding method") +
   ylab("Time in minutes") +
   theme_pubr() +
   theme(
       axis.text.x = element_text(angle = 90), # rotate x-axis labels
       legend.position = "none"                # remove the legend
   );

saveGraph("ASHURE-bench--time");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-8: Graph memory
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Conver memory usage to megabytes (in kilabytes)
graphDT$maxResidentMemory = graphDT$maxResidentMemory / 1000;

graphMem =
    ggplot(
        data = graphDT,
        aes(
            x = test,
            y = maxResidentMemory,
            col = test
        )
);

graphMem =
    graphMem +
    geom_point(
        size = 3,
        alpha = 0.5,
        position = position_jitter(height = 0, width = 0.3)
    ) +
   facet_grid(
       cols = vars(dataset),
       rows = vars(Duplicate)
   ) +
   scale_color_viridis(
        discrete = TRUE,
        option = "D",    # Default pallete
        name = "Pipeline"
    ) +
   xlab("Consensus bulding method") +
   ylab("Maximum resident memory used in megabytes") +
   theme_pubr() +
   theme(
       axis.text.x = element_text(angle = 90), # rotate x-axis labels
       legend.position = "none"                # remove the legend
   );

saveGraph("ASHURE-bench--mem");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-9: Graph extra length
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT$lenAlnDiffI = graphDT$readLength - graphDT$AlignedReadLength;

graphExtraLen =
    ggplot(
        data = graphDT[graphDT$Duplicate != "False +",],
        #data = graphDT,
        aes(
            x = test,
            y = lenAlnDiffI,
            col = test
        )
);

graphExtraLen =
    graphExtraLen +
    geom_point(
        size = 3,
        alpha = 0.5,
        position = position_jitter(height = 0, width = 0.3)
    ) +
   facet_grid(
       cols = vars(dataset),
       #rows = vars(Duplicate),
       scale = "free_y"
   ) +
   scale_color_viridis(
        discrete = TRUE,
        option = "D",    # Default pallete
        name = "False positive"
    ) +
   xlab("Consensus bulding method (False +'s removed)") +
   ylab("Consensus length - Aligned consensus length") +
   theme_pubr() +
   theme(
       axis.text.x = element_text(angle = 90), # rotate x-axis labels
       legend.position = "none"                # remove the legend
   );

saveGraph("ASHURE-bench--extraConLen");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-10: Graph % reads per cluster
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

graphDT[
    ,
    totalReads:=sum(ReadsInCluster),
    by = c("fqFile", "MajorityConsensus", "RaconUsed", "MedakaUsed")
];
graphPercReads =
    ggplot(
        data = graphDT,
        aes(
            x = test,
            y = 100 * (ReadsInCluster / totalReads),
            col = test
        )
);

graphPercReads =
    graphPercReads +
    geom_point(
        size = 3,
        alpha = 0.5,
        position = position_jitter(height = 0, width = 0.3)
    ) +
   geom_hline(aes(yintercept = 1)) + # mark 1% line
   facet_grid(
       cols = vars(dataset),
       rows = vars(Duplicate)
   ) +
   scale_color_viridis(
        discrete = TRUE,
        option = "D",    # Default pallete
        name = "False positive"
    ) +
   xlab("Consensus bulding method") +
   ylab("Percentage of total reads clustered") +
   theme_pubr() +
   theme(
       axis.text.x = element_text(angle = 90), # rotate x-axis labels
       legend.position = "none"                # remove the legend
   );

saveGraph("ASHURE-bench--percReads");


# dataset          
# referenceLength
# readLength
# alignedReadLength
# SNPs
# Insertions
# Deletions
# ReadsInCluster
# DetectedReferences
# DetectableReferences
# TotalReferences
# elapsedTime
# maxResidentMemory
