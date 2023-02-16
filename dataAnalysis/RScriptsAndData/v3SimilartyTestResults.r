library("ggplot2")        # graphing
library("data.table")     # apply functions to dataframe gourps
library("ggpubr")         # theme
#library("RColorBrewer")   # colors
library("viridis")       # alternative color pallete

#*******************************************************************************
# Sec-1 Sub-1:
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

################################################################################
# Name: getDataMed
# Use: get the Median of a data frame column grouped by input columns
# Input: 
#     data: data frame with data to get sum for
#     medCol: String with the column in the data frame to do the Median on
#     catAryStr: String array with the columns to group each row by for the sum
# Output:
#    data frame: with median and columns in catAryStr
#        medCol: Holds the median for each category
# Note: mainly here so takes less space in my code (not really used)
################################################################################
getDataMed = function(data, medColStr, catAryStr)
{ # getDataMed function
    tmpData = setDT(data);
    return(tmpData[, data.frame(medCol = median(get(medColStr), rm.na = TRUE)), 
                     by = catAryStr]);
} # getDataMed function


################################################################################
# Name: plotMedian
# Use: plots the median of each faucet using ggplotstats, stat_summary. The 
#      L
# Input:
#     graph: ggplot graph object to plot the median on
#         Required
#     data: Dataframe with data to graph and get median for, the x-axis column
#           name is in the graph
#     yColStr: Y axis column to get median for (string)
#     catAryStr: Catagories to get medians for (string array or string)
#     colorStr: String or color value with the color to apply to the median 
#         Default: BLACK
#     alphaDbl: Double from 0-1 with the alpha to color
#         Default: 1
#     widthDbl: Double holding the width of the bar markding the medain
#         Default: 1
#     sizeDbl: Height of crossbar (double)
#         Default: 0.5 (geom_crossbar default)
#     showLegend: Boolean to apply the median values to the legend
#         Default: FALSE
# Output: ggplot graph object with the median plotted
# Requires: ggplot2, and ggstatsplot
################################################################################
plotMedian = function(graph, 
                      data,
                      yColStr,
                      catAryStr,
                      colorStr = "BLACK",
                      alphaDbl = 1,
                      widthDbl = 1,
                      sizeDbl = 0.5,
                      showLegend=FALSE)
{ # plot median function

    # Get the medians for all points
    tmpData = getDataMed(data = data,
                         medColStr = yColStr,
                         catAryStr = catAryStr
    ); # get the medians

    # rename median column to y-axis (so ggplot does not complain)
    colnames(tmpData)[which(names(tmpData) == "medCol")] = yColStr;
    graph = graph + 
            geom_crossbar(data = tmpData,
                          aes_string(ymin = (yColStr), ymax = (yColStr)),
                          col = colorStr,
                          alpha = alphaDbl,
                          size = sizeDbl,
                          width = widthDbl, # x-axis size of cross bar
                          show.legend = showLegend # change legend?
    ); # Add the medians as crossbars

    return(graph);
} # plot median function


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: script variables
#   sec-2 sub-1: variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-2 Sub-1: variable declerations
#*******************************************************************************

# file
fileStr = "v3SimilartyTestResults"; # file to work on

# Number of consensuses per test (hardcoded)
underTwoPercInt = 60;
twoPercInt = 80;
threePercInt = 80;

# for data frames or objects
conData = NULL;                                # holds input data
graphData = NULL;                              # holds data going into a graph
graph = NULL;                                  # holds graph from ggpot

# Color variables
extraColStr = c("8DA0CB", "#E78AC3"); # two extra colors for overlaying stats
colorListData = NULL;                 # temporary variable to get color list
colorPalStr = NULL;                   # filled with color hex codes (Sec-3)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3: Read in file & set up for graphing
#   sec-3 sub-1: read in file
#   sec-3 sub-2: Set up color scheme
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-3 Sub-1: variable declerations
#*******************************************************************************

conData = setDT(read.csv(    # set dt converts to data table
              fileStr, 
              header = TRUE, # header of file for datafram names
              sep = "\t"     # deliminator for file
)); # read in my data

#*******************************************************************************
# Sec-3 Sub-2: set up color scheme
#*******************************************************************************


# build up my color palette
colorListData = viridis.map; # get the rgb colors for each palette
colorListData = colorListData[colorListData$opt == "A",]; # magma palette colors

# convert rgb to hex values
colorListData = rgb(colorListData[,1], colorListData[,2], colorListData[,3]);

# grab the colors intrested in (provides about 5 colors
colorPalStr =
   c(
       colorListData[c(185, 60, 210, 150, 115)],  # starting colors
       rgb(0,0,0),                                # black to mark extremes
       colorListData[130]                 # rest of pallete
); # create color pallete

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4: prepare data for graphing
#   sec-4 sub-1: Rename my pipelines
#   sec-4 sub-2: Get & add the percent of major variant reads to pipeline
#   sec-4 sub-3: For indels & mismatches, replace NA's with -2
#   sec-4 sub-4: Get number of misbinned reads
#   sec-4 sub-5: Get number correctly binned reads & detect maj or min variant
#   sec-4 sub-6: Get reference counts
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-4 Sub-1: Rename my pipelines
#*******************************************************************************

conData$pipeline = 
   sub(
       "allReads",
       "None",
        conData$pipeline
); # remove the allReads label

conData$pipeline = 
   sub(
       "newClustPipe",
       "V2-Cluster",
        conData$pipeline
); # remove the allReads label

conData$pipeline = 
   sub(
       "newPipe",
       "V2",
        conData$pipeline
); # remove the allReads label

conData$pipeline = 
   sub(
       "oldPipe",
       "V1",
        conData$pipeline
); # remove the allReads label

conData$pipeline = 
   sub(
       "clair",
       "Clair-",
        conData$pipeline
); # add a dash in the clair label

conData$pipeline = 
   sub(
       "longshot",
       "Longshot-",
        conData$pipeline
); # remove add a dash in the longshot label

conData$pipeline = 
   sub(
       "Ref",
       "",
        conData$pipeline
); # Remove the ref note in pipeline name

conData$pipeline = 
   sub(
       "Dist",
       "Distant",
        conData$pipeline
); # Make using distance more clear

conData$pipeline = 
   sub(
       "v3",
       "V3",
        conData$pipeline
); # Make using distance more clear

conData$Duplicate = 
   sub(
       "TRUE",
       "False Positives",
        conData$Duplicate
); # Make using distance more clear

conData$Duplicate = 
   sub(
       "FALSE",
       "True Positives",
        conData$Duplicate
); # Make using distance more clear

#*******************************************************************************
# Sec-4 Sub-2: Get & add the percent of major variant reads to pipeline
#*******************************************************************************

conData$percMajMin = sapply(
    floor(100 * conData$noMajorRef / (conData$noMajorRef + conData$noMinorRef)),
    FUN =
        (function(x)
             ifelse(x == 79, 80,
                 ifelse(x == 49, 50, x)
        )) # function to convert off roundings to percentage given badread
); # Get percentage of reads given to badread

conData$percDiff = sapply(
    conData$percId,
    FUN =
        (function(x)
             ifelse(x < 1.7, 1.4,
               ifelse(x < 2.4, 2, 3)
        )) # function to convert off roundings to percentage given badread
); # Get percentage of reads given to badread


conData$oldPipe = conData$pipeline; # so can acess later

# replace NA's for indels and mismatches with 0. This is likely what is
#  happening and makes the variant callers look better. My new pipeline and 
#  cluster have no NA's for indels and mismatches, my old pipeline only has
#  NA for indels (only one case (268 entry for the old pipeline)

conData$oldMis = conData$mismatches;

conData$mismatches = 
    tidyr::replace_na(
        data = conData$mismatches,   # only apply to mismatches
        -2                           # so stands out
); # assuming all NA's are 0 mismatches (my pipelines had no NA's)

conData$oldIndel = conData$indels;
conData$indels = 
    tidyr::replace_na(
        data = conData$indels,       # apply to indels
        -2                           # so stands out
); # assuming all NA's are 0 mismatches (my pipelines had no NA's)

#*******************************************************************************
# Sec-4 Sub-3: Change indels > 5 or mismatches > 5 to ">4"
#*******************************************************************************

# get total errors (so can use later)
conData$errorInt = conData$mismatches + conData$indels;

#*******************************************************************************
# Sec-4 Sub-4: Get number of misbinned reads
#*******************************************************************************

conData[
    ,
    misBin:=min(c(noBinMajorRef, noBinMinorRef))/(noBinMajorRef +noBinMinorRef),
    by = c(
           "majorRef",   # major variant reference name
           "minorReft",  # minor variant reference name
           "pipeline",   # pipeline used
           "acutalRef"   # reference consensus was built from
    ) # Make sure every row is done separately
]; # find the number of misbining, & if major or minor ref

#*******************************************************************************
# Sec-4 Sub-5: Get number of correctly binned reads & detect maj or min variant
#*******************************************************************************

conData[
    ,
    corBin:=max(c(noBinMajorRef, noBinMinorRef)),
    by = c(
           "majorRef",   # major variant reference name
           "minorReft",  # minor variant reference name
           "pipeline",   # pipelein used
           "acutalRef"   # reference consensus was built from
    ) # Make sure every row is done separately
]; # find the number of misbining, & if major or minor ref

conData[
    ,
    majOrMin:=ifelse(majorRef == acutalRef, "Major", "Minor"),
    by = c(
           "majorRef",   # major variant reference name
           "minorReft",  # minor variant reference name
           "pipeline",   # pipelein used
           "acutalRef"   # reference consensus was built from
    ) # Make sure every row is done separately
]; # find the number of misbining, & if major or minor ref

#*******************************************************************************
# Sec-4 Sub-5: Get number of correctly binned reads & detect maj or min variant
#*******************************************************************************

conData$pipeline = 
    factor(
        conData$pipeline,
        levels =
            c(
              "None",
              "Longshot-Minor",
              "Longshot-Major",
              "Longshot-Consensus",
              "Longshot-Distant",
              "Clair-Minor",
              "Clair-Major",
              "Clair-Consensus",
              "Clair-Distant",
              "V1",
              "V2",
              "V2-Cluster",
              "V3",
              "V3-sparse"
        ) # levels to use
); # reorder pipelines

conData$percDiff =
    gsub("1.4", "98.5% Similar", conData$percDiff);
conData$percDiff =
    gsub("2", "98% Similar", conData$percDiff);
conData$percDiff =
    gsub("3", "97% Similar", conData$percDiff);

numCon50PercData = 
    data.frame(
        pipeline = c("None", "None", "None"),
        percDiff = c("98.5% Similar", "98% Similar", "97% Similar"),
        percMajMin = c(50, 50, 50),
        V1 = c(underTwoPercInt, twoPercInt, threePercInt)
); # build datafame to hold the max number of valid consensuses

numCon80PercData = 
    data.frame(
        pipeline = c("None", "None", "None"),
        percDiff = c("98.5% Similar", "98% Similar", "97% Similar"),
        percMajMin = c(80, 80, 80),
        V1 = c(underTwoPercInt, twoPercInt, threePercInt)
); # build datafame to hold the max number of valid consensuses


conData$percDiff = 
    factor(
        conData$percDiff,
        levels =
            c(
               "98.5% Similar",
               "98% Similar",
               "97% Similar"
        ) # levels to use
); # reorder percent of reads made

conData$percMajMin = 
    factor(
        conData$percMajMin,
        levels =
            c(
              "50",
              "80"
        ) # levels to use
); # reorder pipelines

#*******************************************************************************
# Sec-4 Sub-7: Get reference counts
#*******************************************************************************

conData[
  ,
  correctCall:=
    length(
      unique(acutalRef) # removes duplicates of refs
  ), # Counts number of unique references
  by = c(
    "pipeline",
    "noMajorRef",
    "noMinorRef",
    "minorReft",
    "majorRef"
  ) # filters to sort data
]; # get number of unique references

conData[
  ,
  offCall:=
    length(memory) - 1, # counts number of references due to noise
  by = c(
    "pipeline",
    "noMajorRef",
    "noMinorRef",
    "minorReft",
    "majorRef",
    "acutalRef"
  ) # filters to sort data
]; # get number of references from nosiy reads

conData[
  ,
  minConError:=min(errorInt),
  by = c(
    "pipeline",
    "noMajorRef",
    "noMinorRef",
    "minorReft",
    "majorRef",
    "acutalRef"
  ) # filters to sort data
]; # get number of references from nosiy reads

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-7: Make graph comparing number of correct consensuses 50-50%
#   sec-7 sub-1: Make graph for comparing correct consensuses
#   sec-7 sub-2: Apply graph settings
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-7 Sub-1: Make graph for comparing correct consensuses
#*******************************************************************************

graphData =
    conData[
        conData$percMajMin == 50 &
        #conData$Duplicate == "FALSE" &   # Make sure to remove duplicates
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph

#graphData$mismatches =
#    tidyr::replace_na(
#        data = graphData$mismatches,
#        "NA"
#); # have unkown level

graphData$mismatches =
    sapply(
        graphData$oldMis,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x)))
); # keep mismatches beneath 4, that way it is easier to compare

graph =
    ggplot(
        data = graphData,
        aes_string(
            x = "pipeline",
            fill = "mismatches" # color bar by mismatches
        )
    ); # graph the mismatches

#*******************************************************************************
# Sec-7 Sub-2: Apply graph settings
#*******************************************************************************

graph =
    graph +
    stat_count(
        position = position_stack(reverse = TRUE),
    ) + # stacked bar chart
    facet_grid(
        cols = vars(percDiff), # make the columns ther perecen difference
        rows = vars(Duplicate) # split by the major-minor ratio
    ) + # split graph into facets
    scale_fill_viridis(                           # for bar chart
        discrete = TRUE,   
        direction = -1,                           # darkest colors = most errors
        option = "D",                            # default virdis color scheme
        name = "Number of False Positive SNPs"
    ) + # apply color
    geom_hline(
        data = numCon50PercData,            # has no co-infection counts
        aes(yintercept = V1)                   # total consensuses tested
    ) +                                      # draw line at predicted
    xlab("Pipeline (50% reads from major strain)") +
    ylab("Number of Consensuses") +
    theme_pubr() +          # theme I like
    theme(
        axis.text.x = element_text(angle = 90), # make x-axis labels 90 degrees
    ); # modify overall look of graph

#saveGraph("../graphs/Similarity--20000--v3--SNP-50-50");
saveGraph("graphs/Similarity--20000--v3--50-50--SNPs");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-8: Make graph comparing number of correct consensuses 80-20%
#   sec-8 sub-1: Make graph for comparing correct consensuses
#   sec-8 sub-2: Apply graph settings
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-8 Sub-1: Make graph for comparing correct consensuses
#*******************************************************************************

graphData =
    conData[
        conData$percMajMin == 80 &
        #conData$Duplicate == "FALSE" &   # Make sure to remove duplicates
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph

graphData$mismatches =
    sapply(
        graphData$oldMis,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x)))
); # keep mismatches beneath 4, that way it is easier to compare

graph =
    ggplot(
        data = graphData,
        aes_string(
            x = "pipeline",
            fill = "mismatches" # color bar by mismatches
        )
    ); # graph the mismatches

#*******************************************************************************
# Sec-8 Sub-2: Apply graph settings
#*******************************************************************************

graph =
    graph +
    stat_count(
        position = position_stack(reverse = TRUE),
    ) + # stacked bar chart
    facet_grid(
        cols = vars(percDiff), # make the columns ther perecen difference
        rows = vars(Duplicate) # split by the major-minor ratio
    ) + # split graph into facets
    scale_fill_viridis(                           # for bar chart
        discrete = TRUE,   
        direction = -1,                           # darkest colors = most errors
        option = "D",                            # default virdis color scheme
        name = "Number of False Positive SNPs"
    ) + # apply color
    geom_hline(
        data = numCon80PercData,            # has no co-infection counts
        aes(yintercept = V1)                   # total consensuses tested
    ) +                                      # draw line at predicted
    xlab("Pipeline (80% reads from major strain)") +
    ylab("Number of Consensuses") +
    theme_pubr() +          # theme I like
    theme(
        axis.text.x = element_text(angle = 90), # make x-axis labels 90 degrees
    ); # modify overall look of graph

#saveGraph("../graphs/Similarity--20000--v3--SNP-80-20");
saveGraph("graphs/Similarity--20000--v3--80-20--SNPs");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-10: Make graph comparing indels in total number of consensuses
#   sec-10 sub-1: Make graph for comparing correct consensuses with indesl
#   sec-10 sub-2: Apply graph settings
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-10 Sub-1: Make graph for comparing correct consensuses with indels
#*******************************************************************************

graphData =
    conData[
        conData$percMajMin == 50 &
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph

graphData$indels =
    sapply(
        graphData$oldIndel,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x)))
); # keep mismatches beneath 4, that way it is easier to compare

graph =
    ggplot(
        data = graphData,
        aes_string(
            x = "pipeline",
            fill = "indels",                # color bar by mismatches
        )
    ); # graph the mismatches

#*******************************************************************************
# Sec-10 Sub-2: Apply graph settings
#*******************************************************************************

graph =
    graph +
    stat_count(
        position = position_stack(reverse = TRUE),
    ) + # stacked bar chart
    facet_grid(
        cols = vars(percDiff), # make the columns ther perecen difference
        rows = vars(Duplicate) # split by the major-minor ratio
    ) + # split graph into facets
    scale_fill_viridis(                           # for bar chart
        discrete = TRUE,   
        direction = -1,                           # darkest colors = most errors
        option = "D",                            # default virdis color scheme
        name = "Number of False Positive Indels", 
    ) + # apply color
    geom_hline(
        data = numCon50PercData,            # has no co-infection counts
        aes(yintercept = V1)                   # total consensuses tested
    ) +                                      # draw line at predicted
    xlab("Pipeline (50% reads from major strain)") +
    ylab("Number of Correct Consensuses") +
    theme_pubr() +          # theme I like
    theme(
        axis.text.x = element_text(angle = 90), # make x-axis labels 90 degrees
    ); # modify overall look of graph

#saveGraph("../graphs/Similarity--20000--v3--50-50--indels");
saveGraph("graphs/Similarity--20000--v3--50-50--indels");

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-11: Make graph comparing indels in total number of consensuses
#   sec-11 sub-1: Make graph for comparing correct consensuses with indesl
#   sec-11 sub-2: Apply graph settings
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*******************************************************************************
# Sec-11 Sub-1: Make graph for comparing correct consensuses with indels
#*******************************************************************************

graphData =
    conData[
        conData$percMajMin == 80 &
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph

graphData$indels =
    sapply(
        graphData$oldIndel,
        FUN = function(x)  # make sure strings are sorted numericly
            ifelse(x > 14, "z >14",
            ifelse(x > 9, paste("A", x), as.character(x)))
); # keep mismatches beneath 4, that way it is easier to compare

graph =
    ggplot(
        data = graphData,
        aes_string(
            x = "pipeline",
            fill = "indels",                # color bar by mismatches
        )
    ); # graph the mismatches

#*******************************************************************************
# Sec-11 Sub-2: Apply graph settings
#*******************************************************************************

graph =
    graph +
    stat_count(
        position = position_stack(reverse = TRUE),
    ) + # stacked bar chart
    facet_grid(
        cols = vars(percDiff), # make the columns ther perecen difference
        rows = vars(Duplicate) # split by the major-minor ratio
    ) + # split graph into facets
    scale_fill_viridis(                           # for bar chart
        discrete = TRUE,   
        direction = -1,                           # darkest colors = most errors
        option = "D",                            # default virdis color scheme
        name = "Number of False Positive Indels", 
    ) + # apply color
    geom_hline(
        data = numCon80PercData,            # has no co-infection counts
        aes(yintercept = V1)                   # total consensuses tested
    ) +                                      # draw line at predicted
    xlab("Pipeline (80% reads from major strain)") +
    ylab("Number of Correct Consensuses") +
    theme_pubr() +          # theme I like
    theme(
        axis.text.x = element_text(angle = 90), # make x-axis labels 90 degrees
    ); # modify overall look of graph

#saveGraph("../graphs/Similarity--20000--v3--80-20--indels");
saveGraph("graphs/Similarity--20000--v3--80-20--indels");

#***********************************************************************
# TIME
#***********************************************************************

graphData =
    conData[
        conData$percMajMin == 50 &
        conData$pipeline != "None" &
        conData$pipeline != "V1" &
        conData$pipeline != "V2-Cluster" &
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph


graph =
  ggplot(
    data = graphData,
    aes(
      x = pipeline,                # number of reads made (roughly)
      y = time / 60,               # graph time in minutes
      color = pipeline,
      shape = pipeline
    ) # aes
); # graph data

#*******************************************************************************
# Sec-5 Sub-2: apply settings to time graph
#*******************************************************************************

graph =
  graph +
  geom_point(
    size = 5,
    alpha = 0.5,                             # so points can overlap
    position =
      position_jitter(
        width = 0.3,
        height = 0
      ) # apply jitters
  ) + # make point chart
  #facet_grid(cols = vars(pipeline)) +
  theme_pubr() +                             # nice clean theme
  scale_shape_discrete() +                    # shape of points
  scale_color_viridis(                          # for point graph
    alpha = 0.5,
    discrete = TRUE,
    option = "D"
  ) + # apply color
  ylab("Time in minutes") +
  xlab("Pipeline with 50% of reads from major strain") +
  ylim(0, NA) +                        # y limits (NA for defalut)
  theme(
    axis.text.x = element_text(angle = 90), # rotate x-axis labels
    #legend.position = "none"                # remove the legend
); # theme, rotate text, remove legend

#saveGraph("../graphs/time");
saveGraph("graphs/Similarity--time--v3");


#*******************************************************************************
# Sec-6 Sub-1: Graph memory data
#*******************************************************************************

graphData = conData;
graphData =
    conData[
        conData$percMajMin == 50 &
        conData$pipeline != "None" &
        conData$pipeline != "V1" &
        conData$pipeline != "Longshot-Major" &
        conData$pipeline != "Longshot-Minor" &
        conData$pipeline != "Longshot-Distant" &
        conData$pipeline != "Clair-Major" &
        conData$pipeline != "Clair-Minor" &
        conData$pipeline != "Clair-Distant"
        ,
]; # remove the extact reference matches from my data to graph


graph =
  ggplot(
    data = graphData,
    aes(
      x = pipeline,                # number of reads made (roughly)
      y = memory / 1000,           # max resident memory in mb
      color = pipeline,
      shape = pipeline
    ) # aes
); # graph data

#*******************************************************************************
# Sec-5 Sub-2: apply settings to time graph
#*******************************************************************************

graph =
  graph +
  geom_point(
    size = 5,
    alpha = 0.5,                             # so points can overlap
    position =
      position_jitter(
        width = 0.3,
        height = 0
      ) # apply jitters
  ) + # make point chart
  #facet_grid(cols = vars(pipeline)) +
  theme_pubr() +                             # nice clean theme
  scale_shape_discrete() +                    # shape of points
  scale_color_viridis(                          # for point graph
    alpha = 0.5,
    discrete = TRUE,
    option = "D"
  ) + # apply color
  ylab("Max resident memory usage in mb") +
  xlab("Pipeline with 50% of reads from major strain") +
  ylim(0, NA) +                        # y limits (NA for defalut)
  theme(
    axis.text.x = element_text(angle = 90), # rotate x-axis labels
    #legend.position = "none"                # remove the legend
); # theme, rotate text, remove legend

#saveGraph("../graphs/memory");
saveGraph("graphs/Similarity--memory--v3");
