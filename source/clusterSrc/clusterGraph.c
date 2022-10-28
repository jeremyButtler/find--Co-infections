/*##############################################################################
# Name: clustGraphGraphFun.c
# Use: Holds the functions for building a graphNode graph
##############################################################################*/

#include "clusterGraph.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clustGraphGraphFun
#    fun-1: insertGraphNode: Inserts a read into a read graph
#    fun-2: parseLine: convert input line to read & ref c-string & get mapq
#    fun-3: assignReadToBridge: assign read to bridge (if needed merge bridges)
#    fun-4: mergeBridges: merges two birdges into one cluster
#    fun-5: printAndFreeGraph: prints & frees all sequencens in a graph
#    fun-6: breakCluster: Breaks nodes in cluster into individual clusters
#    fun-7: breakClusterAndFreeHead: frees head an breaks other cluster nodes 
#                                    into individual clusters
#    fun-8: freeCluster: frees all nodes in a cluster
#    fun-9: printAndFreeCluster: Free and prings all nodes in a cluster
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*##############################################################################
# Output:
#    Returns: the readInfo for the current read on
#             0: If malloc erroed out at some point
#    Modefies: mappedReadsULng to have the read counts for the current read
#              This is reset when a new read starts
#    Modefies: pushes graphNode on clustStack if read maps to new cluster
##############################################################################*/
struct readInfo * insertGraphNode(
    char *lineInCStr,                  /*fromat: readName\treferenceName\tmapq*/
    struct readInfo **readTree, /*readInfo tree to search for read & ref names*/
    struct graphNode **graphRoot, /*Root node of read graph with clusters*/
    struct readInfo *oldReadNode, /*readInfo node from last aligment looked at*/
    struct minValues *minStats,  /*min thresholds to keep alignment or cluster*/
    unsigned long *mappedReadsULng, /*Number of references mapped to this read*/
    struct graphNodeStack **clustStack /*Clusters this read shares edges with*/
) /*Inserts a read into a read graph*/
{ /*insertGraphNode*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # insertGraphNode TOC:
    #    fun-1 sec-1: variable declerations
    #    fun-1 sec-2: Extract & check read name, reference name, & mapq values
    #    fun-1 sec-3: Find the read & determine if new read (set of aligments)
    #    fun-1 sec-4: Check if moveing to new set of aligments (if so prepare)
    #    fun-1 sec-5: Find reference & do final checks for keeping alignment
    #    fun-1 sec-6: Check if the read is a new to graph (if so add)
    #    fun-1 sec-7: Check if the reference can be reasigned to read cluster
    #    fun-1 sec-8: Add reference cluster to list of mapped clusters
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char readCStr[1000],   /*Iterate though the line*/
         refCStr[1000];    /*Start of reference or read name*/

    unsigned char lenReadUChar=0, /*Length of read name*/
                  lenRefUChar=0;  /*Lenght of the reference name*/

    double mapqDbl = 0;           /*Holds mapping quality of alignemnt*/

    struct readInfo *readInfoNode = 0, /*points to read node in readInfo tree*/
                    *refInfoNode = 0;  /*points to reference node in tree */

    struct graphNode *refGraphNode = 0,
                     *readGraphNode = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Extract & check read name, reference name, & mapq values
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    mapqDbl = parseLine(lineInCStr,
                        readCStr,
                        &lenReadUChar,
                        refCStr,
                        &lenRefUChar
    ); /*Get read name, reference name, & mapq from the input line*/

    if(*refCStr == '*' || *readCStr == '*')
        return oldReadNode;       /*No reference/read in aligment (no mapping)*/

    if(mapqDbl < minStats->minMapqDbl)
        return oldReadNode; /*alignment is under min quality to keep*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Find the read & determine if is new read (set of aligments)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(oldReadNode == 0)                 /*1st time*/
    { /*If is teh first time*/
        readInfoNode = findAddNodeToReadTree(readCStr, lenReadUChar, readTree);

        if(readInfoNode == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:92\n");
            return 0;         /*Malloc failed to assign memory for new node*/
        } /*If malloc failed to allocate memory*/

        *mappedReadsULng = 0;
        if(pushReadIntoGraph(graphRoot, readInfoNode) == 0) /*new read*/
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:101\n");
            return 0;         /*Malloc did not assign memory to readGraphNode*/
        } /*If malloc failed to allocate memory*/
    } /*If is the first time*/

    else if(strcmp(oldReadNode->idCStr, readCStr) == 0)
        readInfoNode = oldReadNode; /*have not moved to a new set of aligments*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Check if moveing to new set of aligments (if so prepare)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    else
    { /*else starting on a set of alignments for a new read*/

        readGraphNode = oldReadNode->nodeInGraph;        /*For clarity*/

        if(*mappedReadsULng < minStats->minReadsULng ||
           *mappedReadsULng < minStats->minSharedEdgesUInt 
          )
        { /*If read is not worth keeping*/
            if(breakCluster(&readGraphNode, graphRoot) == 0)
            { /*If malloc failed to allocate memory*/
                fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:125\n");
                 return 0;
            } /*If malloc failed to allocate memory*/
         
            (*oldReadNode).doneChar = 4;  /*Ignoring read for futre*/

            while(*clustStack != 0)
            { /*Resest and free the cluster stacks*/
                (*clustStack)->graphStruct->readCntULng = 0;
                popGraphStackStruct(clustStack);
            } /*Resest and free the cluster stacks*/
        } /*If read is not worth keeping*/

        else                            /*keeping the old read*/
        { /*else need to assign the read to cluster*/
            if(assignReadToCluster(readGraphNode,
                                   clustStack,
                                   minStats,
                                   graphRoot
                                  ) == 0
              )
              { /*If malloc failed to allocate memory*/
                fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:136\n");
                 return 0;
              } /*If malloc failed to allocate memory*/
        } /*else need to assign the read to cluster*/

        /*Move to the new read (set of alignments)*/
        *mappedReadsULng = 0;
        readInfoNode = findAddNodeToReadTree(readCStr, lenReadUChar, readTree);

        if(readInfoNode == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:148\n");
            return 0;         /*Malloc failed to assign memory for new node*/
        } /*If malloc failed to allocate memory*/


        if(readInfoNode->nodeInGraph == 0)
        { /*If need to push the read in the graph*/
            if(pushReadIntoGraph(graphRoot, readInfoNode) == 0) /*Is new read*/
            { /*If malloc failed to allocate memory*/
                fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:157\n");
                return 0;
            } /*If malloc failed to allocate memory*/
        } /*If need to push the read in the graph*/

        if(mvToNewClust(readInfoNode->nodeInGraph, graphRoot) == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:166\n");
            return 0;
        } /*If malloc failed to allocate memory*/
    } /*else starting on a set of alignments for a new read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-5: Find reference & do final checks for keeping alignment
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*readInfoNode).doneChar > 0)
        return readInfoNode;      /*Already looked at alingments for this read*/

    refInfoNode = findAddNodeToReadTree(refCStr, lenRefUChar, readTree);

    if(refInfoNode == 0)
    { /*If malloc failed to allocate memory*/
        fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:180\n");
        return 0;
    } /*If malloc failed to allocate memory*/

    if(refInfoNode->doneChar == 4)
        return readInfoNode;     /*Already marked this ref for deletion*/

     (*mappedReadsULng)++;      /*At this point keeping aligment*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-6: Check if the reference is a new to graph (if so add)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    readGraphNode = readInfoNode->nodeInGraph;   /*clarity in next lines*/

    if((*refInfoNode).nodeInGraph == 0)
    { /*If refence is a new node, add to the read node*/
        if(pushReadIntoGraph(&readGraphNode, refInfoNode) == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:201\n");
            return 0;
        } /*If malloc failed to allocate memory*/

        readGraphNode->clustNode->readCntULng++;   /*Is an new edge*/

        return readInfoNode;
    } /*If refence is a new node add to the read node*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-7: Check if the reference can be reasigned to read cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    refGraphNode = refInfoNode->nodeInGraph;

    if(refGraphNode->clustNode->numNodesULng == 1)
    { /*If the reference is the only node in the cluster (not assigned yet)*/
        readGraphNode->clustNode->readCntULng++;   /*Is an new edge*/

        if(mvToEdge(refGraphNode, readGraphNode, graphRoot) == 0) /*Move ref*/
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:221\n");
            return 0;
        } /*If malloc failed to allocate memory*/

        return readInfoNode;
    } /*If the reference is the only node in the cluster (not assigned yet)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-8: Add reference cluster to list of mapped clusters
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(refGraphNode->clustNode->readCntULng == 0)
    { /*If need to add the references cluster to my list of mapped clusters*/
        if(pushGraphNodeStack(refGraphNode->clustNode, clustStack) == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-1: insertGraphNode: clusterGraph.c:221\n");
            return 0;
        } /*If malloc failed to allocate memory*/
    } /*If need to add the references cluster to my list of mapped clusters*/

    refGraphNode->clustNode->readCntULng++; /*Mark one more mapping to cluster*/
 
    return readInfoNode;
} /*insertGraphNode*/

/*##############################################################################
# Output:
#    Returns: mapq as double
#    Modifies: readCStr/refCStr to hold a copy of read/reference name
#    Modifies: readCStr to hold a copy of read name in lineInCStr
#    Modifies: lenReadUChar & lenRefUChar to hold lengths of read/ref name
##############################################################################*/
double parseLine(
    char *lineInCStr,            /*Line with read name, reference name, & mapq*/
    char *readCStr,              /*Array to hold read name*/
    unsigned char *lenReadUChar, /*Holds length of read name*/
    char *refCStr,              /*Array to hold reference name*/
    unsigned char *lenRefUChar   /*Holds length of reference name*/
) /*Copies read and reference names in input, also returns mapq as double*/
{ /*parseLine*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC: parseLine
    #    fun-2 sec-1: variable declerations
    #    fun-2 sec-2: Copy the read name to readCStr
    #    fun-2 sec-3: Copy the reference name to refCStr
    #    fun-2 sec-4: Convert mapq to C-string, then to double
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: variable declerations
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    char *mapqCStr = 0,
         *tmpCStr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Copy the read name to the readCStr pointer
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    tmpCStr = readCStr;

    while(*lineInCStr != '\t')
    { /*Loop till at end of the read name*/
        *tmpCStr = *lineInCStr;
        lineInCStr++;
        tmpCStr++;
        (*lenReadUChar)++;
    } /*Loop till at end of the read name*/

    *tmpCStr = '\0';    /*make read name seciton into a separate c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Copy the reference name to refCStr
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    lineInCStr++;          /*Move to start of reference name*/
    tmpCStr = refCStr;

    while(*lineInCStr != '\t')
    { /*Loop till at end of the reference name*/
        *tmpCStr = *lineInCStr;
        lineInCStr++;
        tmpCStr++;
        (*lenRefUChar)++;
    } /*Loop till at end of the reference name*/

    *tmpCStr = '\0';    /*make read name seciton into a separate c-string*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-4: Convert mapq to C-string, then to double
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    lineInCStr++;                 /*Move to start of mapq*/
    mapqCStr = lineInCStr;        /*point to start of the mapq entry*/

    while(*lineInCStr != '\n')
        lineInCStr++;

    *lineInCStr = '\0';          /*Turn mapq into c-string*/

    if(*mapqCStr == '*')
        return 0;      /*No mapq value assigned*/

    return strtod(mapqCStr, &tmpCStr); /*convert c-string to double*/
} /*parseLine*/

/*##############################################################################
# Output:
#    Modifies: readToAssigne to be in the final assigned cluster
#    Modifies: If root node cluster is merged, changes root node to new root
#    Frees: clustStack
#    Returns: 1: if succeeded
#             0: if Malloc error at some point
##############################################################################*/
char assignReadToCluster(
    struct graphNode *readToAssign,       /*Read to assign to a cluster*/
    struct graphNodeStack **clustStack,   /*Stack of clusters readToAssign
                                            shares edges with*/
    struct minValues *minStats,  /*min thresholds to keep alignment or cluster*/
    struct graphNode **graphRoot          /*Root node of read graph*/
) /*Assignes a read to a cluster, if mulitple assignments, merges clusters*/
{ /*assignReadToCluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1 TOC: assign read to bridge (if needed merge bridges)
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *clustOn = 0;

    if(readToAssign->readCntULng < minStats->minSharedEdgesUInt)
    { /*If ther are not enough edges to be a new cluster, remove*/
        if(readToAssign->edgesList != 0)
        { /*If need to break the cluster into indivdual clusters*/
          if(breakCluster(&readToAssign, graphRoot) == 0)
          { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-3: assignReadToCluster: clusterGraph.c:360\n");
            return 0;
          } /*If malloc failed to allocate memory*/
        } /*If need to break the cluster into indivdual clusters*/

        readToAssign->numNodesULng = 1;                 /*So always merged*/
    } /*If ther are not enough edges to be a new cluster, remove*/

    else
        readToAssign->readCntULng = 0;           /*So is reset if cluster kept*/

    readToAssign->seqIdNode->doneChar = 1;      /*record looking at this read*/

    while(*clustStack != 0)
    { /*Loop though the stack of clusters*/
        clustOn = (*clustStack)->graphStruct;
        popGraphStackStruct(clustStack);   /*Move to next cluster*/

        if(clustOn->readCntULng > (*minStats).minSharedEdgesUInt)
        { /*If need to merge clusters*/
           clustOn->readCntULng = 0; /*Reset for next set of read alignments*/

           if(mergeClusters(&clustOn, &readToAssign, graphRoot) == 0)
           { /*If malloc failed to allocate memory*/
             fprintf(stderr,"fun-3: assignReadToCluster: clusterGraph.c:386\n");
             return 0;
           } /*If malloc failed to allocate memory*/
        } /*If need to merge clusters*/

        else
            clustOn->readCntULng = 0; /*Reset for next set of read alignments*/
    } /*Loop though the stack of clusters*/

    if(readToAssign->clustNode->numNodesULng < (*minStats).minSharedEdgesUInt)
        readToAssign->seqIdNode->doneChar = 4; /*Mark as ignore*/

    return 1;
} /*assignReadToCluster*/

/*##############################################################################
# Output:
#    Modifies: Merges cluster with fewest reads into the remaining cluster
#    Returns: 1: if sucess
#             0: If malloc failed to assign memory at some point
##############################################################################*/
char mergeClusters(
    struct graphNode **firstClust,   /*First cluster to merge*/
    struct graphNode **secondClust,  /*Second cluster to merge*/
    struct graphNode **graphRoot     /*Root node of the graph. If root node
                                       merged, then reassigns root node*/
) /*Merge two bridges (clusters) into one bridge (cluster)*/
{ /*mergeClusters*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 TOC: mergeClusters
    #    fun-4 sec-1: variable declerations
    #    fun-4 sec-2: Check which cluster to keep
    #    fun-4 sec-3: Reset cluster head nodes for all nodes in cluster to merge
    #    fun-4 sec-4: Complete merge by adding merge cluster kept clusters edges
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-1: Variable decleratins
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode **keptClust = 0,   /*cluster to keep*/
                     **mergeClust = 0,  /*cluster to merge*/
                     *nextNode = 0;

    struct graphNodeStack *graphStack = 0; /*Holds previously visted nodes*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Check which cluster to keep
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*
      This ensures that the smallest bridge set is always removed, which helps
        avoid the worst case of always freeing the larger bridge.
      This would mean transversing the larger merge repeatedly, wich resuts
        in more than n transversals.
    */
    
    if((*firstClust)->clustNode->numNodesULng >=
       (*secondClust)->clustNode->numNodesULng
      )
    { /*If keeping the first bridge*/
        keptClust = firstClust;
        mergeClust = secondClust;
    } /*If keeping the first bridge*/

    else
    { /*else keeping the second bridge*/
        keptClust = secondClust;
        mergeClust = firstClust;
    } /*If keeping the second bridge*/

    *mergeClust = (*mergeClust)->clustNode; /*Make sure head of cluster freeing*/
    *keptClust = (*keptClust)->clustNode;  /*Make sure working at cluster head*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-3: Reset cluster head nodes for all nodes in cluster to merge
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Create my que*/
    if((*mergeClust)->edgesList != 0)
    { /*If need to start a stack (more than one node in cluster*/
        if(pushGraphNodeStack((*mergeClust)->edgesList, &graphStack) == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(stderr, "fun-4: mergeClusters: clusterGraph.c:469\n");
            return 0;
        } /*If malloc failed to allocate memory*/
    } /*If need to start a stack (more than one node in cluster*/

    (*mergeClust)->clustNode = *keptClust; /*Make sure head of cluster freeing*/
    (*keptClust)->numNodesULng += (*mergeClust)->numNodesULng;
    (*mergeClust)->numNodesULng = 0;       /*Make sure head of cluster freeing*/

    while(graphStack != 0)
    { /*Loop through tree till the stack is empty*/
        nextNode = graphStack->graphStruct;     /*Next node to work on*/
        popGraphStackStruct(&graphStack);       /*Move to next node*/

        while(nextNode != 0)
        { /*Loop to the last neighbor node and build up the que*/

             if(nextNode->clustNode != *mergeClust)
                 continue;           /*Not a cluster I am replacing, so ignore*/

             nextNode->clustNode = *keptClust; /*Reset the bridge*/

             if(nextNode->edgesList != 0)
             { /*If have edges to add to stack*/
                if(pushGraphNodeStack( nextNode->edgesList, &graphStack) == 0)
                { /*If malloc failed to allocate memory*/
                    fprintf(
                           stderr,
                           "fun-4: mergeClusters: clusterGraph.c:494\n"
                    ); /*Print out error message*/
                    return 0;
                } /*If malloc failed to allocate memory*/
             } /*If have edges to add to stack*/

             nextNode = nextNode->nextNode;
        } /*Loop to the last neighbor node and build up the que*/
    } /*Loop through tree till the stack is empty*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-4: Complete merge by adding merge cluster kept clusters edges
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*mergeClust == *graphRoot)
    { /*If the nodeing merging was the root of graph*/
        *graphRoot = (*mergeClust)->nextNode;
        (*graphRoot)->prevNode = (*mergeClust)->prevNode;     /*Should be 0*/
    } /*If the nodeing merging was the root of graph*/

    else
    { /*Else need to cut merge clust out*/
        if((*mergeClust)->prevNode != 0)
            (*mergeClust)->prevNode->nextNode = (*mergeClust)->nextNode;

        if((*mergeClust)->nextNode != 0)
            (*mergeClust)->nextNode->prevNode = (*mergeClust)->prevNode;
    } /*Else need to cut merge clust out*/

    (*mergeClust)->prevNode = *keptClust;
    (*mergeClust)->nextNode = (*keptClust)->edgesList;

    if((*keptClust)->edgesList != 0)
        (*keptClust)->edgesList->prevNode = *mergeClust;

    (*keptClust)->edgesList = *mergeClust;

    return 1;
} /*mergeBridges*/

/*##############################################################################
# Output:
#    stdout: prints read name & assigned cluster to stdout (read\tcluster)
#    Returns: 1 if succeded
#             0 if malloc failed to allocate memory
##############################################################################*/
char printAndFreeGraph(
    struct graphNode **graphRoot, /*Root of read graph to print out*/
    struct minValues *minStats,   /*min thresholds to keep alignment or cluster*/
    char printGraphChar          /*1: print each cluster, when meets requirments
                                   0: Do not print*/
) /*prints out the bridges (clusters) and frees a graphNode graph*/
{ /*printAndFreeGraph*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-1 TOC: print read name & cluster number for every read in graph
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long clustOnULng = 0;

    while(*graphRoot != 0)
    { /*While there are clusters to print out*/
        if(printGraphChar == 0 ||
           (*graphRoot)->numNodesULng < minStats->minReadsPerClustULng
          )
        { /*If there are few reads in cluster then required (do not print)*/
            if(freeCluster(graphRoot, graphRoot, clustOnULng) == 0)
            { /*If malloc failed to allocate memory*/
              fprintf(stderr, "fun-5: printAndFreeGraph: clusterGraph.c:556\n");
              return 0;
            } /*If malloc failed to allocate memory*/

            continue;
        } /*If there are few reads in cluster then required (do not print)*/

        /*Print & free the cluster (Function handles swaping graphRoot around)*/
        if(printAndFreeCluster(graphRoot, graphRoot, clustOnULng) == 0)
        { /*If malloc failed to allocate memory in printing cluster*/
            fprintf(stderr, "fun-5: printAndFreeGraph: clusterGraph.c:566\n");
            return 0;
        } /*If malloc failed to allocate memory in printing cluster*/

        clustOnULng++;                         /*Only count assigned clusters*/
    } /*While there are clusters to print out*/

    return 1;
} /*printAndFreeGraph*/

/*##############################################################################
# Output: 
#    Modifies: clustToBeak to have its edges as neighbors (new clusters with
#              one node)
#    Modifies: sets rootGraph to new node when the graphRoot is freeded
#    Frees: if freeHeadNodeChar is 1, frees clustToBreak graphNode
#    Return: 4: If nothing to break up
#            1: if succeeded
#            0: if malloc failed
# Note: This function sets clustToBreak to the head node of the cluster
##############################################################################*/
char breakCluster(
    struct graphNode **clustToBreak, /*Cluster of graph nodes to break up*/
    struct graphNode **graphRoot    /*root of read graph*/
) /*Breaks cluster of graphNode's into single node (self only) clusters*/
{ /*breakCluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 TOC: breakCluster
    #    fun-6 sec-1: varialbe declerations
    #    fun-6 sec-2: Find head node and intialize que
    #    fun-6 sec-3: Break all nodes into separate cluster
    #    fun-6 sec-4: Reset orginal cluster head variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1: varialbe declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNodeStack *graphStack = 0;    /*Itterates though stack*/
    struct graphNode *nextGraphNode = 0,
                     *tmpGraphNode = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-2: Find head node and intialize que
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(clustToBreak == 0)
        return 4; /*Nothing to break up*/

    if((*clustToBreak)->clustNode != *clustToBreak)
        nextGraphNode = (*clustToBreak)->clustNode;
    else
        nextGraphNode = *clustToBreak;

    if(nextGraphNode->edgesList != 0)
    { /*If have edges to break as well*/
        if(pushGraphNodeStack(nextGraphNode->edgesList, &graphStack) == 0)
        { /*If malloc failed to allocate memory in printing cluster*/
            fprintf(stderr, "fun-6: breakCluster: clusterGraph.c:631\n");
            return 0;
        } /*If malloc failed to allocate memory in printing cluster*/
    } /*If have edges to break as well*/
 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-3: Break all nodes into separate cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(graphStack != 0)
    { /*While there are nodes to print out in this cluster*/
       nextGraphNode = graphStack->graphStruct;
       popGraphStackStruct(&graphStack);

       while(nextGraphNode != 0)
       { /*While there is a another neighbor node to print*/
           /*Add edges to stack, free node on, & move onto neihbor node*/
           if(nextGraphNode->edgesList != 0)
           { /*If have edges I need to loop though*/
               if(pushGraphNodeStack(nextGraphNode->edgesList, &graphStack) ==0)
               { /*If malloc failed to allocate memory in printing cluster*/
                   fprintf(stderr, "fun-6: breakCluster: clusterGraph.c:653\n");
                   return 0;
               } /*If malloc failed to allocate memory in printing cluster*/

               nextGraphNode->edgesList = 0;
           } /*If have edges I need to loop though*/

           tmpGraphNode = nextGraphNode->nextNode;

           /*Reset all variables in the new, single read cluster*/

           if((*graphRoot)->nextNode != 0)
               (*graphRoot)->nextNode->prevNode = nextGraphNode;

           nextGraphNode->nextNode = (*graphRoot)->nextNode;
           nextGraphNode->prevNode = *graphRoot;
           nextGraphNode->edgesList = 0;
           nextGraphNode->numNodesULng = 1;
           nextGraphNode->readCntULng = 0;
           nextGraphNode->clustNode = nextGraphNode;
           /*Not using cutNodeFromGraph due to cutNode looping through edges*/

           (*graphRoot)->nextNode = nextGraphNode;   /*Add node as new cluster*/
           nextGraphNode = tmpGraphNode;        /*Move to next node in cluster*/
       } /*While there is a another neighbor node to print*/
    } /*While there are nodes to print out in this cluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-4: Reset orginal cluster head variables
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    (*clustToBreak)->edgesList = 0;
    (*clustToBreak)->clustNode = *clustToBreak;  /*Just being safe*/
    (*clustToBreak)->numNodesULng = 1;
    (*clustToBreak)->readCntULng = 0;

    return 1;
} /*breakCluster*/

/*##############################################################################
# Output: 
#    Frees: frees clustToBreak graphNode
#    Other Output: Check breakCluster function (fun-6)
#    Return: 4: If nothing to break up
#            1: if succeded
#            0: if malloc failed to allocate memory in break cluster
##############################################################################*/
char breakClusterAndFreeHead(
    struct graphNode **clustToBreak, /*Cluster of graph nodes to break up*/
    struct graphNode **graphRoot     /*root of read graph*/
) /*Frees head of cluster & breaks edges into single node (self only) clusters*/
{ /*breakClusterAndFreeHead*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-7 Sec-1 TOC: Break clusters up using brakeCluster & free head node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(clustToBreak == 0)
        return 4; /*Nothing to break up*/

    if(breakCluster(clustToBreak, graphRoot) == 0)
    { /*If malloc failed to allocate memory in printing cluster*/
        fprintf(stderr, "fun-7: breakClusterAndFreeHead: clusterGraph.c:711\n");
        return 0;
    } /*If malloc failed to allocate memory in printing cluster*/

    if(*graphRoot == *clustToBreak)
    { /*If removing the root node*/
        *graphRoot = (*clustToBreak)->nextNode;
        (*graphRoot)->prevNode = 0;
    } /*If removing the root node*/
   
    else
    { /*Else if just need to cut the node from the list before freeing*/
        if((*clustToBreak)->nextNode != 0)
            (*clustToBreak)->nextNode->prevNode = (*clustToBreak)->prevNode;

        if((*clustToBreak)->prevNode != 0)
            (*clustToBreak)->prevNode->nextNode = (*clustToBreak)->nextNode;
    } /*Else if just need to cut the node from the list before freeing*/

    freeGraphNodeStruct(clustToBreak);
    *clustToBreak = 0;

    return 1;
} /*breakClusterAndFreeHead*/

/*##############################################################################
# Output: 
#    Frees: every node in the cluster and sets clustToFree to 0
#    return: 4: If nothing to free
#            1: if succeded
#            0: if malloc failed to make memory for stack
##############################################################################*/
char freeCluster(
    struct graphNode **clustToFree, /*Cluster to free*/
    struct graphNode **graphRoot,   /*Root of graph, reset if cluster is root*/
    unsigned long clustOnULng      /*Id of cluster on (for printing)*/
) /*Frees all nodes in a cluster, may print */
{ /*freeCluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 TOC: Free all nodes in a cluster
    #    fun-8 sec-1: variable declerations
    #    fun-8 sec-2: Get cluster head, cut cluster from graph, set up stack
    #    fun-8 sec-3: print & free the head node of the cluster
    #    fun-8 sec-4: Free all nodes in the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *nextGraphNode = 0,
                     *freeGraphNode = 0;

    struct graphNodeStack *graphStack = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-2: Get cluster head, cut cluster from graph, set up stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*clustToFree == 0)
        return 4; /*Nothing to free*/

    *clustToFree = (*clustToFree)->clustNode;     /*Make sure on cluster head*/
    nextGraphNode = (*clustToFree)->nextNode;

    /*cut cluster from list*/
    if((*clustToFree)->prevNode != 0)
        (*clustToFree)->prevNode->nextNode = (*clustToFree)->nextNode;

    if((*clustToFree)->nextNode != 0)
        (*clustToFree)->nextNode->prevNode = (*clustToFree)->prevNode;

    if((*clustToFree)->edgesList != 0)
    { /*If have edges to push into the stack*/
        if(pushGraphNodeStack((*clustToFree)->edgesList, &graphStack) == 0)
        { /*If malloc failed to allocate memory in printing cluster*/
            fprintf(stderr, "fun-8: freeCluster: clusterGraph.c:783\n");
           return 0;
        } /*If malloc failed to allocate memory in printing cluster*/
    } /*If have edges to push into the stack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-3: print & free the head node of the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    freeGraphNodeStruct(clustToFree);   /*Neighbor nodes are other clusters*/
    *clustToFree = nextGraphNode;       /*Move onto the next cluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-8 Sec-4: Free all nodes in the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(graphStack != 0)
    { /*While there are nodes to print out in this cluster*/
       nextGraphNode = graphStack->graphStruct;
       popGraphStackStruct(&graphStack);

       while(nextGraphNode != 0)
       { /*While there are neighbor nodes to free*/
           if(nextGraphNode->edgesList != 0)
           { /*If have edges to add to the stack*/
              if(pushGraphNodeStack(nextGraphNode->edgesList, &graphStack) == 0)
              { /*If malloc failed to allocate memory in printing cluster*/
                  fprintf(stderr, "fun-8: freeCluster: clusterGraph.c:814\n");
                  return 0;
              } /*If malloc failed to allocate memory in printing cluster*/
           } /*If have edges to add to the stack*/

           freeGraphNode = nextGraphNode;
           nextGraphNode = nextGraphNode->nextNode;
           freeGraphNodeStruct(&freeGraphNode);
       } /*While there are neighbor nodes to free*/
    } /*While there are nodes to print out in this cluster*/

    return 1;
} /*freeCluster*/

/*##############################################################################
# Output: 
#    Frees: every node in the cluster and sets clustToFree to 0
#    stdout: prints each read name & cluster id (name\tid) to stdout
#    return: 4: If nothing to free
#            1: if succeded
#            0: if malloc failed to make memory for stack
##############################################################################*/
char printAndFreeCluster(
    struct graphNode **clustToFree, /*Cluster to free*/
    struct graphNode **graphRoot,   /*Root of graph, reset if cluster is root*/
    unsigned long clustOnULng      /*Id of cluster on (for printing)*/
) /*Frees all nodes in a cluster, may print */
{ /*printAndFreeCluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 TOC: Free and prings all nodes in a cluster
    #    fun-9 sec-1: variable declerations
    #    fun-9 sec-2: Get cluster head, cut cluster from graph, set up stack
    #    fun-9 sec-3: print & free the head node of the cluster
    #    fun-9 sec-4: Free all nodes in the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *nextGraphNode = 0,
                     *freeGraphNode = 0;

    struct graphNodeStack *graphStack = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-2: Get cluster head, cut cluster from graph, set up stack
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*clustToFree == 0)
        return 4; /*Nothing to free*/

    *clustToFree = (*clustToFree)->clustNode;     /*Make sure on cluster head*/
    nextGraphNode = (*clustToFree)->nextNode;

    /*cut cluster from list*/
    if((*clustToFree)->prevNode != 0)
        (*clustToFree)->prevNode->nextNode = (*clustToFree)->nextNode;

    if((*clustToFree)->nextNode != 0)
        (*clustToFree)->nextNode->prevNode = (*clustToFree)->prevNode;

    if((*clustToFree)->edgesList != 0)
    { /*If have edges to push into the stack*/
        if(pushGraphNodeStack((*clustToFree)->edgesList, &graphStack) == 0)
        { /*If malloc failed to allocate memory in printing cluster*/
            fprintf(stderr, "fun-9: freeCluster: clusterGraph.c:783\n");
           return 0;
        } /*If malloc failed to allocate memory in printing cluster*/
    } /*If have edges to push into the stack*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-3: print & free the head node of the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if is a node to print out*/
    if(
       (*clustToFree)->seqIdNode->doneChar != 0 &&
       (*clustToFree)->seqIdNode->doneChar != 4
    ) 
        printf("%s\t%lu\n",(*clustToFree)->seqIdNode->idCStr,clustOnULng);

    freeGraphNodeStruct(clustToFree);   /*Neighbor nodes are other clusters*/
    *clustToFree = nextGraphNode;       /*Move onto the next cluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-9 Sec-4: Free all nodes in the cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(graphStack != 0)
    { /*While there are nodes to print out in this cluster*/
       nextGraphNode = graphStack->graphStruct;
       popGraphStackStruct(&graphStack);

       while(nextGraphNode != 0)
       { /*While there are neighbor nodes to free*/
           if(nextGraphNode->edgesList != 0)
           { /*If have edges to add to the stack*/
              if(pushGraphNodeStack(nextGraphNode->edgesList, &graphStack) == 0)
              { /*If malloc failed to allocate memory in printing cluster*/
                  fprintf(stderr, "fun-9: freeCluster: clusterGraph.c:814\n");
                  return 0;
              } /*If malloc failed to allocate memory in printing cluster*/
           } /*If have edges to add to the stack*/

            /*Check if is a node to print out*/
            if(
               nextGraphNode->seqIdNode->doneChar != 0 &&
               nextGraphNode->seqIdNode->doneChar != 4
            ) 
               printf("%s\t%lu\n",nextGraphNode->seqIdNode->idCStr,clustOnULng);
  
           freeGraphNode = nextGraphNode;
           nextGraphNode = nextGraphNode->nextNode;
           freeGraphNodeStruct(&freeGraphNode);
       } /*While there are neighbor nodes to free*/
    } /*While there are nodes to print out in this cluster*/

    return 1;
} /*printAndFreeCluster*/
