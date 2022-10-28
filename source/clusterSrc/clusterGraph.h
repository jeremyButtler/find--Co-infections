/*##############################################################################
# Name: clustGraphGraphFun.h
# Use: header file for clustGraphGraphFun.c
#      Holds functions needed to manipulate my graphNode graph
##############################################################################*/

#ifndef CLUSTERGRAPH_H
#define CLUSTERGRAPH_H

/*
  <stdlib.h> -> clustGraphReadTree.c -> clustGraphGraphAndTreeStructs.c
   strtod (convert c-string to double)
*/
#include <stdio.h>              /*For print functions*/
#include "clusterStructMinValues.h"
#include "clusterGraphMoveNodes.h"
    /*
      clusterGraphMoveNodes.h includes: clusterReadTree.h
      clusterTree Includes:
       - <string.h> for strcmp
       - clustGraphStackStructs.h
       - clustGraphGraphAndTreeStructs.h (through clustGraphStackStructs)
         - readInfo, graphNode, & bridgeNode structers & functions
    */

struct readInfo * insertGraphNode(
    char *lineInCStr,                  /*fromat: readName\treferenceName\tmapq*/
    struct readInfo **readTree, /*readInfo tree to search for read & ref names*/
    struct graphNode **graphRoot, /*Root node of read graph with clusters*/
    struct readInfo *oldReadNode, /*readInfo node from last aligment looked at*/
    struct minValues *minStats,  /*min thresholds to keep alignment or cluster*/
    unsigned long *mappedReadsULng, /*Number of references mapped to this read*/
    struct graphNodeStack **clustStack /*Clusters this read shares edges with*/
); /*Inserts a read into a read graph*/

double parseLine(
    char *lineInCStr,            /*Line with read name, reference name, & mapq*/
    char *readCStr,              /*Array to hold read name*/
    unsigned char *lenReadUChar, /*Holds length of read name*/
    char *refCStr,              /*Array to hold reference name*/
    unsigned char *lenRefUChar   /*Holds length of reference name*/
); /*Copies read and reference names in input, also returns mapq as double*/

char assignReadToCluster(
    struct graphNode *readToAssign,       /*Read to assign to a cluster*/
    struct graphNodeStack **clustStack,   /*Stack of clusters readToAssign
                                            shares edges with*/
    struct minValues *minStats,  /*min thresholds to keep alignment or cluster*/
    struct graphNode **graphRoot          /*Root node of read graph*/
); /*Assignes a read to a cluster, if mulitple assignments, merges clusters*/

/*##############################################################################
# Output:
#    Modifies: Merges cluster with fewest reads into the remaining cluster
#    Returns: 1: if sucess
#             0: If malloc failed to assign memory at some point
##############################################################################*/
char mergeClusters(
    struct graphNode **firstClust,   /*First cluster/bridge to merge*/
    struct graphNode **secondClust,  /*Second cluster/bridge to merge*/
    struct graphNode **graphRoot      /*Root node of the graph. If root node
                                        merged, then reassigns root node*/
); /*Merge two bridges (clusters) into one bridge (cluster)*/

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
); /*prints out the bridges (clusters) and frees a graphNode graph*/

/*##############################################################################
# Output: 
#    Modifies: clustToBeak to have its edges as neighbors (new clusters with
#              one node)
#    Modifies: sets rootGraph to new node when the graphRoot is freeded
#    Frees: if freeHeadNodeChar is 1, frees clustToBreak graphNode
# Note: This function sets clustToBreak to the head node of the cluster
##############################################################################*/
char breakCluster(
    struct graphNode **clustToBreak, /*Cluster of graph nodes to break up*/
    struct graphNode **graphRoot    /*root of read graph*/
); /*Breaks cluster of graphNode's into single node (self only) clusters*/

/*##############################################################################
# Output: 
#    Frees: frees clustToBreak graphNode
#    Other Output: Check breakCluster function (fun-6)
#    Return: 1: if succeded
#            0: if malloc failed to allocate memory in break cluster
##############################################################################*/
char breakClusterAndFreeHead(
    struct graphNode **clustToBreak, /*Cluster of graph nodes to break up*/
    struct graphNode **graphRoot     /*root of read graph*/
); /*Frees head of cluster & breaks edges into single node (self only) clusters*/

/*##############################################################################
# Output: 
#    Frees: every node in the cluster and sets clustToFree to 0
#    return: 1: if succeded
#            0: if malloc failed to make memory for stack
##############################################################################*/
char freeCluster(
    struct graphNode **clustToFree, /*Cluster to free*/
    struct graphNode **graphRoot,   /*Root of graph, reset if cluster is root*/
    unsigned long clustOnULng      /*Id of cluster on (for printing)*/
); /*Frees all nodes in a cluster, may print */

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
); /*Frees all nodes in a cluster, may print */

/*Used to have freeCluster take char as boolean, but after started ignoreing
  clusters with readInfo->doneChar = 4 it started acting up and not printing 
  out. So duplicated & made dedicated print & free function
*/
#endif
