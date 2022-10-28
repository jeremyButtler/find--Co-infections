/*##############################################################################
# Name: clusterGraphMoveNodes.h
# Use: Holds header for functions to move nodes around in a read graph
##############################################################################*/

#ifndef CLUSTERGRAPHMOVENODES_h
#define CLUSTERGRAPHMOVENODES_h

#include "clusterReadTree.h"

/*##############################################################################
# Output:
#    Modifies: moveToNode to have moveNode as an edge
#        - Does not move moveNode's edges (edges stay in original cluster)
#    Decurments: numNodesULng in moveNode's old cluster
#    Incurments: numNodesULng in moveToNode's cluster
#    Returns: moveNode: If malloc did not error out
#             0: If malloc code not assign memory in cutNodeFromList
##############################################################################*/
struct graphNode * mvToEdge(
    struct graphNode *moveNode,    /*node to move*/
    struct graphNode *moveToNode,  /*Node to assign moveNode to*/
    struct graphNode **graphRoot  /*handle moving root node*/
); /*moves a graphNode to an edge of another graphNode*/

/*##############################################################################
# Output:
#    Modifies: moveToNode to have moveNode as an edge or neigbhor
#    Decurments: numNodesULng in moveNode's old cluster
#    Incurments: numNodesULng in moveToNode's cluster
#    Returns: moveNode: If malloc did not error out
#             0: If malloc code not assign memory in cutNodeFromList
##############################################################################*/
struct graphNode * mvToNeighbor(
    struct graphNode *moveNode,    /*node to move*/
    struct graphNode *moveToNode,  /*Node to assign moveNode to*/
    struct graphNode **graphRoot  /*handle moving root node*/
); /*moves graphNode from edge in a node to a neighbor in another node*/

/*##############################################################################
# Output:
#    Modifies: moveToNode to have moveNode as an edge or neigbhor
#    Decurments: numNodesULng in moveNode's old cluster
#    Incurments: numNodesULng in moveToNode's cluster
#    Returns: moveNode: If malloc did not error out
#             0: If malloc code not assign memory in cutNodeFromList
##############################################################################*/
struct graphNode * mvToNewClust(
    struct graphNode *moveNode,    /*node to move*/
    struct graphNode **graphRoot  /*Each cluster is a neighbor to Root node*/
); /*moves graphNode to be a new cluster & neighbor next to moveToNode*/

/*##############################################################################
# Output:
#    Modifies: graphRoot to be the edges of the input node
#    Modifies: If move nodes is root, cuts move node from list & sets nextNode
#              and prevNode pointers to 0
#    Return: node moved: if suceeded (Assumes you alread did root check)
#            0: if malloc failed to allocate memory
##############################################################################*/
struct graphNode * checkAndMoveRoot(
    struct graphNode *moveNode,  /*Node moving from root (or checking)*/
    struct graphNode **graphRoot /*Root node to check or modify*/
); /*Checks if node to move is root, if so makes another node root*/

/*##############################################################################
# Output:
#    Modifies list to no longer have target graphNode in it
#    Returns: A pointer to the cut out node
##############################################################################*/
struct graphNode * cutNodeFromList(
    struct graphNode *moveNode,  /*Node cutting from list*/
    struct graphNode **graphRoot /*Root node of graph*/
); /*Cuts a node out of its nextnode list*/

/*##############################################################################
# Output:
#    Modifies cluster to point to the new head node
#    Returns: 2: If nothing to do
#             1: if succeded
#             0: If malloc failed to assign memory
##############################################################################*/
char resetCluster(
    struct graphNode *headNode /*node with nodes to reset & new cluster head*/
); /*Resets the cluster of nodes in a cluster to the new head node*/

#endif
