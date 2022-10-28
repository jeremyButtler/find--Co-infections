/*##############################################################################
# Name: clusterGraphMoveNodes.c
# Use: Holds functions to move nodes around in a read graph
##############################################################################*/

#include "clusterGraphMoveNodes.h"

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC: clusterGraphMoveNodes
#    fun-1: mvToEdge: Moves a graph node to an edge
#    fun-2: mvToNeighbor: Moves a graph node to an neighbor
#    fun-3: mvToNewClust: Moves a graph node to a new cluster
#                         Edges stay in original cluster
#    fun-4: checkAndMoveRoot: Check if moving root node, if so prepare to move
#    fun-5: cutNodeFromList: Cuts a node out of its nextnode list
#    fun-6: resetCluster: reset cluster to new head node
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

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
    struct graphNode **graphRoot   /*handle moving root node*/
) /*moves a graphNode to an edge of another graphNode*/
{ /*mvToEdge*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 TOC: mvGraphEdge
    #    fun-1 sec-1: Do checks & cut out move node
    #    fun-1 sec-2: Update and swap clusters
    #    fun-1 sec-3: Move moveNode into moveToNodes edges
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Do checks & cut out move node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *checkMalSuc = 0; /*Checks if malloc errored out*/

    if(moveToNode == 0 || moveNode == 0)
        return moveNode;                     /*Nothing to do*/
    else if(moveToNode == moveNode)
        return moveNode;

    checkMalSuc = cutNodeFromList(moveNode, graphRoot);

    if(checkMalSuc == 0)
    { /*If malloc failed to allocate memory*/
        fprintf(stderr, "fun-1: mvToEdge: clusterGraphMoveNodes.c:55\n");
        return 0;         /*Malloc did not assign memory to readGraphNode*/
    } /*If malloc failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-2: Update and swap clusters
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode->clustNode != 0 && moveNode->clustNode->numNodesULng > 0)
         moveNode->clustNode->numNodesULng--; /*Account for moved node*/

    if(moveToNode->clustNode == 0)                   /*More of a safe guard*/
    { /*If need to make a new cluster*/
        moveToNode->clustNode = moveToNode;
        moveToNode->numNodesULng = 2;
    } /*If need to make a new cluster*/

    else if(moveNode->clustNode != moveToNode->clustNode)
        moveToNode->clustNode->numNodesULng++;         /*Account for moveNode*/

    moveNode->clustNode = moveToNode->clustNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Move moveNode into moveToNodes edges
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    moveNode->prevNode = moveToNode;     /*So can keep parent edge ptr updated*/
    moveNode->nextNode = moveToNode->edgesList;

    if((*moveToNode).edgesList != 0)
        moveToNode->edgesList->prevNode = moveNode;

    moveToNode->edgesList = moveNode;
    return moveNode;
} /*mvToEdge*/

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
    struct graphNode **graphRoot   /*handle moving root node*/
) /*moves graphNode from edge in a node to a neighbor in another node*/
{ /*mvToNeighbor*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 TOC: mvGraphEdge
    #    fun-2 sec-1: Do checks & cut out move node
    #    fun-2 sec-2: Update and swap clusters
    #    fun-2 sec-3: Move node to neihbor
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-1: Do checks & cut out move node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *checkMalSuc = 0; /*Checks if malloc errored out*/

    if(moveToNode == 0 || moveNode == 0)
        return moveNode;                     /*Nothing to do*/
    else if(moveToNode == moveNode)
        return moveNode;

    checkMalSuc = cutNodeFromList(moveNode, graphRoot);

    if(checkMalSuc == 0)
    { /*If malloc failed to allocate memory*/
        fprintf(stderr, "fun-2: mvToNeighbor: clusterGraphMoveNodes.c:128\n");
        return 0;         /*Malloc did not assign memory to readGraphNode*/
    } /*If malloc failed to allocate memory*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-2: Update and swap clusters
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode->clustNode != 0 && moveNode->clustNode->numNodesULng > 0)
         moveNode->clustNode->numNodesULng--; /*Account for moved node*/

    if(moveToNode->clustNode == 0)                   /*More of a safe guard*/
    { /*If need to make a new cluster*/
        moveToNode->clustNode = moveToNode;
        moveToNode->numNodesULng = 2;
    } /*If need to make a new cluster*/

    else if(moveNode->clustNode != moveToNode->clustNode)
        moveToNode->clustNode->numNodesULng++;         /*Account for moveNode*/

    moveNode->clustNode = moveToNode->clustNode;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-2 Sec-3: Move node to neihbor
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    moveNode->prevNode = moveToNode;     /*So can keep parent edge ptr updated*/
    moveNode->nextNode = moveToNode->nextNode;

    if(moveToNode->nextNode != 0)
        moveToNode->nextNode->prevNode = moveNode;

    moveToNode->nextNode = moveNode;

    return moveNode;
} /*mvToNeighbor*/

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
    struct graphNode **graphRoot   /*Each cluster is a neighbor to Root node*/
) /*moves graphNode to be a new cluster & neighbor next to moveToNode*/
{ /*mvToNewClust*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC: mvGraphEdge
    #    fun-3 sec-1: Do checks & cut out move node
    #    fun-3 sec-2: Update old cluster & create new cluster
    #    fun-3 sec-3: Move node to neihbor
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Do checks & cut out move node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *checkMalSuc = 0; /*Checks if malloc errored out*/

    if(*graphRoot == 0 || moveNode == 0)
        return moveNode;                     /*Nothing to do*/

    checkMalSuc = cutNodeFromList(moveNode, graphRoot);

    if(checkMalSuc == 0)
    { /*If malloc failed to allocate memory*/
        fprintf(stderr, "fun-3: mvToNewClusT: clusterGraphMoveNodes.c:199\n");
        return 0;         /*Malloc did not assign memory to readGraphNode*/
    } /*If malloc failed to allocate memory*/

    if(*graphRoot == moveNode)
        return moveNode;          /*Root node only node graph (nothing to do)*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Update old cluster & create new cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode->clustNode != 0 && moveNode->clustNode->numNodesULng > 0)
         moveNode->clustNode->numNodesULng--; /*Account for moved node*/

    moveNode->clustNode = moveNode;
    moveNode->numNodesULng = 1;
    moveNode->readCntULng = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-3: Move node to neihbor
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    moveNode->prevNode = *graphRoot;     /*So can keep parent edge ptr updated*/
    moveNode->nextNode = (*graphRoot)->nextNode;

    if((*graphRoot)->nextNode != 0)
        (*graphRoot)->nextNode->prevNode = moveNode;

    (*graphRoot)->nextNode = moveNode;

    return moveNode;
} /*mvToNewClust*/

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
) /*Checks if node to move is root, if so makes another node root*/
{ /*checkAndMoveRoot*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # TOC: checkAndMoveRoot
    #    fun-4 sec-1: variable declerations
    #    fun-4 sec-2: Check if need to swap root & handle no edge node cases
    #    fun-4 sec-3: Make edges into cluster & see if no edges in new cluster
    #    fun-4 sec-4: Move neighbor nodes in cluster to empty edgesList
    #    fun-4 sec-5: Move neighobr nodes in cluster to full edgesList
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNode *neighList = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-2: Check if need to swap root & handle no edge node cases
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode != *graphRoot)
        return moveNode; /*Not moving root, return*/

    if(moveNode->edgesList == 0)
    { /*If have no edges to swap*/
        if(moveNode->nextNode == 0)
            return moveNode;                    /*Root node only node in graph*/

        *graphRoot = (*graphRoot)->nextNode;     /*No need to save cluster*/
        (*graphRoot)->prevNode = moveNode->prevNode; /*Should be 0*/
        moveNode->nextNode = 0;
        moveNode->prevNode = 0;                      /*Just to be safe*/

        return moveNode;
    } /*If have no edges to swap*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-3: Make edges into cluster & see if no edges in new cluster
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    *graphRoot = moveNode->edgesList;        /*move node is graphRoot*/
    moveNode->edgesList = 0;                 /*Avoid circular lists*/
    (*graphRoot)->prevNode = moveNode->prevNode; /*Always null*/
    (*graphRoot)->clustNode = *graphRoot;
    (*graphRoot)->numNodesULng = 1;
        /*Previously moveNode->edgesList->prevNode was moveNode*/

    if(moveNode->prevNode != 0) /*Should never fire*/
         moveNode->prevNode->nextNode = *graphRoot;

    moveNode->prevNode = 0;        /*Just to make sure moveNode fully cut out*/

    if((*graphRoot)->nextNode == 0)
    { /*If can cut in orginal roots neigbhor nodes*/
        *graphRoot = moveNode->nextNode;    /*Open neighbor node*/
        moveNode->nextNode = 0;

        if((*graphRoot)->nextNode != 0)              /*Will fire*/
            (*graphRoot)->nextNode->prevNode = *graphRoot;

        if(resetCluster(*graphRoot) == 0)
        { /*If malloc failed to allocate memory*/
            fprintf(
                   stderr,
                   "fun-4: checkAndMoveRoot: clusterGraphMoveNodes.c:307\n"
            ); /*Print out the error message*/
            return 0;         /*Malloc did not assign memory to readGraphNode*/
        } /*If malloc failed to allocate memory*/

        return moveNode;
    } /*If can cut in orginal roots neigbhor nodes*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-4: Move neighbor nodes in cluster to empty edgesList
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if((*graphRoot)->edgesList == 0)
    { /*If no edges, only neighbors in same cluster in new cluster*/
        (*graphRoot)->edgesList = (*graphRoot)->nextNode;
        (*graphRoot)->nextNode = moveNode->nextNode;

        if(moveNode->nextNode != 0)
            moveNode->nextNode->prevNode = *graphRoot;

        moveNode->nextNode = 0;

        return moveNode;
    } /*If no edges, only neighbors in same cluster in new cluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-4 Sec-5: Move neighobr nodes in cluster to full edgesList
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    neighList = (*graphRoot)->nextNode;

    while(neighList->nextNode != 0)
        neighList = neighList->nextNode; 

    /*Make space for the old clusters*/
    neighList->nextNode = (*graphRoot)->edgesList;
    (*graphRoot)->edgesList->prevNode = neighList;
    (*graphRoot)->edgesList = (*graphRoot)->nextNode;

    /*Move the old clusters to the new root*/
    (*graphRoot)->nextNode = moveNode->nextNode;
    moveNode->nextNode = 0;

    if((*graphRoot)->nextNode != 0)
        (*graphRoot)->nextNode->prevNode = *graphRoot;

    if(resetCluster(*graphRoot) == 0)/*set clust head to first edge*/
    { /*If malloc failed to allocate memory*/
       fprintf(stderr,  "fun-4: checkAndMoveRoot: clusterGraphMoveNodes.c:357\n");
       return 0;         /*Malloc did not assign memory to readGraphNode*/
    } /*If malloc failed to allocate memory*/

    return moveNode;
} /*checkAndMoveRoot*/

/*##############################################################################
# Output:
#    Modifies list to no longer have target graphNode in it
#    Returns: A pointer to the cut out node
##############################################################################*/
struct graphNode * cutNodeFromList(
    struct graphNode *moveNode,  /*Node cutting from list*/
    struct graphNode **graphRoot /*Root node of graph*/
) /*Cuts a node out of its nextnode list*/
{ /*cutNodeFromList*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 TOC: checkAndMoveRoot
    #    fun-5 sec-1: variable declerations
    #    fun-5 sec-2: Find the tail graphNode of the edgesList in moveNode
    #    fun-5 sec-3: Adjust pointer to the previous graphNode
    #    fun-5 sec-4: Adjust pointer of the next graphNode
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 sec-1: variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char malocSucessChar = 0;               /*Marks if malloc errored out*/
    struct graphNode *edgeListTail = 0;     /*Tail graphNode in list of edges*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-2: Find the tail graphNode of the edgesList in moveNode
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode->clustNode == moveNode)
    { /*If have to change the cluster head*/
        if(moveNode->edgesList != 0)
        { /*If having the old edges keep the cluster*/
            moveNode->edgesList->numNodesULng = moveNode->numNodesULng--;
        } /*If having the old edges keep the cluster*/
    } /*If have to change the cluster head*/

    if(moveNode == *graphRoot)
    { /*If need to change the root node of the graph (cutting it out)*/
        if(checkAndMoveRoot(moveNode, graphRoot) == 0)
        { /*If malloc failed to allocate memory*/
         fprintf(stderr, "fun-5: cutNodeFromList: clusterGraphMoveNodes.c:408\n");
          return 0;         /*Malloc did not assign memory to readGraphNode*/
        } /*If malloc failed to allocate memory*/

        return moveNode;                       /*was root, move node cut out*/
    } /*If need to change the root node of the graph (cutting it out)*/

    if(moveNode->edgesList != 0)
    { /*If have edges I need to move around*/
        edgeListTail = moveNode->edgesList;

        while(edgeListTail->nextNode != 0)
            edgeListTail = edgeListTail->nextNode; 
    } /*If have edges I need to move around*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-3: Adjust pointer to the previous graphNode
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Cut the node to move from the edgesList or neighbor list*/
    if(moveNode->prevNode != 0)
    { /*If move node is in the middle of a neighbor list*/

        if(moveNode == moveNode->prevNode->edgesList)
        { /*If move node is the head of the edges list*/
            if(moveNode->edgesList == 0)
                moveNode->prevNode->edgesList = moveNode->nextNode;
    
            else /*So only the node moves*/
            { /*Else have to add edges in*/
                moveNode->prevNode->edgesList = moveNode->edgesList;
                moveNode->edgesList->prevNode = moveNode->prevNode;
            } /*Else have to add edges in*/

        } /*If move node is the head of the edges list*/

        else
        { /*else move node is another node in the edgesList*/
            if(moveNode->edgesList == 0)
                moveNode->prevNode->nextNode = moveNode->nextNode;
    
            else /*So only the node moves*/
            { /*Else have to add edges in*/
                moveNode->prevNode->nextNode = moveNode->edgesList;
                moveNode->edgesList->prevNode = moveNode->prevNode;
            } /*Else have to add edges in*/
        } /*else move node is another node in the edgesList*/
    } /*If move node is in the middle of a neighbor list*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-5 Sec-4: Adjust pointer of the next graphNode
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(moveNode->nextNode != 0)
    { /*If move node is in the middle of a neighbor list*/
        if(moveNode->edgesList == 0)
            moveNode->nextNode->prevNode = moveNode->prevNode;

        else /*So only the node moves*/
        { /*Else I have to cut the edges into the cluster*/
            moveNode->nextNode->prevNode = edgeListTail;
            edgeListTail->nextNode = moveNode->nextNode;
        } /*Else I have to cut the edges into the cluster*/
    } /*If move node is in the middle of a neighbor list*/

    /*Set cluster to the head edge*/
    malocSucessChar = resetCluster(moveNode->edgesList);
 
    /*Finalize cut*/
    moveNode->prevNode = 0;
    moveNode->nextNode = 0;
    moveNode->edgesList = 0; /*So no longer reference the edges*/

    if(malocSucessChar == 0)/*set clust head to first edge*/
    { /*If malloc failed to allocate memory*/
       fprintf(stderr,  "fun-5: cutNodeFromList: clusterGraphMoveNodes.c:481\n");
       return 0;         /*Malloc did not assign memory to readGraphNode*/
    } /*If malloc failed to allocate memory*/

    return moveNode;
} /*cutNodeFromList*/

/*##############################################################################
# Output:
#    Modifies cluster to point to the new head node
#    Returns: 2: If nothing to do
#             1: if succeded
#             0: If malloc failed to assign memory
##############################################################################*/
char resetCluster(
    struct graphNode *headNode /*node with nodes to reset & new cluster head*/
) /*Resets the cluster of nodes in a cluster to the new head node*/
{ /*resetCluster*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-6 Sec-1 TOC: resetCluster: reset cluster to new head node
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    struct graphNodeStack *graphStack = 0,
                          *checkStackMal = 0;/*Check if malloc assigned memory*/
    struct graphNode *nextGraphNode = 0;
    unsigned long numNodesULng = 0;       /*recount of number nodes in cluster*/

    if(headNode == 0)
        return 2;                           /*Nothing to do*/

    if(headNode->edgesList != 0)
    { /*If have edges to also reset clusters for*/
        if(pushGraphNodeStack(headNode, &graphStack) == 0)
        { /*If malloc failed to allocate memory*/
           fprintf(stderr, "fun-6: resetCluster: clusterGraphMoveNodes.c:516\n");
           return 0;
        } /*If malloc failed to allocate memory*/
    } /*If have edges to also reset clusters for*/

    headNode->clustNode = headNode;

    while(graphStack != 0)
    { /*While there are nodes to print out in this cluster*/
       nextGraphNode = graphStack->graphStruct;
       popGraphStackStruct(&graphStack);

       while(nextGraphNode != 0)
       { /*While there are neighbor nodes to free*/
           if(nextGraphNode->edgesList != 0)
           { /*If have edges to also reset clusters for*/
              checkStackMal = pushGraphNodeStack(nextGraphNode->edgesList,
                                                 &graphStack
              ); /*Push a new node on the stack*/

              if(checkStackMal == 0)
              { /*If malloc failed to assign memory*/
                  while(graphStack != 0) popGraphStackStruct(&graphStack);
                  fprintf(
                      stderr,
                      "fun-6: resetCluster: clusterGraphMoveNodes.c:541\n"
                  ); /*Print out the error*/
                  return 0;
              } /*If malloc failed to assign memory*/
           } /*If have edges to also reset clusters for*/

           nextGraphNode->clustNode = headNode; 
           numNodesULng++;
           nextGraphNode = nextGraphNode->nextNode;
       } /*While there are neighbor nodes to free*/
    } /*While there are nodes to print out in this cluster*/

    headNode->numNodesULng = numNodesULng;
    return 1;
} /*resetCluster*/
