# Cluster:

## Use:

Clusters reads by mapping quality.

## Reason exists:

This is not the best clustering algorithm, so it should not be used. However,
  this is in findCoInfections because it adds some extra flexibility, while
  adding no extra dependencies. It does not assume a certain number of clusters
  to start out with and it does not require n^2 memory usage. In terms of
  accuracy, this clustering algorithm has a tendency to ignore clusters, so it
  is not very good. In terms of speed, it requires n^2, instead of n^2/2
  sequences to be compared, so it will 2x slower. In terms of memory it requires
  O(n), like the single chain hierarchical clustering, which it may have a
  slight similarity to.

There are two ideas behind this. First, clustering of nanopore sequence reads
  to detect co-infections will often detect false positives. So, maybe being a
  bit more conservative is good. Second though memory is cheap, it is not n^2
  cheap. Basically, it takes 10^(2\*6) or one tera byte to process one million
  reads when you in n^2 usage. That being said the time demands alone would 
  likely prevent this large of a comparison. That being said, I think the sparse
  dataset test shows why this should never be used with a single sequence. To
  put it simply, this thing will misbin reads.

To sum it up. It is used here because with filtering using scoreReads and
  mapping qualities, it works. However, I would recommend finding another
  algorithm that has been designed by people how know what they are doing for
  your project. Their algorithm will work better and give better results.

## How works:

For clustering, each alignment represents and edge that connects the query node
  to the reference/mapped read node in a graph
  (for a diagram of clustering Figure-1). For each query we assume that the
  read is a unique cluster and assign any new reads that are edges as a node in
  the query's cluster. If the query was in a previous cluster, before assigning
  we assign it to a new cluster. Each reference (edge) that connects to another
  cluster is counted. After all alignments for the query have been check, all
  clusters that shared on or more edges with the query are counted. The query
  is merged with any cluster that shares over (user input) edges with the query.
  To ensure that no more than n nodes are merged, we always insert the cluster
  with the fewest reads into the merging cluster.

Sorry this is an ascii flow chart. Run this through ditaa with -S -o to build a 
  good png.

```ditaa
[-S -o]
+------------+---------------------+
| Alignment  |   Alignments        |
|            |                     |
|  o 1 -> 2  |    o 1 -> 4         |
|            |    o 1 -> 5         |
|            |                     |
|   +-+      |    +-+              |
|   |1|      |    |1|              |
|   +++      |    +++              |
|    | +-+   |     | +-+ +-+ +-+   |
|    +-+2|   |     +-+2+-+4+-+5|   |
|      +-+   |       +-+ +-+ +-+   |
+------------+---------------------+

+-----------------------+--------------------------+
| Alignment for read 2  | Alignments               |  Cluster two and one share
|                       |                          |  only one edge. So, both
|  o 2 -> 1             |  o 2 -> 3                |  clusters are different.
|                       |  o 2 -> 6                |
|   +-+         +-+     |                          |
|   |1+-=-=-=-->|2|     |   +-+         +-+        |
|   +++         +-+     |   |1+-=-=-=-->|2|        |
|    | +-+ +-+          |   +++         +++        |
|    +-+4+-+5|          |    | +-+ +-+   | +-+ +-+ |
|      +-+ +-+          |    +-+4+-+5|   +-+3+-+6| |
|                       |      +-+ +-+     +-+ +-+ |
+-----------------------+--------------------------+

+-------------------------------+---------------------------------+
| Alignment for read 3          | Alignments                      |
|                               |                                 |
|  o 3 -> 2                     |  o 3 -> 4                       |
|                               |  o 3 -> 5                       |
|                               |  o 3 -> 6                       |
|                               |                                 |
|   +-+         +-+         +-+ |    +-+          +-+         +-+ |
|   |1|         |2+-=-=-=-->|3| |    |1|          |2+-=-=-=-->|3| |
|   +++         +++         +-+ |    +++          +++     +-->+-+ |
|    | +-+ +-+   | +-+          |     | +-+ +-+    | +-+  :   ^^  |
|    +-+4+-+5|   +-+6|          |     +-+4+-+5|    +-+6+--+   ::  |
|      +-+ +-+     +-+          |       +++ +++      +-+      ::  |
|                               |        :   :                ::  |
|                               |        :   +-=-=-=-=-=-=-=--+:  |
|                               |        +-=----=-=-=-=-=-=----+  |
+-------------------------------+---------------------------------+

+-------------------------------+--------------------------------+
| Cluster one and three share   | Read three (old cluster three) |
| more than one edge. Assuming  | and cluster two share more     |
| both clusters are the same.   | than one edge.                 |
|                               |                                |
|   +-+             +-+         |   +-+                          |
|   |1| +-=-=-=-=-=-+2|         |   |1|                          |
|   +++ :           +++         |   +++                          |
|    |  v            |          |    |                           |
|    | +++ +-+ +-+   | +-+      |    | +-+ +-+ +-+ +-+           |
|    +-+3+-+4+-+5|   +-+6|      |    +-+2+-+3+-+4+-+5|           |
|      +-+ +-+ +-+     +++      |      +++ +-+ +-+ +-+           |
|       ^               :       |       |                        |
|       :               :       |      +++                       |
|       +-=-=-=-=-=-=-=-+       |      |6|                       |
|                               |      +-+                       |
+-------------------------------+--------------------------------+

+----------------------------------------------------+
| Move onto alignments for reads 4, 5, 6, ..., and n |
+----------------------------------------------------+

Figure-1: A digram of how the clustering algorthim works. o query -> reference,
  shows a mapping between two reads. Dashed lines indicate that two reads share
  an edge, but are not in the same cluster.
```
