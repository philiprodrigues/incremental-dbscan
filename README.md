# incremental-dbscan

A streaming implementation of [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) clustering for 2D LArTPC hits.

## Motivation

A typical implementation of DBSCAN is given all of its input points in one go. In contrast, in the DUNE DAQ trigger system, hits form a continuous stream in time, and are delivered to the hit-clustering algorithm in order of hit (start) time. `incremental-dbscan` takes one input hit at a time and updates its list of clusters using the usual DBSCAN conditions. Because the input hits are known to be ordered by hit start time, a number of optimizations are possible: for example, to find the neighbours of a newly-added hit, we do not need to look further back in the (time-sorted) list of hits than `eps` (the distance threshold parameter in the DBSCAN algorithm).
