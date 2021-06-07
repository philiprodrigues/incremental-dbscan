#pragma once

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound
#include <set>
#include <list>

#include "Hit.hpp"

//======================================================================
// Find the eps-neighbours of hit q, assuming that the hits vector is sorted by time
int
neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps);

//======================================================================
struct Cluster
{
    // The index of this cluster
    int index{ -1 };
    // A cluster is kComplete if its hits are all kComplete, so no
    // newly-arriving hit could be a neighbour of any hit in the
    // cluster
    Completeness completeness{ Completeness::kIncomplete };
    // The latest time of any hit in the cluster
    float latest_time{ 0 };
    // The latest (largest time) "core" point in the cluster
    Hit* latest_core_point{ nullptr };
    // The hits in this cluster
    HitSet hits;

    // Add hit if it's a neighbour of a hit already in the
    // cluster. Precondition: time of new_hit is >= the time of any
    // hit in the cluster. Returns true if the hit was added
    bool maybe_add_new_hit(Hit* new_hit, float eps, int minPts);

    // Add the hit `h` to this cluster
    void add_hit(Hit* h);

    // Steal all of the hits from cluster `other` and merge them into
    // this cluster
    void steal_hits(Cluster& other);
};

//======================================================================

// A class to hold the state of the DBSCAN algorithm between calls to
// dbscan_partial_add_one
struct DBSCANState
{
    std::vector<Hit*> hits; // All the hits we've seen so far, in time order
    float latest_time{ 0 }; // The latest time of a hit in the vector of hits
    std::list<Cluster> clusters; // All of the currently-active (ie, kIncomplete) clusters
};

//======================================================================
//
// Starting from `seed_hit`, find all the reachable hits and add them
// to `cluster`
void
cluster_reachable(DBSCANState& state,
                  Hit* seed_hit,
                  Cluster& cluster,
                  float eps,
                  unsigned int minPts);

//======================================================================
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement
// that the hits are passed in time order
void
dbscan_partial_add_one(DBSCANState& state,
                       Hit* new_hit,
                       float eps,
                       unsigned int minPts);

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
