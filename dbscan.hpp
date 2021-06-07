#pragma once

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound
#include <set>
#include <list>

#include "Hit.hpp"

//======================================================================
// Find the neighbours of hit q, assuming that the hits vector is sorted by time
int
neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps);

//======================================================================
struct Cluster
{
    int index{-1};
    Completeness completeness{Completeness::kIncomplete};
    float latest_time{0};
    Hit* latest_core_point{nullptr};
    HitSet hits;

    // Add hit if it's a neighbour of a hit already in the
    // cluster. Precondition: time of new_hit is >= the time of any
    // hit in the cluster
    bool maybe_add_new_hit(Hit* new_hit, float eps, int minPts);

    void add_hit(Hit* h);

    void steal_hits(Cluster& other);
};

//======================================================================
struct State
{
    std::vector<Hit*> hits; // All the hits we've seen so far, in time order
    float latest_time{0}; // The latest time of a hit in the vector of hits
    std::list<Cluster> clusters;
};

//======================================================================
void
cluster_reachable(State& state,
                  Hit* seed_hit,
                  Cluster& cluster,
                  float eps,
                  unsigned int minPts);

//======================================================================
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement that the hits are passed in time order
void
dbscan_partial_add_one(State& state,
                       Hit* new_hit,
                       float eps,
                       unsigned int minPts);

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
