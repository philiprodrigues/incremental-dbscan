#pragma once

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound
#include <set>
#include <list>

#include "Hit.hpp"

namespace dbscan {
//======================================================================
// Find the eps-neighbours of hit q, assuming that the hits vector is sorted by
// time
int
neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps, int minPts);

//======================================================================
struct Cluster
{
    Cluster(int index_)
        : index{ index_ }
    {}
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
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement
// that the hits are passed in time order
class IncrementalDBSCAN
{
public:
    IncrementalDBSCAN(float eps, unsigned int minPts, size_t pool_size=100000)
        : m_eps(eps)
        , m_minPts(minPts)
        , m_pool_begin(0)
        , m_pool_end(0)
    {
        for(size_t i=0; i<pool_size; ++i){
            m_hit_pool.emplace_back(0,0);
        }
    }

    void add_point(float time, float channel, std::vector<Cluster>* completed_clusters=nullptr);
    
    // Add a new hit. The hit time *must* be >= the time of all hits
    // previously added
    void add_hit(Hit* new_hit, std::vector<Cluster>* completed_clusters=nullptr);

    void trim_hits();

    std::vector<Hit*> get_hits() const { return m_hits; }

    std::map<int, Cluster> get_clusters() const { return m_clusters; }

private:
    //======================================================================
    //
    // Starting from `seed_hit`, find all the reachable hits and add them
    // to `cluster`
    void cluster_reachable(Hit* seed_hit, Cluster& cluster);

    float m_eps;
    float m_minPts;
    std::vector<Hit> m_hit_pool;
    size_t m_pool_begin, m_pool_end;
    std::vector<Hit*> m_hits; // All the hits we've seen so far, in time order
    float m_latest_time{ 0 }; // The latest time of a hit in the vector of hits
    std::map<int, Cluster>
        m_clusters; // All of the currently-active (ie, kIncomplete) clusters
};

}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
