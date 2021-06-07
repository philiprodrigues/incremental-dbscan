#pragma once

#include <vector>
#include <cmath>
#include <list>

//======================================================================

// Special "cluster numbers" for hits that are not (yet) in a cluster
const int kNoise = -2;
const int kUndefined = -1;

//======================================================================

// Hit classifications in the DBSCAN scheme
enum class Connectedness
{
    // clang-format off
    kUndefined,
    kNoise,     // Fewer than minPts neighbours, not in a cluster
    kCore,      // minPts neighbours or more
    kEdge       // Fewer than minPts neighbours, part of a cluster
    // clang-format on
};

//======================================================================

// As new hits arrive, they push forward the "current" time, and
// eventually a given hit or cluster will know that it cannot be
// modified any further. A hit becomes kComplete when its time is so
// far behind the current time that new hits cannot be neighbours of
// it. A cluster becomes kComplete when its latest hit is complete
enum class Completeness
{
    kIncomplete,
    kComplete,
};

class Hit;

//======================================================================

// An array of unique hits, sorted by time. The actual container
// implementation is a std::vector, which seems to be faster than a
// std::set (needs rechecking)
class HitSet
{
public:
    HitSet();

    // Insert a hit in the set, if not already present. Keeps the
    // array sorted by time
    void insert(Hit* h);

    std::vector<Hit*>::iterator begin() { return hits.begin(); }
    std::vector<Hit*>::iterator end() { return hits.end(); }

    void clear() { hits.clear(); }

    size_t size() { return hits.size(); }

    std::vector<Hit*> hits;
};

//======================================================================
struct Hit
{
    Hit(float _time, int _chan);

    // Add hit `other` to this hit's list of neighbours if they are
    // closer than `eps`. Return true if so
    bool add_potential_neighbour(Hit* other, float eps);

    float time;
    int chan, cluster;
    Connectedness connectedness;
    Completeness completeness;
    HitSet neighbours;
};

//======================================================================
inline float
manhattanDist(const Hit& p, const Hit& q)
{
    return fabs((p.time - q.time)) + fabs(p.chan - q.chan);
}

//======================================================================
inline bool
time_comp_lower(const Hit* hit, const float t)
{
    return hit->time < t;
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
