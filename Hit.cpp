#include "Hit.hpp"

#include <algorithm>

namespace dbscan {
//======================================================================
HitSet::HitSet()
{
    hits.reserve(100);
}

//======================================================================
void
HitSet::insert(Hit* h)
{
    // We're typically inserting hits at or near the end, so do a
    // linear scan instead of full binary search. This turns out to be much faster in our case
    auto it=hits.rbegin();
    while(it!=hits.rend() && (*it)->time > h->time) {
        ++it;
    }
    if(it==hits.rend() || *it!=h) hits.insert(it.base(), h);
}

//======================================================================
Hit::Hit(float _time, int _chan)
    : time(_time)
    , chan(_chan)
    , cluster(kUndefined)
    , connectedness(Connectedness::kUndefined)
{}

//======================================================================

// Return true if hit was indeed a neighbour
bool
Hit::add_potential_neighbour(Hit* other, float eps, int minPts)
{
    if (other != this && euclidean_distance(*this, *other) < eps) {
        neighbours.insert(other);
        if(neighbours.size() + 1 >= minPts){
            connectedness=Connectedness::kCore;
        }
        // Neighbourliness is symmetric
        other->neighbours.insert(this);
        if(other->neighbours.size() + 1 >= minPts){
            other->connectedness=Connectedness::kCore;
        }
        return true;
    }
    return false;
}

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
