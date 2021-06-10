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
    auto time_comp_lower = [](const Hit* hit, const float t) {
        return hit->time < t;
    };
    auto it =
        std::lower_bound(hits.begin(), hits.end(), h->time, time_comp_lower);
    if (it == hits.end() || *it != h) {
        hits.insert(it, h);
    }
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
