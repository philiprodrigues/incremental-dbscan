#include "Hit.h"

#include <algorithm>

HitSet::HitSet()
{
  hits.reserve(100);
}

void HitSet::insert(Hit* h)
{
  auto time_comp_lower=[](const Hit* hit, const float t) { return hit->time < t; };
  auto it=std::lower_bound(hits.begin(), hits.end(), h->time, time_comp_lower);
  if(it==hits.end() || *it!=h) hits.insert(it, h);

  // bool already=false;
  // for(auto& hit: hits){
  //   if(hit==h) already=true;
  // }
  // if(!already) hits.push_back(h);

  // if(std::all_of(hits.begin(), hits.end(), [h](Hit* hit){ return hit != h; })) hits.push_back(h);
}

Hit::Hit(float _time, int _chan)
  : time(_time), chan(_chan), cluster(kUndefined), connectedness(Connectedness::kUndefined), completeness(Completeness::kIncomplete)
{
}

// Return true if hit was indeed a neighbour
bool Hit::add_potential_neighbour(Hit* other, float eps)
{
  if(other!=this && manhattanDist(*this, *other)<eps){
    neighbours.insert(other);
    // Neighbourliness is symmetric
    other->neighbours.insert(this);
    return true;
  }
  return false;
}
