#ifndef HIT_H
#define HIT_H

#include <vector>
#include <cmath>

const int kNoise=-2;
const int kUndefined=-1;

enum class Connectedness
{
    kUndefined,
    kNoise,
    kCore,
    kEdge
};

enum class Completeness
{
    kIncomplete,
    kComplete,
};

class Hit;

class HitSet
{
public:
  HitSet();

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

  // Return true if hit was indeed a neighbour
  bool add_potential_neighbour(Hit* other, float eps);
    
  float time;
  int chan, cluster;
  Connectedness connectedness;
  Completeness completeness;
  HitSet neighbours;
};


//======================================================================
inline float manhattanDist(const Hit& p, const Hit& q)
{
  return fabs((p.time-q.time))+fabs(p.chan-q.chan);
}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
