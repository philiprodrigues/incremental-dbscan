#pragma once

#include <vector>

namespace dbscan {

class Hit;
//======================================================================
//
// Find all of the eps-neighbours of hit `q` in `hits`
std::vector<Hit*>
neighbours(const std::vector<Hit*>& hits, const Hit& q, float eps);

//======================================================================
//
// The original DBSCAN algorithm, transcribed from Wikipedia. Makes no
// assumptions on the sorting or otherwise of the input hits vector
std::vector<std::vector<Hit*>>
dbscan_orig(std::vector<Hit*>& hits, float eps, unsigned int minPts);

}
