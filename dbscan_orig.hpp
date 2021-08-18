#pragma once

#include "dbscan.hpp"

#include "folly/FBVector.h"

namespace dbscan {

class Hit;
//======================================================================
//
// Find all of the eps-neighbours of hit `q` in `hits`
folly::fbvector<Hit*>
neighbours(const folly::fbvector<Hit*>& hits, const Hit& q, float eps);

//======================================================================
//
// The original DBSCAN algorithm, transcribed from Wikipedia. Makes no
// assumptions on the sorting or otherwise of the input hits vector
folly::fbvector<Cluster>
dbscan_orig(folly::fbvector<Hit*>& hits, float eps, unsigned int minPts);

}
