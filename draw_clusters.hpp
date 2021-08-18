#pragma once

#include "dbscan.hpp"
#include "Point.hpp"

#include "folly/FBVector.h"

class TCanvas;

namespace dbscan {
class Hit;

// Draw the clusters in the list. Return the TCanvas in
// which they're drawn
TCanvas*
draw_clusters(const folly::fbvector<Cluster>& clusters, const folly::fbvector<Point>& points);

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
