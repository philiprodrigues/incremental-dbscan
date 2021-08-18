#pragma once

#include "dbscan.hpp"
#include "Point.hpp"

#include <vector>

class TCanvas;

namespace dbscan {
class Hit;

// Draw the clusters in the list. Return the TCanvas in
// which they're drawn
TCanvas*
draw_clusters(const std::vector<Cluster>& clusters, const std::vector<Point>& points);

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
