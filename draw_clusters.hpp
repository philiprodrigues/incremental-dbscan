#pragma once

#include <vector>

class TCanvas;

namespace dbscan {
class Hit;

// Draw the hits in `hits` coloured by cluster. Return the TCanvas in
// which they're drawn
TCanvas*
draw_clusters(const std::vector<Hit*>& hits);

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
