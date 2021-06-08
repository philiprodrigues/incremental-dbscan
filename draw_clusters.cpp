#include "draw_clusters.hpp"

#include "Hit.hpp"

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TColor.h"

#include <vector>
#include <map>

namespace dbscan {
TCanvas*
draw_clusters(const std::vector<Hit*>& hits)
{
    if (hits.empty()){
        return nullptr;
    }
    
    TCanvas* c = new TCanvas;
    std::vector<int> colours{
        TColor::GetColor("#377eb8"),
        TColor::GetColor("#ff7f00"),
        TColor::GetColor("#4daf4a"),
        TColor::GetColor("#f781bf"),
        TColor::GetColor("#a65628"),
        TColor::GetColor("#984ea3"),
        TColor::GetColor("#999999"),
        TColor::GetColor("#e41a1c"),
        TColor::GetColor("#dede00")
    };

    TGraph* grAll = new TGraph;
    std::map<int, TGraph*> grs;
    int colIndex = 0;

    for (auto const& hit : hits) {
        grAll->SetPoint(grAll->GetN(), hit->time, hit->chan);
        if (grs.find(hit->cluster) == grs.end()) {
            TGraph* gr = new TGraph;
            if (hit->cluster == kNoise) {
                gr->SetMarkerColor(kGray);
            } else if (hit->cluster == kUndefined) {
                gr->SetMarkerColor(kGray);
                gr->SetMarkerStyle(2);
            } else {
                gr->SetMarkerColor(colours[(colIndex++) % colours.size()]);
            }
            gr->SetMarkerStyle(kFullSquare);
            grs[hit->cluster] = gr;
        }
        TGraph* gr = grs[hit->cluster];
        gr->SetPoint(gr->GetN(), hit->time, hit->chan);
    }
    grAll->Draw("ap");
    grAll->GetXaxis()->SetTitle("Time");
    grAll->GetYaxis()->SetTitle("Channel");
    for (auto const& gr : grs) {
        gr.second->Draw("p");
    }

    return c;
}

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
