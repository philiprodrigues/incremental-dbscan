#include "draw_clusters.hpp"

#include "Hit.hpp"

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

#include <vector>
#include <map>

TCanvas*
draw_clusters(const std::vector<Hit*>& hits)
{
    if (hits.empty())
        return nullptr;
    TCanvas* c = new TCanvas;
    const int nColours = 6;
    int colours[nColours] = { kRed,         kBlue,       kGreen + 2,
                              kMagenta + 2, kOrange + 2, kCyan + 2 };
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
                gr->SetMarkerColor(colours[(colIndex++) % nColours]);
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

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
