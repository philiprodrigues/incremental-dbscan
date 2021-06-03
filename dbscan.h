#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound

#include "TGraph.h"
#include "TCanvas.h"

const int kNoise=-2;
const int kUndefined=-1;

//======================================================================
struct Hit
{
    Hit(float _time, int _chan)
        : time(_time), chan(_chan), cluster(kUndefined)
        {}

    float time;
    int chan, cluster;
};

//======================================================================
float manhattanDist(const Hit& p, const Hit& q)
{
    return fabs((p.time-q.time))+fabs(p.chan-q.chan);
}

//======================================================================
std::vector<Hit*> neighbours(const vector<Hit*>& hits, const Hit& q, float eps)
{
    std::vector<Hit*> ret;
    for(auto const& hit: hits){
        if(manhattanDist(*hit, q)<=eps){
            ret.push_back(hit);
        }
    }
    return ret;
}

//======================================================================
//
// The original DBSCAN algorithm, transcribed from Wikipedia. Makes no
// assumptions on the sorting or otherwise of the input hits vector
std::vector<std::vector<Hit*> > dbscan_orig(vector<Hit*>& hits, float eps, unsigned int minPts)
{
    std::vector<std::vector<Hit*> > ret;

    int cluster=-1; // The index of the current cluster
    
    for(auto& p: hits){
        if(p->cluster!=kUndefined) continue; // We already did this one
      
        std::vector<Hit*> nbr=neighbours(hits, *p, eps);

        if(nbr.size()<minPts){
            // Not enough neighbours to be a core point. Classify as noise (but we might reclassify later)
            p->cluster=kNoise;
            continue;
        }
        cluster++;
        ret.resize(cluster+1);

        // Assign this core point to the current cluster
        p->cluster=cluster;
        ret[cluster].push_back(p);
        // Seed set is all the neighbours of p except for p
        std::vector<Hit*> seedSet;
        for(auto const& n: nbr){
            if(n!=p){
                seedSet.push_back(n);
            }
        }
        if(seedSet.size()!=nbr.size()-1){
            std::cout << "seedSet.size()=" << seedSet.size() << " but nbr.size()=" << nbr.size() << std::endl;
        }
        assert(seedSet.size()==nbr.size()-1);

        // Loop over all neighbours (and the neighbours of core points, and so on)
        while(!seedSet.empty()){
            Hit* q=seedSet.back();
            seedSet.pop_back();
            // Change noise to a border point
            if(q->cluster==kNoise) q->cluster=cluster;
            if(q->cluster!=kUndefined) continue;
            q->cluster=cluster;
            ret[cluster].push_back(q);
            // Neighbours of q
            std::vector<Hit*> nbrq=neighbours(hits, *q, eps);
            // If q is a core point, add its neighbours to the search list
            if(nbrq.size()>=minPts) seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
        }
    }

    // Put the noise hits in their own clusters
    for(auto& p: hits){
        if(p->cluster==kNoise){
            ret.push_back({p});
        }
    }
    return ret;
}

//======================================================================

// Find the neighbours of hit q, assuming that the hits vector is sorted by time
std::vector<Hit*> neighbours_sorted(const vector<Hit*>& hits, const Hit& q, float eps)
{
    std::vector<Hit*> ret;
    auto time_comp_lower=[](const Hit* hit, const float t) { return hit->time < t; };
    auto time_comp_upper=[](const float t, const Hit* hit) { return hit->time < t; };
    auto start_it=std::lower_bound(hits.begin(), hits.end(), q.time - eps, time_comp_lower);
    auto end_it=std::upper_bound(hits.begin(), hits.end(), q.time + eps, time_comp_upper);
    for(auto hit_it=start_it; hit_it!=end_it; ++hit_it){
        if(manhattanDist(*(*hit_it), q)<=eps){
            ret.push_back(*hit_it);
        }
    }
    return ret;
}

//======================================================================
//
// Modified DBSCAN algorithm that assumes the input in sorted by time
std::vector<std::vector<Hit*> > dbscan_sorted_input(vector<Hit*>& hits, float eps, unsigned int minPts)
{
    std::vector<std::vector<Hit*> > ret;

    int cluster=-1; // The index of the current cluster
    
    for(auto& p: hits){
        if(p->cluster!=kUndefined) continue; // We already did this one
      
        std::vector<Hit*> nbr=neighbours_sorted(hits, *p, eps);

        if(nbr.size()<minPts){
            // Not enough neighbours to be a core point. Classify as noise (but we might reclassify later)
            p->cluster=kNoise;
            continue;
        }
        cluster++;
        ret.resize(cluster+1);

        // Assign this core point to the current cluster
        p->cluster=cluster;
        ret[cluster].push_back(p);
        // Seed set is all the neighbours of p except for p
        std::vector<Hit*> seedSet;
        for(auto const& n: nbr){
            if(n!=p){
                seedSet.push_back(n);
            }
        }
        if(seedSet.size()!=nbr.size()-1){
            std::cout << "seedSet.size()=" << seedSet.size() << " but nbr.size()=" << nbr.size() << std::endl;
        }
        assert(seedSet.size()==nbr.size()-1);

        // Loop over all neighbours (and the neighbours of core points, and so on)
        while(!seedSet.empty()){
            Hit* q=seedSet.back();
            seedSet.pop_back();
            // Change noise to a border point
            if(q->cluster==kNoise) q->cluster=cluster;
            if(q->cluster!=kUndefined) continue;
            q->cluster=cluster;
            ret[cluster].push_back(q);
            // Neighbours of q
            std::vector<Hit*> nbrq=neighbours_sorted(hits, *q, eps);
            // If q is a core point, add its neighbours to the search list
            if(nbrq.size()>=minPts) seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
        }
    }

    // Put the noise hits in their own clusters
    for(auto& p: hits){
        if(p->cluster==kNoise){
            ret.push_back({p});
        }
    }
    return ret;
}

//======================================================================
void draw_clusters(const vector<Hit*>& hits)
{
    if(hits.empty()) return;
    new TCanvas;
    const int nColours=6;
    int colours[nColours]={kRed, kBlue, kGreen+2, kMagenta+2, kOrange+2, kCyan+2};
    TGraph* grAll=new TGraph;
    std::map<int, TGraph*> grs;
    int colIndex=0;

    for(auto const& hit: hits){
        grAll->SetPoint(grAll->GetN(), hit->chan, hit->time);
        if(grs.find(hit->cluster)==grs.end()){
            TGraph* gr=new TGraph;
            if(hit->cluster==kNoise){
                gr->SetMarkerColor(kGray);
            }
            else{
                gr->SetMarkerColor(colours[(colIndex++)%nColours]);
            }
            gr->SetMarkerStyle(kFullSquare);
            grs[hit->cluster]=gr;
        }
        TGraph* gr=grs[hit->cluster];
        gr->SetPoint(gr->GetN(), hit->chan, hit->time);
    }
    grAll->Draw("ap");
    for(auto const& gr: grs){
        gr.second->Draw("p");
    }
}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
