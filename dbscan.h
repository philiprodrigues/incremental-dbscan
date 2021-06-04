#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

const int kNoise=-2;
const int kUndefined=-1;

//======================================================================
struct Hit
{
    Hit(float _time, int _chan)
        : time(_time), chan(_chan), cluster(kUndefined), is_core(false), needs_cluster_reachable(false)
        {}

    float time;
    int chan, cluster;
    bool is_core, needs_cluster_reachable;
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

    auto start_it=std::lower_bound(hits.begin(), hits.end(), q.time - eps, time_comp_lower);
    auto end_it=std::lower_bound(hits.begin(), hits.end(), q.time + eps, time_comp_lower);
    // std::cout << hits.size() << " " << std::distance(hits.begin(), start_it) << " " << std::distance(hits.begin(), end_it) << std::endl;
    for(auto hit_it=start_it; hit_it!=end_it; ++hit_it){
        if(manhattanDist(*(*hit_it), q)<=eps){
            ret.push_back(*hit_it);
        }
    }
    return ret;
}

//======================================================================
//
// Modified DBSCAN algorithm that assumes the input is sorted by time
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

struct State
{
    std::vector<Hit*> hits; // All the hits we've seen so far, in time order
    int cluster{-1}; // The index of the current cluster
    float latest_time{0}; // The latest time of a hit in the vector of hits
    int last_processed_index{0}; // The last index in `hits` that we processed
};

//======================================================================
void cluster_reachable(State& state, std::vector<Hit*> seedSet, int cluster, float eps, unsigned int minPts)
{
    // Loop over all neighbours (and the neighbours of core points, and so on)
    while(!seedSet.empty()){
        Hit* q=seedSet.back();
        seedSet.pop_back();
        if(q->time + eps > state.latest_time){
            q->needs_cluster_reachable=true;
            std::cout << "  hit at " << q->time << ", " << q->chan << " is within eps (" << eps << ") of latest time " << state.latest_time << std::endl;
            continue;
        }

        // Change noise to a border point
        if(q->cluster==kNoise){
            std::cout << "  adding previously-noise hit at " << q->time << ", " << q->chan << " to cluster " << state.cluster << std::endl;
            q->cluster=cluster;
        }
        if(q->cluster!=kUndefined) continue;
        std::cout << "  adding previously-undefined hit at " << q->time << ", " << q->chan << " to cluster " << state.cluster << std::endl;
        q->cluster=cluster;
        // Neighbours of q
        std::vector<Hit*> nbrq=neighbours_sorted(state.hits, *q, eps);
        // If q is a core point, add its neighbours to the search list
        if(nbrq.size()>=minPts){
            q->is_core=true;
            seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
            std::cout << "  hit at " << q->time << ", " << q->chan << " has " << nbrq.size() << " neighbours, so is core" << std::endl;
        }
        else{
            std::cout << "  hit at " << q->time << ", " << q->chan << " has " << nbrq.size() << " neighbours, so is *not* core" << std::endl;
        }
    }
}

//======================================================================
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement that the hits are passed in time order
void dbscan_partial_add_one(State& state, Hit* hit, float eps, unsigned int minPts)
{
    if(hit->cluster != kUndefined){
        throw std::runtime_error("Hit cluster is already set");
    }
    state.latest_time = hit->time;

    for(auto h: state.hits){
        if(h->needs_cluster_reachable && h->time+eps < state.latest_time-1){
            std::cout << "Repeche hit at " << h->time << ", " << h->chan << std::endl;
            std::vector<Hit*> seed{h};
            cluster_reachable(state, seed, h->cluster, eps, minPts);
            h->needs_cluster_reachable=false;
        }
    }

    state.hits.push_back(hit);
    
    // Walk process_until_index forward until it reaches the last hit
    // that is a distance eps or more from the latest hit: for this
    // point (and all those before it), all of the hit's
    // eps-neighbours are already available
    int process_until_index=0;
    while(process_until_index<state.hits.size() && state.hits[process_until_index]->time<state.latest_time - eps) ++process_until_index;

    std::cout << state.hits.size() << ". Got hit with time " << hit->time << ". Processing hits from index " << state.last_processed_index << " to " << process_until_index << std::endl;
    for(int i=state.last_processed_index; i<process_until_index; ++i){
        // Do regular dbscan on the items we haven't processed yet. An item found to be a non-core point here will never be a core point
        Hit* p=state.hits[i];
        if(p->cluster!=kUndefined) continue; // We already did this one
        std::cout << "  hit " << i << " is kUndefined" << std::endl;
        std::vector<Hit*> nbr=neighbours_sorted(state.hits, *p, eps);

        if(nbr.size()<minPts){
            // Not enough neighbours to be a core point. Classify as noise (but we might reclassify later)
            std::cout << "  not enough hits to be a core point" << std::endl;
            p->cluster=kNoise;
            continue;
        }

        state.cluster++;
        // Assign this core point to the current cluster
        p->cluster=state.cluster;
        p->is_core=true;
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

        std::cout << "  " << nbr.size() << " neighbours, which is enough to be a core point, so looping over proper neighbours" << std::endl;
        cluster_reachable(state, seedSet, state.cluster, eps, minPts);
    }

    state.last_processed_index=process_until_index;
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
        grAll->SetPoint(grAll->GetN(), hit->time, hit->chan);
        if(grs.find(hit->cluster)==grs.end()){
            TGraph* gr=new TGraph;
            if(hit->cluster==kNoise){
                gr->SetMarkerColor(kGray);
            }
            if(hit->cluster==kUndefined){
                gr->SetMarkerColor(kGray+2);
                gr->SetMarkerStyle(2);
            }
            else{
                gr->SetMarkerColor(colours[(colIndex++)%nColours]);
            }
            gr->SetMarkerStyle(kFullSquare);
            grs[hit->cluster]=gr;
        }
        TGraph* gr=grs[hit->cluster];
        gr->SetPoint(gr->GetN(), hit->time, hit->chan);
    }
    grAll->Draw("ap");
    grAll->GetXaxis()->SetTitle("Time");
    grAll->GetYaxis()->SetTitle("Channel");
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
