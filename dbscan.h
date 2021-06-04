#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <map>
#include <iostream>
#include <algorithm> // For std::lower_bound
#include <set>

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

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

//======================================================================
struct Hit
{
    Hit(float _time, int _chan)
        : time(_time), chan(_chan), cluster(kUndefined), connectedness(Connectedness::kUndefined), completeness(Completeness::kIncomplete)
        {}

    float time;
    int chan, cluster;
    Connectedness connectedness;
    Completeness completeness;
};

//======================================================================
float manhattanDist(const Hit& p, const Hit& q)
{
    return fabs((p.time-q.time))+fabs(p.chan-q.chan);
}

//======================================================================
std::vector<Hit*> neighbours(const std::vector<Hit*>& hits, const Hit& q, float eps)
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
std::vector<std::vector<Hit*> > dbscan_orig(std::vector<Hit*>& hits, float eps, unsigned int minPts)
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
    // for(auto& p: hits){
    //     if(p->cluster==kNoise){
    //         ret.push_back({p});
    //     }
    // }
    return ret;
}

//======================================================================

// Find the neighbours of hit q, assuming that the hits vector is sorted by time
std::vector<Hit*> neighbours_sorted(const std::vector<Hit*>& hits, const Hit& q, float eps)
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
std::vector<std::vector<Hit*> > dbscan_sorted_input(std::vector<Hit*>& hits, float eps, unsigned int minPts)
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

struct Cluster
{
    int index{-1};
    Completeness completeness{Completeness::kIncomplete};
    float latest_time{0};
    std::set<Hit*> hits;

    void add_hit(Hit* h)
    {
        hits.insert(h);
        h->cluster=index;
        latest_time=std::max(latest_time, h->time);
    }
    
};

struct State
{
    std::vector<Hit*> hits; // All the hits we've seen so far, in time order
    float latest_time{0}; // The latest time of a hit in the vector of hits
    int last_processed_index{0}; // The last index in `hits` that we processed
    std::vector<Cluster> clusters;
};

//======================================================================
void cluster_reachable(State& state, std::vector<Hit*> seedSet, Cluster& cluster, float eps, unsigned int minPts)
{
    // std::cout << "        cluster_reachable with seedSet size " << seedSet.size() << " cluster " << cluster.index << " with initial size " << cluster.hits.size() << std::endl;
    // Loop over all neighbours (and the neighbours of core points, and so on)
    while(!seedSet.empty()){
        Hit* q=seedSet.back();
        seedSet.pop_back();

        // Change noise to a border point
        if(q->connectedness==Connectedness::kNoise){
            // // std::cout << "  adding previously-noise hit at " << q->time << ", " << q->chan << " to cluster " << state.cluster << std::endl;
            cluster.add_hit(q);
        }
        
        if(q->cluster!=kUndefined){
            // // std::cout << "Saw hit in cluster " << q->cluster << " when forming cluster " << cluster.index << std::endl;
            continue;
        }
        
        // std::cout << "        adding previously-undefined hit (" << q->time << ", " << q->chan << ") to cluster " << cluster.index << std::endl;
        cluster.add_hit(q);
        // std::cout << "        cluster " << cluster.index << " now has size " << cluster.hits.size() << std::endl;
        // Neighbours of q
        std::vector<Hit*> nbrq=neighbours_sorted(state.hits, *q, eps);
        // If q is a core point, add its neighbours to the search list
        if(nbrq.size()>=minPts){
            q->connectedness=Connectedness::kCore;
            seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
            // std::cout << "        hit (" << q->time << ", " << q->chan << ") has " << nbrq.size() << " neighbours, so is core" << std::endl;
        }
        else{
            // std::cout << "       hit (" << q->time << ", " << q->chan << ") has " << nbrq.size() << " neighbours, so is *not* core" << std::endl;
        }
    }
}

//======================================================================
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement that the hits are passed in time order
void dbscan_partial_add_one(State& state, Hit* hit, float eps, unsigned int minPts)
{
    static int ncall=0;
    // std::cout << ncall << " Adding hit (" << hit->time << ", " << hit->chan << ")" << std::endl;
    
    state.hits.push_back(hit);
    state.latest_time=hit->time;

    for(Cluster& cluster : state.clusters){
        if(cluster.completeness==Completeness::kComplete){
            // std::cout << "    Cluster " << cluster.index << " is complete. Skipping" << std::endl;
            continue;
        }
        // std::cout << "    Looping over " << cluster.hits.size() << " hits in cluster " << cluster.index << std::endl;
        for(Hit* h : cluster.hits){
            bool completed_this_iteration=h->completeness==Completeness::kIncomplete && h->time < state.latest_time - eps;
            if(completed_this_iteration){
                // std::cout << "        Hit (" << hit->time << ", " << hit->chan << ") is now complete" << std::endl;
                h->completeness=Completeness::kComplete;
            }
            
            if(completed_this_iteration || h->connectedness==Connectedness::kCore){
                // std::cout << "        Calling cluster_reachable for hit (" << hit->time << ", " << hit->chan << ")" << std::endl;
                std::vector<Hit*> nbr=neighbours_sorted(state.hits, *h, eps);

                if(nbr.size()>=minPts){
                    h->connectedness=Connectedness::kCore;
                    
                    // Seed set is all the neighbours of p except for p
                    std::vector<Hit*> seedSet;
                    for(auto const& n: nbr){
                        if(n!=h){
                            seedSet.push_back(n);
                        }
                    }
                    cluster_reachable(state, seedSet, cluster, eps, minPts);
                    // std::cout << "    After cluster_reachable, cluster " << cluster.index << " has size " << cluster.hits.size() << std::endl;
                }
            }
        }
        if(cluster.latest_time < state.latest_time - eps){
            // std::cout << "    cluster " << cluster.index << " is complete. cluster latest_time is " << cluster.latest_time << ", global latest_time is " << state.latest_time << std::endl;
            cluster.completeness=Completeness::kComplete;
        }
    }

    if(hit->connectedness==Connectedness::kUndefined){
        std::vector<Hit*> nbr=neighbours_sorted(state.hits, *hit, eps);
        if(nbr.size()<minPts){
            // hit->connectedness=Connectedness::kNoise;
        }
        else{
            hit->connectedness=Connectedness::kCore;
            Cluster new_cluster;
            new_cluster.index=state.clusters.size();
            // std::cout << ncall << " New cluster idx " << new_cluster.index << " with " << nbr.size() << " points" << std::endl;
            hit->cluster=new_cluster.index;
            new_cluster.completeness=Completeness::kIncomplete;
            new_cluster.add_hit(hit);
            
            state.clusters.push_back(new_cluster);
            
            std::vector<Hit*> seedSet;
            for(auto const& n: nbr){
                if(n!=hit){
                    seedSet.push_back(n);
                }
            }
            cluster_reachable(state, seedSet, state.clusters.back(), eps, minPts);
        }
    }

    ++ncall;
}

//======================================================================
void draw_clusters(const std::vector<Hit*>& hits)
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
            else if(hit->cluster==kUndefined){
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
