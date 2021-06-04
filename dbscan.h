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

struct Hit;
float manhattanDist(const Hit& p, const Hit& q);

//======================================================================
struct Hit
{
    Hit(float _time, int _chan)
        : time(_time), chan(_chan), cluster(kUndefined), connectedness(Connectedness::kUndefined), completeness(Completeness::kIncomplete)
        {
        }

    // Return true if hit was indeed a neighbour
    bool add_potential_neighbour(Hit* other, float eps)
    {
        if(other!=this && manhattanDist(*this, *other)<eps){
            neighbours.insert(other);
            // Neighbourliness is symmetric
            other->neighbours.insert(this);
            return true;
        }
        return false;
    }
    
    float time;
    int chan, cluster;
    Connectedness connectedness;
    Completeness completeness;
    std::set<Hit*> neighbours;
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
            // std::cout << "seedSet.size()=" << seedSet.size() << " but nbr.size()=" << nbr.size() << std::endl;
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
int neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps)
{
    int n=0;
    for(auto hit_it=hits.rbegin(); hit_it!=hits.rend(); ++hit_it){
        if((*hit_it)->time > q.time + eps) continue;
        if((*hit_it)->time < q.time - eps) break;
        
        if(q.add_potential_neighbour(*hit_it, eps)) ++n;
    }
    return n;
}

//======================================================================
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

    void steal_hits(Cluster& other)
    {
        for(auto h: other.hits){
            add_hit(h);
        }
        other.hits.clear();
        other.completeness=Completeness::kComplete;
    }
};

//======================================================================
struct State
{
    std::vector<Hit*> hits; // All the hits we've seen so far, in time order
    std::vector<Hit*> old_hits; // Hits that are part of completed clusters
    float latest_time{0}; // The latest time of a hit in the vector of hits
    int last_processed_index{0}; // The last index in `hits` that we processed
    std::vector<Cluster> clusters;
};

//======================================================================
void cluster_reachable(State& state, Hit* seed_hit, Cluster& cluster, float eps, unsigned int minPts)
{
    // std::cout << "        cluster_reachable with seed_hit->neighbours.size()=" << seed_hit->neighbours.size() << " cluster " << cluster.index << " with initial size " << cluster.hits.size() << std::endl;

    // Loop over all neighbours (and the neighbours of core points, and so on)
    std::vector<Hit*> seedSet(seed_hit->neighbours.begin(), seed_hit->neighbours.end());
    
    while(!seedSet.empty()){
        Hit* q=seedSet.back();
        seedSet.pop_back();

        // Change noise to a border point
        if(q->connectedness==Connectedness::kNoise){
            // std::cout << "        adding previously-noise hit at " << q->time << ", " << q->chan << " to cluster " << cluster.index << std::endl;
            cluster.add_hit(q);
        }
        
        if(q->cluster!=kUndefined){
            // // std::cout << "        Saw hit in cluster " << q->cluster << " when forming cluster " << cluster.index << std::endl;
            continue;
        }
        
        // std::cout << "        adding previously-undefined hit (" << q->time << ", " << q->chan << ") to cluster " << cluster.index << std::endl;
        cluster.add_hit(q);
        // std::cout << "        cluster " << cluster.index << " now has size " << cluster.hits.size() << std::endl;
    
        // If q is a core point, add its neighbours to the search list
        if(q->neighbours.size() + 1 >= minPts){
            q->connectedness=Connectedness::kCore;
            seedSet.insert(seedSet.end(), q->neighbours.begin(), q->neighbours.end());
            // std::cout << "        hit (" << q->time << ", " << q->chan << ") has " << q->neighbours.size() << " neighbours, so is core" << std::endl;
        }
        else{
            // std::cout << "       hit (" << q->time << ", " << q->chan << ") has " << q->neighbours.size() << " neighbours, so is *not* core" << std::endl;
        }
    }
}

//======================================================================
//
// Modified DBSCAN algorithm that takes one hit at a time, with the requirement that the hits are passed in time order
void dbscan_partial_add_one(State& state, Hit* new_hit, float eps, unsigned int minPts)
{
    static int ncall=0;
    
    // std::cout << ncall << " Adding hit (" << new_hit->time << ", " << new_hit->chan << ")" << std::endl;
    
    state.hits.push_back(new_hit);
    state.latest_time=new_hit->time;
    neighbours_sorted(state.hits, *new_hit, eps);

    std::vector<std::pair<int, int>> clusters_to_merge;
    
    for(Cluster& cluster : state.clusters){
        if(cluster.completeness==Completeness::kComplete){
            // std::cout << "    Cluster " << cluster.index << " is complete. Skipping" << std::endl;
            continue;
        }
        // std::cout << "    Looping over " << cluster.hits.size() << " hits in cluster " << cluster.index << std::endl;

        Hit* latest_core_point=nullptr;
        for(Hit* h : cluster.hits){

            // Try adding the new hit to the neighbour list of any incomplete hits
            if(h->completeness==Completeness::kIncomplete){
                if(h->add_potential_neighbour(new_hit, eps)){
                    if(new_hit->cluster!=kUndefined && new_hit->cluster!=cluster.index){
                        // std::cout << "    new_hit already in cluster " << new_hit->cluster << " also neighbours cluster " << cluster.index << ". Should merge" << std::endl;
                        clusters_to_merge.push_back(std::make_pair(new_hit->cluster, cluster.index));
                    }
                    cluster.add_hit(new_hit);
                    if(h->neighbours.size() + 1 >= minPts){
                        h->connectedness = Connectedness::kCore;
                    }
                    else{
                        h->connectedness = Connectedness::kEdge;
                    }
                }
            }
            
            if(h->connectedness == Connectedness::kCore){
                if(!latest_core_point || h->time > latest_core_point->time){
                    latest_core_point=h;
                }
            }

            // Mark any hits that are now completed
            bool completed_this_iteration=h->completeness==Completeness::kIncomplete && h->time < state.latest_time - eps;
            if(completed_this_iteration){
                // std::cout << "        Hit (" << h->time << ", " << h->chan << ") is now complete" << std::endl;
                h->completeness=Completeness::kComplete;
            }
            
            //     // std::cout << "        Calling cluster_reachable for hit (" << hit->time << ", " << hit->chan << ")" << std::endl;
            //     neighbours_sorted(state.hits, *h, eps, nbr);

            //     if(nbr.size()>=minPts){
            //         h->connectedness=Connectedness::kCore;
                    
            //         // Seed set is all the neighbours of p except for p
            //         std::vector<Hit*> seedSet;
            //         for(auto const& n: nbr){
            //             if(n!=h){
            //                 seedSet.push_back(n);
            //             }
            //         }
            //         cluster_reachable(state, seedSet, cluster, eps, minPts, nbr);
            //         // std::cout << "    After cluster_reachable, cluster " << cluster.index << " has size " << cluster.hits.size() << std::endl;
            //     }
            // }
        } // end loop over hits in cluster

        cluster_reachable(state, latest_core_point, cluster, eps, minPts);
        
        if(cluster.latest_time < state.latest_time - eps){
            // std::cout << "    cluster " << cluster.index << " is complete. cluster latest_time is " << cluster.latest_time << ", global latest_time is " << state.latest_time << std::endl;
            cluster.completeness=Completeness::kComplete;
        }
    }

    for(auto cs : clusters_to_merge){
        // std::cout << "   Transferring hits to cluster " << cs.first << " from " << cs.second << std::endl;
        state.clusters[cs.first].steal_hits(state.clusters[cs.second]);
    }
    
    if(new_hit->cluster==kUndefined){
        // If we get here, the new hit wasn't a neighbour of any hit
        // in any cluster. If this hit has enough neighbours, it
        // should seed a new cluster
        
        neighbours_sorted(state.hits, *new_hit, eps);

        if(new_hit->neighbours.size() + 1 >= minPts){
            new_hit->connectedness=Connectedness::kCore;
            int index=state.clusters.size();
            Cluster& new_cluster=state.clusters.emplace_back();
            new_cluster.index=index;
            // std::cout << ncall << " New cluster idx " << new_cluster.index << std::endl;
            new_cluster.completeness=Completeness::kIncomplete;
            new_cluster.add_hit(new_hit);
            for(auto& neighbour: new_hit->neighbours){
                new_cluster.add_hit(neighbour);
            }
            cluster_reachable(state, new_hit, new_cluster, eps, minPts);
            // std::cout << ncall << " After initial cluster_reachable, cluster " << new_cluster.index << " has " << new_cluster.hits.size() << " hits" << std::endl;
        }
    }

    // std::cout << "  dbscan_partial_add_one ends with new_hit->cluster=" << new_hit->cluster << std::endl;
    // // Every so often, remove the hits from completed clusters from consideration, and too-old noise hits
    // if(ncall%1024==0){
    //     size_t norig=state.hits.size();
    //     for(Cluster& cluster : state.clusters){
    //         if(cluster.completeness==Completeness::kComplete){
    //             for(Hit* h : cluster.hits){
    //                 state.old_hits.push_back(h);
    //                 auto time_comp_lower=[](const Hit* hit, const float t) { return hit->time < t; };
    //                 auto start_it=std::lower_bound(state.hits.begin(), state.hits.end(), h->time, time_comp_lower);
    //                 do{
    //                     if(*start_it==h) state.hits.erase(start_it);
    //                 } while((*start_it)->time < h->time);
    //             }
    //         }
    //     }

    //     for (auto it = state.hits.begin(); it != state.hits.end() && (*it)->time < state.latest_time - 100*eps; ) {
    //         if ((*it)->connectedness==Connectedness::kNoise || (*it)->connectedness==Connectedness::kUndefined) {
    //             state.old_hits.push_back(*it);
    //             it = state.hits.erase(it);
    //         } else {
    //             ++it;
    //         }
    //     }

    //     size_t nafter=state.hits.size();
    //     // std::cout << "Cleaned up " << (norig-nafter) << " hits out of " << norig << std::endl;
    // }
    
    ++ncall;
}

//======================================================================
TCanvas* draw_clusters(const std::vector<Hit*>& hits)
{
    if(hits.empty()) return nullptr;
    TCanvas* c=new TCanvas;
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
                gr->SetMarkerColor(kGray);
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

    return c;
}

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
