#include "dbscan_orig.hpp"
#include "dbscan.hpp"

#include "Hit.hpp"
#include <cassert>

namespace dbscan {
std::vector<Hit*>
neighbours(const std::vector<Hit*>& hits, const Hit& q, float eps)
{
    std::vector<Hit*> ret;
    for (auto const& hit : hits) {
        if (euclidean_distance(*hit, q) < eps) {
            ret.push_back(hit);
        }
    }
    return ret;
}
    
std::vector<Cluster>
dbscan_orig(std::vector<Hit*>& hits, float eps, unsigned int minPts)
{
    std::vector<Cluster> ret;

    int clusterIndex = -1; // The index of the current cluster

    for (auto& p : hits) {
        if (p->cluster != kUndefined)
            continue; // We already did this one

        std::vector<Hit*> nbr = neighbours(hits, *p, eps);

        if (nbr.size() < minPts) {
            // Not enough neighbours to be a core point. Classify as noise (but
            // we might reclassify later)
            p->cluster = kNoise;
            continue;
        }
        clusterIndex++;
        ret.emplace_back(clusterIndex);
        Cluster& current_cluster = ret.back();

        // Assign this core point to the current cluster
        p->cluster = clusterIndex;
        current_cluster.add_hit(p);
        // Seed set is all the neighbours of p except for p
        std::vector<Hit*> seedSet;
        for (auto const& n : nbr) {
            if (n != p) {
                seedSet.push_back(n);
            }
        }
        if (seedSet.size() != nbr.size() - 1) {
            // std::cout << "seedSet.size()=" << seedSet.size() << " but
            // nbr.size()=" << nbr.size() << std::endl;
        }
        assert(seedSet.size() == nbr.size() - 1);

        // Loop over all neighbours (and the neighbours of core points, and so
        // on)
        while (!seedSet.empty()) {
            Hit* q = seedSet.back();
            seedSet.pop_back();
            // Change noise to a border point
            if (q->cluster == kNoise)
                q->cluster = clusterIndex;
            if (q->cluster != kUndefined)
                continue;
            q->cluster = clusterIndex;
            current_cluster.add_hit(q);
            // Neighbours of q
            std::vector<Hit*> nbrq = neighbours(hits, *q, eps);
            // If q is a core point, add its neighbours to the search list
            if (nbrq.size() >= minPts)
                seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
        }
    }
    return ret;
}

}
