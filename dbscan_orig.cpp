#include "dbscan_orig.hpp"

#include "Hit.hpp"
#include <cassert>

std::vector<Hit*>
neighbours(const std::vector<Hit*>& hits, const Hit& q, float eps)
{
    std::vector<Hit*> ret;
    for (auto const& hit : hits) {
        if (manhattanDist(*hit, q) <= eps) {
            ret.push_back(hit);
        }
    }
    return ret;
}
std::vector<std::vector<Hit*>>
dbscan_orig(std::vector<Hit*>& hits, float eps, unsigned int minPts)
{
    std::vector<std::vector<Hit*>> ret;

    int cluster = -1; // The index of the current cluster

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
        cluster++;
        ret.resize(cluster + 1);

        // Assign this core point to the current cluster
        p->cluster = cluster;
        ret[cluster].push_back(p);
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
                q->cluster = cluster;
            if (q->cluster != kUndefined)
                continue;
            q->cluster = cluster;
            ret[cluster].push_back(q);
            // Neighbours of q
            std::vector<Hit*> nbrq = neighbours(hits, *q, eps);
            // If q is a core point, add its neighbours to the search list
            if (nbrq.size() >= minPts)
                seedSet.insert(seedSet.end(), nbrq.begin(), nbrq.end());
        }
    }
    return ret;
}
