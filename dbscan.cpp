#include "dbscan.hpp"

#include <cassert>

int
neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps)
{
    int n = 0;
    // Loop over the hits starting from the latest hit, since we will
    // ~always be adding a hit at recent times
    for (auto hit_it = hits.rbegin(); hit_it != hits.rend(); ++hit_it) {
        if ((*hit_it)->time > q.time + eps)
            continue;
        if ((*hit_it)->time < q.time - eps)
            break;

        if (q.add_potential_neighbour(*hit_it, eps))
            ++n;
    }
    return n;
}

bool
Cluster::maybe_add_new_hit(Hit* new_hit, float eps, int minPts)
{
    bool do_add = false;

    auto begin_it = std::lower_bound(
        hits.begin(), hits.end(), new_hit->time - eps, time_comp_lower);

    for (auto it = begin_it; it != hits.end(); ++it) {
        Hit* h = *it;
        if (h->add_potential_neighbour(new_hit, eps)) {
            do_add = true;
            if (h->neighbours.size() + 1 >= minPts) {
                h->connectedness = Connectedness::kCore;
            } else {
                h->connectedness = Connectedness::kEdge;
            }
        }

        // Mark any hits that are now completed
        bool completed_this_iteration =
            h->completeness == Completeness::kIncomplete &&
            h->time < new_hit->time - eps;
        if (completed_this_iteration) {
            // std::cout << "        Hit (" << h->time << ", " << h->chan << ")
            // is now complete" << std::endl;
            h->completeness = Completeness::kComplete;
        }
    } // end loop over hits in cluster

    if (do_add) {
        add_hit(new_hit);
    }

    return do_add;
}
void
Cluster::add_hit(Hit* h)
{
    hits.insert(h);
    h->cluster = index;
    latest_time = std::max(latest_time, h->time);
    if (h->connectedness == Connectedness::kCore &&
        (!latest_core_point || h->time > latest_core_point->time))
        latest_core_point = h;
}
void
Cluster::steal_hits(Cluster& other)
{
    for (auto h : other.hits) {
        assert(h);
        add_hit(h);
    }
    other.hits.clear();
    other.completeness = Completeness::kComplete;
}
void
cluster_reachable(State& state,
                  Hit* seed_hit,
                  Cluster& cluster,
                  float eps,
                  unsigned int minPts)
{
    // std::cout << "        cluster_reachable with
    // seed_hit->neighbours.size()=" << seed_hit->neighbours.size() << " cluster
    // " << cluster.index << " with initial size " << cluster.hits.size() <<
    // std::endl;

    // Loop over all neighbours (and the neighbours of core points, and so on)
    std::vector<Hit*> seedSet(seed_hit->neighbours.begin(),
                              seed_hit->neighbours.end());

    while (!seedSet.empty()) {
        Hit* q = seedSet.back();
        seedSet.pop_back();

        // Change noise to a border point
        if (q->connectedness == Connectedness::kNoise) {
            // std::cout << "        adding previously-noise hit at " << q->time
            // << ", " << q->chan << " to cluster " << cluster.index <<
            // std::endl;
            cluster.add_hit(q);
        }

        if (q->cluster != kUndefined) {
            // // std::cout << "        Saw hit in cluster " << q->cluster << "
            // when forming cluster " << cluster.index << std::endl;
            continue;
        }

        // std::cout << "        adding previously-undefined hit (" << q->time
        // << ", " << q->chan << ") to cluster " << cluster.index << std::endl;
        cluster.add_hit(q);
        // std::cout << "        cluster " << cluster.index << " now has size "
        // << cluster.hits.size() << std::endl;

        // If q is a core point, add its neighbours to the search list
        if (q->neighbours.size() + 1 >= minPts) {
            q->connectedness = Connectedness::kCore;
            seedSet.insert(
                seedSet.end(), q->neighbours.begin(), q->neighbours.end());
            // std::cout << "        hit (" << q->time << ", " << q->chan << ")
            // has " << q->neighbours.size() << " neighbours, so is core" <<
            // std::endl;
        } else {
            // std::cout << "       hit (" << q->time << ", " << q->chan << ")
            // has " << q->neighbours.size() << " neighbours, so is *not* core"
            // << std::endl;
        }
    }
}
void
dbscan_partial_add_one(State& state,
                       Hit* new_hit,
                       float eps,
                       unsigned int minPts)
{
    static int ncall = 0;
    static int next_cluster_index = 0;

    // std::cout << ncall << " Adding hit (" << new_hit->time << ", " <<
    // new_hit->chan << ")" << std::endl;

    state.hits.push_back(new_hit);
    state.latest_time = new_hit->time;
    neighbours_sorted(state.hits, *new_hit, eps);

    // All the clusters that this hit neighboured. If there are
    // multiple clusters neighbouring this hit, we'll merge them at
    // the end
    using cluster_iterator = std::list<Cluster>::iterator;
    std::vector<cluster_iterator> clusters_to_merge;

    auto clust_it = state.clusters.begin();

    while (clust_it != state.clusters.end()) {
        Cluster& cluster = *clust_it;
        if (cluster.completeness == Completeness::kComplete) {
            // if(ncall>5500) std::cout << "Erasing cluster " << cluster.index
            // << std::endl;
            clust_it = state.clusters.erase(clust_it);
            continue;
        }

        if (cluster.maybe_add_new_hit(new_hit, eps, minPts)) {
            clusters_to_merge.push_back(clust_it);
        }

        cluster_reachable(
            state, cluster.latest_core_point, cluster, eps, minPts);

        if (cluster.latest_time < state.latest_time - eps) {
            cluster.completeness = Completeness::kComplete;
        }
        ++clust_it;
    } // end loop over clusters

    // Merge any clusters that need merging
    if (clusters_to_merge.size() >= 2) {
        cluster_iterator into = clusters_to_merge.front();
        for (size_t i = 1; i < clusters_to_merge.size(); ++i) {
            into->steal_hits(*clusters_to_merge[i]);
        }
    }

    if (new_hit->cluster == kUndefined) {
        // If we get here, the new hit wasn't a neighbour of any hit
        // in any cluster. If this hit has enough neighbours, it
        // should seed a new cluster

        neighbours_sorted(state.hits, *new_hit, eps);

        if (new_hit->neighbours.size() + 1 >= minPts) {
            new_hit->connectedness = Connectedness::kCore;
            Cluster& new_cluster = state.clusters.emplace_back();
            new_cluster.index = next_cluster_index++;
            // std::cout << ncall << " New cluster idx " << new_cluster.index <<
            // std::endl;
            new_cluster.completeness = Completeness::kIncomplete;
            new_cluster.add_hit(new_hit);
            for (auto& neighbour : new_hit->neighbours) {
                new_cluster.add_hit(neighbour);
            }
            cluster_reachable(state, new_hit, new_cluster, eps, minPts);
            // std::cout << ncall << " After initial cluster_reachable, cluster
            // " << new_cluster.index << " has " << new_cluster.hits.size() << "
            // hits" << std::endl;
        }
    }

    ++ncall;
}
