#include "dbscan.hpp"

#include <cassert>

namespace dbscan {

//======================================================================
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

//======================================================================
bool
Cluster::maybe_add_new_hit(Hit* new_hit, float eps, int minPts)
{
    // Should we add this hit?
    bool do_add = false;

    // Hits earlier than new_hit time minus `eps` can't possibly be
    // neighbours, so start the search there in the sorted list of hits in this cluster
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
            h->completeness = Completeness::kComplete;
        }
    } // end loop over hits in cluster

    if (do_add) {
        add_hit(new_hit);
    }

    return do_add;
}

//======================================================================
void
Cluster::add_hit(Hit* h)
{
    hits.insert(h);
    h->cluster = index;
    latest_time = std::max(latest_time, h->time);
    if (h->connectedness == Connectedness::kCore &&
        (!latest_core_point || h->time > latest_core_point->time)){
        latest_core_point = h;
    }
}

//======================================================================
void
Cluster::steal_hits(Cluster& other)
{
    // TODO: it might be faster to do some sort of explicit "merge" of the hits, eg:
    //
    // this->hits.insert(other hits); // Inserts at end
    // std::inplace_merge(...)
    //
    // This might save some reallocations of the vector
    for (auto h : other.hits) {
        assert(h);
        add_hit(h);
    }
    other.hits.clear();
    other.completeness = Completeness::kComplete;
}

//======================================================================
void
IncrementalDBSCAN::cluster_reachable(Hit* seed_hit,
                                     Cluster& cluster)
{
    // Loop over all neighbours (and the neighbours of core points, and so on)
    std::vector<Hit*> seedSet(seed_hit->neighbours.begin(),
                              seed_hit->neighbours.end());

    while (!seedSet.empty()) {
        Hit* q = seedSet.back();
        seedSet.pop_back();

        // Change noise to a border point
        if (q->connectedness == Connectedness::kNoise) {
            cluster.add_hit(q);
        }

        if (q->cluster != kUndefined) {
            continue;
        }

        cluster.add_hit(q);

        // If q is a core point, add its neighbours to the search list
        if (q->neighbours.size() + 1 >= m_minPts) {
            q->connectedness = Connectedness::kCore;
            seedSet.insert(
                seedSet.end(), q->neighbours.begin(), q->neighbours.end());
        }
    }
}

//======================================================================
void
IncrementalDBSCAN::add_hit(Hit* new_hit)
{
    static int ncall = 0;
    static int next_cluster_index = 0;

    m_hits.push_back(new_hit);
    m_latest_time = new_hit->time;
    // TODO: is it necessary to do this here? do we need the list of
    // neighbours before looping over the clusters?
    neighbours_sorted(m_hits, *new_hit, m_eps);


    // All the clusters that this hit neighboured. If there are
    // multiple clusters neighbouring this hit, we'll merge them at
    // the end
    using cluster_iterator = std::list<Cluster>::iterator;
    std::vector<cluster_iterator> clusters_to_merge;

    auto clust_it = m_clusters.begin();

    while (clust_it != m_clusters.end()) {
        Cluster& cluster = *clust_it;

        // If this cluster was already marked as complete, the new hit
        // can't be part of it. Skip it and remove from the list
        if (cluster.completeness == Completeness::kComplete) {
            clust_it = m_clusters.erase(clust_it);
            continue;
        }

        // Try adding the new hit to this cluster
        if (cluster.maybe_add_new_hit(new_hit, m_eps, m_minPts)) {
            clusters_to_merge.push_back(clust_it);
        }

        // TODO: should we only be doing this if we actually added the hit to this cluster?
        cluster_reachable(cluster.latest_core_point, cluster);

        // Mark the cluster complete if appropriate
        if (cluster.latest_time < m_latest_time - m_eps) {
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

        neighbours_sorted(m_hits, *new_hit, m_eps);

        if (new_hit->neighbours.size() + 1 >= m_minPts) {
            new_hit->connectedness = Connectedness::kCore;
            Cluster& new_cluster = m_clusters.emplace_back();
            new_cluster.index = next_cluster_index++;
            new_cluster.completeness = Completeness::kIncomplete;
            new_cluster.add_hit(new_hit);
            for (auto& neighbour : new_hit->neighbours) {
                new_cluster.add_hit(neighbour);
            }
            cluster_reachable(new_hit, new_cluster);
        }
    }

    ++ncall;
}

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
