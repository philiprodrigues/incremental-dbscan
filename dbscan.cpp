#include "dbscan.hpp"
#include "Hit.hpp"

#include <cassert>
#include <limits>

namespace dbscan {

//======================================================================
int
neighbours_sorted(const std::vector<Hit*>& hits, Hit& q, float eps, int minPts)
{
    int n = 0;
    // Loop over the hits starting from the latest hit, since we will
    // ~always be adding a hit at recent times
    for (auto hit_it = hits.rbegin(); hit_it != hits.rend(); ++hit_it) {
        if ((*hit_it)->time > q.time + eps)
            continue;
        if ((*hit_it)->time < q.time - eps)
            break;

        if (q.add_potential_neighbour(*hit_it, eps, minPts))
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
    // neighbours, so start the search there in the sorted list of hits in this
    // cluster
    auto begin_it = std::lower_bound(
        hits.begin(), hits.end(), new_hit->time - eps, time_comp_lower);

    for (auto it = begin_it; it != hits.end(); ++it) {
        Hit* h = *it;
        if (h->add_potential_neighbour(new_hit, eps, minPts)) {
            do_add = true;
            if (h->neighbours.size() + 1 >= minPts) {
                h->connectedness = Connectedness::kCore;
            } else {
                h->connectedness = Connectedness::kEdge;
            }
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
        (!latest_core_point || h->time > latest_core_point->time)) {
        latest_core_point = h;
    }
}

//======================================================================
void
Cluster::steal_hits(Cluster& other)
{
    // TODO: it might be faster to do some sort of explicit "merge" of the hits,
    // eg:
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
IncrementalDBSCAN::cluster_reachable(Hit* seed_hit, Cluster& cluster)
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
void IncrementalDBSCAN::add_point(float time, float channel)
{
    Hit& new_hit=m_hit_pool[m_pool_end];
    new_hit.reset(time, channel);
    ++m_pool_end;
    if(m_pool_end==m_hit_pool.size()) m_pool_end=0;
    add_hit(&new_hit);
}
    
//======================================================================
void
IncrementalDBSCAN::add_hit(Hit* new_hit)
{
    // TODO: this should be a member variable, not a static, in case
    // there are multiple IncrementalDBSCAN instances
    static int next_cluster_index = 0;

    m_hits.push_back(new_hit);
    m_latest_time = new_hit->time;

    // All the clusters that this hit neighboured. If there are
    // multiple clusters neighbouring this hit, we'll merge them at
    // the end
    std::set<int> clusters_neighbouring_hit;

    // Find all the hit's neighbours
    neighbours_sorted(m_hits, *new_hit, m_eps, m_minPts);

    for (auto neighbour : new_hit->neighbours) {
        if (neighbour->cluster != kUndefined && neighbour->cluster != kNoise &&
            neighbour->neighbours.size() + 1 >= m_minPts) {
            // This neighbour is a core point in a cluster. Add the cluster to the list of
            // clusters that will contain this hit
            clusters_neighbouring_hit.insert(neighbour->cluster);
        }
    }

    if (clusters_neighbouring_hit.empty()) {
        // This hit didn't match any existing cluster. See if we can
        // make a new cluster out of it. Otherwise mark it as noise

        if (new_hit->neighbours.size() + 1 >= m_minPts) {
            // std::cout << "New cluster starting at hit time " << new_hit->time << " with " << new_hit->neighbours.size() << " neighbours" << std::endl;
            new_hit->connectedness = Connectedness::kCore;
            auto new_it = m_clusters.emplace_hint(
                m_clusters.end(), next_cluster_index, next_cluster_index);
            Cluster& new_cluster = new_it->second;
            new_cluster.completeness = Completeness::kIncomplete;
            new_cluster.add_hit(new_hit);
            next_cluster_index++;
            cluster_reachable(new_hit, new_cluster);
        }
        else{
            // std::cout << "New hit time " << new_hit->time << " with " << new_hit->neighbours.size() << " neighbours is noise" << std::endl;
        }
    } else {
        // This hit neighboured at least one cluster. Add the hit and
        // its noise neighbours to the first cluster in the list, then
        // merge the rest of the clusters into it

        auto index_it = clusters_neighbouring_hit.begin();

        auto it = m_clusters.find(*index_it);
        assert(it != m_clusters.end());
        Cluster& cluster = it->second;
        // std::cout << "Adding hit time " << new_hit->time << " with " << new_hit->neighbours.size() << " neighbours to existing cluster" << std::endl;
        cluster.add_hit(new_hit);

        // TODO: this seems wrong: we're adding this hit's neighbours
        // to the cluster even if this hit isn't a core point, but if
        // I wrap the whole thing in "if(new_hit is core)" then the
        // results differ from classic DBSCAN
        for (auto q : new_hit->neighbours) {
            if (q->cluster == kUndefined || q->cluster == kNoise) {
                // std::cout << "  Adding hit time " << q->time << " to existing cluster" << std::endl;
                cluster.add_hit(q);
            }
            // If the neighbouring hit q has exactly m_minPts
            // neighbours, it must have become a core point by the
            // addition of new_hit. Add q's neighbours to the cluster
            if(q->neighbours.size() + 1 == m_minPts){
                for (auto r : q->neighbours) {
                    cluster.add_hit(r);
                }
            }
        }


        ++index_it;

        for (; index_it != clusters_neighbouring_hit.end(); ++index_it) {

            auto other_it = m_clusters.find(*index_it);
            assert(other_it != m_clusters.end());
            Cluster& other_cluster = other_it->second;
            cluster.steal_hits(other_cluster);
        }
    }

    // Last case: new_hit and its neighbour are both noise, but the
    // addition of new_hit makes the neighbour a core point. So we
    // start a new cluster at the neighbour, and walk out from there
    for (auto& neighbour : new_hit->neighbours) {
        if(neighbour->neighbours.size() + 1 >= m_minPts){
            // std::cout << "new_hit's neighbour at " << neighbour->time << " has " << neighbour->neighbours.size() << " neighbours, so is core" << std::endl;
            if(neighbour->cluster==kNoise || neighbour->cluster==kUndefined){
                if(new_hit->cluster==kNoise || new_hit->cluster==kUndefined){
                    auto new_it = m_clusters.emplace_hint(
                                                          m_clusters.end(), next_cluster_index, next_cluster_index);
                    Cluster& new_cluster = new_it->second;
                    new_cluster.completeness = Completeness::kIncomplete;
                    new_cluster.add_hit(neighbour);
                    next_cluster_index++;
                    cluster_reachable(neighbour, new_cluster);
                }
            }
        }
        else {
            // std::cout << "new_hit's neighbour at " << neighbour->time << " has " << neighbour->neighbours.size() << " neighbours, so is NOT core" << std::endl;
        }
    }


    // Delete any completed clusters from the list
    auto clust_it = m_clusters.begin();
    while (clust_it != m_clusters.end()) {
        Cluster& cluster = clust_it->second;

        if (cluster.latest_time < m_latest_time - m_eps) {
            cluster.completeness = Completeness::kComplete;
        }

        if (cluster.completeness == Completeness::kComplete) {
            clust_it = m_clusters.erase(clust_it);
            continue;
        } else {
            ++clust_it;
        }
    }
}

void
IncrementalDBSCAN::trim_hits()
{
    // Find the earliest time of a hit in any cluster in the list (active or
    // not)
    float earliest_time = std::numeric_limits<float>::max();

    for (auto& cluster : m_clusters) {
        earliest_time =
            std::min(earliest_time, (*cluster.second.hits.begin())->time);
    }

    // If there were no clusters, set the earliest_time to the latest time
    // (otherwise it would still be FLOAT_MAX)
    if (m_clusters.empty()) {
        earliest_time = m_latest_time;
    }
    auto last_it = std::lower_bound(m_hits.begin(),
                                    m_hits.end(),
                                    earliest_time - 10 * m_eps,
                                    time_comp_lower);

    m_hits.erase(m_hits.begin(), last_it);
}

}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
