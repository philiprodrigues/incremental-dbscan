#include "dbscan.h"
#include "TStopwatch.h"

#include <thread>
#include <chrono>

void test_dbscan()
{
    std::vector<Hit*> hits;
    hits.push_back(new Hit(8.3, 101));
    hits.push_back(new Hit(2.6, 103));
    hits.push_back(new Hit(5.3, 104));
    hits.push_back(new Hit(6.1, 105));
    hits.push_back(new Hit(6.8, 106));
    hits.push_back(new Hit(7.3, 107));
    hits.push_back(new Hit(7.9, 108));
    hits.push_back(new Hit(8.0, 109));
    hits.push_back(new Hit(8.7, 110));


    hits.push_back(new Hit(16.1, 105));
    hits.push_back(new Hit(16.8, 106));
    hits.push_back(new Hit(17.3, 107));
    hits.push_back(new Hit(17.9, 108));
    hits.push_back(new Hit(18.0, 109));
    hits.push_back(new Hit(18.7, 110));


    TStopwatch ts;
    dbscan_orig(hits, 5, 2);
    ts.Stop();
    ts.Print();
    draw_clusters(hits);

    std::vector<Hit*> hits_sorted(hits);
    std::sort(hits.begin(), hits.end(), [](Hit* a, Hit* b) { return a->time < b->time; });
    dbscan_sorted_input(hits_sorted, 5, 2);
    draw_clusters(hits_sorted);

    
}
