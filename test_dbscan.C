#include "dbscan.h"
#include "TStopwatch.h"

#include <thread>
#include <chrono>

void test_dbscan()
{
    std::vector<Hit*> hits;
    hits.push_back(new Hit(83, 101, 1, 1));
    hits.push_back(new Hit(26, 103, 1, 1));
    hits.push_back(new Hit(53, 104, 1, 1));
    hits.push_back(new Hit(61, 105, 1, 1));
    hits.push_back(new Hit(68, 106, 1, 1));
    hits.push_back(new Hit(73, 107, 1, 1));
    hits.push_back(new Hit(79, 108, 1, 1));
    hits.push_back(new Hit(80, 109, 1, 1));
    hits.push_back(new Hit(87, 110, 1, 1));


    hits.push_back(new Hit(161, 105, 1, 1));
    hits.push_back(new Hit(168, 106, 1, 1));
    hits.push_back(new Hit(173, 107, 1, 1));
    hits.push_back(new Hit(179, 108, 1, 1));
    hits.push_back(new Hit(180, 109, 1, 1));
    hits.push_back(new Hit(187, 110, 1, 1));


    TStopwatch ts;
    dbscan_orig(hits, 5, 2);

    ts.Stop();
    ts.Print();
    draw_clusters(hits);
}
