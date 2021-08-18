#include "Hit.hpp"
#include "Point.hpp"

#include "dbscan.hpp"
#include "dbscan_orig.hpp"
#include "draw_clusters.hpp"

#include "TStopwatch.h"
#include "TRint.h"
#include "TCanvas.h"

#include <thread>
#include <chrono>
#include <fstream>
#include <string>
#include <cassert>

#ifdef HAVE_PROFILER
#include "gperftools/profiler.h"
#endif

#include "CLI11.hpp"

//======================================================================
std::vector<Point>
get_points(std::string name, int nhits, int nskip)
{
    std::vector<Point> points;

    std::ifstream fin(name);
    uint64_t timestamp, first_timestamp{ 0 };
    int channel;
    int i = 0;
    while (fin >> channel >> timestamp) {
        if (first_timestamp == 0)
            first_timestamp = timestamp;
        if (i++ < nskip)
            continue;
        if (nhits > 0 && i > nskip + nhits)
            break;

        points.push_back({channel, (timestamp - first_timestamp) / 100});
    }

    return points;
}

std::vector<dbscan::Hit*>
points_to_hits(const std::vector<Point>& points)
{
    std::vector<dbscan::Hit*> ret;
    for(auto const& p: points){
        ret.push_back(new dbscan::Hit(p.time, p.chan));
    }
    return ret;
}

//======================================================================
bool
cluster_has_hit(const dbscan::Cluster& cluster, const dbscan::Hit* test_hit)
{
    for(auto const& hit : cluster.hits){
        if(hit->time == test_hit->time &&
           hit->chan == test_hit->chan) {
            return true;
        }
    }
    return false;
}

//======================================================================
void print_cluster_hits(const dbscan::Cluster& cluster)
{
    for(auto const& hit : cluster.hits){
        std::cout << std::hex << hit << std::dec << " " << hit->time << ", " << hit->chan << std::endl;
    }
}

//======================================================================
bool
compare_clusters(std::vector<dbscan::Cluster>& clusters1, std::vector<dbscan::Cluster>& clusters2)
{
    bool ok=true;
    
    if(clusters1.size() != clusters2.size()){
        std::cout << "clusters1 has " << clusters1.size() << " clusters but clusters2 has " << clusters2.size() << " clusters" << std::endl;
        ok=false;
    }

    for(auto const& cluster1 : clusters1){
        // First, find the cluster in the other list that contains the first hit from cluster1
        const dbscan::Cluster* other_cluster=nullptr;
        dbscan::Hit* hit1=cluster1.hits.hits[0];
        for(auto const& cluster2 : clusters2){
            if(cluster_has_hit(cluster2, hit1)){
                other_cluster=&cluster2;
                break;
            }
        }
        
        if(!other_cluster){
            std::cout << "(" << hit1->time << ", " << hit1->chan
                      << ") has cluster " << hit1->cluster << " but is not present in clusters2" << std::endl;
            ok=false;
        }

        if(cluster1.hits.size() != other_cluster->hits.size()){
            std::cout << "cluster1 has " << cluster1.hits.size() << " hits but other_cluster has " << other_cluster->hits.size() << " hits" << std::endl;
            std::cout << "cluster1 hits:" << std::endl;
            print_cluster_hits(cluster1);
            std::cout << "other_cluster hits:" << std::endl;
            print_cluster_hits(*other_cluster);
            ok=false;
        }
        
        for(auto const& hit : cluster1.hits){
            if(!cluster_has_hit(*other_cluster, hit)){
                std::cout << "Hit (" << hit1->time << ", " << hit1->chan << ") is present in cluster1 but not other_cluster" << std::endl;
                ok=false;
            }
        }
    }

    return ok;
    // bool all_same = true;
    // // Map from cluster index in v1 to cluster index in v2
    // std::map<int, int> index_map;

    // // Noise must map to noise
    // index_map[-1]=-1;

    // for (size_t i = 0; i < v1.size(); ++i) {
    //     for (size_t j = 0; j < v1[i].hits.size(); ++j) {
    //     dbscan::Hit* hit1 = v1[i];
    //     dbscan::Hit* hit2 = v2[i];
    //     if (hit1->time != hit2->time || hit1->chan != hit2->chan) {
    //         std::cout << "Mismatched input vectors" << std::endl;
    //         all_same = false;
    //         break;
    //     }

    //     int index1 = hit1->cluster;
    //     if (index1 < 0)
    //         index1 = -1;
    //     int index2 = hit2->cluster;
    //     if (index2 < 0)
    //         index2 = -1;

    //     bool differ = false;
    //     if((index1 < 0 && index2 >= 0) ||
    //        (index1 >= 0 && index2 < 0)){
    //         // One is noise, other is in a cluster
    //         differ = true;
    //     }
    //     else if (index_map.count(index1)) {
    //         if (index2 != index_map[index1]) {
    //             differ = true;
    //         }
    //     } else {
    //         index_map[index1] = index2;
    //     }

    //     if(differ){
    //         std::cout << "(" << hit1->time << ", " << hit1->chan
    //                   << ") has cluster " << hit1->cluster << " but ("
    //                   << hit2->time << ", " << hit2->chan
    //                   << ") has cluster " << hit2->cluster << std::endl;
    //         all_same=false;
    //     }

    // }
    // return all_same;
}

//======================================================================
void
test_dbscan(std::string filename,
            int nhits,
            int nskip,
            bool test,
            bool plot,
            std::string profile_filename,
            int minPts,
            float eps)
{
    std::cout << "Reading hits" << std::endl;
    auto points = get_points(filename, nhits, nskip);
    std::cout << "Sorting hits" << std::endl;
    // Sort the hits by time for the incremental DBSCAN, which
    // requires it. We'll also give regular DBSCAN the sorted hits,
    // which will make later comparisons easier
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.time < b.time;
    });

    std::vector<dbscan::Cluster> clusters_orig;
    if (test) {
        // Run the naive DBSCAN implementation for comparison with the
        // incremental one
        auto hits=points_to_hits(points);
        std::cout << "Running dbscan_orig" << std::endl;
        auto clusters=dbscan::dbscan_orig(hits, eps, minPts);
        clusters_orig=clusters;
        if(plot){
            TCanvas* c = draw_clusters(clusters, points);
            c->Print("dbscan-orig.png");
        }
    }

    // // We make a copy so we can compare the output of dbscan_orig and
    // // IncrementalDBSCAN
    // std::cout << "Copying hit vector for incremental dbscan" << std::endl;
    // std::vector<dbscan::Hit*> hits_inc;
    // for (auto h : hits) {
    //     hits_inc.push_back(new dbscan::Hit(h->time, h->chan));
    // }

#ifdef HAVE_PROFILER
    // Start the profiler here, so the profile only measures the
    // incremental DBSCAN, not the hit reading and the original DBSCAN
    if (profile_filename != "")
        ProfilerStart(profile_filename.c_str());
#else
    if (profile_filename != "")
        std::cerr << "profile filename specified, but run_dbscan built without "
                     "profiler support"
                  << std::endl;
#endif

    std::cout << "Running incremental dbscan" << std::endl;
    dbscan::IncrementalDBSCAN dbscanner(eps, minPts);
    TStopwatch ts;
    int i = 0;
    double last_real_time = 0;
    std::vector<dbscan::Cluster> clusters;
    for (auto p : points) {
        dbscanner.add_point(p.time, p.chan, &clusters);
        if (++i % 100000 == 0) {
            double real_time = ts.RealTime();
            ts.Continue();
            std::cout << "100k hits took " << (real_time - last_real_time)
                      << "s" << std::endl;
            last_real_time = real_time;
        }
        dbscanner.trim_hits();
    }

    // Give it a far-future hit so it goes through all of the hits
    Point future_point{110, 10000000};
    dbscanner.add_point(future_point.time, future_point.chan, &clusters);
    ts.Stop();

#ifdef HAVE_PROFILER
    if (profile_filename != "")
        ProfilerStop();
#endif

    // Clock is 50 MHz, but we divided the time by 100 when we read in the hits
    double data_time = (points.back().time - points.front().time) / 50e4;
    double processing_time = ts.RealTime();
    std::cout << "Found " << clusters.size() << " clusters total" << std::endl;
    std::cout << "Processed " << points.size() << " hits representing "
              << data_time << "s of data in " << processing_time
              << "s. Ratio=" << (data_time / processing_time) << std::endl;

    if (plot) {
        TCanvas* c = dbscan::draw_clusters(clusters, points);
        c->Print("dbscan-incremental.png");
    }

    if (test) {
        bool same = compare_clusters(clusters_orig, clusters);
        if (same) {
            std::cout << "dbscan_orig and incremental results matched"
                      << std::endl;
        } else {
            std::cout << "dbscan_orig and incremental results differed"
                      << std::endl;
        }
    }
}

//======================================================================
int
main(int argc, char** argv)
{
    CLI::App cliapp{ "Run incremental DBSCAN" };

    std::string filename;
    ;
    cliapp.add_option("-f,--file", filename, "Input file of hits");
    bool test = false;
    cliapp.add_flag(
        "-t,--test", test, "Test mode (compare to original dbscan)");
    bool plot = false;
    cliapp.add_flag("--plot", plot, "Plot results");
    std::string profile;
    cliapp.add_option(
        "-p,--profile", profile, "Run perftools profiler with output to file");
    int nskip = 0;
    cliapp.add_option(
        "-s,--nskip", nskip, "Number of hits at start of file to skip");
    int nhits = -1;
    cliapp.add_option(
        "-n,--nhits", nhits, "Maximum number of hits to read from file");
    int minPts = 2;
    cliapp.add_option(
        "-m,--minpts", minPts, "Minimum number of hits to form a cluster");
    float eps=10;
    cliapp.add_option(
        "-d,--distance", eps, "Distance threshold for points to be neighbours");

    CLI11_PARSE(cliapp, argc, argv);

#ifndef HAVE_PROFILER
    if (profile != "") {
        std::cerr << "Profile filename specified but run_dbscan built without "
                     "profiler support"
                  << std::endl;
        exit(1);
    }
#endif

    int dummy_argc = 1;
    const char* dummy_argv[] = { "foo" };
    // TRint is here to start up the ROOT event loop so we can display the
    // canvases on screen
    TRint* app = nullptr;
    if (plot)
        app = new TRint("foo", &dummy_argc, const_cast<char**>(dummy_argv));

    test_dbscan(filename, nhits, nskip, test, plot, profile, minPts, eps);
    if (plot)
        app->Run();
    delete app;
    return 0;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
