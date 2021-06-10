#include "Hit.hpp"
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
std::vector<dbscan::Hit*>
get_hits(std::string name, int nhits, int nskip)
{
    std::vector<dbscan::Hit*> hits;
    if (name == "simple") {
        hits.push_back(new dbscan::Hit(8.3, 101));
        hits.push_back(new dbscan::Hit(2.6, 103));
        hits.push_back(new dbscan::Hit(5.3, 104));
        hits.push_back(new dbscan::Hit(6.1, 105));
        hits.push_back(new dbscan::Hit(6.8, 106));
        hits.push_back(new dbscan::Hit(7.3, 107));
        hits.push_back(new dbscan::Hit(7.9, 108));
        hits.push_back(new dbscan::Hit(8.0, 109));
        hits.push_back(new dbscan::Hit(8.7, 110));
        hits.push_back(new dbscan::Hit(16.1, 105));
        hits.push_back(new dbscan::Hit(16.8, 106));
        hits.push_back(new dbscan::Hit(17.3, 107));
        hits.push_back(new dbscan::Hit(17.9, 108));
        hits.push_back(new dbscan::Hit(18.0, 109));
        hits.push_back(new dbscan::Hit(18.7, 110));
    } else {
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

            hits.push_back(
                new dbscan::Hit((timestamp - first_timestamp) / 100, channel));
        }
    }
    return hits;
}

//======================================================================
bool compare_clusters(std::vector<dbscan::Hit*>& v1,
                      std::vector<dbscan::Hit*>& v2)
{
    assert(v1.size()==v2.size());

    bool same=true;
    // Map from cluster index in v1 to cluster index in v2
    std::map<int, int> index_map;

    for(size_t i=0; i<v1.size(); ++i){
        dbscan::Hit* hit1=v1[i];
        dbscan::Hit* hit2=v2[i];
        if(hit1->time!=hit2->time || hit1->chan!=hit2->chan){
            std::cout << "Mismatched input vectors" << std::endl;
            same=false;
            break;
        }

        int index1=hit1->cluster;
        if(index1<0) index1=-1;
        int index2=hit2->cluster;
        if(index2<0) index2=-1;

        if(index_map.count(index1)){
            if(index2 != index_map[index1]){
                std::cout << "(" << hit1->time << ", " << hit1->chan << ") has cluster " << hit1->cluster << " but (" << hit2->time << ", " << hit2->chan << ") has cluster " << hit2->cluster << std::endl;
                same=false;
            }
        }
        else{
            index_map[index1]=index2;
        }
    }
    return same;
}

//======================================================================
void
test_dbscan(std::string filename,
            int nhits,
            int nskip,
            bool test,
            std::string profile_filename)
{
    const int minPts = 2;
    const float eps = 10;

    std::cout << "Reading hits" << std::endl;
    auto hits = get_hits(filename, nhits, nskip);
    std::cout << "Sorting hits" << std::endl;
    // Sort the hits by time for the incremental DBSCAN, which
    // requires it. We'll also give regular DBSCAN the sorted hits,
    // which will make later comparisons easier
    std::sort(hits.begin(), hits.end(), [](dbscan::Hit* a, dbscan::Hit* b) {
        return a->time < b->time;
    });

    if (test) {
        // Run the naive DBSCAN implementation for comparison with the
        // incremental one
        std::cout << "Running dbscan_orig" << std::endl;
        dbscan::dbscan_orig(hits, eps, minPts);
        TCanvas* c = draw_clusters(hits);
        c->Print("dbscan-orig.png");
    }

    // We make a copy so we can compare the output of dbscan_orig and IncrementalDBSCAN
    std::cout << "Copying hit vector for incremental dbscan" << std::endl;
    std::vector<dbscan::Hit*> hits_inc;
    for(auto h: hits) {
        hits_inc.push_back(new dbscan::Hit(h->time, h->chan));
    }

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
    for (auto h : hits_inc) {
        dbscanner.add_hit(h);
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
    dbscan::Hit future_hit(10000, 110);
    dbscanner.add_hit(&future_hit);
    ts.Stop();

#ifdef HAVE_PROFILER
    if (profile_filename != "")
        ProfilerStop();
#endif

    // Clock is 50 MHz, but we divided the time by 100 when we read in the hits
    double data_time =
        (hits_inc.back()->time - hits_inc.front()->time) / 50e4;
    double processing_time = ts.RealTime();
    std::cout << "Processed " << hits_inc.size() << " hits representing "
              << data_time << "s of data in " << processing_time
              << "s. Ratio=" << (data_time / processing_time) << std::endl;
    if (test) {
        TCanvas* c = dbscan::draw_clusters(hits_inc);
        c->Print("dbscan-incremental.png");

        bool same=compare_clusters(hits, hits_inc);
        if(same){
            std::cout << "dbscan_orig and incremental results matched" << std::endl;
        }
        else{
            std::cout << "dbscan_orig and incremental results differed" << std::endl;
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
        "-t,--test", test, "Test mode (show event display with clusters)");
    std::string profile;
    cliapp.add_option(
        "-p,--profile", profile, "Run perftools profiler with output to file");
    int nskip = 0;
    cliapp.add_option(
        "-s,--nskip", nskip, "Number of hits at start of file to skip");
    int nhits = -1;
    cliapp.add_option(
        "-n,--nhits", nhits, "Maximum number of hits to read from file");

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
    if (test)
        app = new TRint("foo", &dummy_argc, const_cast<char**>(dummy_argv));

    test_dbscan(filename, nhits, nskip, test, profile);
    if (test)
        app->Run();
    delete app;
    return 0;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
