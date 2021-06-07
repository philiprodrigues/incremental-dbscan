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

#ifdef HAVE_PROFILER
#include "gperftools/profiler.h"
#endif

#include "CLI11.hpp"

//======================================================================
std::vector<Hit*> get_hits(std::string name, int nhits, int nskip)
{
    std::vector<Hit*> hits;
    if(name=="simple"){
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
    }
    else{
        std::ifstream fin(name);
        uint64_t timestamp, first_timestamp{0};
        int channel;
        int i=0;
        while(fin >> channel >> timestamp){
            if(first_timestamp==0) first_timestamp=timestamp;
            if(i++<nskip) continue;
            if(nhits>0 && i>nskip+nhits) break;

            hits.push_back(new Hit((timestamp-first_timestamp)/100, channel));
        }
    }
    return hits;
}

//======================================================================
void test_dbscan(std::string filename, int nhits, int nskip, bool test, std::string profile_filename)
{
    const int minPts=2;
    const float eps=10;

    auto hits=get_hits(filename, nhits, nskip);

    if(test){
        // Run the naive DBSCAN implementation for comparison with the
        // incremental one
        dbscan_orig(hits, eps, minPts);
        TCanvas* c=draw_clusters(hits);
        c->Print("dbscan-orig.png");
    }

    // Sort the hits by time for the incremental DBSCAN, which requires it
    std::vector<Hit*> hits_sorted(hits);
    std::sort(hits.begin(), hits.end(), [](Hit* a, Hit* b) { return a->time < b->time; });
    // If test==true, we've already set cluster indexes for the hits,
    // so we have to reset them here
    for(auto h: hits_sorted) h->cluster=kUndefined;

#ifdef HAVE_PROFILER
    // Start the profiler here, so the profile only measures the
    // incremental DBSCAN, not the hit reading and the original DBSCAN
    if(profile_filename!="") ProfilerStart(profile_filename.c_str());
#else
    if(profile_filename!="") std::cerr << "profile filename specified, but run_dbscan built without profiler support" << std::endl;
#endif

    DBSCANState state;
    TStopwatch ts;
    int i=0;
    double last_real_time=0;
    for(auto h: hits_sorted){
        dbscan_partial_add_one(state, h, eps, minPts);
        if(++i % 100000 == 0){
            double real_time=ts.RealTime();
            ts.Continue();
            std::cout << "100k hits took " << (real_time-last_real_time) << "s" << std::endl;
            last_real_time=real_time;
        }
    }

    // Give it a far-future hit so it goes through all of the hits
    Hit future_hit(10000, 110);
    dbscan_partial_add_one(state, &future_hit, eps, minPts);
    ts.Stop();

#ifdef HAVE_PROFILER
    if(profile_filename!="") ProfilerStop();
#endif

    // Clock is 50 MHz, but we divided the time by 100 when we read in the hits
    double data_time=(hits_sorted.back()->time - hits_sorted.front()->time)/50e4;
    double processing_time=ts.RealTime();
    std::cout << "Processed " << hits_sorted.size() << " hits representing " << data_time << "s of data in " << processing_time << "s. Ratio=" << (data_time/processing_time) << std::endl;
    if(test){
        TCanvas* c=draw_clusters(hits_sorted);
        c->Print("dbscan-incremental.png");
    }
}

//======================================================================
int main(int argc, char** argv)
{
    CLI::App cliapp{"Run incremental DBSCAN"};

    std::string filename;;
    cliapp.add_option("-f,--file", filename, "Input file of hits");
    bool test=false;
    cliapp.add_flag("-t,--test", test, "Test mode (show event display with clusters)");
    std::string profile;
    cliapp.add_option("-p,--profile", profile, "Run perftools profiler with output to file");
    int nskip=0;
    cliapp.add_option("-s,--nskip", nskip, "Number of hits at start of file to skip");
    int nhits=-1;
    cliapp.add_option("-n,--nhits", nhits, "Maximum number of hits to read from file");

    CLI11_PARSE(cliapp, argc, argv);

#ifndef HAVE_PROFILER
    if(profile!=""){
        std::cerr << "Profile filename specified but run_dbscan built without profiler support" << std::endl;
        exit(1);
    }
#endif
    
    int dummy_argc=1;
    const char* dummy_argv[]={"foo"};
    // TRint is here to start up the ROOT event loop so we can display the canvases on screen
    TRint* app=nullptr;
    if(test) app=new TRint("foo", &dummy_argc, const_cast<char**>(dummy_argv));

    test_dbscan(filename, nhits, nskip, test, profile);
    if(test) app->Run();
    delete app;
    return 0;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
