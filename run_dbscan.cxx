#include "dbscan.h"
#include "TStopwatch.h"

#include <thread>
#include <chrono>
#include <fstream>

std::vector<Hit*> get_hits(std::string name)
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
        while(fin >> channel >> timestamp){
            if(first_timestamp==0) first_timestamp=timestamp;
            hits.push_back(new Hit((timestamp-first_timestamp)/100, channel));
        }
    }
    return hits;
}

void test_dbscan(const char* filename)
{
    auto hits=get_hits(filename);

    TStopwatch ts;
    // dbscan_orig(hits, 5, 2);
    // ts.Stop();
    // ts.Print();
    // draw_clusters(hits);

    std::vector<Hit*> hits_sorted(hits);
    std::sort(hits.begin(), hits.end(), [](Hit* a, Hit* b) { return a->time < b->time; });
    // for(auto h: hits_sorted) h->cluster=kUndefined;
    // dbscan_sorted_input(hits_sorted, 5, 2);
    // draw_clusters(hits_sorted);

    for(auto h: hits_sorted) h->cluster=kUndefined;
    
    State state;
    ts.Reset();
    ts.Start();
    for(auto h: hits_sorted){
        dbscan_partial_add_one(state, h, 5, 2);
    }

    // Give it a far-future hit so it goes through all of the hits
    Hit future_hit(10000, 110);
    dbscan_partial_add_one(state, &future_hit, 5, 2);
    ts.Print();
    // draw_clusters(hits_sorted);
    
}

int main(int argc, char** argv)
{
    if(argc!=2){
        std::cout << "Usage: run_dbscan input_file" << std::endl;
        return 1;
    }
    test_dbscan(argv[1]);
    return 0;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// c-file-style: "linux"
// End:
