// Walks that start and end at the rooted, two-ended interface.
// The state is a physical noncrossing matching with a distinguished gap.
// Edge construction changes its actual array of partners; there is no
// reflection quotient, reachability pruning, or spectral hypothesis.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

using Count=std::uint64_t;
using Id=std::uint32_t;
using Matching=std::array<unsigned char,28>;
using Row=std::array<Id,4>;
struct Offsets {Id even,odd;};
struct Graph {
    std::vector<Row> rows[2];
    Id root;
};

// Bracket flags are used as a lookup key only. Assigning a transition
// always joins/inserts/retains actual pairs in physical endpoint order.
Id key(const Matching& mate,int width) {
    Id flags=0;
    for(int i=0;i<width;++i) if(i<mate[i]) flags|=Id(1)<<i;
    return flags;
}
std::vector<std::vector<Matching>> matchings(int bound) {
    std::vector<std::vector<Matching>> all(bound+1);
    all[0].push_back(Matching{});
    for(int w=2;w<=bound;w+=2) for(int end=1;end<w;end+=2)
        for(const auto& inside:all[end-1]) for(const auto& after:all[w-end-1]) {
            Matching m{};m[0]=end;m[end]=0;
            for(int j=1;j<end;++j) m[j]=inside[j-1]+1;
            for(int j=end+1;j<w;++j)m[j]=after[j-end-1]+end+1;
            all[w].push_back(m);
        }
    return all;
}
Graph graph(int bound) {
    const auto all=matchings(bound);
    std::vector<std::unordered_map<Id,Offsets>> starts(bound+1);
    Id sizes[2]={1,1}; // Missing successors refer to zero, with an all-zero row.
    for(int w=2;w<=bound;w+=2) {
        starts[w].reserve(all[w].size());
        for(const auto& m:all[w]){
            auto inserted=starts[w].emplace(key(m,w),Offsets{sizes[0],sizes[1]});
            assert(inserted.second);
            sizes[0]+=w/2+1;sizes[1]+=w/2;
        }
    }
    Graph g;g.rows[0].resize(sizes[0]);g.rows[1].resize(sizes[1]);
    auto target=[&](int w,const Matching& m,int split)->Id {
        auto start=starts[w].at(key(m,w));
        return (split%2?start.odd:start.even)+split/2;
    };
    for(int w=2;w<=bound;w+=2) for(const auto& m:all[w]){
        for(int split=0;split<=w;++split){
            Row& r=g.rows[split%2][target(w,m,split)];
            if(split>0)r[0]=target(w,m,split-1);
            if(split<w)r[1]=target(w,m,split+1);
            if(w<bound){
                Matching child{};
                auto idx=[&](int j){return j+(j>=split?2:0);};
                for(int j=0;j<w;++j)child[idx(j)]=idx(m[j]);
                child[split]=split+1;child[split+1]=split;
                r[2]=target(w+2,child,split+1);
            }
            if(split>0 && split<w && m[split-1]!=split) {
                Matching child{};const int a=m[split-1],b=m[split];
                auto idx=[&](int j){return j-(j>split?2:0);};
                for(int j=0;j<w;++j)if(j!=split-1 && j!=split){
                    const int partner=(j==a?b:(j==b?a:m[j]));
                    child[idx(j)]=idx(partner);
                }
                r[3]=target(w-2,child,split-1);
            }
        }
    }
    g.root=target(2,all[2][0],1);
    return g;
}

inline Count sum(const Row& r,const std::vector<Count>& x) {
    return x[r[0]]+x[r[1]]+x[r[2]]+x[r[3]];
}
template<int P,int Q> void process(const Graph& g,int slots,int every) {
    constexpr Count seed=Count(1)<<45;
    constexpr Count cap=(Count(1)<<50)-1;
    static_assert(Q<=1000 && (16*cap)<=std::numeric_limits<Count>::max()/Q);
    std::vector<Count> current(g.rows[1].size()),half(g.rows[0].size());
    current[g.root]=1;
    for(int i=0;i<=9;++i){
        std::cout<<"small "<<i+1<<" "<<current[g.root]<<"\n";
        for(Id j=1;j<half.size();++j)half[j]=sum(g.rows[0][j],current);
        for(Id j=1;j<current.size();++j)current[j]=sum(g.rows[1][j],half);
    }
    std::fill(current.begin(),current.end(),0);current[g.root]=seed;
    Count max_current=seed;int exponent=0;
    std::uint64_t caps=0;
    std::cout<<"seed "<<seed<<" cap "<<cap<<"\n";
    for(int m=1;m<=slots;++m){
        // Keeping max <= 2^48 leaves a factor of sixteen in the two-slot sum.
        // The shift multiplies the rigorous total scale by 2^shift.
        const int shift=max_current>=(Count(1)<<46)?4:3;
        #pragma omp parallel for schedule(static)
        for(Id j=1;j<half.size();++j)half[j]=sum(g.rows[0][j],current);
        Count maximum=0,saturations=0;
        #pragma omp parallel for schedule(static) reduction(max:maximum) reduction(+:saturations)
        for(Id j=1;j<current.size();++j){
            const Count v=sum(g.rows[1][j],half);
            Count t;
            if constexpr(P==0)t=v>>shift;
            else {t=(v*Q)/P;if(t>cap){++saturations;t=cap;}}
            current[j]=t;maximum=std::max(maximum,t);
        }
        max_current=maximum;caps+=saturations;
        if constexpr(P==0){exponent+=shift;assert(max_current<=(Count(1)<<47));}
        const bool passed=P && current[g.root]>seed;
        if(m==1 || m%every==0 || m==slots || passed)
            std::cout<<"return "<<m<<" "<<current[g.root]<<" "<<exponent<<" "<<maximum<<" "<<caps<<"\n"<<std::flush;
        if(passed)break;
    }
}
int main(int argc,char** argv){
    if(argc!=6){std::cerr<<"return WIDTH SLOTS MODE THREADS INTERVAL\n";return 2;}
    const int w=std::stoi(argv[1]),m=std::stoi(argv[2]),threads=std::stoi(argv[4]),interval=std::stoi(argv[5]);
    if(w<2 || w>28 || w%2 || m<1 || m>1000000 || threads<1 || interval<1)return 2;
    omp_set_num_threads(threads);auto g=graph(w);
    std::cout<<"graph "<<w<<" "<<g.rows[0].size()-1<<" "<<g.rows[1].size()-1<<" "<<g.root<<"\n"<<std::flush;
    const std::string mode=argv[3];
    if(mode=="binary")process<0,1>(g,m,interval);
    else if(mode=="11.4")process<57,5>(g,m,interval);
    else if(mode=="11.381")process<11381,1000>(g,m,interval);
    else return 2;
}
