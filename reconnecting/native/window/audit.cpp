// Separate checker: join equal component labels on a padded list of open ends.
// No transition, state numbering or power iterations are imported from the generator.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
using U=std::uint64_t;
using Wide=__uint128_t;
using Labels=std::array<int,38>;

std::vector<U> read(const std::string& filename) {
    std::ifstream f(filename,std::ios::binary|std::ios::ate);auto size=f.tellg();
    if(size<8 || size%8) throw std::runtime_error("unexpected size");
    std::vector<U> a(size/8);f.seekg(0);f.read(reinterpret_cast<char*>(a.data()),size);
    if(!f) throw std::runtime_error("read error");
    return a;
}

Labels decode(U code,int m) {
    Labels labels{},open{};int depth=0;
    for(int i=0;i<m+2;++i) labels[i]=i; // the first and last labels represent fresh padding
    for(int i=1;i<=m;++i) {
        switch(code&3) {
          case 0: break;
          case 1: open[depth++]=i;break;
          case 2: if(!depth) throw std::runtime_error("bad close");labels[i]=open[--depth];break;
          default: throw std::runtime_error("bad symbol");
        }
        code>>=2;
    }
    if(depth||code) throw std::runtime_error("bad state");
    return labels;
}

U encode(const Labels& labels,int first,int size) {
    Labels seen;seen.fill(-1);U code=0;
    for(int j=0;j<size;++j) {
        auto &prior=seen[labels[first+j]];
        if(prior>=0) code|=(U(1)<<(2*prior))|(U(2)<<(2*j));
        else prior=j;
    }
    return code;
}

// Immutable key catalogue. Indices are stored in power-of-two buckets with
// triangular probing. A matching hash always requires a full word comparison.
class Catalogue {
    static constexpr unsigned empty=std::numeric_limits<unsigned>::max();
    const std::vector<U>& words;
    std::vector<unsigned> buckets;
    size_t position(U word) const {
        U hash=word;
        hash^=hash>>33;hash*=0xff51afd7ed558ccdULL;
        hash^=hash>>33;hash*=0xc4ceb9fe1a85ec53ULL;hash^=hash>>33;
        size_t mask=buckets.size()-1,slot=hash&mask,step=0;
        while(buckets[slot]!=empty && words[buckets[slot]]!=word) slot=(slot+(++step))&mask;
        return slot;
    }
public:
    explicit Catalogue(const std::vector<U>& keys):words(keys),
        buckets(std::bit_ceil(keys.size()+keys.size()/2+1),empty) {
        if(keys.size()>=empty) throw std::invalid_argument("too many states");
        for(unsigned i=0;i<keys.size();++i) {
            auto& cell=buckets[position(keys[i])];
            if(cell!=empty) throw std::runtime_error("duplicate state");
            cell=i;
        }
    }
    unsigned at(U word) const {
        unsigned i=buckets[position(word)];
        if(i==empty) throw std::runtime_error("missing possible successor");
        return i;
    }
};

int main(int argc,char**argv) {
    if(argc<7||argc>8) throw std::invalid_argument("prefix k p q numerator denominator [--fold]");
    bool fold=argc==8&&std::string(argv[7])=="--fold";
    if(argc==8&&!fold) throw std::invalid_argument("unrecognised option");
    auto codes=read(std::string(argv[1])+".states"),v=read(std::string(argv[1])+".vector");
    const int k=std::stoi(argv[2]),m=2*k;
    const U p=std::stoull(argv[3]),q=std::stoull(argv[4]),hi=std::stoull(argv[5]),den=std::stoull(argv[6]);
    if(k<1||k>16||p<q||q<1||p>200||!hi||!den||den>1000000000000ULL||
       Wide(hi)*p*q>std::numeric_limits<U>::max()||v.size()!=codes.size()||codes[0]!=0)
        throw std::runtime_error("input outside exact arithmetic bounds");
    const std::array<U,4> weight{p*p,p*q,p*q,q*q};
    if(*std::min_element(v.begin(),v.end())==0) throw std::runtime_error("zero entry");
    Catalogue id(codes);
    U rejects=0,edges=0;
    for(size_t i=0;i<v.size();++i) {
        auto l=decode(codes[i],m);Wide sum=0;
        for(int op=0;op<4;++op) {
            auto work=l;int first=0;
            if(op==0) {
                // Insert twice the same unused label; trim two extreme slots at either end.
                std::move_backward(work.begin()+k+1,work.begin()+m+2,work.begin()+m+4);
                work[k+1]=work[k+2]=m+3;first=2;
            } else if(op==2) first=2;
            else if(op==3) {
                int a=l[k],b=l[k+1];
                if(a==b) {++rejects;continue;}
                for(int j=0;j<m+2;++j) if(work[j]==b) work[j]=a;
                std::move(work.begin()+k+2,work.begin()+m+2,work.begin()+k);
            }
            U next=encode(work,first,m);
            if(fold) {
                std::reverse(work.begin()+first,work.begin()+first+m);
                next=std::min(next,encode(work,first,m));
            }
            sum+=Wide(weight[op])*v[id.at(next)];++edges;
        }
        if(sum*den>Wide(hi)*p*q*v[i]) throw std::runtime_error("failed upper comparison");
    }
    std::cout<<"{\"states\":"<<v.size()<<",\"admitted_edges\":"<<edges<<",\"rejected_closures\":"<<rejects
             <<",\"all_rows_pass\":true,\"v_min\":"<<*std::min_element(v.begin(),v.end())
             <<",\"v_empty\":"<<v.front()<<"}\n";
}
