constexpr unsigned gap_weight[]={611,624,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,626,624,611};
// Check one saved integer vector against the two physical slit moves, insertion
// and splice. Standalone: this reader neither imports the generator nor uses
// its rank, reflection, edge tables, parity sweeps or fraction comparisons.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using U = std::uint64_t;
using Word = std::uint32_t;
using Wide = unsigned __int128;
constexpr U scale = U(1) << 52;
struct Ratio {
    U n, d;
    bool below(Ratio other) const { return Wide(n)*other.d < Wide(other.n)*d; }
};

struct Layer {
    int k = 0;
    std::array<U,2> first{};
    Word mask = 0;
    std::vector<Word> words, prefix, suffix;
    U index(int parity, Word r, int u) const {
        return first[parity] + U(r)*((k+2-parity)/2) + u/2;
    }
    Word rank(Word p) const { return prefix[p & mask]+suffix[p >> k]; }
    void enumerate(int pairs, std::array<U,2>& counts) {
        k = pairs;
        first = counts;
        mask = (Word(1)<<k)-1;
        std::vector<std::vector<Word>> tails(k+1);
        std::vector<int> height(mask+1,-1);
        prefix.resize(mask+1);
        suffix.resize(mask+1);
        for (Word part=0; part<=mask; ++part) {
            int h=0, minimum=0;
            for (int pos=0; pos<k; ++pos) {
                h += ((part>>pos)&1) ? 1 : -1;
                minimum=std::min(minimum,h);
            }
            if (!minimum) height[part]=h; // valid first half
            if (h==minimum) {            // valid second half from -h down to zero
                suffix[part]=tails[-h].size();
                tails[-h].push_back(part);
            }
        }
        for (Word head=0; head<=mask; ++head) if (height[head]>=0) {
            prefix[head]=words.size();
            for (Word tail:tails[height[head]]) words.push_back(head | (tail<<k));
        }
        for (int p=0; p<2; ++p) counts[p]+=U(words.size())*((k+2-p)/2);
    }
};

struct Check {
    std::vector<Layer> layers;
    std::array<std::vector<U>,2> values;
    int cap;
    explicit Check(int half):layers(half+1),cap(half) {
        std::array<U,2> counts{};
        for (int k=1; k<=cap; ++k) layers[k].enumerate(k,counts);
        for (int p=0; p<2; ++p) values[p].resize(counts[p]);
    }
    void sweep(int p, const std::vector<U>* original=nullptr) {
        auto& dest=values[p];
        const auto& src=values[1-p];
        std::array<Ratio,2> lower{{{~U(0),1},{~U(0),1}}}, upper{{{0,1},{0,1}}};
        std::array<bool,2> positive{true,true};
        for (int k=1; k<=cap; ++k) {
            const Layer& l=layers[k];
            #pragma omp parallel
            {
                std::array<Ratio,2> lo{{{~U(0),1},{~U(0),1}}}, hi{{{0,1},{0,1}}};
                std::array<bool,2> valid{true,true};
                #pragma omp for schedule(static)
                for (Word r=0; r<l.words.size(); ++r) {
                    Word word=l.words[r], reflected=0;
                    std::array<int,32> mate{}, stack{};
                    int depth=0;
                    for (int i=0; i<2*k; ++i) {
                        if ((word>>i)&1) stack[depth++]=i;
                        else {
                            int j=stack[--depth];
                            mate[j]=i; mate[i]=j;
                            reflected |= Word(1)<<(2*k-1-i);
                        }
                    }
                    for (int u=p; u<=k; u+=2) {
                        U z=0;
                        // Both shifts are labelled; equality of targets keeps both terms.
                        U left = u ? src[l.index(1-p,r,u-1)] : 0;
                        U right = u<k ? src[l.index(1-p,r,u+1)]
                                      : src[l.index(1-p,l.rank(reflected),u-1)];
                        z = left + right;
                        if (k == cap) {
                            Wide add = Wide(gap_weight[u]) * right;
                            if (u) add += Wide(gap_weight[u-1]) * left;
                            z += U((add + 4095) >> 12);
                        }
                        if (k<cap) {
                            Word low=(Word(1)<<u)-1;
                            Word wider=(word&low) | (Word(1)<<u) | ((word>>u)<<(u+2));
                            const auto& t=layers[k+1];
                            z+=src[t.index(1-p,t.rank(wider),u+1)];
                        }
                        if (u && mate[u-1]!=u) {
                            int a=mate[u-1], b=mate[u];
                            Word spliced=word|(Word(1)<<std::min(a,b));
                            spliced &= ~(Word(1)<<std::max(a,b));
                            Word below=(Word(1)<<(u-1))-1;
                            Word smaller=(spliced&below) | ((spliced>>(u+1))<<(u-1));
                            const auto& t=layers[k-1];
                            z+=src[t.index(1-p,t.rank(smaller),u-1)];
                        }
                        U index=l.index(p,r,u);
                        if (p) {
                            const int candidates=original ? 2 : 1;
                            for (int t=0; t<candidates; ++t) {
                                U n=z, d=dest[index];
                                if (original) {
                                    U a=t==0 ? 2 : 5, b=t==0 ? 23 : 58;
                                    U previous=(*original)[index];
                                    if (a*z<=b*d || a*d<=b*previous) {valid[t]=false;continue;}
                                    n=a*z-b*d;
                                    d=a*d-b*previous;
                                }
                                Ratio next{n,d};
                                if (next.below(lo[t])) lo[t]=next;
                                if (hi[t].below(next)) hi[t]=next;
                            }
                        }
                        dest[index]=z;
                    }
                }
                #pragma omp critical
                {
                    for (int t=0; t<2; ++t) {
                        if (lo[t].below(lower[t])) lower[t]=lo[t];
                        if (upper[t].below(hi[t])) upper[t]=hi[t];
                        positive[t]=positive[t]&&valid[t];
                    }
                }
            }
        }
        if (p) for (int t=0; t<(original ? 2 : 1); ++t) {
            if (original) std::cout<<"shift "<<(t==0 ? "23/2" : "58/5")<<" ";
            if (positive[t]) std::cout<<"min "<<lower[t].n<<"/"<<lower[t].d
                                      <<" max "<<upper[t].n<<"/"<<upper[t].d<<"\n";
            else std::cout<<"unavailable\n";
        }
    }
};

int main(int argc,char** argv) {
    if (argc != 2) return 2;
    std::ifstream file(argv[1],std::ios::binary);
    std::array<U,5> hdr{};
    file.read(reinterpret_cast<char*>(hdr.data()),sizeof(hdr));
    if (!file || hdr[0]!=0x5354524950373301ULL || hdr[1]!=16 || hdr[3]!=540 || hdr[4]>8)
        throw std::runtime_error("unsupported little-endian checkpoint header");
    Check check(hdr[1]);
    auto& v=check.values[1];
    if (v.size()!=hdr[2]) throw std::runtime_error("size mismatch");
    file.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(U));
    if (!file || file.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("payload size");
    U largest=0;
    for (U& entry:v) {
        entry >>= hdr[4];
        if (entry==0 || entry>scale) throw std::runtime_error("inadmissible coordinate");
        largest=std::max(largest,entry);
    }
    std::cout<<"width "<<2*hdr[1]<<" round "<<hdr[3]<<" odd "<<v.size()
             <<" even "<<check.values[0].size()<<" root "<<v.front()<<" largest "<<largest<<"\n";
    std::vector<U> original;
    if (argc==3) original=v;
    check.sweep(0);
    check.sweep(1);
    if (argc==3) {
        check.sweep(0);
        check.sweep(1,&original);
    }
}
