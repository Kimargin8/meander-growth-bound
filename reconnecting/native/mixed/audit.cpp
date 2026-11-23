constexpr int WANTED=32, DENOM=4096, even_preserved=5;
constexpr unsigned one[]={582,593,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,595,593,582};
constexpr unsigned two[]={63,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,63};
constexpr unsigned three[]={83,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,83};
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
    std::vector<U> first,bare;

    void compare_previous(const std::vector<U>& normalized){
        const char* name=std::getenv("LRWR_GROWTH77_PREVIOUS");
        if(!name)throw std::runtime_error("missing old witness");
        std::ifstream input(name,std::ios::binary);
        std::array<U,5> h;input.read(reinterpret_cast<char*>(h.data()),sizeof(h));
        if(!input||h[0]!=0x5354524950373301ULL||h[1]!=16||h[2]!=normalized.size()||h[3]!=500||h[4]>8)
            throw std::runtime_error("old-witness header");
        Ratio lower{~U(0),1};U at=0;
        std::vector<U> buffer(1<<19);
        for(U i=0;i<normalized.size();i+=buffer.size()){
            U count=std::min<U>(buffer.size(),normalized.size()-i);
            input.read(reinterpret_cast<char*>(buffer.data()),count*sizeof(U));
            if(!input)throw std::runtime_error("old-witness data");
            for(U j=0;j<count;++j){
                U old=buffer[j]>>h[4];
                if(!old||old>scale)throw std::runtime_error("old-witness value");
                Ratio q{old,normalized[i+j]};
                if(q.below(lower)){lower=q;at=i+j;}
            }
        }
        if(input.peek()!=std::char_traits<char>::eof())throw std::runtime_error("old-witness size");
        std::cout<<"previous over current "<<lower.n<<"/"<<lower.d<<" at "<<at<<"\n";
    }

    Word reversed(Word x) const {
        x=~x;
        x=((x&0x55555555U)<<1)|((x>>1)&0x55555555U);
        x=((x&0x33333333U)<<2)|((x>>2)&0x33333333U);
        x=((x&0x0f0f0f0fU)<<4)|((x>>4)&0x0f0f0f0fU);
        return __builtin_bswap32(x)>>(32-2*cap);
    }

    void partners(Word word,std::array<int,32>& mate) const {
        int stack[16],top=0;
        for(int i=0;i<2*cap;++i){
            if((word>>i)&1)stack[top++]=i;
            else{int j=stack[--top];mate[i]=j;mate[j]=i;}
        }
    }

    // Insert the two ends physically at width 2*cap+2, splice there, and
    // only then return to the bank-folded index of the smaller word.
    Word outside_target(Word word,const std::array<int,32>& mate,int u,int d) const {
        const int j=u+1+d, left=j-1, right=j;
        auto old=[&](int i){return i<u?i:i-2;};
        auto moved=[&](int i){return i<u?i:i+2;};
        const int a=moved(mate[old(left)]), b=moved(mate[old(right)]);
        if(a==right)return std::numeric_limits<Word>::max();
        U lo=(U(1)<<u)-1;
        U wide=(word&lo)|(U(1)<<u)|((U(word)>>u)<<(u+2));
        wide=(wide|(U(1)<<std::min(a,b)))&~(U(1)<<std::max(a,b));
        U mask=(U(1)<<left)-1;
        Word reduced=Word((wide&mask)|((wide>>(right+1))<<left));
        if(u+d>cap)reduced=reversed(reduced);
        return layers[cap].rank(reduced);
    }

    void audit_top(int p,Ratio& minimum,Ratio& maximum){
        const auto& l=layers[cap];const auto& small=layers[cap-1];
        const int width=2*cap;const std::size_t stride=cap/2+1;
        auto& dst=values[p];const auto& src=values[1-p];
        first.resize(l.words.size()*stride);
        bare.resize(l.words.size()*stride);
        const auto none=std::numeric_limits<Word>::max();
        #pragma omp parallel for schedule(static)
        for(Word r=0;r<l.words.size();++r){
            const Word word=l.words[r],mirror=l.rank(reversed(word));
            std::array<int,32> mate;partners(word,mate);
            for(int u=p;u<=cap;u+=2){
                U left=u?src[l.index(1-p,r,u-1)]:0;
                U right=u<cap?src[l.index(1-p,r,u+1)]:src[l.index(1-p,mirror,u-1)];
                Wide extra=Wide(right)*one[u]+(u?Wide(left)*one[u-1]:0);
                U plain=left+right;
                if(u && mate[u-1]!=u){
                    int a=mate[u-1],b=mate[u];
                    Word edit=(word|Word(1)<<std::min(a,b))&~(Word(1)<<std::max(a,b));
                    Word mask=(Word(1)<<(u-1))-1;
                    Word out=(edit&mask)|((edit>>(u+1))<<(u-1));
                    plain+=src[small.index(1-p,small.rank(out),u-1)];
                }
                for(int d:{-3,3}){
                    int g=u+d;if(g<0 || g>width)continue;
                    Word dest=outside_target(word,mate,u,d);
                    if(dest!=none)extra+=Wide(three[std::min(g,u)])*src[l.index(1-p,dest,std::min(g,width-g))];
                }
                auto address=std::size_t(r)*stride+u/2;
                bare[address]=plain;
                first[address]=plain+U(extra>>12);
            }
        }
        Ratio domination{1,1};
        #pragma omp parallel
        {
            Ratio lo{~U(0),1},hi{0,1},ordinary{1,1};
            #pragma omp for schedule(static)
            for(Word r=0;r<l.words.size();++r){
                Word word=l.words[r],mirror=l.rank(reversed(word));
                std::array<int,32> mate;partners(word,mate);
                for(int u=p;u<=cap;u+=2){
                    Wide extra=0;
                    for(int d:{-2,2}){
                        int g=u+d;if(g<0 || g>width)continue;
                        auto slot=std::min(g,width-g)/2;
                        extra+=Wide(even_preserved)*first[std::size_t(g>cap?mirror:r)*stride+slot];
                        Word dest=outside_target(word,mate,u,d);
                        if(dest!=none)extra+=Wide(two[std::min(g,u)])*first[std::size_t(dest)*stride+slot];
                    }
                    auto idx=l.index(p,r,u);
                    auto address=std::size_t(r)*stride+u/2;
                    U val=first[address]+U(extra>>12);
                    Ratio retained{bare[address],val};
                    if(retained.below(ordinary))ordinary=retained;
                    if(p){
                        Ratio ratio{val,dst[idx]};
                        if(ratio.below(lo))lo=ratio;
                        if(hi.below(ratio))hi=ratio;
                    }
                    dst[idx]=val;
                }
            }
            #pragma omp critical
            {
                if(lo.below(minimum))minimum=lo;
                if(maximum.below(hi))maximum=hi;
                if(ordinary.below(domination))domination=ordinary;
            }
        }
        std::cout<<"ordinary fraction parity "<<p<<" "<<domination.n<<"/"<<domination.d<<"\n";
    }

    void sweep(int p, const std::vector<U>* original=nullptr) {
        auto& dest=values[p];
        const auto& src=values[1-p];
        std::array<Ratio,2> lower{{{~U(0),1},{~U(0),1}}}, upper{{{0,1},{0,1}}};
        std::array<bool,2> positive{true,true};
        for (int k=1; k<=cap; ++k) {
            if(k==cap){audit_top(p,lower[0],upper[0]);continue;}
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
                        if (u) z+=src[l.index(1-p,r,u-1)];
                        if (u<k) z+=src[l.index(1-p,r,u+1)];
                        else z+=src[l.index(1-p,l.rank(reflected),u-1)];
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
    if (argc!=2) return 2;
    std::ifstream file(argv[1],std::ios::binary);
    std::array<U,5> hdr{};
    file.read(reinterpret_cast<char*>(hdr.data()),sizeof(hdr));
    if (!file || hdr[0]!=0x5354524950373301ULL || 2*hdr[1]!=WANTED || hdr[4]>8)
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
    if(WANTED==32)check.compare_previous(v);
    std::vector<U> original;
    if (argc==3) original=v;
    check.sweep(0);
    check.sweep(1);
    if (argc==3) {
        check.sweep(0);
        check.sweep(1,&original);
    }
}
