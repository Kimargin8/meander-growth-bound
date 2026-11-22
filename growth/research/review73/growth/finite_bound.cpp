// Quantitative finite-n bound from one saved positive witness, checked against
// directed lower approximations to every coordinate of A^r e_root.
#define main original_confined_main
#include "../../review57/strip/confined.cpp"
#undef main
#include <fstream>
using U=std::uint64_t;
using D=unsigned __int128;
struct Ratio{U a=~U(0),b=1;Index index=0;};
static bool before(Ratio l,Ratio r){return D(l.a)*r.b<D(r.a)*l.b;}
struct Interval{
    Ratio lo,hi{0,1,0};
    U min_v=~U(0),max_v=0;
    void add(Ratio r){
        if(before(r,lo))lo=r;
        if(before(hi,r))hi=r;
        min_v=std::min(min_v,r.b);max_v=std::max(max_v,r.b);
    }
    void merge(const Interval& other){
        if(before(other.lo,lo))lo=other.lo;
        if(before(hi,other.hi))hi=other.hi;
        min_v=std::min(min_v,other.min_v);max_v=std::max(max_v,other.max_v);
    }
};

int main(int argc,char**argv){
    if(argc!=3)return 1;
    std::ifstream file(argv[1],std::ios::binary);
    int end=std::stoi(argv[2]);
    std::array<U,5> h{};
    file.read(reinterpret_cast<char*>(h.data()),sizeof(h));
    if(!file || h[0]!=0x5354524950373301ULL || h[1]<1 || h[1]>16 ||
       h[3]>100000 || h[4]>8 || end<1 || end>120)return 2;
    const int cap=h[1];
    constexpr U scale=U(1)<<52;
    std::vector<U> witness(h[2]);
    file.read(reinterpret_cast<char*>(witness.data()),std::streamsize(witness.size()*sizeof(U)));
    if(!file||file.peek()!=std::char_traits<char>::eof())throw std::runtime_error("checkpoint size");
    U wmin=~U(0),wmax=0;
    #pragma omp parallel for schedule(static) reduction(min:wmin) reduction(max:wmax)
    for(std::size_t i=0;i<witness.size();++i){
        witness[i]>>=h[4];
        wmin=std::min(wmin,witness[i]);wmax=std::max(wmax,witness[i]);
    }
    if(!wmin||wmax>scale)throw std::runtime_error("witness range");
    std::cout<<"witness width "<<2*cap<<" round "<<h[3]<<" odd "<<h[2]<<" shift "<<h[4]
             <<" root "<<witness.front()<<" minimum "<<wmin<<" maximum "<<wmax<<std::endl;
    std::vector<Layer> ls(cap+1);
    std::vector<U> offsets(cap+1);
    U cells=0,odd=0;
    for(int k=1;k<=cap;++k){
        ls[k].generate(k,cells);offsets[k]=odd;
        cells+=ls[k].words.size()*(k+1);odd+=ls[k].words.size()*((k+1)/2);
    }
    if(odd!=witness.size()||cells>=std::numeric_limits<Index>::max())return 3;
    for(int k=1;k<=cap;++k)make_edges(ls,k);
    std::vector<U> values(cells);
    const Index root=ls[1].cell(0,1);
    values[root]=1;
    auto sweep=[&](int parity,int shift){
        U max_value=0;
        for(int k=1;k<=cap;++k){
            const Layer& layer=ls[k];
            #pragma omp parallel for schedule(static) reduction(max:max_value)
            for(Index r=0;r<layer.words.size();++r){
                auto i=layer.cell(r,0);
                auto p=layer.words[r];
                for(int u=parity;u<=k;u+=2){
                    U s=u?(values[i+u-1]>>shift):0;
                    s+=values[u<k?i+u+1:layer.cell(layer.reflection[r],u-1)]>>shift;
                    if(k<cap){
                        Word mask=(1u<<u)-1;
                        Word inserted=(p&mask)|(1u<<u)|((p>>u)<<(u+2));
                        s+=values[ls[k+1].cell(ls[k+1].rank(inserted),u+1)]>>shift;
                    }
                    if(u){
                        Index d=layer.down[std::size_t(r)*k+u-1].get();
                        if(d!=missing)s+=values[ls[k-1].cell(d,u-1)]>>shift;
                    }
                    values[i+u]=s;max_value=std::max(max_value,s);
                }
            }
        }
        return max_value;
    };
    // Exactly two shifted candidates, fixed before seeing the final data.
    // The small-width settings exercise the same integer comparison in tests.
    const std::array<U,3> mult{1,2,5};
    const std::array<U,3> subtract=cap==16?std::array<U,3>{0,23,58}:std::array<U,3>{0,1,8};
    std::vector<U> image(odd);
    auto copy_odd=[&](bool into_board){
        for(int k=1;k<=cap;++k){
            const Layer& l=ls[k];
            #pragma omp parallel for schedule(static)
            for(Index r=0;r<l.words.size();++r)
                for(int u=1;u<=k;u+=2){
                    auto j=offsets[k]+std::size_t(r)*((k+1)/2)+u/2;
                    if(into_board)values[l.cell(r,u)]=witness[j];
                    else image[j]=values[l.cell(r,u)];
                }
        }
    };
    copy_odd(true);sweep(0,0);sweep(1,0);copy_odd(false);sweep(0,0);sweep(1,0);
    std::array<Interval,3> intervals;
    auto den=[&](int h,U i){
        if(!h)return witness[i];
        U a=mult[h]*image[i],b=subtract[h]*witness[i];
        if(a<=b)throw std::runtime_error("nonpositive shifted witness");
        return a-b;
    };
    for(int k=1;k<=cap;++k){
        const auto& l=ls[k];
        #pragma omp parallel
        {
            std::array<Interval,3> tmp;
            #pragma omp for schedule(static)
            for(Index r=0;r<l.words.size();++r)for(int u=1;u<=k;u+=2){
                auto j=offsets[k]+std::size_t(r)*((k+1)/2)+u/2;
                auto id=l.cell(r,u);
                for(int h=0;h<3;++h){
                    U top=image[j];
                    if(h){
                        U a=mult[h]*values[id],b=subtract[h]*image[j];
                        if(a<=b)std::abort();
                        top=a-b;
                    }
                    tmp[h].add({top,den(h,j),id});
                }
            }
            #pragma omp critical
            {for(int t=0;t<3;++t)intervals[t].merge(tmp[t]);}
        }
    }
    for(int h=0;h<3;++h){
        auto& t=intervals[h];
        std::cout<<"candidate "<<mult[h]<<" "<<subtract[h]<<" lower "<<t.lo.a<<"/"<<t.lo.b
                 <<" upper "<<t.hi.a<<"/"<<t.hi.b<<" root "<<den(h,0)
                 <<" minimum "<<t.min_v<<" maximum "<<t.max_v<<std::endl;
    }
    std::fill(values.begin(),values.end(),0);values[root]=1;
    U max_value=1;int exponent=0;
    for(int m=0;m<=end;++m){
        std::cout<<"level "<<m<<" exponent "<<exponent<<" root "<<values[root]
                 <<" maximum "<<max_value<<'\n';
        if(m==cap-1 || m==end){
            std::array<Ratio,3> best;
            for(int k=1;k<=cap;++k){
                const auto& l=ls[k];
                #pragma omp parallel
                {
                    std::array<Ratio,3> local;
                    #pragma omp for schedule(static)
                    for(Index r=0;r<l.words.size();++r){
                        for(int u=1;u<=k;u+=2){
                            Index j=l.cell(r,u);
                            auto pos=offsets[k]+std::size_t(r)*((k+1)/2)+u/2;
                            for(int h=0;h<3;++h){
                                Ratio q{values[j],den(h,pos),j};
                                if(before(q,local[h]))local[h]=q;
                            }
                        }
                    }
                    #pragma omp critical
                    {for(int t=0;t<3;++t)if(before(local[t],best[t]))best[t]=local[t];}
                }
            }
            for(int h=0;h<3;++h){
                const auto& t=best[h];
                if(!t.a)throw std::runtime_error("zero lower coordinate after return horizon");
                std::cout<<"domination "<<m<<" "<<mult[h]<<" "<<subtract[h]<<" "<<t.a<<"/"<<t.b<<std::endl;
            }
        }
        if(m==end)break;
        // Before the first 15 jobs, at most 16^15 paths in a coordinate.
        // Thereafter input <= 2^52 and each whole step adds at most 16 terms.
        int shift=0;
        if(m>=std::min(15,cap-1))while((max_value>>shift)>scale)++shift;
        exponent+=shift;
        sweep(0,shift);max_value=sweep(1,0);
    }
    return 0;
}
