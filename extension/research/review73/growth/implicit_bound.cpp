// Single-board bank-folded certificate, without storing a four-edge row per cell.
// An optional private checkpoint contains only the odd slit coordinates.
#define main original_confined_main
#include "../../review57/strip/confined.cpp"
#undef main
#include <cstdio>
#include <fstream>
using U=std::uint64_t;
using D=unsigned __int128;
constexpr U scale=U(1)<<52;
struct Ratio{
    U a,b;Index row;
    bool less(Ratio r)const{return D(a)*r.b<D(r.a)*b;}
};
struct Stats{
    Ratio lo{~U(0),1,0},hi{0,1,0};
    U top=0,positive=0;
    void merge(const Stats& x){
        if(x.lo.less(lo))lo=x.lo;
        if(hi.less(x.hi))hi=x.hi;
        top=std::max(top,x.top);positive+=x.positive;
    }
};
struct Board{
    int cap;
    std::vector<Layer> layers;
    std::vector<U> v;
    U odd=0;
    Index root;
    explicit Board(int W):cap(W/2),layers(cap+1){
        U total=0;
        for(int k=1;k<=cap;++k){
            layers[k].generate(k,total);
            total+=layers[k].words.size()*(k+1);
            odd+=layers[k].words.size()*((k+1)/2);
            if(total>=std::numeric_limits<Index>::max())throw std::overflow_error("cells");
        }
        for(int k=1;k<=cap;++k)make_edges(layers,k);
        v.resize(total);
        root=layers[1].cell(0,1);
        std::cout<<"width "<<W<<" cells "<<total<<" root "<<root<<" odd "<<odd
                 <<" even "<<total-odd<<std::endl;
    }
    template<bool boolean,bool measured>
    Stats sweep(int parity,int readshift,int oldshift){
        Stats all;
        for(int k=1;k<=cap;++k){
            const auto& layer=layers[k];const auto& small=layers[k-1];
            const auto& big=layers[std::min(k+1,cap)];
            #pragma omp parallel
            {
                Stats local;
                #pragma omp for schedule(static)
                for(Index r=0;r<layer.words.size();++r){
                    Word word=layer.words[r];Index start=layer.cell(r,0);
                    for(int u=parity;u<=k;u+=2){
                        U z=u?(v[start+u-1]>>readshift):0;
                        auto side=u<k?start+u+1:layer.cell(layer.reflection[r],u-1);
                        z+=v[side]>>readshift;
                        if(k<cap){
                            Word mask=(1u<<u)-1;
                            Word opened=(word&mask)|(1u<<u)|((word>>u)<<(u+2));
                            z+=v[big.cell(big.rank(opened),u+1)]>>readshift;
                        }
                        if(u){
                            Index d=layer.down[std::size_t(r)*k+u-1].get();
                            if(d!=missing)z+=v[small.cell(d,u-1)]>>readshift;
                        }
                        const Index id=start+u;
                        if constexpr(measured){
                            U old=v[id]>>oldshift;
                            if(!z || !old)std::abort();
                            Ratio ratio{z,old,id};
                            if(ratio.less(local.lo))local.lo=ratio;
                            if(local.hi.less(ratio))local.hi=ratio;
                        }
                        v[id]=boolean?(z>0):z;
                        if constexpr(boolean)local.positive+=z>0;
                        else local.top=std::max(local.top,z);
                    }
                }
                #pragma omp critical
                all.merge(local);
            }
        }
        return all;
    }
    void checkpoint(const std::string& path,bool read,int& round,int& shift){
        if(path.empty())return;
        std::fstream io(read?path:path+".pending",std::ios::binary|
                        (read?std::ios::in:(std::ios::out|std::ios::trunc)));
        if(!io)throw std::runtime_error("checkpoint open");
        std::array<U,5> header{0x5354524950373301ULL,U(cap),odd,U(round),U(shift)};
        if(read){
            std::array<U,5> seen{};io.read(reinterpret_cast<char*>(seen.data()),sizeof(seen));
            if(!io || !std::equal(header.begin(),header.begin()+3,seen.begin())
               || seen[3]>100000 || seen[4]>8)throw std::runtime_error("checkpoint header");
            round=seen[3];shift=seen[4];
        }else io.write(reinterpret_cast<const char*>(header.data()),sizeof(header));
        std::vector<U> buffer(1<<16);
        std::size_t pos=0;
        for(int k=1;k<=cap;++k){
            const auto& l=layers[k];
            for(Index r=0;r<l.words.size();++r)for(int u=1;u<=k;u+=2){
                U& value=v[l.cell(r,u)];
                auto j=pos%buffer.size();
                if(read){
                    if(!j){
                        auto count=std::min<U>(odd-pos,buffer.size());
                        io.read(reinterpret_cast<char*>(buffer.data()),sizeof(U)*count);
                        if(!io)throw std::runtime_error("short checkpoint");
                    }
                    value=buffer[j];
                    if(!value || value>16*scale)throw std::runtime_error("invalid checkpoint");
                }else{
                    buffer[j]=value;
                    if(j+1==buffer.size() || pos+1==odd)
                        io.write(reinterpret_cast<const char*>(buffer.data()),sizeof(U)*(j+1));
                }
                ++pos;
            }
        }
        if(pos!=odd || !io)throw std::runtime_error("checkpoint mismatch");
        io.close();
        if(!read && std::rename((path+".pending").c_str(),path.c_str()))throw std::runtime_error("rename");
    }
};
int main(int argc,char** argv){
    if(argc<3||argc>5)return 1;
    int width=std::stoi(argv[1]),end=std::stoi(argv[2]);
    if(width<2||width>32||width%2||end<0)return 2;
    std::string path=argc>=4?argv[3]:"";
    bool resume=argc==5 && std::string(argv[4])=="resume";
    Board b(width);
    int round=0,shift=0;
    if(resume)b.checkpoint(path,true,round,shift);
    else{
        b.v[b.root]=1;
        for(int n=1;n<=12;++n){
            std::cout<<"count "<<n<<" "<<b.v[b.root]<<'\n';
            b.sweep<false,false>(0,0,0);b.sweep<false,false>(1,0,0);
        }
        std::fill(b.v.begin(),b.v.end(),0);b.v[b.root]=1;
        for(int p=1;;++p){
            b.sweep<true,false>(0,0,0);
            auto ret=b.sweep<true,false>(1,0,0);
            if(ret.positive==b.odd){
                std::cout<<"all_odd_return_in "<<2*p<<" moves"<<std::endl;break;
            }
            if(p>2*width)throw std::runtime_error("return");
        }
        std::fill(b.v.begin(),b.v.end(),scale);
    }
    for(;round<=end;){
        b.sweep<false,false>(0,shift,0);
        bool check=round%10==0 || round==end;
        auto a=check?b.sweep<false,true>(1,0,shift):b.sweep<false,false>(1,0,0);
        if(check)std::cout<<"round "<<round<<" min "<<a.lo.a<<"/"<<a.lo.b<<" row "
                          <<a.lo.row<<" max "<<a.hi.a<<"/"<<a.hi.b<<" row "<<a.hi.row<<std::endl;
        // v_new is floor(raw / 2^shift), with every entry at most 2^52.
        shift=0;while((a.top>>shift)>scale)++shift;
        ++round;
        if(!path.empty() && round%100==0)b.checkpoint(path,false,round,shift);
    }
    return 0;
}
