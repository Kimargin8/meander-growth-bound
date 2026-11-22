// Independent bank-folded construction. Uses the validated rank and edge data
// of the confined count. Compile with g++ -O3 -march=native -fopenmp.
#define main original_confined_main
#include "../../review57/strip/confined.cpp"
#undef main
using U=std::uint64_t;
using D=unsigned __int128;
using Row=std::array<Index,5>;

struct Graph{
    std::array<std::vector<Row>,2> rows;
    Index root,size;
    explicit Graph(int width){
        const int cap=width/2;
        std::vector<Layer> layers(cap+1);
        U total=0;
        for(int k=1;k<=cap;++k){
            layers[k].generate(k,total);
            total+=layers[k].words.size()*(k+1);
            if(total>=std::numeric_limits<Index>::max())throw std::overflow_error("cells");
        }
        size=total;
        root=layers[1].cell(0,1);
        for(int k=1;k<=cap;++k)make_edges(layers,k);
        std::size_t counts[2]={};
        for(int k=1;k<=cap;++k)for(int u=0;u<=k;++u)counts[u%2]+=layers[k].words.size();
        for(int p:{0,1})rows[p].reserve(counts[p]);
        for(int k=1;k<=cap;++k){
            const auto& layer=layers[k];
            for(Index r=0;r<layer.words.size();++r){
                const auto pattern=layer.words[r];
                for(int u=0;u<=k;++u){
                    Row row{layer.cell(r,u),size,size,size,size};
                    if(u)row[1]=layer.cell(r,u-1);
                    row[2]=u<k?layer.cell(r,u+1):layer.cell(layer.reflection[r],u-1);
                    if(k<cap){
                        Word mask=(1u<<u)-1;
                        Word opened=(pattern&mask)|(1u<<u)|((pattern>>u)<<(u+2));
                        const auto& wider=layers[k+1];
                        row[3]=wider.cell(wider.rank(opened),u+1);
                    }
                    if(u){
                        Index j=layer.down[std::size_t(r)*k+u-1].get();
                        if(j!=missing)row[4]=layers[k-1].cell(j,u-1);
                    }
                    rows[u%2].push_back(row);
                }
            }
        }
    }
};

struct Ratio{
    U n,d;Index row;
    bool less(const Ratio& r)const{return D(n)*r.d<D(r.n)*d;}
};
struct Stats{
    Ratio lo{~U(0),1,0},hi{0,1,0};
    U biggest=0,positive=0;
    void record(Ratio r){
        if(r.less(lo))lo=r;
        if(hi.less(r))hi=r;
    }
};
Stats sweep(const std::vector<Row>& rows,std::vector<U>& x,bool boolean,bool check){
    Stats answer;
    #pragma omp parallel
    {
        Stats local;
        #pragma omp for schedule(static)
        for(std::size_t i=0;i<rows.size();++i){
            const auto& r=rows[i];
            U sum=x[r[1]]+x[r[2]]+x[r[3]]+x[r[4]];
            if(check){
                if(!x[r[0]] || !sum)std::abort();
                local.record({sum,x[r[0]],r[0]});
            }
            x[r[0]]=boolean?(sum>0):sum;
            local.positive+=sum>0;
            local.biggest=std::max(local.biggest,sum);
        }
        #pragma omp critical
        {
            if(local.lo.less(answer.lo))answer.lo=local.lo;
            if(answer.hi.less(local.hi))answer.hi=local.hi;
            answer.positive+=local.positive;
            answer.biggest=std::max(answer.biggest,local.biggest);
        }
    }
    return answer;
}

int main(int argc,char** argv){
    if(argc!=3)return 1;
    const int W=std::stoi(argv[1]),rounds=std::stoi(argv[2]);
    if(W<2||W>30||W%2||rounds<0)return 2;
    Graph g(W);
    std::vector<U> x(std::size_t(g.size)+1);
    std::cout<<"width "<<W<<" cells "<<g.size<<" root "<<g.root
             <<" odd "<<g.rows[1].size()<<" even "<<g.rows[0].size()<<std::endl;
    x[g.root]=1;
    for(int n=1;n<=12;++n){
        std::cout<<"count "<<n<<" "<<x[g.root]<<'\n';
        sweep(g.rows[0],x,false,false);sweep(g.rows[1],x,false,false);
    }
    std::fill(x.begin(),x.end(),0);x[g.root]=1;
    for(int m=1;;++m){
        sweep(g.rows[0],x,true,false);
        Stats ret=sweep(g.rows[1],x,true,false);
        if(ret.positive==g.rows[1].size()){
            std::cout<<"all_odd_return_in "<<2*m<<" moves"<<std::endl;break;
        }
        if(m>2*W)throw std::runtime_error("could not certify return paths");
    }
    constexpr U scale=U(1)<<52;
    for(const auto& row:g.rows[1])x[row[0]]=scale;
    for(int round=0;round<=rounds;++round){
        sweep(g.rows[0],x,false,false);
        bool check=round%10==0 || round==rounds;
        Stats a=sweep(g.rows[1],x,false,check);
        if(check)std::cout<<"round "<<round<<" min "<<a.lo.n<<"/"<<a.lo.d<<" row "<<a.lo.row
                          <<" max "<<a.hi.n<<"/"<<a.hi.d<<" row "<<a.hi.row<<std::endl;
        if(round==rounds)break;
        U gain=(D(scale)<<32)/a.biggest;
        #pragma omp parallel for schedule(static)
        for(std::size_t i=0;i<g.rows[1].size();++i){
            Index j=g.rows[1][i][0];x[j]=(D(x[j])*gain)>>32;
        }
    }
    return 0;
}
