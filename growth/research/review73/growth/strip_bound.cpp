// A finite, fully labelled Jensen transfer and an integer Collatz certificate.
// Each nonempty bank split is retained, without quotient by bank exchange.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

using U=std::uint64_t;
using D=unsigned __int128;
using I=std::uint32_t;
using Row=std::array<I,4>;
constexpr I none=std::numeric_limits<I>::max();

static void dycks(int width,int i,int height,I word,std::vector<I>& out){
    if(i==width){out.push_back(word);return;}
    if(width-i>height)dycks(width,i+1,height+1,word|(I(1)<<i),out);
    if(height)dycks(width,i+1,height-1,word,out);
}

struct Transfer{
    std::vector<Row> graph;
    std::vector<I> active;
    I root;
    explicit Transfer(int W){
        std::vector<std::vector<I>> patterns(W+1);
        std::vector<std::unordered_map<I,I>> index(W+1);
        U size=0;
        for(int w=2;w<=W;w+=2){
            dycks(w,0,0,0,patterns[w]);
            for(I p:patterns[w]){
                index[w].emplace(p,size);
                size+=w+1;
            }
        }
        if(size>=none)throw std::overflow_error("state index");
        graph.resize(size);
        root=index[2].at(1)+1;
        for(int w=2;w<=W;w+=2)for(I p:patterns[w]){
            int partner[32],stack[16],top=0;
            for(int i=0;i<w;++i){
                if((p>>i)&1)stack[top++]=i;
                else{int j=stack[--top];partner[j]=i;partner[i]=j;}
            }
            const I base=index[w].at(p);
            for(int u=0;u<=w;++u){
                const I low=(I(1)<<u)-1;
                Row edges{none,none,none,none};
                if(u>0)edges[0]=base+u-1;
                if(u<w)edges[1]=base+u+1;
                if(w<W){
                    I q=(p&low)|(I(1)<<u)|((p&~low)<<2);
                    edges[2]=index[w+2].at(q)+u+1;
                }
                if(u>0 && u<w && partner[u-1]!=u){
                    int a=partner[u-1],b=partner[u];
                    I q=(p|(I(1)<<std::min(a,b)))&~(I(1)<<std::max(a,b));
                    I lower=(I(1)<<(u-1))-1;
                    q=(q&lower)|((q>>(u+1))<<(u-1));
                    edges[3]=index[w-2].at(q)+u-1;
                }
                graph[base+u]=edges;
            }
        }
        // Two breadth-first searches keep precisely the root's communicating class.
        std::vector<char> forward(size),backward(size);
        std::vector<I> queue{root};
        forward[root]=true;
        for(std::size_t i=0;i<queue.size();++i){
            for(I j:graph[queue[i]])if(j!=none && !forward[j]){
                forward[j]=true;queue.push_back(j);
            }
        }
        const auto reachable=queue.size();
        std::vector<I> offsets(size+1);
        for(auto& row:graph)for(I j:row)if(j!=none)++offsets[j+1];
        for(I i=0;i<size;++i)offsets[i+1]+=offsets[i];
        std::vector<I> incoming(offsets.back()),cursor=offsets;
        for(I i=0;i<size;++i)for(I j:graph[i])if(j!=none)incoming[cursor[j]++]=i;
        queue.assign(1,root);backward[root]=true;
        for(std::size_t i=0;i<queue.size();++i){
            I j=queue[i];
            for(I h=offsets[j];h<offsets[j+1];++h){
                I k=incoming[h];
                if(!backward[k]){backward[k]=true;queue.push_back(k);}
            }
        }
        for(I i=0;i<size;++i)if(forward[i]&&backward[i])active.push_back(i);
        for(I i:active)for(I& j:graph[i])if(j!=none && !(forward[j]&&backward[j]))j=none;
        std::cout<<"W "<<W<<" states "<<size<<" forward "<<reachable<<" backward "
                 <<queue.size()<<" active "<<active.size()<<'\n';
    }
    void apply(const std::vector<U>& v,std::vector<U>& out)const{
        // All callers keep v<=2^52, v after one step <=4*2^52, so sums fit.
        #pragma omp parallel for schedule(static)
        for(std::size_t h=0;h<active.size();++h){
            I i=active[h];U sum=0;
            for(I j:graph[i])if(j!=none)sum+=v[j];
            out[i]=sum;
        }
    }
};

int main(int argc,char** argv){
    if(argc<3||argc>5)return 1;
    int width=std::stoi(argv[1]),iterations=std::stoi(argv[2]);
    if(width<2||width>26||width%2||iterations<0)return 2;
    Transfer T(width);
    const I N=T.graph.size();
    std::vector<U> x(N),y(N),z(N);
    x[T.root]=1;
    for(int n=1;n<=12;++n){
        std::cout<<"count "<<n<<" "<<x[T.root]<<'\n';
        T.apply(x,y);T.apply(y,z);x.swap(z);
    }
    constexpr U scale=U(1)<<52;
    for(I i:T.active)x[i]=scale;
    for(int round=0;round<=iterations;++round){
        T.apply(x,y);T.apply(y,z);
        if(round%10==0 || round==iterations){
            U mn=std::numeric_limits<U>::max(),md=1,mx=0,xd=1;
            I imin=0,imax=0;
            for(I i:T.active){
                if(x[i]==0)throw std::runtime_error("zero witness entry");
                if(D(z[i])*md<D(mn)*x[i]){mn=z[i];md=x[i];imin=i;}
                if(D(z[i])*xd>D(mx)*x[i]){mx=z[i];xd=x[i];imax=i;}
            }
            std::cout<<"round "<<round<<" min "<<mn<<"/"<<md<<" node "<<imin
                     <<" max "<<mx<<"/"<<xd<<" node "<<imax<<std::endl;
        }
        if(round==iterations)break;
        U largest=*std::max_element(z.begin(),z.end());
        const U gain=(D(scale)<<32)/largest;
        for(I i:T.active)x[i]=(D(z[i])*gain)>>32;
    }
    if(argc>=4){
        std::ofstream out(argv[3],std::ios::binary);
        out.write(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(U));
        if(!out)throw std::runtime_error("witness write");
    }
}
