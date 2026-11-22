// Exhaust every folded state and decode its outgoing edges as literal labelled
// involutions. This verifies the stored down, reflection, rank and open targets
// used by the scalar implicit/folded transfers, including the top width.
#define main original_confined_main
#include "../../review57/strip/confined.cpp"
#undef main

using U=std::uint64_t;

static void require(bool test){if(!test)std::abort();}

int main(int argc,char** argv){
    int wmax=argc==2?std::stoi(argv[1]):0;
    if(wmax<2||wmax>32||wmax%2)return 2;
    int cap=wmax/2;
    std::vector<Layer> layers(cap+1);
    U total=0,catalan=1;
    for(int k=1;k<=cap;++k){
        layers[k].generate(k,total);
        catalan=catalan*2*(2*k-1)/(k+1);
        require(layers[k].words.size()==catalan);
        total+=catalan*(k+1);
        require(total<std::numeric_limits<Index>::max());
    }
    for(int k=1;k<=cap;++k)make_edges(layers,k);
    U states=0,ranks=0,opens=0,closes=0,forbidden=0,shifts=0;
    for(int k=1;k<=cap;++k){
        const Layer& layer=layers[k];
        U valid=0,blocked=0,added=0,slid=0;
        #pragma omp parallel for schedule(static) reduction(+:valid,blocked,added,slid)
        for(Index r=0;r<layer.words.size();++r){
            Word p=layer.words[r];
            require(layer.rank(p)==r);
            int partner[32],stack[16],h=0;
            for(int i=0;i<2*k;++i){
                if((p>>i)&1){
                    require(h<k);stack[h++]=i;
                }else{
                    require(h>0);
                    int earlier=stack[--h];partner[i]=earlier;partner[earlier]=i;
                }
            }
            require(!h);
            Word reflected=0;
            for(int old=0;old<2*k;++old){
                int new_site=2*k-1-old,other=2*k-1-partner[old];
                if(new_site<other)reflected|=Word(1)<<new_site;
            }
            require(layer.reflection[r]<layer.words.size());
            require(layer.words[layer.reflection[r]]==reflected);
            for(int u=0;u<=k;++u){
                slid+=1+(u>0);
                if(k<cap){
                    const Layer& wider=layers[k+1];
                    Word expect=0;
                    for(int new_site=0;new_site<2*k+2;++new_site){
                        int other;
                        if(new_site==u)other=u+1;
                        else if(new_site==u+1)other=u;
                        else{
                            int old_site=new_site-(new_site>u+1?2:0);
                            int old_other=partner[old_site];
                            other=old_other+(old_other>=u?2:0);
                        }
                        if(new_site<other)expect|=Word(1)<<new_site;
                    }
                    Word mask=(Word(1)<<u)-1;
                    Word encoded=(p&mask)|(Word(1)<<u)|((p>>u)<<(u+2));
                    Index dest=wider.rank(encoded);
                    require(dest<wider.words.size()&&wider.words[dest]==expect);
                    ++added;
                }
                if(!u)continue;
                Index d=layer.down[std::size_t(r)*k+u-1].get();
                int left=u-1,right=u;
                if(partner[left]==right){require(d==missing);++blocked;continue;}
                const Layer& smaller=layers[k-1];
                require(d!=missing&&d<smaller.words.size());
                Word expected=0;
                for(int new_site=0;new_site<2*k-2;++new_site){
                    int old_site=new_site+(new_site>=left?2:0);
                    int old_other=partner[old_site];
                    if(old_other==left)old_other=partner[right];
                    else if(old_other==right)old_other=partner[left];
                    if(old_site<old_other)expected|=Word(1)<<new_site;
                }
                require(smaller.words[d]==expected);++valid;
            }
        }
        states+=layer.words.size()*(k+1);ranks+=layer.words.size();
        closes+=valid;forbidden+=blocked;opens+=added;shifts+=slid;
        std::cout<<"width "<<2*k<<" matchings "<<layer.words.size()<<" shifts "<<slid
                 <<" opens "<<added<<" closes "<<valid<<" early_closures "<<blocked<<std::endl;
    }
    require(states==total);
    std::cout<<"TOTAL states "<<states<<" matchings "<<ranks<<" shifts "<<shifts
             <<" opens "<<opens<<" closes "<<closes<<" early_closures "<<forbidden<<std::endl;
}
