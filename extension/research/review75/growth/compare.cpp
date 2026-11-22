// Coordinatewise comparison between the two normalized stored odd witnesses.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using U=std::uint64_t;
using Wide=unsigned __int128;
int main(int argc,char** argv){
    if(argc!=3)return 2;
    std::ifstream a(argv[1],std::ios::binary), b(argv[2],std::ios::binary);
    std::array<U,5> ha{},hb{};
    a.read(reinterpret_cast<char*>(ha.data()),sizeof(ha));
    b.read(reinterpret_cast<char*>(hb.data()),sizeof(hb));
    if(ha[0]!=0x5354524950373301ULL || ha[1]!=16 || ha[2]!=386058577 || ha[3]!=300
       || hb[3]!=440 || !std::equal(ha.begin(),ha.begin()+3,hb.begin()) || ha[4]>8 || hb[4]>8)
        throw std::runtime_error("unexpected witnesses");
    std::vector<U> old(1<<17), fresh(old.size());
    U numerator=~U(0),denominator=1, at=0;
    for(U offset=0;offset<ha[2];offset+=old.size()){
        auto count=std::min<U>(ha[2]-offset,old.size());
        a.read(reinterpret_cast<char*>(old.data()),sizeof(U)*count);
        b.read(reinterpret_cast<char*>(fresh.data()),sizeof(U)*count);
        if(!a || !b)throw std::runtime_error("short witnesses");
        for(U i=0;i<count;++i){
            U o=old[i]>>ha[4], f=fresh[i]>>hb[4];
            if(o==0 || f==0 || std::max(o,f)>(U(1)<<52))throw std::runtime_error("invalid entries");
            if(Wide(o)*denominator<Wide(f)*numerator){numerator=o;denominator=f;at=offset+i;}
        }
        if(offset==0)std::cout<<"roots "<<(old[0]>>ha[4])<<" "<<(fresh[0]>>hb[4])<<"\n";
    }
    if(a.peek()!=std::char_traits<char>::eof() || b.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("trailing payload");
    std::cout<<"rows "<<ha[2]<<" min "<<numerator<<"/"<<denominator<<" index "<<at<<"\n";
}
