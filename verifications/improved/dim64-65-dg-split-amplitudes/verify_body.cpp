#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using U=uint64_t;

U mon(unsigned u){
  U z=0;
  for(unsigned i=0;i<64;++i) if((i&u)==u) z|=U(1)<<i;
  return z;
}

int main(int argc,char**argv){
  assert(argc==2);
  std::ifstream in(argv[1]);
  std::vector<U> quad, cub, C;
  for(unsigned u=0;u<64;++u){
    if(std::popcount(u)<=2) quad.push_back(mon(u));
    if(std::popcount(u)==3) cub.push_back(mon(u));
  }
  assert(quad.size()==22&&cub.size()==20);
  for(unsigned j;in>>j;){
    assert(j<(1u<<20));
    U w=0;
    for(unsigned i=0;i<20;++i) if((j>>i)&1) w^=cub[i];
    C.push_back(w);
  }
  assert(in.eof()&&C.size()==80&&std::set(C.begin(),C.end()).size()==80);
  std::vector<U> Q(1u<<22);
  for(size_t i=1;i<Q.size();++i)
    Q[i]=Q[i-1]^quad[std::countr_zero(i)];
  std::set<U> D;
  for(U a:C)for(U b:C)D.insert(a^b);
  // An exhaustive minimum-distance proof for the actual cubic palette.
  for(U d:D){
    int m=64;
    for(size_t i=(d==0);i<Q.size();++i)m=std::min(m,std::popcount(d^Q[i]));
    assert(m>=16);
  }
  std::clog<<"PASS parent minimum 16 (checked "<<D.size()<<" cubic differences)\n";
  std::vector<U>B;
  for(U c:C)for(U q:Q)if(std::popcount(c^q)==16)B.push_back(c^q);
  std::sort(B.begin(),B.end());
  assert(B.size()==37996&&std::adjacent_find(B.begin(),B.end())==B.end());
  for(U z:B)std::cout<<z<<"\n";
}
