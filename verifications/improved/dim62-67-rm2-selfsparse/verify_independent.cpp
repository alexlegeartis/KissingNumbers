#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using W=uint64_t;
int main(int argc,char** argv){
 std::ifstream f(argc>1?argv[1]:"witness.txt");
 std::vector<unsigned> lab; unsigned a;
 while(f>>a)lab.push_back(a);
 assert(f.eof()&&lab.size()==80&&std::set<unsigned>(lab.begin(),lab.end()).size()==80&&*std::max_element(lab.begin(),lab.end())<(1U<<20));
 std::vector<uint32_t> low(1,0),aff(1,0);
 for(int m=0;m<32;m++){
  unsigned w=__builtin_popcount(m);if(w>2)continue;
  uint32_t s=0;for(int x=0;x<32;x++)if((x&m)==m)s|=uint32_t(1)<<x;
  auto add=[&](std::vector<uint32_t>& z){size_t n=z.size();for(size_t i=0;i<n;i++)z.push_back(z[i]^s);};
  add(low);if(w<=1)add(aff);
 }
 assert(low.size()==65536&&aff.size()==64);
 assert(std::set<uint32_t>(low.begin(),low.end()).size()==65536);
 assert(std::set<uint32_t>(aff.begin(),aff.end()).size()==64);
 std::vector<W> q;
 for(auto v:aff)for(auto u:low)q.push_back(W(u)|(W(u^v)<<32));
 std::set<unsigned> D;
 for(auto u:lab)for(auto v:lab)D.insert(u^v);
 for(auto d:D){
  unsigned masks[20],i=0;for(int m=0;m<64;m++)if(__builtin_popcount(m)==3)masks[i++]=m;
  W x=0;for(unsigned m=0;m<64;m++){
   int v=0;
   for(unsigned b=0;b<20;b++)if(((d>>b)&1)&&((m&masks[b])==masks[b]))v^=1;
   x|=W(v)<<m;
  }
  int best=64;
  for(W w:q){ if(!d&&!w)continue;best=std::min(best,__builtin_popcountll(x^w)); }
  assert(best==16);
 }
 std::cout<<"PASS independent Plotkin-split RM recursion, "<<D.size()<<" cubic differences, min16\n";
}
