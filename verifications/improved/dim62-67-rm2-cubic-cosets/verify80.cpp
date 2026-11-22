#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

using u64=uint64_t;
using u32=uint32_t;
u64 monomial(int mask){
  u64 out=0;
  for(int x=0;x<64;x++)if((x&mask)==mask)out|=1ull<<x;
  return out;
}
int main(int argc,char**argv){
  std::ifstream file(argc==2?argv[1]:"witness.txt");
  std::vector<u32> labels;
  for(u32 v;file>>v;)labels.push_back(v);
  assert(file.eof()&&labels.size()>=80);
  assert(std::set<u32>(labels.begin(),labels.end()).size()==labels.size());
  std::vector<u64> cube,base;
  for(int i=0;i<64;i++) {
    int w=__builtin_popcount(i);
    if(w<=2)base.push_back(monomial(i));
    if(w==3)cube.push_back(monomial(i));
  }
  assert(base.size()==22&&cube.size()==20);
  // Full rank of the degree<=2 monomials by exact binary elimination.
  std::array<u64,64> pivot{};
  for(u64 x:base){
    while(x&&pivot[63-__builtin_clzll(x)])x^=pivot[63-__builtin_clzll(x)];
    assert(x); pivot[63-__builtin_clzll(x)]=x;
  }
  std::set<u32> differences{0};
  for(u32 a:labels){assert(a<(1u<<20));for(u32 b:labels)differences.insert(a^b);}
  std::vector<u64> code(1u<<22);
  for(u32 i=1;i<code.size();i++)code[i]=code[i-1]^base[__builtin_ctz(i)];
  std::array<unsigned,65> hist{};
  for(u32 d:differences){
    u64 x=0;
    for(u32 y=d;y;y&=y-1)x^=cube[__builtin_ctz(y)];
    int md=65;
    // Zero difference is also checked against every nonzero base codeword.
    for(u32 i=(d==0);i<code.size();i++)md=std::min(md,__builtin_popcountll(code[i]^x));
    assert(md>=16);++hist[md];
  }
  std::cout<<"verified RM(2,6) rank "<<base.size()<<", cubic representatives "<<labels.size()
           <<", pair-difference types incl self "<<differences.size()<<"\n";
  for(int i=16;i<=64;i++)if(hist[i])std::cout<<"minimum distance "<<i<<": "<<hist[i]<<" difference types\n";
  std::cout<<"dense 64D vectors "<<(u64(labels.size())<<22)<<"\n";
}
