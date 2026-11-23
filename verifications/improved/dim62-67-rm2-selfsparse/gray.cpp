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
  std::ifstream file(argc>=2?argv[1]:"witness.txt");
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
  std::vector<u64> supports;
  for(u32 d:labels){
    u64 x=0;for(u32 y=d;y;y&=y-1)x^=cube[__builtin_ctz(y)];
    for(u64 v:code)if(__builtin_popcountll(v^x)==16)supports.push_back(v^x);
  }
  std::sort(supports.begin(),supports.end());
  assert(supports.size()==37996&&std::adjacent_find(supports.begin(),supports.end())==supports.end());
  // The same exhaustive minimum-distance check proves distinct supports intersect in <=8.
  int incidence[64]={},pair[64][64]={},qb[64]={},qp[64][64]={};
  auto add=[&](u64 u,int* I,int (*P)[64]){
    for(u64 s=u;s;s&=s-1){int i=__builtin_ctzll(s);++I[i];
      for(u64 t=s&(s-1);t;t&=t-1)++P[i][__builtin_ctzll(t)];}
  };
  for(u64 s:supports)add(s,incidence,pair); 
  std::set<u64> blocks;
  for(int a=0;a<64;a++)for(int b=a+1;b<64;b++)for(int c=b+1;c<64;c++)
    blocks.insert((1ull<<a)|(1ull<<b)|(1ull<<c)|(1ull<<(a^b^c)));
  assert(blocks.size()==10416);
  for(u64 s:blocks)add(s,qb,qp);
  std::array<int,32> dual{}; int parity[5]={65535,0,0,0,0};
  for(int i=0;i<16;i++)for(int j=0;j<4;j++)if(i>>j&1)parity[j+1]|=1<<i;
  for(int i=1;i<32;i++)dual[i]=dual[i-1]^parity[__builtin_ctz(i)];
  std::set<int> ds(dual.begin(),dual.end());
  assert(ds.size()==32);
  for(int w:ds){int d=__builtin_popcount(w);assert(d==0||d==8||d==16);}
  int signs=0;
  for(int w=0;w<65536;w++){
    bool legal=true;
    for(int h:parity)if(__builtin_popcount(w&h)%2){legal=false;break;}
    if(legal){++signs;assert(w==0||__builtin_popcount(w)>=4);}
  }
  assert(signs==2048);
  int max63=0,I=0,J=0,U=0,V=0;
  for(int i=0;i<64;i++)for(int j=i+1;j<64;j++){
    int both=pair[i][j],neither=int(supports.size())-incidence[i]-incidence[j]+both;
    int ub=qp[i][j],un=int(blocks.size())-qb[i]-qb[j]+ub;
    int n=(labels.size()<<21)+2048*neither+1024*both+16*un+8*ub+124;
    if(n>max63){max63=n;I=i;J=j;U=neither;V=both;}
  }
  int n64=(labels.size()<<22)+signs*supports.size()+16*blocks.size()+128;
  assert(n64==413526912);
  std::cout<<"weight-16 supports "<<supports.size()<<"\n";
  std::cout<<"best 63D pair "<<I<<" "<<J<<", neither "<<U<<", both "<<V<<", lower bound "<<max63<<"\n";
  for(int d=64;d<=68;d++)std::cout<<"dimension "<<d<<" lower bound "<<n64+2*(d-64)<<"\n";
  if(argc==3){ std::ofstream f(argv[2],std::ios::binary);f.write((const char*)supports.data(),sizeof(u64)*supports.size()); }

}
