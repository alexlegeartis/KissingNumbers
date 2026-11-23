#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using U=uint64_t;
using W=uint32_t;
U mono(int m,int n){U v=0;for(int x=0;x<n;x++)if((x&m)==m)v|=U(1)<<x;return v;}
U eval(unsigned x,const std::vector<U>& b){
  assert(x<(1u<<b.size()));
  U v=0;for(;x;x&=x-1)v^=b[std::countr_zero(x)];return v;
}
std::vector<U> span(const std::vector<U>& b){
  std::vector<U> a={0};
  for(auto x:b){size_t n=a.size();for(size_t j=0;j<n;j++)a.push_back(a[j]^x);}
  return a;
}
int main(int argc,char** argv){
  assert(argc>=2);
  std::vector<U> C,Q,C5,Q5;
  for(int m=0;m<64;m++){
   if(std::popcount(unsigned(m))==3)C.push_back(mono(m,64));
   if(std::popcount(unsigned(m))<=2)Q.push_back(mono(m,64));
   if(m<32)(std::popcount(unsigned(m))<=2?Q5:C5).push_back(mono(m,32));
  }
  assert(C.size()==20&&Q.size()==22&&C5.size()==16&&Q5.size()==16);
  auto A=span(Q),B=span(Q5);assert(A.size()==(1<<22)&&B.size()==(1<<16));
  std::ifstream parents("cubic80.txt");
  std::vector<unsigned> P;
  for(unsigned c;parents>>c;)P.push_back(c);
  assert(parents.eof()&&P.size()==80);
  for(int ia=1;ia<argc;ia++)for(int j=1;j<ia;j++){
    unsigned u=std::stoul(argv[ia]),v=std::stoul(argv[j]);assert(u!=v);
    U w=eval(u^v,C);int best=64;
    for(U q:A)best=std::min(best,std::popcount(w^q));
    std::cout<<"PAIR "<<u<<' '<<v<<' '<<best<<'\n';
  }
  for(int ia=1;ia<argc;ia++){
   unsigned s=std::stoul(argv[ia]);U donor=eval(s,C);
   int glob=64;
   for(auto p:P){
    assert(p!=s);int best=64;
    U v=eval(p^s,C);
    for(U q:A)best=std::min(best,std::popcount(v^q));
    assert(best>=14);glob=std::min(glob,best);
   }
   std::cout<<"PARENT "<<s<<' '<<glob<<' '<<P.size()*A.size()<<'\n';
   std::ifstream hh("hyperplane_heads.txt");
   int n=0;
   for(unsigned m,k,g,a;hh>>m>>k>>g>>a;n++){
    assert(m>0&&m<64&&k<=1&&a==15);
    W v=0;int c=0;
    for(unsigned x=0;x<64;x++)if((std::popcount(x&m)&1u)==k){
     v|=W((donor>>x)&1u)<<c++;
    }
    assert(c==32);
    v^=W(eval(g,C5));
    int best=32;
    for(U q:B)best=std::min(best,std::popcount(W(q)^v));
    std::cout<<"HEAD "<<s<<' '<<n<<' '<<best<<'\n';
   }
   assert(hh.eof()&&n==126);
  }
}
