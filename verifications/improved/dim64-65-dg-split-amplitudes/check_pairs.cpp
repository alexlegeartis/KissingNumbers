#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using U=uint64_t;
int main(int argc,char**argv){
  assert(argc==2);
  std::ifstream f(argv[1]);
  unsigned n,m;
  f>>n>>m;
  assert(n==40414&&m==14723);
  std::vector<U>B(n);
  for(U&u:B){f>>u;assert(std::popcount(u)==16);}
  assert(std::set(B.begin(),B.end()).size()==n);
  std::vector<std::vector<unsigned>> E(n);
  for(unsigned i=0,a,b;i<m;++i){f>>a>>b;assert(n>a&&a>b);E[a].push_back(b);}
  unsigned extra;
  assert(f && !(f>>extra) && f.eof());
  unsigned checked=0;
  for(unsigned i=0;i<n;++i){
    std::sort(E[i].begin(),E[i].end());
    unsigned k=0;
    for(unsigned j=0;j<i;++j){
      unsigned d=std::popcount(B[i]&B[j]);
      if(d>8){assert(d<=10&&k<E[i].size()&&E[i][k]==j);++k;}
    }
    assert(k==E[i].size());
    checked+=k;
  }
  assert(checked==m);
  std::cout<<"PASS every one of "<<U(n)*(n-1)/2<<" support pairs, "<<m<<" guards\n";
}
