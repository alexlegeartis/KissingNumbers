#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
int main(int argc,char**argv){
 if(argc!=5)return 2;int n=std::stoi(argv[1]),limit=std::stoi(argv[4]);
 std::ifstream in(argv[2]);std::vector<uint64_t> supports;uint64_t m;
 while(in>>m)supports.push_back(m);std::ofstream out(argv[3]);long long count=0;
 for(int a=0;a<n;a++){
  std::vector<int> containing;for(int i=0;i<(int)supports.size();i++)if(supports[i]>>a&1)containing.push_back(i);
  for(int b=0;b<n;b++)if(b!=a)for(int c=b+1;c<n;c++)if(c!=a)
  for(int d=c+1;d<n;d++)if(d!=a)for(int e=d+1;e<n;e++)if(e!=a){
   uint64_t leaves=(1ULL<<b)|(1ULL<<c)|(1ULL<<d)|(1ULL<<e);std::vector<int> bad;
   for(int i:containing)if(__builtin_popcountll(supports[i]&leaves)>=3){bad.push_back(i);if((int)bad.size()>limit)break;}
   if((int)bad.size()>limit)continue;out<<a<<' '<<leaves;
   for(int i:bad)out<<' '<<i;out<<'\n';count++;
  }
 }
 std::cerr<<"near-star supports "<<count<<'\n';
}
