#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=4){cerr<<"usage: enumerate_stars dimension supports output\n";return 2;} int n=stoi(argv[1]);ifstream in(argv[2]);vector<uint64_t> full;uint64_t m;while(in>>m)full.push_back(m);ofstream out(argv[3]);long long count=0;
 for(int a=0;a<n;a++)for(int b=0;b<n;b++)if(b!=a)for(int c=b+1;c<n;c++)if(c!=a)for(int d=c+1;d<n;d++)if(d!=a)for(int e=d+1;e<n;e++)if(e!=a){
  uint64_t leaves=(1ULL<<b)|(1ULL<<c)|(1ULL<<d)|(1ULL<<e);bool ok=true;
  for(uint64_t s:full)if(s>>a&1 && __builtin_popcountll(s&leaves)>=3){ok=false;break;}
  if(ok){out<<a<<" "<<leaves<<"\n";count++;}
 }
 cerr<<"n="<<n<<" compatible stars="<<count<<"\n";
}
