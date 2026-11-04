#include <fstream>
#include <vector>
#include <cstdint>
#include <iostream>
using namespace std;
int main(int argc,char**argv){
 if(argc!=4)return 2;
 ifstream in(argv[1],ios::binary);uint32_t n,g;in.read((char*)&n,4);in.read((char*)&g,4);
 vector<uint32_t> gen(n*g);in.read((char*)gen.data(),4*gen.size());
 for(int sym=0;sym<2;sym++){
  vector<uint32_t> orb(n*n,UINT32_MAX),reps,queue;queue.reserve(n*n);
  for(uint32_t start=0;start<n*n;start++){
   if(orb[start]!=UINT32_MAX)continue;
   uint32_t id=reps.size();reps.push_back(start);queue.clear();queue.push_back(start);orb[start]=id;
   for(size_t i=0;i<queue.size();i++){
    uint32_t a=queue[i]/n,b=queue[i]%n;
    auto add=[&](uint32_t x){if(orb[x]==UINT32_MAX){orb[x]=id;queue.push_back(x);}};
    for(uint32_t j=0;j<g;j++)add(gen[j*n+a]*n+gen[j*n+b]);
    if(sym){add(b*n+a);add(a*n+(a^b));}
   }
  }
  ofstream out(argv[2+sym],ios::binary);uint32_t k=reps.size();
  out.write((char*)&k,4);out.write((char*)orb.data(),4*orb.size());out.write((char*)reps.data(),4*reps.size());
  cout<<(sym?"unordered triple":"ordered pair")<<" orbits "<<k<<endl;
 }
}
