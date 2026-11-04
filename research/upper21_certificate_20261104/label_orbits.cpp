#include <fstream>
#include <vector>
#include <cstdint>
#include <iostream>
using namespace std;
int main(int argc,char**argv){
 if(argc!=3)return 2;
 ifstream in(argv[1],ios::binary);uint32_t n,g;in.read((char*)&n,4);in.read((char*)&g,4);
 vector<uint32_t>gen(n*g),p(21);in.read((char*)gen.data(),4*gen.size());in.read((char*)p.data(),84);
 vector<uint32_t>lg(g*21);for(uint32_t j=0;j<g;j++)for(int a=0;a<21;a++)for(int b=0;b<21;b++)if(gen[j*n+p[a]]==p[b])lg[j*21+a]=b;
 vector<uint32_t>orb(n*441,UINT32_MAX),reps,queue;
 for(uint32_t start=0;start<orb.size();start++){
  if(orb[start]!=UINT32_MAX)continue;
  uint32_t id=reps.size();reps.push_back(start);queue.clear();queue.push_back(start);orb[start]=id;
  for(size_t i=0;i<queue.size();i++){
   uint32_t d=queue[i]/441,a=queue[i]%441/21,b=queue[i]%21;
   auto add=[&](uint32_t x){if(orb[x]==UINT32_MAX){orb[x]=id;queue.push_back(x);}};
   for(uint32_t j=0;j<g;j++)add((gen[j*n+d]*21+lg[j*21+a])*21+lg[j*21+b]);
   add((d*21+b)*21+a);
  }
 }
 ofstream out(argv[2],ios::binary);uint32_t k=reps.size();out.write((char*)&k,4);out.write((char*)orb.data(),4*orb.size());out.write((char*)reps.data(),4*reps.size());
 cout<<"label pair orbits "<<k<<endl;
}
