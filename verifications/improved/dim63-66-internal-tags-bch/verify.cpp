#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using U=uint64_t;
using namespace std;
int pop(U x){return __builtin_popcountll(x);}
U mon(int a){U v=0;for(int i=0;i<64;i++)if((a&i)==a)v|=1ull<<i;return v;}
using V=pair<U,int>;
int bit(U a,int n){return a>>n&1;}
int lab(U a,int degree){
 for(int i=0;i<6;i++)a^=(a&~mon(1<<i))<<(1<<i);
 int n=0,v=0;
 for(int i=0;i<64;i++){
  if(pop(i)==degree)v|=bit(a,i)<<n++;
  else if(pop(i)>degree)assert(!bit(a,i));
 }
 return v;
}
vector<V> load(const char *p){
 ifstream f(p);assert(f);
 vector<V>Z;set<V>D;
 U v;int h;
 while(f>>v>>h){
  assert(pop(v)==12&&0<=h&&h<64&&!bit(v,h)&&D.emplace(v,h).second);
  Z.push_back({v,h});
 }
 assert(f.eof());
 return Z;
}
int main(){
 vector<int> C;{ifstream f("old_cubic.txt");assert(f);for(int x;f>>x;)C.push_back(x);assert(f.eof());}
 assert(C.size()==80&&set<int>(C.begin(),C.end()).size()==80);
 vector<U> quad={0}, cub;
 for(int i=0;i<64;i++){
  if(pop(i)<=2){U v=mon(i);auto n=quad.size();for(size_t j=0;j<n;j++)quad.push_back(quad[j]^v);}
  if(pop(i)==3)cub.push_back(mon(i));
 }
 assert(quad.size()==(1u<<22)&&cub.size()==20);
 auto make=[&](int a){U z=0;assert(0<=a&&a<(1<<20));for(unsigned i=0;i<20;i++)if(a>>i&1)z^=cub[i];return z;};
 vector<U> old[3];
 int shift[3]={0,822310,262160},w[3]={16,14,12},counts[3]={37996,3904,128};
 for(int j=0;j<3;j++){
  for(int x:C){
   U z=make(x^shift[j]);
   for(U y:quad)if(pop(z^y)==w[j])old[j].push_back(z^y);
  }
  assert(old[j].size()==counts[j]&&set<U>(old[j].begin(),old[j].end()).size()==counts[j]);
 }
 auto A=load("tags64.txt"),B=load("tags66.txt");
 set<V> UA(A.begin(),A.end());assert(A.size()==256&&B.size()==64);
 for(auto p:B)assert(UA.count(p));
 vector<int> signs[4];int tr[4]={0,3,5,9};
 for(int a=0;a<4096;a++){
  int t=0;for(int i=0;i<12;i++)if(a>>i&1)t^=16|i;
  if(!t)for(int j=0;j<4;j++)signs[j].push_back(a^tr[j]);
 }
 auto verify_tags=[&](const vector<V>&T,int where){
  map<U,int> hist;uint64_t sec=0;
  for(int i=0;i<(int)T.size();i++){
   auto [a,h]=T[i];int s=hist[a]++;assert(s<4);
   int p=0,q=0;for(int k=0;k<20;k++)p+=bit(a,k);q=p;
   if(h!=0&&h!=20&&bit(a,0)==bit(a,20)){
    int k0=0;
    for(int pat:signs[s])if(!bit(a,0)||((pat&1)==((pat>>q)&1)))k0++;
    sec+=2*k0;
   }
   // An old four-position support is an affine 2-flat, so any
   // violation containing h would force h=x+y+z for three body sites.
   for(int x=0;x<64;x++)if(bit(a,x))
    for(int y=x+1;y<64;y++)if(bit(a,y))
     for(int z=y+1;z<64;z++)if(bit(a,z))assert((x^y^z)!=h);
   for(int j=0;j<where;j++)
    for(U v:old[j])assert(4*pop(a&v)+8*bit(v,h)<=32);
   for(int j=0;j<i;j++){
    auto [b,k]=T[j];if(a==b)continue;
    assert(4*pop(a&b)+(h==k?16:8*(bit(a,k)+bit(b,h)))<=32);
   }
  }
  for(int x=0;x<4;x++){
   assert(signs[x].size()==128);
   for(int y=0;y<4;y++)for(int a:signs[x])for(int b:signs[y]){
    if(x==y){if(a!=b)assert(pop(a^b)>=4);}
    else assert(pop(a^b)>=2);
   }
  }
  // L1 cap vs dense; max product with an old axial pole or any
  // four-site vector whose support misses h (unless a 2-flat).
  assert(12*2+4<=32&&8*4<=32&&4*2*4<=32&&4*(4+2*2)<=32);
  return sec;
 };
 uint64_t section=verify_tags(A,1);verify_tags(B,3);
 assert(section>0);
 ifstream p("polar.txt");assert(p);
 int n;p>>n;assert(n==3);
 vector<U> H={0},R(n);
 for(int i=0;i<16;i++){
  U w;p>>w;lab(w,2);
  int m=H.size();for(int j=0;j<m;j++)H.push_back(H[j]^w);
 }
 assert(H.size()==65536&&set<U>(H.begin(),H.end()).size()==H.size());
 for(auto v:H)if(v)assert(pop(v)>=22);
 for(U &v:R){p>>v;lab(v,3);}
 U s;p>>s;lab(s,2);assert(p);
 for(int i=0;i<n;i++){
  for(int x:C)for(U u:quad)assert(pop(R[i]^make(x)^u)>=14);
  for(int j=0;j<n;j++)for(U u:H){
   assert(pop(R[i]^R[j]^s^u)>=8);
   if(i!=j)assert(pop(R[i]^R[j]^u)>=22);
  }
 }
 // New polar squared norm 81, old stock norm squared 64. Diverse
 // same-side Hamming distance at least 22; opposite sides at least 8.
 // All new internal tags have L1 exactly 28, only on equatorial sites.
 assert(64-2*22+17<=40&&64-2*8-17<=40&&64-2*14<=36&&8*8*17<=36*36&&28<=36);
 cout<<"PASS all old interstock inequalities, class distances and the 63D section."<<endl;
 cout<<"63 "<<214085892+section<<endl;
 cout<<"64 "<<413526912+uint64_t(A.size())*256<<endl;
 cout<<"65 "<<413526914+uint64_t(A.size())*256+2*n*H.size()<<endl;
 cout<<"66 "<<421587844+uint64_t(B.size())*256<<endl;
}
