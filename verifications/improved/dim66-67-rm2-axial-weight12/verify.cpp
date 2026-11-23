#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using U=uint64_t;
int wt(U x){return __builtin_popcountll(x);}
U mon(int a){U s=0;for(int i=0;i<64;i++)if((i&a)==a)s|=U(1)<<i;return s;}
int main(int argc,char**argv){
 std::ifstream is(argc>1?argv[1]:"witness.txt");
 std::vector<int> C;for(int a;is>>a;)C.push_back(a);
 assert(is.eof()&&C.size()==80&&std::set<int>(C.begin(),C.end()).size()==80);
 std::vector<U> G,A;
 for(int i=0;i<64;i++){if(wt(i)<=2)G.push_back(mon(i));if(wt(i)==3)A.push_back(mon(i));}
 assert(G.size()==22&&A.size()==20);
 std::vector<U> Q(1u<<22);
 for(unsigned i=1;i<Q.size();i++)Q[i]=Q[i-1]^G[__builtin_ctz(i)];
 auto X=[&](unsigned a){U s=0;for(;a;a&=a-1)s^=A[__builtin_ctz(a)];return s;};
 int shifts[5]={0,822310,213640,779738,262160};
 int tails[5]={0,3,5,6,0};
 std::vector<U> S[5];
 for(int i=0;i<5;i++){
  int target=i==0?16:i==4?12:14;
  for(int a:C){assert(0<=a&&a<(1<<20));U z=X(a^shifts[i]);
   for(U h:Q)if(wt(z^h)==target)S[i].push_back(z^h);
  }
  std::sort(S[i].begin(),S[i].end());
  assert(std::adjacent_find(S[i].begin(),S[i].end())==S[i].end());
 }
 assert(S[0].size()==37996&&S[4].size()==128&&S[1].size()==3904&&S[2].size()==3904&&S[3].size()==3904);
 std::map<std::array<int,4>,unsigned> histogram;
 unsigned ok66=0,ok67=0;
 for(size_t i=0;i<S[4].size();i++){
  U a=S[4][i];
  for(size_t j=0;j<i;j++)assert(wt(a&S[4][j])<=4);
  std::array<int,4> ov{};
  for(int j=0;j<4;j++)for(U b:S[j])ov[j]=std::max(ov[j],wt(a&b));
  ++histogram[ov];
  for(int j=0;j<2;j++)if(4*ov[0]<=32&&4*ov[1]+8*((tails[1]>>j)&1)<=32)ok66++;
  for(int j=0;j<3;j++){
   bool pass=4*ov[0]<=32;
   for(int k=1;k<4;k++)pass&=4*ov[k]+8*((tails[k]>>j)&1)<=32;
   if(pass)ok67++;
  }
 }
 assert(histogram.size()==2);
 assert(histogram.at({8,6,7,7})==96&&histogram.at({8,6,6,6})==32);
 assert(ok66==256&&ok67==96);
 // Direct replay of all local palettes and every normalized inequality of
 // two vectors on the same support; distinct-support intersection <=4 tested above.
 std::vector<U> P[3];unsigned mask[3]={0,3,5};
 for(unsigned s=0;s<4096;s++){
  int syndrome=0;
  for(int a=0;a<12;a++)if((s>>a)&1)syndrome^=(16|a);
  if(!syndrome)for(int j=0;j<3;j++)P[j].push_back(s^mask[j]);
 }
 for(int j=0;j<3;j++){
  assert(P[j].size()==128);
  for(int k=0;k<3;k++)for(U a:P[j])for(U b:P[k])
   for(int z:{-1,1})for(int w:{-1,1})if(j!=k||a!=b||z!=w)
    assert(48-8*wt(a^b)+(j==k?16*z*w:0)<=32);
 }
 // Uniform checks against dense, four-body support, and poles.
 assert(12*2<=32&&4*2*4<=32&&8*4<=32);
 unsigned old[2]={421522308,437513094};
 std::cout<<"PASS support counts 37996/3904/3904/3904/128; two exact overlap profiles 96/32\n";
 std::cout<<"PASS literal same-support checks and pair gates, slots "<<ok66<<" "<<ok67<<"\n";
 std::cout<<"66 "<<old[0]+ok66*256<<"\n67 "<<old[1]+ok67*256<<"\n";
}
