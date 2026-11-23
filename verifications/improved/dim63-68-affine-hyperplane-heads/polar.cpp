#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
int wt(U x){return popcount(x);}
U mon(int a){U s=0;for(int i=0;i<64;i++)if((i&a)==a)s|=U(1)<<i;return s;}
vector<U> span(const vector<U>&bs) {
 vector<U>A={0};for(U b:bs){size_t n=A.size();for(size_t i=0;i<n;i++) A.push_back(A[i]^b);}
 assert(set<U>(A.begin(),A.end()).size()==A.size());return A;
}
int main(){
 vector<U>qb,cb;
 for(int i=0;i<64;i++){if(wt(i)<=2)qb.push_back(mon(i));if(wt(i)==3)cb.push_back(mon(i));}
 auto Q=span(qb);
 ifstream cf("cubic80.txt");vector<unsigned>C;for(unsigned c;cf>>c;)C.push_back(c);
 assert(cf.eof()&&C.size()==80&&set<unsigned>(C.begin(),C.end()).size()==80);
 auto word=[&](unsigned j){U a=0;for(;j;j&=j-1){auto k=countr_zero(j);assert(k<cb.size());a^=cb[k];}return a;};
 array<vector<U>,5>B;
 array<unsigned,5>S={0,822310,213640,779738,262160};
 array<int,5> target={16,14,14,14,12};
 for(int i=0;i<5;i++)for(auto c:C){U e=word(c^S[i]);for(U q:Q)if(wt(e^q)==target[i])B[i].push_back(e^q);}
 for(int i=0;i<5;i++){sort(B[i].begin(),B[i].end());assert(adjacent_find(B[i].begin(),B[i].end())==B[i].end());}
 assert(B[0].size()==37996&&B[1].size()==3904&&B[2].size()==3904&&B[3].size()==3904&&B[4].size()==128);
 ifstream wf("signed_witness.txt");vector<U> E;
 for(int i,j,a,b;wf>>i>>j>>a>>b;)if(j==-1){U e=uint32_t(a)|(U(uint32_t(b))<<32);assert(wt(e)==16);E.push_back(e);}
 assert(wf.eof()&&E.size()==2354);
 unsigned checks=0,n66=0,n67=0;
 array<int,4> tails={0,3,5,6};map<array<int,4>,unsigned> census;
 for(U a:B[4]){
   array<int,4>ov={};for(int i=0;i<4;i++)for(U b:B[i])ov[i]=max(ov[i],wt(a&b));
   ++census[ov];
   for(int k=0;k<2;k++)n66+=ov[0]<=8 && ov[1]+2*((tails[1]>>k)&1)<=8;
   for(int k=0;k<3;k++){bool ok=ov[0]<=8;for(int i=1;i<4;i++)ok&=ov[i]+2*((tails[i]>>k)&1)<=8;n67+=ok;}
   for(U b:E){assert(wt(a&b)<=8);++checks;}
   for(U b:B[4])if(a!=b)assert(wt(a&b)<=4);
 }
 assert(checks==301312&&n66==256&&n67==96&&census.size()==2);
 assert(census.at({8,6,7,7})==96&&census.at({8,6,6,6})==32);
 array<vector<U>,3> P;array<int,3> shift={0,3,5};
 for(unsigned i=0;i<4096;i++){
   unsigned s=0;for(int k=0;k<12;k++)if(i>>k&1)s^=16|k;
   if(!s)for(int j=0;j<3;j++)P[j].push_back(i^shift[j]);
 }
 for(int j=0;j<3;j++){
   assert(P[j].size()==128);
   for(int k=0;k<3;k++)for(U a:P[j])for(U b:P[k])
      for(int s:{-1,1})for(int t:{-1,1})if(j!=k||a!=b||s!=t)
       assert(48-8*wt(a^b)+(j==k?16*s*t:0)<=32);
 }
 assert(12*2<=32&&4*2*4<=32&&8*4<=32);
 cout<<"AXIAL exact new intersections "<<checks<<", max<=8; additions 66 "<<256*n66<<" 67 "<<256*n67<<endl;

 ifstream hf("tall.txt");
 vector<U> h(16);for(U &x:h)hf>>x;auto H=span(h);
 vector<U> A=span(vector<U>(h.begin(),h.begin()+7));
 assert(H.size()==65536&&A.size()==128);
 for(U a:H)if(a)assert(wt(a)>=22);
 for(U a:A)if(a)assert(wt(a)>=26);
 array<U,3> low;for(U&x:low)hf>>x;U shift2,high_master;hf>>shift2>>high_master;
 vector<U> high(8);for(U&x:high)hf>>x;
 U extra; assert(hf&&!(hf>>extra));
 assert(set<U>(high.begin(),high.end()).size()==8&&high[0]==high_master);
 for(U d:high)assert(find(H.begin(),H.end(),d^high_master)!=H.end());
 for(int i=0;i<3;i++){
   for(unsigned c:C){U z=low[i]^word(c);for(U q:Q)assert(wt(z^q)>=14);}
   for(int j=0;j<3;j++)for(U a:H){
     if(i!=j)assert(wt(low[i]^low[j]^a)>=22);
     assert(wt(low[i]^low[j]^a^shift2)>=8);
   }
   for(U a:H){assert(wt(low[i]^high_master^a)>=22);assert(wt(low[i]^high_master^a^shift2)>=13);}
 }
 for(unsigned c:C){U z=high_master^word(c);for(U q:Q)assert(wt(z^q)>=12);}
 for(int i=0;i<8;i++)for(int j=0;j<8;j++)if(i!=j)for(U a:A)
   assert(wt(high[i]^high[j]^a)>=26);
 // sqrt(17) / 6 layers, old shell and mixed radii. Low/tall worst
 // opposite-side distance 13 or better checked above.
 assert(64-2*22+17<=40&&64-2*8-17<=40&&64-2*14<=36);
 assert(64-2*26+36<=50&&64-36<=50&&64-2*12<=40);
 assert(36*36<=45*45&&36*17<25*25);
 assert(64*17<=36*36&&32<=36&&32<=40);
 cout<<"POLAR 65 393216 low replacement points and 2048 tall points; remove 2 old poles"<<endl;
 cout<<"65 "<<418741122+2046<<"\n66 "<<426343300+256*n66<<"\n67 "<<442334086+256*n67<<endl;
}
