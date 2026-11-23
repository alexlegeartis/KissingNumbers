#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using U=uint64_t;
static int B(U x){return __builtin_popcountll(x);}
static int F(U x){return __builtin_ctzll(x);}
struct Node{
 uint32_t group;
 uint16_t ray,b0,b1,b2;
 uint32_t row;
};
vector<unsigned char> chosen;
vector<vector<int>> A;
vector<Node> W;
vector<vector<unsigned char>> Conflict;
// Exact maximum clique on the complement for small connected components;
// a deterministic node cap on larger ones leaves only a lower bound.
struct Tiny {
 U C[64],best=0;
 int b=0,nodes=0,cap=0;
 void rec(U P,U S,int q){
   if(++nodes>cap)return;
   if(!P){ if(q>b){b=q;best=S;} return; }
   int ord[64],col[64],t=0,h=0;
   U left=P;
   while(left){
     ++h;U avail=left;
     while(avail){int i=F(avail);U x=U(1)<<i;
       ord[t]=i;col[t++]=h;left^=x;avail&=~(x|C[i]);
     }
   }
   for(int j=t-1;j>=0;--j){
     if(q+col[j]<=b)return;
     int i=ord[j];U x=U(1)<<i;
     rec(P&C[i],S|x,q+1);
     if(nodes>cap)return;
     P^=x;
   }
 }
};
static unsigned token(unsigned z){
 z^=z>>16;z*=0x7feb352du;z^=z>>15;z*=0x846ca68bu;return z^(z>>16);
}
void improve(const vector<int>& pos,int off,vector<int>& mp){
 const int n=(int)pos.size();
 if(n<=1){if(n==1)chosen[W[off+pos[0]].row]=1;return;}
 vector<int> best;
 for(int j:pos)if(chosen[W[off+j].row])best.push_back(j);
 vector<int> cand,ord=pos;
 for(unsigned it=0;it<5;++it){
  sort(ord.begin(),ord.end(),[&](int x,int y){
   if(it<3&&A[x].size()!=A[y].size())return A[x].size()<A[y].size();
   return token(unsigned(x)^it*0x9e3779b9u)<token(unsigned(y)^it*0x9e3779b9u);
  });
  for(int j:pos)mp[j]=0;
  cand.clear();
  for(int j:ord)if(!mp[j]){cand.push_back(j);for(int k:A[j])mp[k]=1;}
  if(cand.size()>best.size())best=cand;
 }
 if(n<=64 && n>=3){
   const U full=n==64?~U(0):(U(1)<<n)-1;
   Tiny s;s.cap=n<=16?5000:1400;
   for(int j=0;j<n;++j)mp[pos[j]]=j;
   for(int j=0;j<n;++j){
     U a=U(1)<<j;for(int k:A[pos[j]])a|=U(1)<<mp[k];
     s.C[j]=full^a;
   }
   for(int j:best)s.best|=U(1)<<mp[j];
   s.b=B(s.best);
   s.rec(full,0,0);
   if(s.b>(int)best.size()){
      best.clear();while(s.best){int j=F(s.best);best.push_back(pos[j]);s.best&=s.best-1;}
   }
 }
 for(int j:pos)chosen[W[off+j].row]=0;
 for(int j:best)chosen[W[off+j].row]=1;
}
int main(){
 ios::sync_with_stdio(false);
 uint32_t r,n;cin.read((char*)&r,4);cin.read((char*)&n,4);
 assert(r<=453&&n==13824*r);
 Conflict.assign(r,vector<unsigned char>(r));
 for(auto &p:Conflict)cin.read((char*)p.data(),r);
 W.resize(n);chosen.resize(n);
 cin.read((char*)W.data(),sizeof(Node)*n);static_assert(sizeof(Node)==16);
 cin.read((char*)chosen.data(),n);assert(cin.good());
 sort(W.begin(),W.end(),[](const Node&a,const Node&b){return a.group!=b.group?a.group<b.group:a.row<b.row;});
 uint64_t old=accumulate(chosen.begin(),chosen.end(),uint64_t(0)),edges=0;
 uint32_t groups=0,maxg=0,maxc=0;int comps=0;
 for(uint32_t i=0;i<n;){
   uint32_t j=i+1;while(j<n&&W[j].group==W[i].group)++j;
   int k=j-i;groups++;maxg=max(maxg,(uint32_t)k);
   A.assign(k,{});
   for(int q=0;q<k;q++)for(int p=0;p<q;p++){
    const Node&x=W[i+q];const Node&y=W[i+p];
    if(!Conflict[x.ray][y.ray])continue;
    int same=(x.b0==y.b0)+(x.b1==y.b1)+(x.b2==y.b2);
    const unsigned mode=Conflict[x.ray][y.ray];
    if(mode==1)assert(same<3);
    if(same!=mode+1)continue;
    assert(!(chosen[x.row]&&chosen[y.row]));
    A[p].push_back(q);A[q].push_back(p);++edges;
   }
   vector<unsigned char>seen(k,0);vector<int>mp(k),pos;
   for(int p=0;p<k;p++)if(!seen[p]){
    seen[p]=1;pos.assign(1,p);
    for(size_t x=0;x<pos.size();x++)for(int y:A[pos[x]])if(!seen[y]){seen[y]=1;pos.push_back(y);}
    maxc=max(maxc,(uint32_t)pos.size());comps++;
    improve(pos,i,mp);
   }
   for(int p=0;p<k;++p) if(chosen[W[i+p].row])for(int q:A[p])assert(!chosen[W[i+q].row]);
   i=j;
 }
 uint64_t nw=accumulate(chosen.begin(),chosen.end(),uint64_t(0));assert(nw>=old);
 cerr<<"components "<<comps<<" groups "<<groups<<" maxgroup "<<maxg<<" maxcomponent "<<maxc<<" edges "<<edges<<" seed "<<old<<" result "<<nw<<" gained "<<nw-old<<"\n";
 cout.write((char*)chosen.data(),n);
}
