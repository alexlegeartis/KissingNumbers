#include<bits/stdc++.h>
using namespace std;using U=uint64_t;
int main(int argc,char**argv){
 assert(argc==2); ifstream f(argv[1]);vector<int>L;
 for(int x;f>>x;)L.push_back(x);
 assert(L.size()==80);vector<U>M,Q;vector<U>A{0},B{0};
 auto mon=[](int k,int p){U s=0;for(int x=0;x<1<<p;x++)if((x&k)==k)s|=U(1)<<x;return s;};
 for(int k=0;k<64;k++)if(__builtin_popcount(k)==3)M.push_back(mon(k,6));
 for(int k=0;k<32;k++){
  int w=__builtin_popcount(k); if(w>2)continue; U s=mon(k,5);
  size_t n=A.size();for(size_t i=0;i<n;i++)A.push_back(A[i]^s);
  if(w<=1){size_t n=B.size();for(size_t i=0;i<n;i++)B.push_back(B[i]^s);}
 }
 for(U d:B)for(U a:A)Q.push_back(a|((a^d)<<32));
 assert(Q.size()==(1u<<22));
 auto tab=[&](int d){U w=0;for(;d;d&=d-1)w^=M[__builtin_ctz(d)];return w;};
 const int N=4;int shifts[N]={0,822310,213640,779738};
 int tail[N]={0,3,5,6};vector<U> supp[N];
 set<int> DD;
 for(int a:L)for(int b:L)DD.insert(a^b);
 set<int> D;
 for(int a=0;a<N;a++)for(int b=0;b<a;b++)for(int z:DD)D.insert(z^shifts[a]^shifts[b]);
 assert(!D.count(0));
 for(int z:D){
  U t=tab(z); int md=64;
  for(U q:Q){int w=__builtin_popcountll(q^t);if(w<md)md=w;}
  assert(md>=14);
 }
 for(int a=0;a<N;a++){
  int wt=16-__builtin_popcount(tail[a]);
  for(int d:L){U t=tab(d^shifts[a]);for(U q:Q)if(__builtin_popcountll(q^t)==wt)supp[a].push_back(q^t);}
  sort(supp[a].begin(),supp[a].end());
  assert(adjacent_find(supp[a].begin(),supp[a].end())==supp[a].end());
  assert(supp[a].size()==(a?3904:37996));
 }
 // Independent literal check, which verifies the actual augmented supports
 // irrespective of the uniform minimum used to design the shifts.
 for(int a=0;a<N;a++)for(int b=0;b<=a;b++){
  int overlap=__builtin_popcount(tail[a]&tail[b]);
  for(size_t i=0;i<supp[a].size();i++){
   size_t end=a==b?i:supp[b].size();
   for(size_t j=0;j<end;j++)assert(__builtin_popcountll(supp[a][i]&supp[b][j])+overlap<=8);
  }
 }
 cout<<"PASS "<<D.size()<<" translated pair difference types, all distance >=14; exact pairwise 16-support overlaps <=8\n";
 cout<<"tail support counts";for(auto &v:supp)cout<<" "<<v.size();cout<<endl;
}
