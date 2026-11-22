#include <array>
#include <cassert>
#include <cstdint>
#include <set>
#include <iostream>
using U=uint32_t; using LL=uint64_t;
constexpr int n=24;
constexpr U mask=(1u<<24)-1;
int G[n][2*n];
struct T{U a,b;};
T operator+(T x,T y){
 U z=~(x.a|x.b), w=~(y.a|y.b);
 return {(x.a&w)|(z&y.a)|(x.b&y.b),(x.b&w)|(z&y.b)|(x.a&y.a)};
}
std::array<T,2> A[n];
U fix[2],ax[2];
LL nodes[8],h[8],H[5]; int k,side;
int wt(U x){return __builtin_popcount(x);}
void dfs(int i,int m,T y,U support){
 ++nodes[m];
 if(wt(y.a|y.b)==15-m && ((y.a|y.b)&fix[1-side])==0){
  ++h[m];
  int a=wt(support&ax[side])+wt((y.a|y.b)&ax[1-side]);
  assert(a>=0&&a<=k);
  ++H[a];
 }
 if(m>=7)return;
 for(int j=i;j<24;j++)if(!(fix[side]>>j&1)){
  dfs(j+1,m+1,y+A[j][0],support|(1u<<j));
  dfs(j+1,m+1,y+A[j][1],support|(1u<<j));
 }
}
int main(){
 auto& f=std::cin;
 for(int i=0;i<24;i++)for(int j=0;j<48;j++){
  f>>G[i][j];assert(G[i][j]>=0&&G[i][j]<3);
  if(j<24)assert(G[i][j]==(i==j));
 }
 for(int i=0;i<24;i++)for(int j=0;j<24;j++){
  int x=0;
  for(int t=0;t<48;t++)x+=G[i][t]*G[j][t];
  assert(x%3==0);
 }
 fix[0]=(1u<<0)|(1u<<3);
 std::set<int> used{0,3};
 for(int i=0;i<4;i++){int t;f>>t;assert(0<=t&&t<48&&used.insert(t).second);ax[t/24]|=1u<<(t%24);}
 k=4;
 for(side=0;side<2;side++){
  for(int i=0;i<n;i++){
   U a=0,b=0;
   for(int j=0;j<n;j++){
    int v=side? (3-G[j][i+n])%3:G[i][j+n];
    if(v==1)a|=1u<<j;
    if(v==2)b|=1u<<j;
   }
   A[i][0]={a,b}; A[i][1]={b,a};
  }
  dfs(0,0,{},0);
  std::cout<<"side "<<side<<" nodes";for(LL x:nodes)std::cout<<" "<<x;
  std::cout<<" outer";for(LL x:h)std::cout<<" "<<x;std::cout<<"\n";
 }
 LL t=0, g=0;
 LL w[5]={7812,6561,2187,729,729};
 std::cout<<"histogram";
 for(int i=0;i<=k;i++){std::cout<<" "<<H[i];t+=H[i];g+=H[i]*w[i];}
 std::cout<<" count "<<t<<" sparse "<<g<<"\n";
}
