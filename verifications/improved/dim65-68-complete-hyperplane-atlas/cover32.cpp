#include<bits/stdc++.h>
#include<immintrin.h>
using namespace std;using U=uint32_t;using V=uint64_t;
U mono(int k){U w=0;for(int x=0;x<32;x++)if((x&k)==k)w|=1u<<x;return w;}
unsigned gf(unsigned a,unsigned b){unsigned c=0;for(;b;b>>=1){if(b&1)c^=a;a<<=1;if(a&32)a^=37;}return c;}
vector<U> genwords(const vector<U>&G){vector<U>S={0};for(U u:G){auto n=S.size();for(size_t i=0;i<n;i++)S.push_back(S[i]^u);}return S;}
int main(){
 vector<U>K={~U(0)}, Qb,Cb;
 for(int i=0;i<32;i++)(popcount(unsigned(i))<=2?Qb:Cb).push_back(mono(i));
 for(int i=0;i<5;i++)K.push_back(mono(1<<i));
 for(int i=0;i<5;i++){
   U a=0;
   for(unsigned x=0;x<32;x++){
    unsigned z=gf(gf(gf(x,x),x),1u<<i),t=0;for(int r=0;r<5;r++){t^=z;z=gf(z,z);}
    assert(t<2);a|=t<<x;
   }
   K.push_back(a);
 }
 U e[32]={},select=~U(0);
 for(U u:K){
   while(u){int x=31-countl_zero(u);if(!e[x]){e[x]=u;select^=1u<<x;break;}u^=e[x];}
   assert(u);
 }
 assert(popcount(select)==21);
 auto label=[&](U u){for(int i=31;i>=0;i--)if(e[i]&&(u>>i&1))u^=e[i];return _pext_u32(u,select);};
 vector<U> Qdom;for(U u:genwords(Qb))Qdom.push_back(label(u));
 sort(Qdom.begin(),Qdom.end());Qdom.erase(unique(Qdom.begin(),Qdom.end()),Qdom.end());assert(Qdom.size()==32);
 vector<unsigned> D(1<<21,255),W(1<<21);W[0]=0;D[0]=0;size_t len=1;
 U units[32];for(int j=0;j<32;j++)units[j]=label(1u<<j);
 for(size_t i=0;i<len;i++)for(U b:units){
  U w=W[i]^b;
  if(D[w]==255){W[len++]=w;D[w]=D[W[i]]+1;}
 }
 assert(len==1u<<21);
 array<int,33>Krad={};for(int r:D)Krad[r]++;
 cout<<"All quotient distances / radius";for(int d:Krad)if(d)cout<<" "<<d;cout<<endl;
 auto XC=genwords(Cb);array<int,33>RM2rad={};map<array<int,33>,int> shapes;
 int all5=0,all5_high=0,maxgood=0,minLow=32;
 for(unsigned i=0;i<XC.size();i++){
  auto p=label(XC[i]);array<int,33> H={};int m=32,M=0;
  for(unsigned y:Qdom){int v=D[p^y];m=min(m,v);M=max(M,v);H[v]++;}
  RM2rad[m]++;
  if(m>=5){++all5;if(M>=9){++all5_high;minLow=min(minLow,H[5]+H[6]);++shapes[H];}}
  if(M>=9)maxgood=max(maxgood,32-H[0]-H[1]-H[2]-H[3]-H[4]-H[5]-H[6]-H[7]);
 }
 cout<<"RM2 covering distribution";for(int d:RM2rad)if(d)cout<<" "<<d;cout<<endl;
 cout<<"coset d>=5 "<<all5<<" & max nine "<<all5_high<<" min count <=6 "<<minLow<<" max count >=8 "<<maxgood<<endl;
 for(auto [m,n]:shapes){cout<<"shape "<<n;for(int d:m)if(d)cout<<" "<<(&d?"":"")<<d;cout<<endl;}
}
