#include<bits/stdc++.h>
using namespace std;using U=uint64_t;
int main(){
 const int n=1<<22;
 vector<U> Q(n);vector<int> f(n),g(n);
 vector<U> b;
 for(int m=0;m<64;m++)if(popcount((unsigned)m)<=2){
  U w=0;for(int j=0;j<64;j++)if((j&m)==m)w|=1ull<<j;
  b.push_back(w);
 }
 for(int i=1;i<n;i++){Q[i]=Q[i&(i-1)]^b[countr_zero((unsigned)i)];int w=popcount(Q[i]);f[i]=(w==16);if(f[i])g[i]=w;}
 const int degree=count(f.begin(),f.end(),1);
 for(int d=1;d<n;d*=2)for(int i=0;i<n;i+=d*2)for(int j=0;j<d;j++){int x=f[i+j],y=f[i+j+d];f[i+j]=x+y;f[i+j+d]=x-y;}
 map<int,int> h;for(int z:f)h[z]++;
 for(auto [k,v]:h)cout<<k<<" "<<v<<"\n";
 int e=-h.begin()->first;cout<<"degree "<<degree<<" cap "<<(uint64_t(n)*e)<<"/"<<degree+e<<" = "<<double(n)*e/(degree+e)<<"\n";
}
