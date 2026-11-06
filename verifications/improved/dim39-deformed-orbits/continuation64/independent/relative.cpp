#include <climits>
#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <limits>
using namespace std;
int main(){int n,d,m; long long L; cin>>n>>d>>m>>L; vector<uint64_t> gs(m);for(auto &g:gs)cin>>g;vector<array<long long,39>> a(n);for(auto &v:a)for(auto &x:v)cin>>x;long long mx=LLONG_MIN;uint64_t count=0;int wi=-1,wj=-1,wg=-1;for(int g=0;g<m;g++){for(int i=0;i<n;i++){array<long long,39> v;for(int k=0;k<39;k++)v[k]=((gs[g]>>k)&1)?-a[i][k]:a[i][k];for(int j=i;j<n;j++){if(i==j&&gs[g]==0)continue;long long s=0;for(int k=0;k<39;k++)s+=v[k]*a[j][k];++count;if(s>mx){mx=s;wi=i;wj=j;wg=g;}if(2*s>L*L){cerr<<"VIOLATION "<<i<<" "<<j<<" "<<g<<" "<<s<<endl;return 1;}}}cerr<<"group "<<g+1<<"/"<<m<<endl;}cout<<count<<" "<<mx<<" "<<wi<<" "<<wj<<" "<<wg<<endl;}
