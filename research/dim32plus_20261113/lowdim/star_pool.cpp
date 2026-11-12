// Complete search pool for replacing every octad containing one coordinate.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=atoi(argv[1]);string prefix=argv[3];vector<uint64_t>seed;uint64_t x;ifstream f(argv[2]);while(f>>x)seed.push_back(x);
 uint32_t C[64][6]={};for(int i=0;i<64;i++){C[i][0]=1;for(int j=1;j<6;j++)C[i][j]=i?C[i-1][j]+C[i-1][j-1]:0;}
 vector<array<int,5>>sub;for(int a=0;a<4;a++)for(int b=a+1;b<5;b++)for(int c=b+1;c<6;c++)for(int d=c+1;d<7;d++)for(int e=d+1;e<8;e++)sub.push_back({a,b,c,d,e});
 vector<int>owner(C[n][5],-1);for(int id=0;id<(int)seed.size();id++){int p[8],k=0;for(x=seed[id];x;x&=x-1)p[k++]=__builtin_ctzll(x);for(auto s:sub){int t=0;for(int j=0;j<5;j++)t+=C[p[s[j]]][j+1];if(owner[t]>=0)abort();owner[t]=id;}}
 vector<vector<uint64_t>>out(n);uint64_t all=(1ULL<<n)-1,tested=0;
 for(int a=0;a<n-7;a++)for(int b=a+1;b<n-6;b++)for(int c=b+1;c<n-5;c++)for(int d=c+1;d<n-4;d++)for(int e=d+1;e<n-3;e++)for(int f=e+1;f<n-2;f++)for(int g=f+1;g<n-1;g++)for(int h=g+1;h<n;h++){
  int p[8]={a,b,c,d,e,f,g,h};uint64_t common=all;
  for(auto q:sub){int t=0;for(int j=0;j<5;j++)t+=C[p[q[j]]][j+1];int o=owner[t];if(o>=0){common&=seed[o];if(!common)break;}}
  tested++;if(!common)continue;x=0;for(int v:p)x|=1ULL<<v;
  for(uint64_t q=common;q;q&=q-1)out[__builtin_ctzll(q)].push_back(x);
 }
 for(int j=0;j<n;j++){ofstream o(prefix+"_"+to_string(j)+".txt");for(auto x:out[j])o<<x<<'\n';cout<<j<<" "<<out[j].size()<<endl;}
 cout<<"tested "<<tested<<endl;
}
