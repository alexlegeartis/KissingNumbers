#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
int main(int argc,char**argv){if(argc<3){cerr<<"usage: bridge_search seed.txt output.txt [seconds]\n";return 2;}ifstream f(argv[1]);vector<uint64_t>s;uint64_t x;while(f>>x)s.push_back(x);mt19937_64 rng(2026102804);vector<vector<uint64_t>> bucket(s.size());auto t=chrono::steady_clock::now();uint64_t trials=0,swaps=0;int best=s.size();
auto save=[&](){ofstream o(argv[2]);for(auto m:s)o<<m<<"\n";};
while(chrono::duration<double>(chrono::steady_clock::now()-t).count()<(argc>3?stod(argv[3]):150)){uint64_t m=0;
if(trials%4){m=s[rng()%s.size()];int changes=1+rng()%3;for(int k=0;k<changes;k++){int a=rng()%61,b=rng()%61;if((m>>a&1)&&!(m>>b&1))m^=(uint64_t(1)<<a)|(uint64_t(1)<<b);}}
else{while(__builtin_popcountll(m)<8)m|=uint64_t(1)<<(rng()%61);}
int rows=0;for(int r=0;r<15;r++)rows+=bool(m&(uint64_t(15)<<(4*r)));if(rows>5)continue;
int nc=0,u=-1;for(int j=0;j<(int)s.size();j++)if(__builtin_popcountll(m&s[j])>4){u=j;if(++nc>1)break;}
trials++;
if(!nc){s.push_back(m);bucket.emplace_back();}
else if(nc==1){bool done=false;for(auto b:bucket[u]){if(__builtin_popcountll(m&b)>4)continue;bool ok=true;for(int j=0;j<(int)s.size();j++)if(j!=u&&__builtin_popcountll(b&s[j])>4){ok=false;break;}if(ok){s[u]=m;s.push_back(b);bucket[u].clear();bucket.emplace_back();done=true;break;}}
if(!done){if(bucket[u].size()<12)bucket[u].push_back(m);else bucket[u][rng()%12]=m;if(trials%200==0){s[u]=m;bucket[u].clear();swaps++;}}}
if((int)s.size()>best){best=s.size();save();cout<<"supports "<<best<<" trials "<<trials<<" swaps "<<swaps<<endl;}}
save();cout<<"DONE supports "<<s.size()<<" trials "<<trials<<" swaps "<<swaps<<endl;}
