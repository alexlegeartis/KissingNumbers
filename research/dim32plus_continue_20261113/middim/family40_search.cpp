#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc<4)return 1;ifstream f(argv[1]);int n;f>>n;vector<int>mu(6561),pw(8);vector<vector<int>>dig(6561,vector<int>(8));
 for(auto&x:pw)f>>x;for(auto&x:mu)f>>x;for(auto&row:dig)for(auto&x:row)f>>x;
 vector<int>cl(n),r(n),weight(n);for(int i=0;i<n;i++)f>>cl[i]>>r[i]>>weight[i];
 int nw=(n+63)/64;vector<vector<uint64_t>>adj(n,vector<uint64_t>(nw));
 for(int i=0;i<n;i++)for(int j=0;j<i;j++){
  int delta=0;for(int k=0;k<8;k++)delta+=((dig[cl[i]][k]-dig[cl[j]][k]+3)%3)*pw[k];
  int b=r[i]+r[j]-18+mu[delta];
  if(b<=0||b*b<=r[i]*r[j]){adj[i][j/64]|=1ULL<<(j%64);adj[j][i/64]|=1ULL<<(i%64);}
 }
 mt19937_64 gen(20261113);uniform_real_distribution<double>U(0,1);long long best=-1;vector<int>bestset;vector<double>score(n);vector<uint64_t>avail(nw);auto begin=chrono::steady_clock::now();
 for(int it=0;;it++){
  fill(avail.begin(),avail.end(),~0ULL);if(n%64)avail.back()=(1ULL<<(n%64))-1;
  double noise=(it%5==0?.1:it%5==1?.5:it%5==2?1.2:it%5==3?2.5:5.);
  for(int i=0;i<n;i++)score[i]=weight[i]*exp(noise*(U(gen)-.5));
  vector<int>cur;long long value=0;
  if(it && it%3){auto seed=bestset;shuffle(seed.begin(),seed.end(),gen);int keep=max(0,int(seed.size())-1-int(U(gen)*seed.size()));
   for(int a=0;a<keep;a++){int x=seed[a];cur.push_back(x);value+=weight[x];for(int w=0;w<nw;w++)avail[w]&=adj[x][w];}}
  while(true){int bi=-1;double bs=-1;for(int w=0;w<nw;w++){auto bits=avail[w];while(bits){int j=__builtin_ctzll(bits)+64*w;bits&=bits-1;if(score[j]>bs){bs=score[j];bi=j;}}}if(bi<0)break;
   cur.push_back(bi);value+=weight[bi];for(int w=0;w<nw;w++)avail[w]&=adj[bi][w];}
  if(value>best){best=value;bestset=cur;ofstream out(argv[2]);out<<"{\"value\":"<<value<<",\"iteration\":"<<it<<",\"vertices\":[";for(size_t j=0;j<cur.size();j++){if(j)out<<",";out<<cur[j];}out<<"]}\n";cerr<<"BEST "<<best<<" families "<<cur.size()<<"\n";}
  if(chrono::duration<double>(chrono::steady_clock::now()-begin).count()>stod(argv[3]))break;
 }
}
