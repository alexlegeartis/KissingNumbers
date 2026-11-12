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
struct Vertex {int p; double alpha,w; vector<double> h,z;};
int main(int argc,char**argv){
 if(argc<4)return 1;ifstream f(argv[1]);int n,k,r;f>>n>>k>>r;vector<Vertex>v(n);
 for(auto &a:v){f>>a.p>>a.alpha>>a.w;a.h.resize(r);for(auto &x:a.h)f>>x;a.z.resize(k);for(auto &x:a.z)f>>x;}
 int nw=(n+63)/64;vector<vector<uint64_t>>adj(n,vector<uint64_t>(nw));vector<int>deg(n);
 for(int i=0;i<n;i++)for(int j=0;j<i;j++){double dot=0;for(int l=0;l<k;l++)dot+=v[i].z[l]*v[j].z[l];double hh=0;for(int l=0;l<r;l++)hh+=v[i].h[l]*v[j].h[l];double p=v[i].alpha*v[j].alpha*((v[i].p==v[j].p?18:9)-dot)+hh;
 if(p<=9+1e-10){adj[i][j/64]|=1ULL<<(j%64);adj[j][i/64]|=1ULL<<(i%64);deg[i]++;deg[j]++;}}
 mt19937_64 gen(2026111301);uniform_real_distribution<double>U(0,1);double best=-1;vector<int>bestset;int maxiter=stoi(argv[3]);
 vector<uint64_t>avail(nw);vector<double>score(n);auto begin=chrono::steady_clock::now();
 for(int it=0;it<maxiter;it++){
  fill(avail.begin(),avail.end(),~0ULL);if(n%64)avail.back()=(1ULL<<(n%64))-1;
  vector<int>cur;double value=0;double noise=(it%4==0?.08:it%4==1?.4:it%4==2?1.0:2.5);
  for(int i=0;i<n;i++)score[i]=v[i].w*exp(noise*(U(gen)-.5));
  if(it>0 && it%3){auto seed=bestset;shuffle(seed.begin(),seed.end(),gen);int keep=max(0,(int)seed.size()-1-(int)(U(gen)*seed.size()));
   for(int a=0;a<keep;a++){int x=seed[a];cur.push_back(x);value+=v[x].w;for(int w=0;w<nw;w++)avail[w]&=adj[x][w];}}
  while(true){int bi=-1;double bs=-1;
   for(int w=0;w<nw;w++){uint64_t bits=avail[w];while(bits){int j=__builtin_ctzll(bits)+64*w;bits&=bits-1;if(score[j]>bs){bs=score[j];bi=j;}}}
   if(bi<0)break;cur.push_back(bi);value+=v[bi].w;for(int w=0;w<nw;w++)avail[w]&=adj[bi][w];
  }
  if(value>best+.1){best=value;bestset=cur;ofstream out(argv[2]);out.precision(17);out<<"{\"estimate\":"<<value<<",\"iteration\":"<<it<<",\"vertices\":[";for(int j=0;j<(int)cur.size();j++){if(j)out<<",";out<<cur[j];}out<<"]}\n";cerr<<"it "<<it<<" best "<<llround(best)<<" slices "<<cur.size()<<"\n";}
  if(chrono::duration<double>(chrono::steady_clock::now()-begin).count()> (argc>4?stod(argv[4]):45))break;
 }
}
