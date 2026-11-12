// Break the simultaneous order-three phase symmetry: select individual points.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
using namespace std;
using O=array<int,5>;int8_t T[253][253][3];
uint64_t key(O v){return (((v[0]*253ULL+v[1])*253+v[2])*3+v[3])*3+v[4];}
struct Set{vector<int>data,pos;Set(int n):pos(n,-1){}void erase(int i){int p=pos[i];if(p<0)return;int j=data.back();data[p]=j;pos[j]=p;data.pop_back();pos[i]=-1;}void insert(int i){if(pos[i]>=0)return;pos[i]=data.size();data.push_back(i);}};
int main(int argc,char**argv){
 if(argc!=4){cerr<<"usage: refill_points selected_orbits.bin seconds rng_seed\n";return 2;}
 ifstream tab("T.bin",ios::binary);tab.read((char*)T,sizeof(T));if(!tab)return 3;
 vector<O>vs;O v;ifstream pool("small.bin",ios::binary);while(pool.read((char*)v.data(),sizeof(v)))vs.push_back(v);if(vs.empty())return 4;
 int K=vs.size(),N=3*K;unordered_map<uint64_t,int>lookup;for(int i=0;i<K;i++)lookup[key(vs[i])]=i;
 vector<vector<int>>adj(N);vector<int>rows[3][253];long long edges=0;
 for(int i=0;i<K;i++){
  const O&u=vs[i];
  for(int k=0;k<3;k++){
   for(int j:rows[k][u[k]]){
    const O&w=vs[j];if(k>0&&u[0]==w[0])continue;if(k>1&&u[1]==w[1])continue;
    for(int t=0;t<3;t++)if(T[u[0]][w[0]][t]+T[u[1]][w[1]][(w[3]-u[3]+t+3)%3]+T[u[2]][w[2]][(w[4]-u[4]+t+3)%3]>6){
     for(int r=0;r<3;r++){int a=3*i+r,b=3*j+(r+t)%3;adj[a].push_back(b);adj[b].push_back(a);edges++;}
    }
   }
   rows[k][u[k]].push_back(i);
  }
 }
 cerr<<"vertices "<<N<<" edges "<<edges<<"\n";
 vector<char>chosen(N);vector<int>count(N),best;Set zeros(N),ones(N);int size=0;
 auto refresh=[&](int i){zeros.erase(i);ones.erase(i);if(!chosen[i]){if(count[i]==0)zeros.insert(i);if(count[i]==1)ones.insert(i);}};
 auto add=[&](int i){if(chosen[i])abort();chosen[i]=true;size++;refresh(i);for(int j:adj[i]){count[j]++;refresh(j);}};
 auto del=[&](int i){if(!chosen[i])abort();chosen[i]=false;size--;refresh(i);for(int j:adj[i]){count[j]--;refresh(j);}};
 ifstream seed(argv[1],ios::binary);while(seed.read((char*)v.data(),sizeof(v))){int i=lookup.at(key(v));for(int r=0;r<3;r++)add(3*i+r);}
 for(int i=0;i<N;i++)refresh(i);
 mt19937 rng(stoul(argv[3]));auto start=chrono::steady_clock::now();double seconds=stod(argv[2]);long long steps=0;
 auto save=[&](){if(size>(int)best.size()){best.clear();ofstream f("best_points.bin",ios::binary);for(int i=0;i<N;i++)if(chosen[i]){best.push_back(i);array<int,6>point;copy(vs[i/3].begin(),vs[i/3].end(),point.begin());point[5]=i%3;f.write((char*)point.data(),sizeof(point));}cerr<<"best "<<size<<" steps "<<steps<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";}};
 auto conflict=[&](int a,int b){const auto&u=vs[a/3];const auto&w=vs[b/3];int t=(b%3-a%3+3)%3;return T[u[0]][w[0]][t]+T[u[1]][w[1]][(w[3]-u[3]+t+3)%3]+T[u[2]][w[2]][(w[4]-u[4]+t+3)%3]>6;};
 while(true){
  if((steps&1023)==0&&chrono::duration<double>(chrono::steady_clock::now()-start).count()>=seconds)break;
  while(!zeros.data.empty()){
   int i=zeros.data[rng()%zeros.data.size()];for(int t=0;t<3;t++){int j=zeros.data[rng()%zeros.data.size()];if(adj[j].size()<adj[i].size())i=j;}add(i);
  }
  save();
  if(steps%50==0){unordered_map<int,vector<int>>groups;bool gain=false;
   for(int i:ones.data){int owner=-1;for(int j:adj[i])if(chosen[j]){owner=j;break;}auto&g=groups[owner];for(int k:g)if(!conflict(i,k)){del(owner);add(i);add(k);gain=true;break;}if(gain)break;g.push_back(i);}if(gain){steps++;continue;}}
  if(!ones.data.empty()&&steps%4000!=3999){
   int i=ones.data[rng()%ones.data.size()];int owner=-1;for(int j:adj[i])if(chosen[j]){owner=j;break;}del(owner);add(i);
  }else{
   for(int i=0;i<N;i++)if(chosen[i])del(i);for(int i:best)add(i);
   for(int t=0;t<1+(int)(rng()%4);t++){int i=rng()%N;if(chosen[i])continue;for(int j:adj[i])if(chosen[j])del(j);add(i);}
  }
  steps++;
 }
 cerr<<"finished steps="<<steps<<" best="<<best.size()<<"\n";
}
