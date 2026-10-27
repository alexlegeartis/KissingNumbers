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
int maximum_dot=6; using O=array<int,5>;int8_t T[253][253][3];
bool conflict(O u,O v){
 if(u[0]!=v[0]&&u[1]!=v[1]&&u[2]!=v[2])return false;
 for(int t=0;t<3;t++)if(T[u[0]][v[0]][t]+T[u[1]][v[1]][(v[3]-u[3]+t+3)%3]+T[u[2]][v[2]][(v[4]-u[4]+t+3)%3]>maximum_dot)return true;
 return false;
}
long long key(O v){return (((v[0]*253LL+v[1])*253+v[2])*3+v[3])*3+v[4];}
struct Set {vector<int> data,pos;Set(int n):pos(n,-1){};void erase(int i){int p=pos[i];if(p<0)return;int j=data.back();data[p]=j;pos[j]=p;data.pop_back();pos[i]=-1;}void insert(int i){if(pos[i]>=0)return;pos[i]=data.size();data.push_back(i);}};
int main(int argc,char**argv){
 if(argc!=4){cerr<<"usage: refill steps seed maximum_dot\n";return 2;} maximum_dot=stoi(argv[3]);
 ifstream table("T.bin",ios::binary);table.read((char*)T,sizeof(T));
 vector<O> vs;ifstream pool("small.bin",ios::binary);O v;while(pool.read((char*)v.data(),sizeof(v)))vs.push_back(v);int N=vs.size();
 unordered_map<long long,int> lookup;for(int i=0;i<N;i++)lookup[key(vs[i])]=i;
 vector<vector<int>> adj(N);vector<int> rows[3][253];long long E=0;
 for(int i=0;i<N;i++){
  auto u=vs[i];for(int k=0;k<3;k++){
   for(int j:rows[k][u[k]]){auto w=vs[j];if(k>0&&u[0]==w[0])continue;if(k>1&&u[1]==w[1])continue;
    if(conflict(u,w)){adj[i].push_back(j);adj[j].push_back(i);E++;}}
   rows[k][u[k]].push_back(i);
  }
 }
 cerr<<"vertices "<<N<<" edges "<<E<<"\n";
 vector<char> chosen(N);vector<int> count(N),best;Set zeros(N),ones(N);int size=0;
 auto refresh=[&](int i){zeros.erase(i);ones.erase(i);if(!chosen[i]){if(count[i]==0)zeros.insert(i);if(count[i]==1)ones.insert(i);}};
 auto add=[&](int i){chosen[i]=true;size++;refresh(i);for(int j:adj[i]){count[j]++;refresh(j);}};
 auto del=[&](int i){chosen[i]=false;size--;refresh(i);for(int j:adj[i]){count[j]--;refresh(j);}};
 ifstream seed("selected.bin",ios::binary);while(seed.read((char*)v.data(),sizeof(v)))add(lookup.at(key(v)));
 for(int i=0;i<N;i++)refresh(i);
 mt19937 rng(argc>2?stoi(argv[2]):361028);auto start=chrono::steady_clock::now();long long max_steps=stoll(argv[1]);int steps=0;
 auto save=[&](){if(size>(int)best.size()){best.clear();for(int i=0;i<N;i++)if(chosen[i])best.push_back(i);ofstream out("best.bin",ios::binary);for(int i:best)out.write((char*)vs[i].data(),sizeof(O));cerr<<"best "<<size<<" steps "<<steps<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";}};
 while(steps<max_steps){
  while(!zeros.data.empty()){
   int i=zeros.data[rng()%zeros.data.size()];for(int t=0;t<3;t++){int j=zeros.data[rng()%zeros.data.size()];if(adj[j].size()<adj[i].size())i=j;}add(i);
  }
  save();
  // Explicit (1,2) improvements.
  if(steps%50==0){
   unordered_map<int,vector<int>> owners;bool gain=false;
   for(int i:ones.data){int owner=-1;for(int j:adj[i])if(chosen[j]){owner=j;break;}auto &group=owners[owner];for(int k:group)if(!conflict(vs[i],vs[k])){del(owner);add(i);add(k);gain=true;break;}if(gain)break;group.push_back(i);}
   if(gain){steps++;continue;}
  }
  if(!ones.data.empty() && steps%4000!=3999){
   int i=ones.data[rng()%ones.data.size()];int owner=-1;for(int j:adj[i])if(chosen[j]){owner=j;break;}del(owner);add(i);
  }else{
   // Restore the incumbent, then force one to three outside orbits.
   for(int i=0;i<N;i++)if(chosen[i])del(i);for(int i:best)add(i);
   for(int t=0;t<1+(int)(rng()%3);t++){int i=rng()%N;if(chosen[i])continue;for(int j:adj[i])if(chosen[j])del(j);add(i);}
  }
  steps++;
 }
 cerr<<"finished steps="<<steps<<" best="<<best.size()<<"\n";
}
