// Search simultaneous exchanges through shared blocker sets of size at most four.
// All accepted states use the exact relative-phase table. A separate coordinate
// checker, not this heuristic, certifies the resulting kissing configuration.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
using O=array<int,5>;
int8_t T[253][253][3];
bool conflict(const O&u,const O&v){
 if(u[0]!=v[0]&&u[1]!=v[1]&&u[2]!=v[2])return false;
 for(int t=0;t<3;t++)if(T[u[0]][v[0]][t]+T[u[1]][v[1]][(v[3]-u[3]+t+3)%3]+T[u[2]][v[2]][(v[4]-u[4]+t+3)%3]>6)return true;
 return false;
}
uint64_t key(O v){return (((v[0]*253ULL+v[1])*253+v[2])*3+v[3])*3+v[4];}
struct VH{size_t operator()(const vector<int>&v)const{size_t h=0;for(int i:v)h=(h*1000003)^size_t(i+1);return h;}};
int main(int argc,char**argv){
 if(argc!=3){cerr<<"usage: shared_trades seed.bin seconds\n";return 2;}
 double limit=stod(argv[2]);auto started=chrono::steady_clock::now();
 auto elapsed=[&](){return chrono::duration<double>(chrono::steady_clock::now()-started).count();};
 ifstream tab("T.bin",ios::binary);tab.read((char*)T,sizeof(T));if(!tab)return 3;
 vector<O>vs;O v;ifstream pool("small.bin",ios::binary);while(pool.read((char*)v.data(),sizeof(v)))vs.push_back(v);if(vs.empty())return 4;
 unordered_map<uint64_t,int>lookup;for(int i=0;i<(int)vs.size();i++)lookup[key(vs[i])]=i;
 vector<char>chosen(vs.size());ifstream seed(argv[1],ios::binary);while(seed.read((char*)v.data(),sizeof(v)))chosen[lookup.at(key(v))]=true;
 mt19937 rng(361113);
 auto save=[&](){ofstream f("best_trades.bin",ios::binary);int n=0;for(int i=0;i<(int)vs.size();i++)if(chosen[i]){f.write((char*)vs[i].data(),sizeof(O));n++;}cerr<<"best "<<n<<" seconds "<<elapsed()<<"\n";};
 save();int round=0;
 while(elapsed()<limit){
  vector<int>bylabel[3][253];
  for(int i=0;i<(int)vs.size();i++)if(chosen[i])for(int k=0;k<3;k++)bylabel[k][vs[i][k]].push_back(i);
  unordered_map<vector<int>,vector<int>,VH> groups;
  vector<int>zeros;
  for(int i=0;i<(int)vs.size();i++)if(!chosen[i]){
   vector<int>b;const auto&u=vs[i];
   for(int k=0;k<3&&b.size()<=4;k++)for(int j:bylabel[k][u[k]]){
    const auto&w=vs[j];if(k>0&&u[0]==w[0])continue;if(k>1&&u[1]==w[1])continue;
    if(conflict(u,w))b.push_back(j);if(b.size()>4)break;
   }
   if(b.size()>4)continue;
   sort(b.begin(),b.end());if(b.empty())zeros.push_back(i);else groups[b].push_back(i);
  }
  if(!zeros.empty()){
   for(int i:zeros){bool okay=true;for(int j:zeros)if(chosen[j]&&conflict(vs[i],vs[j])){okay=false;break;}if(okay)chosen[i]=true;}
   save();continue;
  }
  unordered_set<vector<int>,VH> keys;
  unordered_map<int,vector<int>> pair_neighbours;
  for(const auto&g:groups){keys.insert(g.first);if(g.first.size()==2){int a=g.first[0],b=g.first[1];pair_neighbours[a].push_back(b);pair_neighbours[b].push_back(a);}}
  for(auto&g:pair_neighbours){auto&z=g.second;sort(z.begin(),z.end());z.erase(unique(z.begin(),z.end()),z.end());
   for(int a=0;a<(int)z.size();a++)for(int b=a+1;b<(int)z.size();b++){vector<int>w={g.first,z[a],z[b]};sort(w.begin(),w.end());keys.insert(w);}}
  vector<vector<int>> triples;for(const auto&k:keys)if(k.size()==3)triples.push_back(k);
  for(const auto&k:triples)for(int a:k)for(int b:pair_neighbours[a]){
   if(binary_search(k.begin(),k.end(),b))continue;auto w=k;w.push_back(b);sort(w.begin(),w.end());keys.insert(w);
  }
  vector<vector<int>>work(keys.begin(),keys.end());shuffle(work.begin(),work.end(),rng);
  stable_sort(work.begin(),work.end(),[](const auto&a,const auto&b){return a.size()<b.size();});
  cerr<<"round "<<round++<<" blocker groups "<<groups.size()<<" regions "<<work.size()<<" seconds "<<elapsed()<<"\n";
  bool improved=false;int regions=0;
  for(const auto&k:work){
   if((regions++&1023)==0&&elapsed()>=limit)break;
   vector<int>local;
   for(int mask=1;mask<(1<<(int)k.size());mask++){
    vector<int>b;for(int j=0;j<(int)k.size();j++)if(mask>>j&1)b.push_back(k[j]);
    auto it=groups.find(b);if(it!=groups.end())local.insert(local.end(),it->second.begin(),it->second.end());
   }
   int target=k.size()+1;if((int)local.size()<target)continue;
   if(local.size()>128)shuffle(local.begin(),local.end(),rng),local.resize(128);
   int n=local.size();vector<vector<char>>bad(n,vector<char>(n));vector<int>degree(n);
   for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(conflict(vs[local[i]],vs[local[j]])){bad[i][j]=bad[j][i]=true;degree[i]++;degree[j]++;}
   vector<int>avail(n);iota(avail.begin(),avail.end(),0);stable_sort(avail.begin(),avail.end(),[&](int a,int b){return degree[a]<degree[b];});
   vector<int>answer,partial;int nodes=0;
   function<bool(const vector<int>&,int)>dfs=[&](const vector<int>&a,int need){
    if(need==0){answer=partial;return true;}if((int)a.size()<need||++nodes>50000)return false;
    for(int pos=0;pos<=(int)a.size()-need;pos++){
     int v=a[pos];vector<int>next;for(int j=pos+1;j<(int)a.size();j++)if(!bad[v][a[j]])next.push_back(a[j]);
     partial.push_back(v);if(dfs(next,need-1))return true;partial.pop_back();
    }return false;
   };
   if(!dfs(avail,target))continue;
   vector<int>selected;for(int i:answer)selected.push_back(local[i]);
   for(int i:local)if(find(selected.begin(),selected.end(),i)==selected.end()){
    bool okay=true;for(int j:selected)if(conflict(vs[i],vs[j])){okay=false;break;}if(okay)selected.push_back(i);
   }
   for(int i:k){if(!chosen[i])return 5;chosen[i]=false;}for(int i:selected)chosen[i]=true;
   cerr<<"accepted "<<k.size()<<" -> "<<selected.size()<<" checked regions "<<regions<<"\n";save();improved=true;break;
  }
  if(!improved){cerr<<"no improving searched region; seconds "<<elapsed()<<"\n";break;}
 }
 save();
}
