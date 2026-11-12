// Exact finite support-pool generation followed by heuristic large-neighborhood search.
// A search failure proves no upper bound. All candidate acceptances use integer incidence.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_set>
#include <vector>
#include <omp.h>
using namespace std;
uint32_t bc[64][9];
vector<array<int,5>> subidx;
array<uint32_t,56> subs(uint64_t v){
  int p[8],k=0;while(v){p[k++]=__builtin_ctzll(v);v&=v-1;}
  if(k!=8)abort();array<uint32_t,56> ans;
  for(int j=0;j<56;j++){auto a=subidx[j];ans[j]=bc[p[a[0]]][1]+bc[p[a[1]]][2]+bc[p[a[2]]][3]+bc[p[a[3]]][4]+bc[p[a[4]]][5];}
  return ans;
}
int main(int argc,char**argv){
  if(argc<6){cerr<<"dimension seed-file output-prefix maximum-initial-blockers seconds [randomseed]\n";return 2;}
  int n=atoi(argv[1]),limit=atoi(argv[4]);double seconds=atof(argv[5]);string prefix=argv[3];
  for(int i=0;i<64;i++){bc[i][0]=1;for(int j=1;j<9;j++)bc[i][j]=(i?bc[i-1][j]+bc[i-1][j-1]:0);}
  for(int a=0;a<4;a++)for(int b=a+1;b<5;b++)for(int c=b+1;c<6;c++)for(int d=c+1;d<7;d++)for(int e=d+1;e<8;e++)subidx.push_back({a,b,c,d,e});
  vector<uint64_t> seed;uint64_t v;ifstream in(argv[2]);while(in>>v)seed.push_back(v);
  vector<int> owner(bc[n][5],-1);for(int i=0;i<(int)seed.size();i++)for(auto t:subs(seed[i])){if(owner[t]!=-1)abort();owner[t]=i;}
  vector<uint64_t> pool=seed;
  auto started=chrono::steady_clock::now();
  string cache=prefix+".pool.bin";ifstream cached(cache,ios::binary);
  if(cached){pool.clear();while(cached.read((char*)&v,8))pool.push_back(v);cout<<"cached pool "<<pool.size()<<endl;}
  else{
    vector<uint64_t> extra;
    #pragma omp parallel
    {
      vector<uint64_t> local;
      #pragma omp for schedule(dynamic,1)
      for(int a=0;a<n-7;a++){
        for(int b=a+1;b<n-6;b++)for(int c=b+1;c<n-5;c++)for(int d=c+1;d<n-4;d++)for(int e=d+1;e<n-3;e++)for(int f=e+1;f<n-2;f++)for(int g=f+1;g<n-1;g++)for(int h=g+1;h<n;h++){
          int p[8]={a,b,c,d,e,f,g,h};int blockers[64],nb=0;
          for(auto q:subidx){int t=bc[p[q[0]]][1]+bc[p[q[1]]][2]+bc[p[q[2]]][3]+bc[p[q[3]]][4]+bc[p[q[4]]][5];int o=owner[t];if(o<0)continue;bool found=false;for(int j=0;j<nb;j++)if(blockers[j]==o){found=true;break;}if(!found){blockers[nb++]=o;if(nb>limit)break;}}
          if(nb<=limit){uint64_t x=0;for(int z:p)x|=1ULL<<z;local.push_back(x);}
        }
      }
      #pragma omp critical
      extra.insert(extra.end(),local.begin(),local.end());
    }
    sort(extra.begin(),extra.end());unordered_set<uint64_t> old(seed.begin(),seed.end());for(auto x:extra)if(!old.count(x))pool.push_back(x);
    ofstream save(cache,ios::binary);save.write((char*)pool.data(),pool.size()*8);
    cout<<"generated pool "<<pool.size()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;
  }
  int m=pool.size();vector<array<uint32_t,56>> ss(m);vector<int> sizes(bc[n][5]);
  for(int i=0;i<m;i++){ss[i]=subs(pool[i]);for(auto x:ss[i])sizes[x]++;}
  vector<int> offsets(sizes.size()+1);for(int i=0;i<(int)sizes.size();i++)offsets[i+1]=offsets[i]+sizes[i];
  vector<int> members(offsets.back()),head=offsets;
  for(int i=0;i<m;i++)for(auto x:ss[i])members[head[x]++]=i;
  cout<<"incidence entries "<<members.size()<<" memory_bytes "<<members.size()*4+ss.size()*sizeof(ss[0])<<endl;
  vector<int> cnt(m),chosen,position(m,-1),tabu(m),best;vector<uint8_t>sel(m);fill(owner.begin(),owner.end(),-1);
  auto add=[&](int id){if(cnt[id])abort();sel[id]=1;position[id]=chosen.size();chosen.push_back(id);for(auto x:ss[id]){if(owner[x]!=-1)abort();owner[x]=id;for(int a=offsets[x];a<offsets[x+1];a++)cnt[members[a]]++;}};
  auto remove=[&](int id,vector<int>&zeros,int step){if(!sel[id])abort();sel[id]=0;tabu[id]=step+5;int p=position[id];int last=chosen.back();chosen[p]=last;position[last]=p;chosen.pop_back();position[id]=-1;for(auto x:ss[id]){if(owner[x]!=id)abort();owner[x]=-1;for(int a=offsets[x];a<offsets[x+1];a++){int j=members[a];if(--cnt[j]==0)zeros.push_back(j);}}};
  for(int i=0;i<(int)seed.size();i++)add(i);best=chosen;
  mt19937 rng(argc>6?atoi(argv[6]):20261113);int bestn=chosen.size(),lastsave=0;started=chrono::steady_clock::now();
  auto save=[&](){ofstream o(prefix+".best.txt");for(int i:best)o<<pool[i]<<'\n';};save();
  for(int step=1;chrono::duration<double>(chrono::steady_clock::now()-started).count()<seconds;step++){
    vector<int> zero,removed;
    int target=-1,score=INT32_MAX;
    for(int t=0;t<96;t++){
      int j=rng()%m;if(sel[j]||tabu[j]>step)continue;
      // Distinct blockers, rather than five-subset incidence multiplicity,
      // are the actual deletion price. A seven-coordinate overlap is one
      // blocker, not twenty-one blockers.
      int ids[56],k=0;for(auto x:ss[j]){int id=owner[x];if(id<0)continue;bool seen=false;for(int u=0;u<k;u++)if(ids[u]==id){seen=true;break;}if(!seen)ids[k++]=id;}
      int val=10*k+(rng()%20);if(val<score){target=j;score=val;}
    }
    if(target<0)continue;
    for(auto x:ss[target]){int id=owner[x];if(id>=0&&find(removed.begin(),removed.end(),id)==removed.end())removed.push_back(id);}
    if(step%25==0){int num=10+rng()%60;while((int)removed.size()<num){int id=chosen[rng()%chosen.size()];if(find(removed.begin(),removed.end(),id)==removed.end())removed.push_back(id);}}
    for(int id:removed)remove(id,zero,step);add(target);
    shuffle(zero.begin(),zero.end(),rng);
    for(int id:zero)if(!sel[id]&&!cnt[id]&&tabu[id]<=step)add(id);
    for(int id:removed)if(!sel[id]&&!cnt[id])add(id);
    if((int)chosen.size()>bestn){bestn=chosen.size();best=chosen;save();lastsave=step;cout<<"BEST "<<bestn<<" step "<<step<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;}
    if(step%1000==0)cout<<"progress step "<<step<<" current "<<chosen.size()<<" best "<<bestn<<endl;
    if((int)chosen.size()<bestn-8 || step%5000==0){vector<int> tmp=chosen;for(int id:tmp)remove(id,zero,step);fill(tabu.begin(),tabu.end(),0);for(int id:best)add(id);}
  }
  cout<<"FINAL best "<<bestn<<" source "<<seed.size()<<endl;
}
