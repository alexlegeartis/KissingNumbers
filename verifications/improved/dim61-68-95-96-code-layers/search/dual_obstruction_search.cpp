// Finite fill/correct search for supports in thirteen or fourteen rows.
// Search output is untrusted until both independent certificate replays pass.
// argv: start, dual pairs, output, trials, neutral period, single duals.
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
struct S{uint64_t lo,hi,dual,dual2=0;};
int pop(uint64_t x){return __builtin_popcountll(x);}
int cell(S x,int i){return i<16?(x.lo>>(4*i))&15:(x.hi>>(4*(i-16)))&15;}
void put(S &x,int i,int v){if(i<16)x.lo=(x.lo&~(15ULL<<(4*i)))|(uint64_t(v)<<(4*i));else x.hi=(x.hi&~(15ULL<<(4*(i-16))))|(uint64_t(v)<<(4*(i-16)));}
int dot(S a,S b){return pop(a.lo&b.lo)+pop(a.hi&b.hi);}
int score(S x){const int t[]={16384,0,8192,0,6144,0,5120,0,8960};return t[pop(x.lo&255)];}
void save(const std::vector<S>&f,const char *name){std::ofstream o(name);o<<f.size()<<'\n';for(auto x:f)o<<x.lo<<' '<<x.hi<<' '<<x.dual<<' '<<x.dual2<<'\n';}
int main(int argc,char **argv){
 if(argc!=7)return 2;std::ifstream in(argv[1]),du(argv[2]);int n,m;in>>n;du>>m;std::vector<S>f(n);for(auto&x:f)in>>x.lo>>x.hi>>x.dual>>x.dual2;
 std::vector<std::pair<uint64_t,uint64_t>>ds(m);for(auto &d:ds)du>>d.first>>d.second;std::ifstream d13(argv[6]);int n13;d13>>n13;std::vector<uint64_t>ds13(n13);for(auto &d:ds13)d13>>d;
 std::mt19937_64 rng(281313);uint64_t trials=std::stoull(argv[4]);unsigned period=std::stoul(argv[5]);
 std::vector<std::vector<S>> buckets(f.size());uint64_t adds=0,trades=0,neutral=0,one=0,zero=0;long long total=0;for(auto x:f)total+=score(x);
 auto started=std::chrono::steady_clock::now();
 for(uint64_t trial=0;trial<trials;trial++){
  bool large=rng()%2;auto [d,e]=ds[rng()%ds.size()];if(!large){d=ds13[rng()%ds13.size()];e=0;}int occupied=large?14:13,single=large?4:2;S x{0,0,d,e};std::vector<int>occ;
  for(int i=0;i<24;i++)if(((d|e)>>(2*i))&3)occ.push_back(i);
  while(occ.size()<occupied){int i=rng()%24;if(std::find(occ.begin(),occ.end(),i)==occ.end())occ.push_back(i);}
  std::shuffle(occ.begin(),occ.end(),rng);std::stable_partition(occ.begin(),occ.end(),[&](int i){return (((d|e)>>(2*i))&3)==1;});int parity[2]={0,0};
  for(int pos=0;pos<occupied;pos++){
   int i=occ[pos],a=((d|e)>>(2*i))&3,v=0,k=((e>>(2*i))&3)?1:0;
   if(pos<single){int c=rng()%4;v=1<<c;if(a==1)parity[k]^=(c>>1);if(a==2)parity[k]^=((c&1)^((c>>1)&1));if(a==3)parity[k]^=(c&1);}
   else if(a){int b=rng()%2;v=a==2?(b?6:9):(b?10:5);parity[k]^=b;}
   else{const int p[]={5,6,9,10};v=p[rng()%4];}
   put(x,i,v);
  }
  for(int k=0;k<(large?2:1);k++)if(!parity[k]){uint64_t dk=k?e:d;for(int pos=0;pos<occupied;pos++){int i=occ[pos],a=(dk>>(2*i))&3;if(a){int v=cell(x,i);if(pos>=single)put(x,i,v^15);else{int c=__builtin_ctz(v);c^=(a==1?2:1);put(x,i,1<<c);}break;}}}
  int blocker=-1;bool reject=false;
  for(int j=0;j<int(f.size());j++)if(dot(x,f[j])>12){if(blocker>=0){reject=true;break;}blocker=j;}
  if(reject)continue;
  if(blocker<0){zero++;f.push_back(x);buckets.emplace_back();total+=score(x);adds++;save(f,argv[3]);std::cout<<"add "<<trial<<" size "<<f.size()<<" score "<<total<<std::endl;}
  else{
   one++;bool traded=false;
   for(auto y:buckets[blocker]){
    if(dot(x,y)>12||score(x)+score(y)<score(f[blocker]))continue;
    bool ok=true;for(int j=0;j<int(f.size());j++)if(j!=blocker&&dot(y,f[j])>12){ok=false;break;}
    if(!ok)continue;
    total+=score(x)+score(y)-score(f[blocker]);f[blocker]=x;f.push_back(y);buckets[blocker].clear();buckets.emplace_back();trades++;traded=true;save(f,argv[3]);
    std::cout<<"trade "<<trial<<" size "<<f.size()<<" score "<<total<<std::endl;break;
   }
   if(!traded){auto &b=buckets[blocker];if(b.size()<64)b.push_back(x);else b[rng()%b.size()]=x;
    if(rng()%period==0&&score(x)>=score(f[blocker])){total+=score(x)-score(f[blocker]);f[blocker]=x;b.clear();neutral++;}}
  }
  if(trial%100000==0)std::cout<<"trial "<<trial<<" n "<<f.size()<<" singles "<<one<<" free "<<zero<<std::endl;
 }
 save(f,argv[3]);std::cout<<"done n "<<f.size()<<" adds "<<adds<<" trades "<<trades<<" neutral "<<neutral<<" singles "<<one<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;
}
