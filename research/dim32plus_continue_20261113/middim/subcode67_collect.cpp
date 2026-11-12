#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;using U=uint64_t;
struct Rec{U a;array<uint32_t,8>old;};
struct Fixed{U s=0,phase=0;vector<U>checks;};
struct Eq{unsigned mask,rhs;};
struct Free{int kind,index;U intersection,chosen;};
unsigned local(U s,U v){unsigned x=0;int j=0;for(int i=0;i<64;i++)if(s>>i&1){if(v>>i&1)x|=1u<<j;j++;}return x;}
U expand(U s,unsigned x){U v=0;int j=0;for(int i=0;i<64;i++)if(s>>i&1){if(x>>j&1)v|=U(1)<<i;j++;}return v;}
struct Basis{
 array<unsigned,16>b{};array<unsigned,16>rhs{};int rank=0;
 bool add(unsigned x,unsigned y=0){for(int j=15;j>=0;j--)if(x>>j&1){if(b[j]){x^=b[j];y^=rhs[j];}else{b[j]=x;rhs[j]=y;rank++;return true;}}return y==0;}
 bool contains(unsigned x)const{for(int j=15;j>=0;j--)if(x>>j&1){if(!b[j])return false;x^=b[j];}return true;}
 unsigned solution()const{unsigned s=0;for(int j=0;j<16;j++)if(b[j])if((__builtin_popcount(b[j]&s)&1)!=rhs[j])s|=1u<<j;return s;}
};
vector<unsigned> eights(unsigned mask){vector<int>p;for(int i=0;i<16;i++)if(mask>>i&1)p.push_back(i);vector<unsigned>out;if(p.size()==8)out.push_back(mask);else if(p.size()==9)for(int i:p)out.push_back(mask^(1u<<i));else if(p.size()==10)for(size_t i=0;i<p.size();i++)for(size_t j=0;j<i;j++)out.push_back(mask^(1u<<p[i])^(1u<<p[j]));return out;}
bool constrain(U a,const Fixed&b,Basis&eq,vector<Eq>&rows){U z=a&b.s;int n=__builtin_popcountll(z);if(n<=8)return true;if(n>10)return false;U h=0;for(U x:b.checks)if((x&z)==x){h=x;break;}if(!h)return false;unsigned r=local(a,h),v=1u^(__builtin_popcountll(b.phase&h)&1);rows.push_back({r,v});return eq.add(r,v);}
int main(int argc,char**argv){
 if(argc<6)return 1;ofstream feasible(argv[5]);ifstream in(argv[1]);int no,nn,ne;in>>no>>nn>>ne;vector<Fixed>old(no),fixed(nn);vector<U>extra(ne);
 auto readfixed=[&](Fixed&b){int n;in>>b.s>>b.phase>>n;b.checks.resize(n);for(U&x:b.checks)in>>x;};
 for(auto&b:old)readfixed(b);for(auto&b:fixed)readfixed(b);for(U&x:extra)in>>x;
 ifstream f(argv[2],ios::binary|ios::ate);size_t bytes=f.tellg();f.seekg(0);vector<Rec>pool(bytes/sizeof(Rec));f.read((char*)pool.data(),bytes);
 mt19937_64 gen(2026111317);shuffle(pool.begin(),pool.end(),gen);int best=0,valid=0,tested=0;auto start=chrono::steady_clock::now();
 for(auto&r:pool){tested++;Basis eq;vector<Eq>forced;vector<Free>free;bool good=true;
  for(unsigned oi:r.old){const auto&b=old[oi];if(b.checks.empty())free.push_back({0,int(oi),r.a&b.s,0});else if(!constrain(r.a,b,eq,forced)){good=false;break;}}
  if(!good)continue;for(auto&b:fixed)if(!constrain(r.a,b,eq,forced)){good=false;break;}if(!good)continue;
  for(int j=0;j<ne;j++){U z=r.a&extra[j];int n=__builtin_popcountll(z);if(n>10){good=false;break;}if(n>8)free.push_back({1,j,z,0});}if(!good)continue;
  valid++;best=0;unsigned phase=eq.solution();vector<vector<unsigned>>choices;for(auto&b:free)choices.push_back(eights(local(r.a,b.intersection)));
  for(int trial=0;trial<32;trial++){
   Basis dual;vector<unsigned>checks;auto add=[&](unsigned x){if(!dual.contains(x)){dual.add(x);checks.push_back(x);}};
   for(auto e:forced)add(e.mask);array<int,16>order;for(int j=0;j<16;j++)order[j]=j;
   shuffle(order.begin(),order.end(),gen);
   if(!forced.empty()&&trial%2==0)stable_sort(order.begin(),order.end(),[&](int a,int b){return bool(forced[0].mask>>a&1)>bool(forced[0].mask>>b&1);});
   add(65535);for(int bit=0;bit<4;bit++){unsigned h=0;for(int j=0;j<16;j++)if(j>>bit&1)h|=1u<<order[j];add(h);}
   auto selected=free;vector<int>ids(free.size());for(int j=0;j<int(ids.size());j++)ids[j]=j;shuffle(ids.begin(),ids.end(),gen);
   for(int id:ids){unsigned h=0;for(auto x:choices[id])if(dual.contains(x)){h=x;break;}
    if(!h){h=choices[id][gen()%choices[id].size()];add(h);}selected[id].chosen=expand(r.a,h);}
   int count=1<<(16-dual.rank);if(count>best){best=count;ofstream out(argv[3]);out<<"{\"support\":"<<r.a<<",\"phase\":"<<phase<<",\"dimension\":"<<(16-dual.rank)<<",\"count\":"<<count<<",\"checks\":[";
    for(size_t j=0;j<checks.size();j++){if(j)out<<",";out<<checks[j];}out<<"],\"forced\":[";
    for(size_t j=0;j<forced.size();j++){if(j)out<<",";out<<"["<<forced[j].mask<<","<<forced[j].rhs<<"]";}out<<"],\"free\":[";
    for(size_t j=0;j<selected.size();j++){if(j)out<<",";auto b=selected[j];out<<"["<<b.kind<<","<<b.index<<","<<b.chosen<<"]";}out<<"],\"tested\":"<<tested<<",\"phase_feasible\":"<<valid<<"}\n";
    cerr<<"BEST count "<<best<<" tested "<<tested<<" phase feasible "<<valid<<" forced "<<forced.size()<<" free "<<free.size()<<"\n";
   }
  }
  {ifstream tmp(argv[3]);string line;getline(tmp,line);feasible<<line<<endl;}
  if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>stod(argv[4]))break;
 }
 cerr<<"DONE best "<<best<<" tested "<<tested<<" feasible "<<valid<<"\n";
}
