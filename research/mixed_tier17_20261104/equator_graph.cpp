#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;using V=array<int8_t,16>;using M=array<int32_t,4>;using N=array<int32_t,3>;using R=array<int32_t,6>;
template<class T>vector<T> read(const char*p){vector<T>v;T x;ifstream f(p,ios::binary);while(f.read((char*)&x,sizeof(x)))v.push_back(x);return v;}
int main(){string r="research/mixed_tier17_20261104/";auto C=read<V>((r+"mixed17_six_all.bin").c_str()),E=read<V>((r+"mixed17_eq_orbits.bin").c_str()),T=read<V>((r+"mixed17_seed_tier.bin").c_str());auto em=read<M>((r+"mixed17_eq_meta.bin").c_str());auto cm=read<N>((r+"mixed17_six_meta.bin").c_str());auto rows=read<R>((r+"mixed17_equator_frontier.bin").c_str());auto H=read<uint16_t>((r+"mixed17_rank4_group.bin").c_str()),P=read<uint16_t>((r+"mixed17_seed_flatplus.bin").c_str()),F=read<uint16_t>((r+"mixed17_seed_flatminus.bin").c_str());if(C.size()!=33600||E.size()!=480||rows.empty())return 1;
 auto maxdot=[&](const V&a,const V&b){int p[16],s=0,best=-1000;for(int i=0;i<16;i++){p[i]=a[i]*b[i];s+=p[i];}for(auto h:H){int d=s;for(int k=0;k<16;k++)if(h>>k&1)d-=2*p[k];if(d>best)best=d;}return best;};
 auto flatmax=[&](const V&a,uint16_t w){int best=-1000;for(auto h:H){int d=0;auto z=w^h;for(int k=0;k<16;k++)d+=a[k]*((z>>k&1)?-1:1);if(d>best)best=d;}return best;};
 vector<int32_t>weights;for(auto row:rows)weights.push_back(2*cm[row[0]][2]);for(auto e:em)weights.push_back(e[2]);for(int i=0;i<32;i++)weights.push_back(32);for(int i=0;i<24;i++)weights.push_back(16);ofstream w((r+"mixed17_switch_weights.bin").c_str(),ios::binary);w.write((char*)weights.data(),weights.size()*4);
 ofstream out((r+"mixed17_switch_edges.bin").c_str(),ios::binary);long long edges=0;auto add=[&](int a,int b){out.write((char*)&a,4);out.write((char*)&b,4);edges++;};int n=rows.size();
 for(int i=0;i<n;i++){auto v=C[rows[i][0]];int count=0;for(int j=0;j<(int)E.size();j++)if(maxdot(v,E[j])>16){add(i,n+j);count+=em[j][2];}int nt=0;for(int j=0;j<32;j++)if(maxdot(v,T[j])>8){add(i,n+480+j);nt++;}int np=0,nm=0;for(int j=0;j<12;j++)if(flatmax(v,P[j])>8){add(i,n+512+j);np++;}for(int j=0;j<12;j++)if(flatmax(v,F[j])>8){add(i,n+524+j);nm++;}if(count!=rows[i][2]||nt!=rows[i][3]||np!=rows[i][4]||nm!=rows[i][5])return 2;
  for(int j=0;j<i;j++)if(maxdot(v,C[rows[j][0]])>8)add(i,j);
 }
 cout<<"candidate_orbits "<<n<<" seed_atoms 536 edges "<<edges<<endl;
}
