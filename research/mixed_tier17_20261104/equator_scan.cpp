#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;using V=array<int8_t,16>;using M=array<int32_t,4>;using N=array<int32_t,3>;
template<class T>vector<T> read(const char*p){vector<T>v;T x;ifstream f(p,ios::binary);while(f.read((char*)&x,sizeof(x)))v.push_back(x);return v;}
int main(){const char*R="research/mixed_tier17_20261104/";auto C=read<V>((string(R)+"mixed17_six_all.bin").c_str()),E=read<V>((string(R)+"mixed17_eq_orbits.bin").c_str()),T=read<V>((string(R)+"mixed17_seed_tier.bin").c_str());auto em=read<M>((string(R)+"mixed17_eq_meta.bin").c_str());auto cm=read<N>((string(R)+"mixed17_six_meta.bin").c_str());auto H=read<uint16_t>((string(R)+"mixed17_rank4_group.bin").c_str()),P=read<uint16_t>((string(R)+"mixed17_seed_flatplus.bin").c_str()),M=read<uint16_t>((string(R)+"mixed17_seed_flatminus.bin").c_str());if(C.size()!=33600||E.size()!=480||H.size()!=16||T.size()!=32||P.size()!=12||M.size()!=12)return 1;
 auto maxdot=[&](const V&a,const V&b){int p[16],sum=0,best=-1000;for(int i=0;i<16;i++){p[i]=a[i]*b[i];sum+=p[i];}for(auto h:H){int d=sum;for(int k=0;k<16;k++)if(h>>k&1)d-=2*p[k];if(d>best)best=d;}return best;};
 auto flatmax=[&](const V&a,uint16_t w){int best=-1000;for(auto h:H){int d=0;uint16_t z=w^h;for(int k=0;k<16;k++)d+=a[k]*((z>>k&1)?-1:1);if(d>best)best=d;}return best;};
 ofstream out((string(R)+"mixed17_equator_frontier.bin").c_str(),ios::binary);map<int,int>hist;int valid=0,best=-100000;
 for(int i=0;i<(int)C.size();i++){
  int self=0;bool good=true;for(auto h:H){int d=0;bool changed=false;for(int k=0;k<16;k++){int z=C[i][k];if(z&&(h>>k&1)){z=-z;changed=true;}d+=C[i][k]*z;}if(changed&&d>8){good=false;break;}}if(!good)continue;valid++;
  int eqcost=0;for(int j=0;j<(int)E.size();j++)if(maxdot(C[i],E[j])>16)eqcost+=em[j][2];
  int tier=0;bool existing=false;for(auto&v:T){int d=maxdot(C[i],v);if(d==24){existing=true;break;}if(d>8)tier++;}if(existing)continue;
  int fp=0,fm=0;for(auto w:P)if(flatmax(C[i],w)>8)fp++;for(auto w:M)if(flatmax(C[i],w)>8)fm++;
  int benefit=2*cm[i][2],cost=eqcost+2*16*tier+16*(fp+fm);int net=benefit-cost;hist[net]++;if(net>best){best=net;cout<<"candidate "<<i<<" net "<<net<<" eqcost "<<eqcost<<" tier_each "<<tier<<" flat "<<fp<<","<<fm<<" orbit "<<cm[i][2]<<endl;}
  if(net>=-256){array<int32_t,6>row{i,net,eqcost,tier,fp,fm};out.write((char*)row.data(),sizeof(row));}
 }
 cout<<"valid "<<valid<<" best "<<best<<" histogram_top";int n=0;for(auto it=hist.rbegin();it!=hist.rend()&&n<15;++it,++n)cout<<" "<<it->first<<":"<<it->second;cout<<endl;
}
