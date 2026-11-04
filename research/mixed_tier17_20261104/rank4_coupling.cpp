#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(){vector<array<int8_t,16>> V;array<int8_t,16>v;ifstream f("research/mixed_tier17_20261104/mixed17_rank4_orbits.bin",ios::binary);while(f.read((char*)v.data(),16))V.push_back(v);vector<uint16_t> H,F;uint16_t h;ifstream g("research/mixed_tier17_20261104/mixed17_rank4_group.bin",ios::binary),p("research/mixed_tier17_20261104/mixed17_rank4_flat_reps.bin",ios::binary);while(g.read((char*)&h,2))H.push_back(h);while(p.read((char*)&h,2))F.push_back(h);if(V.size()!=3840||H.size()!=16||F.size()!=128)return 1;
 vector<uint8_t> same(V.size()*128),opposite(V.size()*128);for(int i=0;i<(int)V.size();i++)for(int j=0;j<128;j++){int maxd=-100;for(auto h:H){uint16_t w=F[j]^h;int d=0;for(int k=0;k<16;k++)d+=V[i][k]*((w>>k&1)?-1:1);maxd=max(maxd,d);}same[128*i+j]=maxd>8;opposite[128*i+j]=maxd>16;}
 vector<uint8_t> flat(128*128);for(int i=0;i<128;i++)for(int j=0;j<128;j++){int d=16;for(auto h:H)d=min(d,__builtin_popcount((unsigned)(F[i]^F[j]^h)));if(i!=j)flat[128*i+j]=d<6;else if(d!=0)return 2;}
 ofstream a("research/mixed_tier17_20261104/mixed17_rank4_same.bin",ios::binary);a.write((char*)same.data(),same.size());ofstream b("research/mixed_tier17_20261104/mixed17_rank4_opposite.bin",ios::binary);b.write((char*)opposite.data(),opposite.size());ofstream c("research/mixed_tier17_20261104/mixed17_rank4_flat_graph.bin",ios::binary);c.write((char*)flat.data(),flat.size());long x=0,y=0,z=0;for(auto q:same)x+=q;for(auto q:opposite)y+=q;for(auto q:flat)z+=q;cout<<"same conflicts "<<x<<" opposite conflicts "<<y<<" flat edges "<<z/2<<endl;
}
