#include <array>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
using U=std::uint64_t;
using W=std::uint32_t;
int main(int argc,char**argv){
 assert(argc==4);
 auto mono=[](int m,int n){
  U v=0;for(int i=0;i<n;i++) if((i&m)==m)v|=(U(1)<<i);
  return v;
 };
 auto unpack=[&](unsigned mask, const std::vector<U>& base){
  U out=0;for(;mask;mask&=mask-1){auto j=std::countr_zero(mask);assert(j<base.size());out^=base[j];}
  return out;
 };
 std::vector<U> c6,c5,qb;
 for(int i=0;i<64;i++)if(std::popcount(unsigned(i))==3)c6.push_back(mono(i,64));
 for(int i=0;i<32;i++)if(std::popcount(unsigned(i))<=2)qb.push_back(mono(i,32));else if(std::popcount(unsigned(i))==3)c5.push_back(mono(i,32));
 assert(c6.size()==20&&c5.size()==10&&qb.size()==16);
 // The published head selector spans all degree>=3 terms (including quartic
 // and quintic), in increasing monomial order. Its surviving actual cosets
 // are tested below for weight at least four and even, not assumed even.
 c5.clear();
 for (int i=0;i<32;i++)if(std::popcount(unsigned(i))>2)c5.push_back(mono(i,32));
 assert(c5.size()==16);
 std::vector<W> Q(1<<16);
 for(unsigned i=1;i<Q.size();i++)Q[i]=W(unpack(i,qb));
 assert(std::set<W>(Q.begin(),Q.end()).size()==Q.size());
 std::vector<U> P;
 std::ifstream in(argv[1]);
 for(unsigned f;in>>f;)P.push_back(unpack(f,c6));
 assert(in.eof()&&P.size()==80&&std::set<U>(P.begin(),P.end()).size()==80);
 std::set<int> heads{};
 std::array<U,80>M{};
 std::ifstream spec(argv[2]);
 int g;spec>>g;assert(g==21);
 U unionmarks=0;
 for(int j=0;j<21;j++){
  int c,x,y,z;spec>>c>>x>>y>>z;
  assert(spec&&c>=0&&c<80&&!M[c]&&x>=0&&x<64&&y>=0&&y<64&&z>=0&&z<64&&x!=y&&y!=z&&z!=x);
  M[c]=(U(1)<<x)|(U(1)<<y)|(U(1)<<z);
  assert(!(unionmarks&M[c]));unionmarks|=M[c];
 }
 assert(std::popcount(unionmarks)==63);
 for(int h;spec>>h;) assert(h>=0&&h<126&&heads.insert(h).second);
 assert(spec.eof()&&heads.size()==44);
 std::ifstream hdata(argv[3]);
 int d=0,n=0;std::size_t tests=0,w4=0;std::set<U> used;
 for(int a,b,h,u;hdata>>a>>b>>h>>u;d++){
  if(!heads.contains(d))continue;
  assert(a>0&&a<64&&(b==0||b==1)&&h>=0&&h<65536&&u==15);
  std::vector<int>E;
  U s=0;
  for(int i=0;i<64;i++)if((std::popcount(unsigned(a&i))&1)==b){E.push_back(i);s|=U(1)<<i;}
  assert(E.size()==32&&used.insert(s).second);
  int pivot=std::countr_zero(unsigned(a));
  for(int i=0;i<32;i++){int x=E[i];assert(((x&((1<<pivot)-1))|((x>>(pivot+1))<<pivot))==i);}
  auto proj=[&](U v){W w=0;for(unsigned j=0;j<E.size();j++)w|=W((v>>E[j])&1)<<j;return w;};
  W B=W(unpack(h,c5));
  for(int c=0;c<80;c++) if(M[c]){
   W mm=proj(M[c]),v=proj(P[c])^B;
   for(W q:Q){
    W w=q^v;int k=std::popcount(w);
    assert(k>=4&&(k%2==0));
    if(k==4){assert(!(w&mm));w4++;}
    tests++;
   }
  }
  n++;
 }
 assert(hdata.eof()&&d==126&&n==44&&tests==44ULL*21*65536);
 std::cout<<"PASS "<<tests<<" triple/atlas differences, "<<w4<<" weight-four masks\n";
}
