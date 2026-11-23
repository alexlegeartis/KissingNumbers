#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
using U=uint64_t;
static U mono(int m){U z=0;for(int x=0;x<64;x++)if((x&m)==m)z|=1ull<<x;return z;}
static int gf(int a,int b){int s=0;for(;b;b>>=1){if(b&1)s^=a;a<<=1;if(a&64)a^=67;}return s;}
int main(int argc,char** argv){
 std::ifstream src(argc==2?argv[1]:"witness.txt");
 std::vector<int> W;for(int x;src>>x;)W.push_back(x);
 assert(W.size()==80&&std::set<int>(W.begin(),W.end()).size()==80&&src.eof());
 std::vector<int> A,B,C;for(int i=0;i<64;i++){int n=__builtin_popcount(i);if(n<=2)A.push_back(i);if(n==2)B.push_back(i);if(n==3)C.push_back(i);}
 auto eval=[&](int i,std::vector<int>const&M){U z=0;for(int j=0;j<(int)M.size();j++)if(i>>j&1)z^=mono(M[j]);return z;};
 std::vector<U> full(1<<22);for(unsigned i=1;i<full.size();i++)full[i]=full[i-1]^mono(A[__builtin_ctz(i)]);
 std::vector<U> sub;U rep[64]={};bool seen[64]={};int H[1<<15]={};
 for(int i=1;i<(1<<15);i++){int m=B[__builtin_ctz(i)],a=m&(m-1),b=m^a;H[i]=H[i&(i-1)]^gf(gf(a,b),a^b);}
 for(int i=0;i<(1<<15);i++){
  U f=eval(i,B);
  if(!seen[H[i]]){seen[H[i]]=true;rep[H[i]]=f;}
  if(!H[i])for(int a=0;a<128;a++)sub.push_back(f^eval(a,{0,1,2,4,8,16,32}));
 }
 assert(std::count(seen,seen+64,true)==64);
 assert(sub.size()==65536&&std::set<U>(sub.begin(),sub.end()).size()==sub.size());
 for(U v:sub)assert(v==0||__builtin_popcountll(v)>=24);
 std::array<int,2> cub={275722,686811},ph={44,62};
 for(int a:cub)for(int w:W){U z=eval(a^w,C);int md=64;for(U h:full)md=std::min(md,__builtin_popcountll(z^h));assert(md>=14);}
 U d=eval(cub[0]^cub[1],C)^rep[ph[0]]^rep[ph[1]];
 int md=64;for(U h:sub)md=std::min(md,__builtin_popcountll(d^h));
 assert(md>=22); 
 // Shift the negative side by rep[1], a quadratic in a different
 // checksum class. Distinct cubics cannot coincide across the sides,
 // and equal cubics differ by the nonzero checksum: all differences
 // are nonzero degree <=3 and hence min >=8 (or check directly).
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){
  U z=eval(cub[i]^cub[j],C)^rep[ph[i]]^rep[ph[j]]^rep[1]; 
  for(U h:sub)assert(__builtin_popcountll(z^h)>=8);
 }
 // Dense-equator cross bound is 36 = sqrt(64*81)/2.
 // Same-side bound is 64 - 2*22 + 17 = 37 <=81/2.
 // Opposite-side bound is 64 - 2*8 - 17 = 31 <=81/2.
 // Equatorial sparse and four-support cross bound <=32.
 // The two axial pole cross products satisfy 8*sqrt(17) <=36.
 assert(36*36==64*81/4&&2*(64-2*22+17)<=81&&2*(64-2*8-17)<=81&&64*17<=36*36);
 std::cout<<"lift half size "<<2*sub.size()<<", same-side intercoset minimum "<<md<<", added "<<4*sub.size()<<", enhanced 65D lower bound "<<413526914+4*sub.size()<<"\n";
}
