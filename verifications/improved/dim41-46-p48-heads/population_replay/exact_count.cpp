// Exact shifted enumeration in a rigorously enclosing integral LDL ellipsoid.
// A coordinate row w is constrained modulo the supplied discriminant exponent.
// Physical squared norm is w G w^t / modulus^2.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
using I=__int128_t; using U=__uint128_t; constexpr int N=40;
long long G[N][N],M[N][N],D[N],w[N],c[N],P,target,modulus;
uint64_t nodes=0,leaves=0,count=0;
I parse128(const std::string&s){I v=0;for(char x:s){assert(x>='0'&&x<='9');assert(v<(I(1)<<123));v=v*10+x-'0';}return v;}
long long floor_div(long long a,long long b){long long q=a/b,r=a%b;return q-(r<0);}
long long ceil_div(long long a,long long b){return -floor_div(-a,b);}
uint64_t isqrt(U n){uint64_t k=(uint64_t)sqrtl((long double)n);while((U)k*k>n)--k;while((U)(k+1)*(k+1)<=n)++k;return k;}
void visit(int k,I rem){
 ++nodes;
 if(k<0){++leaves;I q=0;for(int i=0;i<N;i++)for(int j=0;j<N;j++)q+=(I)w[i]*G[i][j]*w[j];if(q==target)++count;return;}
 I offset=0;for(int j=k+1;j<N;j++)offset+=(I)M[j][k]*w[j];assert(offset>-(I(1)<<60)&&offset<(I(1)<<60));long long off=(long long)offset;
 assert(rem>=0);uint64_t b=isqrt((U)(rem/D[k]));assert(b<(1ULL<<60));
 long long lo=ceil_div(-(long long)b-off,P),hi=floor_div((long long)b-off,P);
 lo+=((c[k]-lo)%modulus+modulus)%modulus;
 assert(lo>-(1LL<<20)&&hi<(1LL<<20));
 for(long long v=lo;v<=hi;v+=modulus){w[k]=v;I t=(I)P*v+off;assert(t>-(I(1)<<44)&&t<(I(1)<<44));I used=(I)D[k]*t*t;assert(used<=rem);visit(k-1,rem-used);}
}
int main(int argc,char**argv){
 if(argc!=2&&argc!=3)return 2;int skip=argc==3?std::stoi(argv[2]):0;std::ifstream in(argv[1]);int n,number;in>>n>>P>>modulus>>number;assert(n==N&&P>0&&P<(1LL<<40)&&modulus>0&&modulus<(1LL<<20));
 for(auto&x:D){in>>x;assert(x>0&&x<(1LL<<36));}
 for(auto&r:M)for(auto&x:r){in>>x;assert(x>-(1LL<<40)&&x<(1LL<<40));}
 for(auto&r:G)for(auto&x:r){in>>x;assert(x>-(1LL<<20)&&x<(1LL<<20));}
 for(int job=0;job<number;job++){
  int cid;std::string rs;in>>cid>>target>>rs;I radius=parse128(rs);assert(radius>0&&radius<(I(1)<<120));assert(target>0&&target%modulus==0);
  for(auto&x:c)in>>x;assert(in.good());if(job<skip)continue;nodes=leaves=count=0;auto start=std::chrono::steady_clock::now();visit(N-1,radius);
  double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::cout<<"{\"class_id\":"<<cid<<",\"body_norm_numerator\":"<<target/modulus<<",\"body_norm_denominator\":"<<modulus<<",\"count\":"<<count<<",\"nodes\":"<<nodes<<",\"leaves\":"<<leaves<<",\"seconds\":"<<sec<<"}\n"<<std::flush;
 }
}
