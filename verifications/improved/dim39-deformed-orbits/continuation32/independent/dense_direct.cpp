#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <chrono>
#include <omp.h>
#ifdef __AVX2__
#include <immintrin.h>
#endif
int64_t dot40(const int32_t*a,const int32_t*b){
#ifdef __AVX2__
 __m256i acc=_mm256_setzero_si256();
 for(int k=0;k<40;k+=8){auto x=_mm256_loadu_si256((const __m256i*)(a+k));auto y=_mm256_loadu_si256((const __m256i*)(b+k));acc=_mm256_add_epi64(acc,_mm256_mul_epi32(x,y));acc=_mm256_add_epi64(acc,_mm256_mul_epi32(_mm256_srli_epi64(x,32),_mm256_srli_epi64(y,32)));}
 alignas(32) int64_t r[4];_mm256_store_si256((__m256i*)r,acc);return r[0]+r[1]+r[2]+r[3];
#else
 int64_t s=0;for(int k=0;k<40;k++)s+=(int64_t)a[k]*b[k];return s;
#endif
}
int main(int argc,char**argv){std::ifstream in(argv[1],std::ios::binary);int64_t N,d,order,lo,hi;in.read((char*)&N,8);in.read((char*)&d,8);in.read((char*)&order,8);in.read((char*)&lo,8);in.read((char*)&hi,8);if(N!=10256||d!=39||order!=32||lo!=99999996||hi!=100000004)return 2;std::vector<int32_t>V(N*40),W(N*40);std::vector<uint64_t>group(order);in.read((char*)V.data(),V.size()*4);in.read((char*)group.data(),order*8);if(!in)return 3;int64_t mn=INT64_MAX,mx=0;
 for(int i=0;i<N;i++){int64_t n=0;for(int k=0;k<40;k++){int64_t x=V[i*40+k];if(std::abs(x)>100000000|| (k==39&&x))return 4;n+=x*x;}if(n<lo*lo||n>hi*hi)return 5;mn=std::min(mn,n);mx=std::max(mx,n);}
 int64_t maximum=INT64_MIN,wi=-1,wj=-1,wg=-1;uint64_t count=0;int64_t allowed=lo*lo/2;auto start=std::chrono::steady_clock::now();
 for(int g=0;g<order;g++){
  for(int i=0;i<N;i++)for(int k=0;k<40;k++)W[i*40+k]=(k<39&&((group[g]>>k)&1))?-V[i*40+k]:V[i*40+k];
  for(int t=0;t<256;t++){int i=(g*103+t*31)%N,j=(g*307+t*67)%N;int64_t scalar=0;for(int k=0;k<40;k++)scalar+=(int64_t)V[i*40+k]*W[j*40+k];if(scalar!=dot40(&V[i*40],&W[j*40]))return 7;}
  int64_t gmax=INT64_MIN,gi=-1,gj=-1;uint64_t gc=0;
  #pragma omp parallel
  {int64_t local=INT64_MIN,li=-1,lj=-1;uint64_t lc=0;
   #pragma omp for schedule(dynamic,16)
   for(int i=0;i<N;i++)for(int j=i;j<N;j++){if(i==j&&group[g]==0)continue;int64_t s=dot40(&V[i*40],&W[j*40]);lc++;if(s>local){local=s;li=i;lj=j;}}
   #pragma omp critical
   {gc+=lc;if(local>gmax){gmax=local;gi=li;gj=lj;}}
  }
  count+=gc;if(gmax>maximum){maximum=gmax;wi=gi;wj=gj;wg=g;}
  std::cerr<<"group "<<g+1<<" / "<<order<<" max "<<gmax<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
  if(gmax>allowed){std::cerr<<"VIOLATION "<<gi<<" "<<gj<<" "<<group[g]<<"\n";return 6;}
 }
 std::cout<<"{\"passed\":true,\"representatives\":"<<N<<",\"group_order\":"<<order<<",\"products_checked\":"<<count<<",\"minimum_squared_norm\":"<<mn<<",\"maximum_squared_norm\":"<<mx<<",\"maximum_integer_dot\":"<<maximum<<",\"allowed_integer_dot\":"<<allowed<<",\"witness\":["<<wi<<","<<wj<<","<<group[wg]<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"threads\":"<<omp_get_max_threads()<<"}\n";
}
