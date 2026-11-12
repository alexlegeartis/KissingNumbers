#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
int main(int argc,char**argv) {
 if(argc!=2)return 2;
 std::ifstream f(argv[1],std::ios::binary);std::uint64_t n;
 f.read(reinterpret_cast<char*>(&n),8);if(!f||n<1||n>100000)return 2;
 std::vector<std::int64_t> G(48*48),X(n*48),XG(n*48,0);
 f.read(reinterpret_cast<char*>(G.data()),G.size()*8);
 f.read(reinterpret_cast<char*>(X.data()),X.size()*8);if(!f||f.peek()!=EOF)return 2;
 for(auto x:G)if(x>1000||x<-1000)return 3;
 for(auto x:X)if(x>1000||x<-1000)return 3;
 // Every XG partial sum is <=48*1000^2; every final dot <=48^2*1000^3<2^63.
 for(std::uint64_t i=0;i<n;i++)for(int j=0;j<48;j++)for(int k=0;k<48;k++)XG[i*48+j]+=X[i*48+k]*G[k*48+j];
 std::uint64_t pairs=0;std::int64_t maxdot=-1000000;
 for(std::uint64_t i=0;i<n;i++){
  std::int64_t norm=0;for(int k=0;k<48;k++)norm+=XG[i*48+k]*X[i*48+k];
  if(norm!=6){std::cerr<<"bad norm "<<i<<" "<<norm<<"\n";return 4;}
  for(std::uint64_t j=0;j<i;j++){
   std::int64_t dot=0;for(int k=0;k<48;k++)dot+=XG[i*48+k]*X[j*48+k];
   if(dot>2){std::cerr<<"bad pair "<<i<<" "<<j<<" "<<dot<<"\n";return 5;}
   maxdot=std::max(maxdot,dot);pairs++;
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"signed_points\":"<<n<<",\"unordered_pairs\":"<<pairs<<",\"maximum_product\":"<<maxdot<<",\"arithmetic\":\"bounded int64 scalar accumulation\"}\n";
}
