// Compress a conflict graph into exactly equivalent at-most-one constraints.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc,char** argv){
 if(argc!=3)return 2;std::ifstream in(argv[1],std::ios::binary);uint32_t n,magic;in.read((char*)&n,4);in.read((char*)&magic,4);if(magic!=0x4b4e4731)return 3;
 std::vector<uint64_t> off(n+1);in.read((char*)off.data(),8*(n+1));std::vector<uint32_t> nb(off[n]);in.read((char*)nb.data(),4*nb.size());if(!in)return 4;
 const uint32_t width=(n+63)/64;std::vector<uint64_t> adj(size_t(n)*width),uncovered;
 for(uint32_t i=0;i<n;++i)for(uint64_t e=off[i];e<off[i+1];++e)adj[size_t(i)*width+nb[e]/64]|=1ull<<(nb[e]%64);
 uncovered=adj;uint64_t remaining=nb.size()/2,cliques=0,literals=0;uint32_t maxsize=0;std::mt19937_64 rng(17020);
 std::string temporary=std::string(argv[2])+".tmp";std::ofstream out(temporary);
 std::vector<uint64_t> cand(width);std::vector<uint32_t> candidates,clique;
 for(uint32_t first=0;first<n;++first){
  for(;;){
   uint32_t second=n;
   for(uint32_t k=0;k<width;++k)if(uncovered[size_t(first)*width+k]){second=64*k+__builtin_ctzll(uncovered[size_t(first)*width+k]);break;}
   if(second==n)break;
   clique={first,second};
   for(uint32_t k=0;k<width;++k)cand[k]=adj[size_t(first)*width+k]&adj[size_t(second)*width+k];
   for(;;){
    candidates.clear();for(uint32_t k=0;k<width;++k){uint64_t bits=cand[k];while(bits){int b=__builtin_ctzll(bits);candidates.push_back(64*k+b);bits&=bits-1;}}
    if(candidates.empty())break;
    uint32_t best=candidates[rng()%candidates.size()];int bestscore=-1;
    for(uint32_t trial=0;trial<std::min(size_t(24),candidates.size());++trial){
     uint32_t v=candidates[rng()%candidates.size()];int gain=0,degree=0;
     for(auto c:clique)gain+=(uncovered[size_t(v)*width+c/64]>>(c%64))&1;
     for(uint32_t k=0;k<width;++k)degree+=__builtin_popcountll(cand[k]&adj[size_t(v)*width+k]);
     int score=10000*gain+degree;if(score>bestscore){bestscore=score;best=v;}
    }
    clique.push_back(best);for(uint32_t k=0;k<width;++k)cand[k]&=adj[size_t(best)*width+k];
   }
   for(size_t a=0;a<clique.size();++a)for(size_t b=a+1;b<clique.size();++b){
    uint32_t x=clique[a],y=clique[b];if(!(adj[size_t(x)*width+y/64]>>(y%64)&1))return 5;
    if(uncovered[size_t(x)*width+y/64]>>(y%64)&1){--remaining;uncovered[size_t(x)*width+y/64]&=~(1ull<<(y%64));uncovered[size_t(y)*width+x/64]&=~(1ull<<(x%64));}
   }
   out<<clique.size();for(auto v:clique)out<<' '<<v;out<<'\n';++cliques;literals+=clique.size();maxsize=std::max(maxsize,uint32_t(clique.size()));
  }
 }
 if(remaining)return 6;
 out.close();if(std::rename(temporary.c_str(),argv[2]))return 7;
 std::cout<<"{\"vertices\":"<<n<<",\"edges\":"<<nb.size()/2<<",\"cliques\":"<<cliques<<",\"literals\":"<<literals<<",\"max_clique_size\":"<<maxsize<<",\"uncovered_edges\":0}"<<std::endl;
}
