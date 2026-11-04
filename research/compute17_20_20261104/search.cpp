// Independent-set ruin/recreate search. Graph is shared read-only across workers.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <random>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

struct Graph {
 uint32_t n; const uint64_t* off; const uint32_t* nb;
 explicit Graph(const char* path) {
  int fd=open(path,O_RDONLY); if(fd<0)std::exit(2); struct stat st;fstat(fd,&st);
  auto p=(const uint8_t*)mmap(nullptr,st.st_size,PROT_READ,MAP_SHARED,fd,0);close(fd);
  n=*(const uint32_t*)p; if(*(const uint32_t*)(p+4)!=0x4b4e4731)std::exit(3);
  off=(const uint64_t*)(p+8);nb=(const uint32_t*)(p+8+8*(n+1));
 }
};

struct State {
 const Graph& g;const std::vector<int>& w;std::mt19937_64 rng;
 std::vector<int> cnt,pos,apos;std::vector<uint32_t> chosen,avail;
 std::vector<double> penalty;int score=0;
 State(const Graph& gg,const std::vector<int>& ww,uint64_t s):g(gg),w(ww),rng(s),cnt(g.n),pos(g.n,-1),apos(g.n,-1),penalty(g.n){}
 void az(uint32_t v){ if(apos[v]<0){apos[v]=avail.size();avail.push_back(v);} }
 void dz(uint32_t v){int p=apos[v];if(p>=0){apos[avail.back()]=p;avail[p]=avail.back();avail.pop_back();apos[v]=-1;}}
 void add(uint32_t v){
  if(cnt[v]||pos[v]>=0)std::abort(); dz(v);pos[v]=chosen.size();chosen.push_back(v);score+=w[v];
  for(uint64_t e=g.off[v];e<g.off[v+1];++e){uint32_t b=g.nb[e];if(cnt[b]++==0)dz(b);}
 }
 void remove(uint32_t v){
  int p=pos[v];if(p<0)std::abort();pos[chosen.back()]=p;chosen[p]=chosen.back();chosen.pop_back();pos[v]=-1;score-=w[v];az(v);
  for(uint64_t e=g.off[v];e<g.off[v+1];++e){uint32_t b=g.nb[e];if(--cnt[b]==0&&pos[b]<0)az(b);}
 }
 void reset(const std::vector<uint32_t>& s){
  std::fill(cnt.begin(),cnt.end(),0);std::fill(pos.begin(),pos.end(),-1);std::fill(apos.begin(),apos.end(),-1);chosen.clear();avail.clear();score=0;
  for(auto v:s){if(cnt[v])std::abort();pos[v]=chosen.size();chosen.push_back(v);score+=w[v];for(uint64_t e=g.off[v];e<g.off[v+1];++e)++cnt[g.nb[e]];}
  for(uint32_t v=0;v<g.n;++v)if(!cnt[v]&&pos[v]<0)az(v);
 }
 void fill(int samples,double alpha,double pen){
  while(!avail.empty()){
   uint32_t best=avail[rng()%avail.size()];double bs=-1;
   for(int k=0;k<samples;++k){auto v=avail[rng()%avail.size()];double s=w[v]/std::pow(1.+double(g.off[v+1]-g.off[v]),alpha)/(1.+pen*penalty[v]);
    s*=.9+.2*double(rng()%10000)/10000.;if(s>bs){bs=s;best=v;}}
   add(best);
  }
 }
 int cost(uint32_t v)const{int c=0;for(uint64_t e=g.off[v];e<g.off[v+1];++e){auto b=g.nb[e];if(pos[b]>=0)c+=w[b];}return c;}
 void force(uint32_t v){
  std::vector<uint32_t> gone;for(uint64_t e=g.off[v];e<g.off[v+1];++e){auto b=g.nb[e];if(pos[b]>=0)gone.push_back(b);}
  for(auto b:gone)remove(b);if(pos[v]<0)add(v);
 }
 bool valid()const{for(auto v:chosen)for(uint64_t e=g.off[v];e<g.off[v+1];++e)if(pos[g.nb[e]]>=0)return false;return true;}
};

int main(int argc,char** argv){
 if(argc!=8){std::cerr<<"usage: search graph meta output seconds threads random_seed stop_at\n";return 2;}
 Graph g(argv[1]);std::ifstream meta(argv[2]);int n,fixed,dim;meta>>n>>fixed>>dim;if(n!=int(g.n))return 4;
 std::vector<int>w(n);for(auto& x:w)meta>>x;int k;meta>>k;std::vector<uint32_t> seed(k);for(auto& x:seed)meta>>x;if(!meta)return 5;
 int base=0;for(auto v:seed)base+=w[v];std::atomic<int> best(base);std::atomic<bool> stop(false);std::atomic<uint64_t> iterations(0);
 double start=omp_get_wtime(),seconds=std::atof(argv[4]);int threads=std::atoi(argv[5]);uint64_t random_seed=std::stoull(argv[6]);int stop_at=std::atoi(argv[7]);
 omp_set_num_threads(threads);std::string out(argv[3]);
 #pragma omp parallel
 {
  int worker=omp_get_thread_num();State s(g,w,random_seed+104729*worker);s.reset(seed);if(!s.valid())std::abort();
  int localbest=base;std::vector<uint32_t> local=seed;uint64_t rounds=0;double last=start;
  while(!stop&&omp_get_wtime()-start<seconds){
   auto old=s.chosen;int oldscore=s.score;
   int mode=s.rng()%5;
   if(mode<=1){
    // Large neighbourhoods are essential: old searches only freed <=200 vertices.
    double fraction=std::vector<double>{.12,.25,.40,.60,.80}[s.rng()%5];
    int remove=std::max(1,int(s.chosen.size()*fraction));
    while(remove--&&!s.chosen.empty())s.remove(s.chosen[s.rng()%s.chosen.size()]);
   }else{
    int forces=1+s.rng()%((worker%3==0)?16:4);
    while(forces--){
     uint32_t b=0;double bc=1e30;bool found=false;
     for(int z=0;z<64;++z){uint32_t v=s.rng()%g.n;if(s.pos[v]>=0)continue;
      double c=double(s.cost(v))/w[v]+.15*s.penalty[v]+double(s.rng()%10000)/20000.;
      if(c<bc){bc=c;b=v;found=true;}}
     if(found)s.force(b);
    }
   }
   s.fill(2+(worker%4)*6,(worker%4)*.3,.015);
   ++rounds;++iterations;
   if(s.score>localbest){localbest=s.score;local=s.chosen;}
   if(s.score>best.load()){
    #pragma omp critical(checkpoint)
    {
     if(s.score>best.load()){
      if(!s.valid())std::abort();best=s.score;
      std::ofstream f(out+"_best.txt");for(auto v:s.chosen)f<<v<<'\n';
      std::cout<<"{\"event\":\"improvement\",\"dimension\":"<<dim<<",\"count\":"<<fixed+s.score<<",\"worker\":"<<worker<<",\"seconds\":"<<omp_get_wtime()-start<<"}"<<std::endl;
      if(fixed+s.score>=stop_at)stop=true;
     }
    }
   }
   for(auto v:s.chosen)s.penalty[v]+=1.;
   // Occasional downhill steps avoid repeatedly restoring the same maximal set.
   double temperature=(worker%4==0?4.0:(worker%4==1?1.0:.1))*(*std::min_element(w.begin(),w.end()));
   bool accept=s.score>=oldscore||double(s.rng()%1000000)/1000000.<std::exp(double(s.score-oldscore)/temperature);
   if(!accept)s.reset(old);
   if(rounds%200==0&&s.score<localbest)s.reset(local);
   double now=omp_get_wtime();if(worker==0&&now-last>30){
    #pragma omp critical(checkpoint)
    std::cout<<"{\"event\":\"progress\",\"dimension\":"<<dim<<",\"best_count\":"<<fixed+best.load()<<",\"iterations\":"<<iterations.load()<<",\"seconds\":"<<now-start<<"}"<<std::endl;
    last=now;
   }
  }
 }
 std::ofstream report(out+"_result.json");
 report<<"{\"dimension\":"<<dim<<",\"baseline\":"<<fixed+base<<",\"best_count\":"<<fixed+best.load()<<",\"iterations\":"<<iterations.load()<<",\"seconds\":"<<omp_get_wtime()-start<<",\"workers\":"<<threads<<",\"seed\":"<<random_seed<<"}\n";
 std::cout<<"{\"event\":\"complete\",\"best_count\":"<<fixed+best.load()<<",\"iterations\":"<<iterations.load()<<",\"seconds\":"<<omp_get_wtime()-start<<"}"<<std::endl;
}
