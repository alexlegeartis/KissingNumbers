// Fixed-cardinality individual-ray exchanges on an exactly constructed graph.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <omp.h>
template<class T> void read(const std::string& path,std::vector<T>& v){
  std::ifstream f(path,std::ios::binary);f.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(T));
  if(!f)throw std::runtime_error(path);
}
int main(int argc,char**argv){
  if(argc<5){std::cerr<<"trade_tabu POOL OUT SECONDS THREADS [SEED]\n";return 2;}
  std::string dir=argv[1],out=argv[2];double seconds=std::stod(argv[3]);int threads=std::stoi(argv[4]);
  uint64_t seed=argc>5?std::stoull(argv[5]):202611050842ULL;
  std::filesystem::create_directories(out);std::vector<uint32_t> h(4);read(dir+"/header.bin",h);
  int n=h[0],ns=h[1],base=h[2],K=h[3];
  if(K<=base || K>=n || threads<1 || threads>127)throw std::runtime_error("header or threads");
  std::vector<uint8_t> graph(size_t(n)*n);read(dir+"/weights.bin",graph);
  std::vector<uint16_t> seeds(size_t(ns)*base);read(dir+"/seeds.bin",seeds);
  auto start=std::chrono::steady_clock::now();std::atomic<bool> solved(false);
  std::atomic<int> global(100000),saved(0);std::atomic<uint64_t> moves(0);
  omp_set_num_threads(threads);
  #pragma omp parallel
  {
    int tid=omp_get_thread_num();std::mt19937_64 rng(seed+tid);
    std::vector<int> chosen(K),counts(n);std::vector<uint64_t> tabu(n);std::vector<uint8_t> in(n);
    uint64_t it=0,last=0;int cost=0,best=100000,ever=100000,restarts=0;
    auto reset=[&](){
      std::fill(counts.begin(),counts.end(),0);std::fill(in.begin(),in.end(),0);std::fill(tabu.begin(),tabu.end(),0);
      int s=rng()%ns;std::vector<int> order(base);for(int j=0;j<base;++j)order[j]=j;
      std::shuffle(order.begin(),order.end(),rng);int keep_options[7]={base,base-4,base-16,base-32,base-64,base-128,base/2};
      int keep=keep_options[rng()%7];cost=0;
      for(int j=0;j<K;++j){int v=-1,b=100000,ties=0;
        if(j<keep)v=seeds[size_t(s)*base+order[j]];
        else for(int a=0;a<n;++a)if(!in[a]){
          if(counts[a]<b){b=counts[a];v=a;ties=1;}
          else if(counts[a]==b && rng()%++ties==0)v=a;
        }
        if(v<0)std::abort();chosen[j]=v;in[v]=1;cost+=counts[v];
        for(int a=0;a<n;++a)counts[a]+=graph[size_t(v)*n+a];
      }
      best=cost;last=it;++restarts;
    };
    reset();
    while(!solved.load(std::memory_order_relaxed)){
      if((it&511)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>seconds)break;
      if((it&8191)==0){int exact=0;for(int a=0;a<K;++a)for(int b=0;b<a;++b)exact+=graph[size_t(chosen[a])*n+chosen[b]];
        if(exact!=cost)std::abort();
        for(int test=0;test<8;++test){int v=rng()%n,sum=0;for(int u:chosen)sum+=graph[size_t(u)*n+v];if(sum!=counts[v])std::abort();}
      }
      if(cost<ever){ever=cost;
        std::ofstream f(out+"/thread_"+std::to_string(tid)+"_cost_"+std::to_string(cost)+".txt");
        f<<cost<<"\n";for(int v:chosen)f<<v<<" ";f<<"\n";++saved;
      }
      if(cost<global.load()){
        #pragma omp critical(report)
        {if(cost<global.load()){global=cost;std::cout<<"best_conflicts "<<cost<<" thread "<<tid<<" iteration "<<it<<std::endl;}}
      }
      if(!cost){solved=true;break;}
      ++it;int pos=-1,worst=-1,ties=0;
      if(rng()%12==0)pos=rng()%K;
      else for(int j=0;j<K;++j){int c=counts[chosen[j]];
        if(c>worst){worst=c;pos=j;ties=1;}else if(c==worst && rng()%++ties==0)pos=j;}
      int old=chosen[pos],oldcost=counts[old],replacement=-1,next=100000;ties=0;
      auto oldrow=graph.data()+size_t(old)*n;
      for(int v=0;v<n;++v)if(!in[v]){int c=counts[v]-oldrow[v];
        if(tabu[v]>it && cost-oldcost+c>=best)continue;
        if(c<next){next=c;replacement=v;ties=1;}else if(c==next && rng()%++ties==0)replacement=v;}
      if(replacement<0){reset();continue;}auto newrow=graph.data()+size_t(replacement)*n;
      for(int v=0;v<n;++v)counts[v]+=int(newrow[v])-int(oldrow[v]);
      chosen[pos]=replacement;in[old]=0;in[replacement]=1;cost+=next-oldcost;
      tabu[old]=it+15+rng()%65;
      if(cost<best){best=cost;last=it;}
      if(it-last>12000)reset();
    }
    moves+=it;
  }
  std::ofstream f(out+"/summary.json");
  f<<"{\"vertices\":"<<n<<",\"target\":"<<K<<",\"threads\":"<<threads
   <<",\"seconds_requested\":"<<seconds<<",\"moves\":"<<moves<<",\"best_conflicts\":"<<global
   <<",\"saved\":"<<saved<<",\"solved\":"<<(solved?"true":"false")<<",\"seed\":"<<seed<<"}\n";
  std::cout<<"finished best "<<global<<" moves "<<moves<<std::endl;
}
