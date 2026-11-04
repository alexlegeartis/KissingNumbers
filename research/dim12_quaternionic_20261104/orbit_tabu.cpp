// Fixed-cardinality tabu search on the exact eight-point-orbit conflict graph.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
#include <omp.h>

template<class T> void read_file(const std::string& p,std::vector<T>& v){
  std::ifstream f(p,std::ios::binary);f.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(T));
  if(!f)throw std::runtime_error("Cannot read "+p);
}
int main(int argc,char** argv){
  if(argc<5){std::cerr<<"orbit_tabu POOL_DIR OUTPUT_DIR SECONDS THREADS\n";return 2;}
  std::string dir=argv[1],out=argv[2];double seconds=std::stod(argv[3]);int threads=std::stoi(argv[4]);
  std::filesystem::create_directories(out);
  std::vector<uint32_t> header(2);read_file(dir+"/header.bin",header);
  const int n=header[0],ns=header[1],K=97;
  if(n<=K || ns<1)throw std::runtime_error("Invalid pool header");
  std::vector<uint8_t> graph(size_t(n)*n);read_file(dir+"/weights.bin",graph);
  std::vector<uint16_t> seeds(size_t(ns)*96);read_file(dir+"/seeds.bin",seeds);
  std::vector<double> xyz(size_t(n)*12);read_file(dir+"/representatives.bin",xyz);
  auto begin=std::chrono::steady_clock::now();std::atomic<int> global_best(100000),saved(0);
  std::atomic<bool> solved(false);std::atomic<uint64_t> total_iterations(0);
  omp_set_num_threads(threads);
  #pragma omp parallel
  {
    int tid=omp_get_thread_num();std::mt19937_64 rng(842202611040ULL+tid);
    std::vector<int> chosen(K),count(n),tabu(n);std::vector<uint8_t> selected(n);
    uint64_t iteration=0,last_improvement=0;int cost=0,local_best=100000,restarts=0,local_saved=0,best_saved=100000;
    std::unordered_set<uint64_t> seen;
    auto edge=[&](int a,int b){return graph[size_t(a)*n+b];};
    auto reset=[&](){
      std::fill(count.begin(),count.end(),0);std::fill(selected.begin(),selected.end(),0);
      std::fill(tabu.begin(),tabu.end(),0);int s=rng()%ns;
      std::vector<int> order(96);for(int j=0;j<96;++j)order[j]=j;std::shuffle(order.begin(),order.end(),rng);
      int keep_options[6]={32,48,64,80,88,96};int keep=keep_options[rng()%6];
      for(int j=0;j<keep;++j){int v=seeds[size_t(s)*96+order[j]];chosen[j]=v;selected[v]=1;
        const uint8_t* row=graph.data()+size_t(v)*n;
        for(int a=0;a<n;++a)count[a]+=row[a];}
      cost=0;
      for(int j=keep;j<K;++j){
        int best=100000,v=-1,ties=0;
        for(int a=0;a<n;++a)if(!selected[a]){
          if(count[a]<best){best=count[a];v=a;ties=1;}
          else if(count[a]==best && rng()%++ties==0)v=a;
        }
        if(v<0)std::abort();chosen[j]=v;selected[v]=1;cost+=count[v];
        const uint8_t* row=graph.data()+size_t(v)*n;
        for(int a=0;a<n;++a)count[a]+=row[a];
      }
      ++restarts;local_best=cost;last_improvement=iteration;
    };
    reset();
    while(!solved.load(std::memory_order_relaxed)){
      if((iteration&255)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()>seconds)break;
      if((iteration&4095)==0){
        int exact=0;for(int a=0;a<K;++a)for(int b=0;b<a;++b)exact+=edge(chosen[a],chosen[b]);
        if(exact!=cost)std::abort();
        for(int test=0;test<8;++test){int v=rng()%n,c=0;for(int s:chosen)c+=edge(v,s);if(c!=count[v])std::abort();}
      }
      if(cost<global_best.load()){
        #pragma omp critical(progress)
        {if(cost<global_best.load()){global_best=cost;std::cout<<"best_conflicts "<<cost<<" thread "<<tid<<" iteration "<<iteration<<std::endl;}}
      }
      if(cost<=600 && (cost<best_saved || (local_saved<24 && iteration%2048==0))){
        auto sorted=chosen;std::sort(sorted.begin(),sorted.end());uint64_t hash=1469598103934665603ULL;
        for(int v:sorted){hash^=v;hash*=1099511628211ULL;}
        if(seen.insert(hash).second){
          int a=-1,b=-1,edges=0;double peak=0;
          for(int j=0;j<K;++j)for(int k=0;k<j;++k)if(edge(chosen[j],chosen[k])){
            int u=chosen[j],v=chosen[k];double value=0;++edges;
            for(int block=0;block<3;++block){double dot=0;
              for(int d=0;d<4;++d)dot+=xyz[size_t(u)*12+4*block+d]*xyz[size_t(v)*12+4*block+d];value+=std::abs(dot);}
            if(value>peak){peak=value;a=u;b=v;}
          }
          if(cost==0 || (peak<.87 && edges<24)){
            int number=saved.fetch_add(1);++local_saved;best_saved=std::min(best_saved,cost);
            std::ofstream f(out+"/candidate_"+std::to_string(number)+".txt");
            if(!f)std::abort();
            f<<cost<<" "<<a<<" "<<b<<" "<<peak<<" "<<edges<<"\n";
            for(int v:chosen)f<<v<<" ";f<<"\n";
            if(cost==0){solved=true;std::cout<<"FOUND_97_ORBITS candidate "<<number<<std::endl;break;}
          }
        }
      }
      ++iteration;
      // Remove a highly conflicting vertex; occasional random removal diversifies.
      int pos=-1,score=-1,ties=0;
      if(rng()%20==0)pos=rng()%K;
      else for(int j=0;j<K;++j){int c=count[chosen[j]];
        if(c>score){score=c;pos=j;ties=1;}else if(c==score && rng()%++ties==0)pos=j;}
      int removed=chosen[pos],old=count[removed],best=100000,replacement=-1;ties=0;
      const uint8_t* oldrow=graph.data()+size_t(removed)*n;
      for(int v=0;v<n;++v)if(!selected[v]){
        int next=count[v]-oldrow[v];
        if(uint64_t(tabu[v])>iteration && cost-old+next>=local_best)continue;
        if(next<best){best=next;replacement=v;ties=1;}
        else if(next==best && rng()%++ties==0)replacement=v;
      }
      if(replacement<0){reset();continue;}
      const uint8_t* newrow=graph.data()+size_t(replacement)*n;
      for(int v=0;v<n;++v)count[v]+=int(newrow[v])-int(oldrow[v]);
      cost+=best-old;selected[removed]=0;selected[replacement]=1;chosen[pos]=replacement;
      tabu[removed]=iteration+10+rng()%30;
      if(cost<local_best){local_best=cost;last_improvement=iteration;}
      if(iteration-last_improvement>12000)reset();
    }
    total_iterations.fetch_add(iteration);
    #pragma omp critical(progress)
    {std::cout<<"thread_done "<<tid<<" iterations "<<iteration<<" restarts "<<restarts<<std::endl;}
  }
  std::ofstream f(out+"/summary.json");
  f<<"{\n  \"threads\": "<<threads<<",\n  \"seconds_requested\": "<<seconds
   <<",\n  \"iterations\": "<<total_iterations<<",\n  \"best_weighted_conflict_cost\": "<<global_best
   <<",\n  \"saved_candidates\": "<<saved<<",\n  \"found_97_compatible_orbits\": "<<(solved?"true":"false")<<"\n}\n";
  return 0;
}
