// Research prototype: simultaneous hole/line assignments on a Cayley component.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
int main(int argc,char**argv){
 if(argc!=6)return 2;
 std::ifstream in(argv[1]);int n,l;in>>n>>l;
 std::vector<int>d(l),base(n),best;for(auto&x:d)in>>x;for(auto&x:base)in>>x;
 int restarts=std::stoi(argv[3]),steps=std::stoi(argv[4]);std::mt19937_64 rng(std::stoull(argv[5]));
 int bestcount=std::count_if(base.begin(),base.end(),[](int x){return x>=0;});best=base;
 auto bad=[](int a,int b,int p){return a>=0&&b>=0&&(a==p||b==p||a==b);};
 for(int r=0;r<restarts;r++){
  auto state=(r%3==0?base:best);int count=bestcount,edges=0;
  count=std::count_if(state.begin(),state.end(),[](int x){return x>=0;});
  for(int t=0;t<steps;t++){
   int v=rng()%n,old=state[v],a=int(rng()%(l+1))-1;if(a==old)continue;
   int delta=0;for(int p=0;p<l;p++){int k=state[v^d[p]];delta+=bad(a,k,p)-bad(old,k,p);}
   int dn=(a>=0)-(old>=0);double phase=double(t)/steps;
   double temp=(0.8+0.3*(r%4))*std::pow(0.015,phase);
   double penalty=0.55+1.45*phase;
   double energy=penalty*delta-dn;
   if(energy<=0||std::generate_canonical<double,53>(rng)<std::exp(-energy/temp)){
    state[v]=a;edges+=delta;count+=dn;
    if(edges==0&&count>bestcount){bestcount=count;best=state;std::cout<<"gain "<<r<<" "<<t<<" "<<bestcount<<std::endl;}
   }
  }
  std::cout<<"restart "<<r<<" count "<<count<<" conflicts "<<edges<<" best "<<bestcount<<std::endl;
 }
 std::ofstream out(argv[2]);for(int x:best)out<<x<<' ';out<<'\n';
}
