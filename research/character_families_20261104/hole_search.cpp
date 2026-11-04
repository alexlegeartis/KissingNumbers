// Generic Cayley dominating-set proposal. A valid labelling is still required.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
int main(int argc,char**argv){
 if(argc!=6)return 2;
 std::ifstream in(argv[1]);int n,l;in>>n>>l;std::vector<int>d(l),seed(n);for(int&x:d)in>>x;for(int&x:seed)in>>x;d.push_back(0);
 int restarts=std::stoi(argv[3]),steps=std::stoi(argv[4]);std::mt19937_64 rng(std::stoull(argv[5]));
 std::vector<int>original;for(int i=0;i<n;i++)if(seed[i]<0)original.push_back(i);
 int best=n;std::vector<int>bestholes;
 for(int r=0;r<restarts;r++){
  auto holes=original;holes.erase(holes.begin()+rng()%holes.size());std::vector<int>cover(n),is_hole(n),mark(n);int tick=0;
  for(int h:holes){is_hole[h]=1;for(int x:d)cover[h^x]++;}
  int missing=std::count(cover.begin(),cover.end(),0);
  for(int t=0;t<steps;t++){
   if(missing<best){best=missing;bestholes=holes;std::cout<<"best "<<best<<" restart "<<r<<" step "<<t<<std::endl;}
   if(!missing)break;
   int u=rng()%n;if(rng()%10){while(cover[u])u=rng()%n;}
   int a=u^d[rng()%d.size()];if(is_hole[a])continue;
   int bi=rng()%holes.size(),b=holes[bi];++tick;int delta=0;double concentration=0;
   for(int x:d){int v=a^x;mark[v]=tick;if(!cover[v])--delta;concentration+=2*cover[v]+1;}
   for(int x:d){int v=b^x;if(cover[v]==1&&mark[v]!=tick)++delta;int c=cover[v]+(mark[v]==tick);concentration+=1-2*c;}
   double phase=double(t%100000)/100000.;double temperature=.45*std::pow(.03,phase);
   double energy=delta+.015*concentration;
   if(energy<=0||std::generate_canonical<double,53>(rng)<std::exp(-energy/temperature)){
    for(int x:d){cover[a^x]++;cover[b^x]--;}
    is_hole[a]=1;is_hole[b]=0;holes[bi]=a;missing+=delta;
   }
  }
  std::cout<<"restart "<<r<<" uncovered "<<missing<<" best "<<best<<std::endl;
  if(!best)break;
 }
 std::ofstream out(argv[2]);out<<best<<'\n';for(int h:bestholes)out<<h<<' ';out<<'\n';
}
