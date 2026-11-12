#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace std;
mt19937 rng(1232026111);
double now(){return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();}
struct Profile{array<int,31>a;int weight,tailtotal;};
vector<Profile>P;array<array<int,3>,21>B;array<array<int,2>,21>E;
array<int,24>pc,bpc;array<int,7>pt,bpt;
vector<vector<unsigned char>>blocked(21);vector<int> countblocked;
int cost=0,best=100000000;long long moves=0;
bool bad(int r,int z){auto&p=P[z];return p.a[pc[B[r][0]]]+p.a[pc[B[r][1]]]+p.a[pc[B[r][2]]]+p.tailtotal-p.a[24+pt[E[r][0]]]-p.a[24+pt[E[r][1]]]>16;}
void remove(const vector<int>&rs){for(int r:rs)for(int z=0;z<(int)P.size();z++)if(blocked[r][z]&&!--countblocked[z])cost-=P[z].weight;}
void insert(const vector<int>&rs){for(int r:rs)for(int z=0;z<(int)P.size();z++){bool b=bad(r,z);blocked[r][z]=b;if(b&&countblocked[z]++==0)cost+=P[z].weight;}}
void init(){countblocked.assign(P.size(),0);cost=0;vector<int>all(21);iota(all.begin(),all.end(),0);insert(all);}
void save(string dir){ofstream o(dir+"/layer_3plus5_best.txt");o<<best<<"\n";for(int x:bpc)o<<x<<" ";o<<"\n";for(int x:bpt)o<<x<<" ";o<<"\n";}
int main(int argc,char**argv){string dir=argv[1];double budget=argc>2?stod(argv[2]):60;
 ifstream f(dir+"/absolute_profiles31.txt");int m;f>>m;P.resize(m);for(auto&p:P){f>>p.weight;p.tailtotal=0;for(int j=0;j<31;j++){f>>p.a[j];if(j>=24)p.tailtotal+=p.a[j];}}
 ifstream g(dir+"/layer_3plus5_design.txt");for(int r=0;r<21;r++){g>>E[r][0]>>E[r][1];for(int j=0;j<3;j++)g>>B[r][j];}
 for(auto&v:blocked)v.resize(P.size());iota(pc.begin(),pc.end(),0);iota(pt.begin(),pt.end(),0);init();cout<<"initial "<<cost<<endl;
 double start=now(),epoch=0;int restarts=0;
 while(now()-start<budget){
  if(moves%25000==0){if(restarts%3==0){shuffle(pc.begin(),pc.end(),rng);shuffle(pt.begin(),pt.end(),rng);}else{pc=bpc;pt=bpt;}init();epoch=now();restarts++;}
  bool tail=(rng()%4==0);int n=tail?7:24,a=rng()%n,b=rng()%n;if(a==b)continue;
  vector<int>rs;for(int r=0;r<21;r++){bool affected=false;if(tail){affected=(E[r][0]==a||E[r][0]==b||E[r][1]==a||E[r][1]==b);}else for(int j=0;j<3;j++)if(B[r][j]==a||B[r][j]==b)affected=true;if(affected)rs.push_back(r);}
  int prev=cost;remove(rs);if(tail)swap(pt[a],pt[b]);else swap(pc[a],pc[b]);insert(rs);
  double temp=500*pow(.001,double(moves%25000)/25000);
  if(cost>prev&&uniform_real_distribution<double>(0,1)(rng)>exp(-(cost-prev)/temp)){remove(rs);if(tail)swap(pt[a],pt[b]);else swap(pc[a],pc[b]);insert(rs);}
  moves++;
  if(cost<best){best=cost;bpc=pc;bpt=pt;save(dir);cout<<"best "<<best<<" net "<<2688-best<<" moves "<<moves<<" sec "<<now()-start<<endl;}
 }
 cout<<"final "<<best<<" net "<<2688-best<<" moves "<<moves<<" restarts "<<restarts<<" sec "<<now()-start<<endl;
}
