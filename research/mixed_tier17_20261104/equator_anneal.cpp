#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
int main(){string R="research/mixed_tier17_20261104/";vector<int>w;int x;ifstream f(R+"mixed17_switch_weights.bin",ios::binary);while(f.read((char*)&x,4))w.push_back(x);int n=w.size()-536;if(n<=0)return 1;vector<vector<int>>adj(n),blocks(n);ifstream e(R+"mixed17_switch_edges.bin",ios::binary);int a,b;while(e.read((char*)&a,4)&&e.read((char*)&b,4)){if(a>=n)return 2;if(b>=n)blocks[a].push_back(b-n);else{adj[a].push_back(b);adj[b].push_back(a);}}vector<int>count(536);vector<char>chosen(n);vector<int>best;long score=0,bestscore=0;
 auto add=[&](int i){if(chosen[i])exit(3);chosen[i]=1;score+=w[i];for(int j:blocks[i])if(count[j]++==0)score-=w[n+j];};
 auto del=[&](int i){if(!chosen[i])exit(4);chosen[i]=0;score-=w[i];for(int j:blocks[i])if(--count[j]==0)score+=w[n+j];};
 mt19937 rng(170106);uniform_real_distribution<double>U(0,1);cout<<"candidates "<<n<<endl;
 for(int step=0;step<10000000;step++){
  int i=rng()%n;if(chosen[i])continue;long before=score;vector<int>removed;for(int j:adj[i])if(chosen[j]){removed.push_back(j);del(j);}add(i);long delta=score-before;double temp=8.0+64.0*(1.0-(step%250000)/250000.0);
  if(delta<0&&U(rng)>exp(double(delta)/temp)){del(i);for(int j:removed)add(j);}
  if(score>bestscore){bestscore=score;best.clear();for(int j=0;j<n;j++)if(chosen[j])best.push_back(j);ofstream o(R+"mixed17_switch_anneal_best.bin",ios::binary);for(int j:best)o.write((char*)&j,4);cout<<"step "<<step<<" gain "<<bestscore<<" selected "<<best.size()<<endl;}
  if(step%1000000==0)cout<<"step "<<step<<" current "<<score<<" best "<<bestscore<<endl;
 }
 cout<<"finished best gain "<<bestscore<<endl;
}
