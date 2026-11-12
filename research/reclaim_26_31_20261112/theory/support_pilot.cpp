#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>
using namespace std;
using U=uint32_t;
mt19937 rng(20261112);
long long bc[33][9];
const int N=31,K=7;
U tailmask=((1u<<7)-1)<<24;
vector<U> pool, golay, parent;
set<U> golayset;
vector<int> owner;
vector<U> selected;
double now(){return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();}
template<class F> bool fives(U x,F fn){
 int a[8],p=0;while(x){int j=countr_zero(x);x&=x-1;a[p++]=j;}
 for(int i=0;i<4;i++)for(int j=i+1;j<5;j++)for(int k=j+1;k<6;k++)for(int l=k+1;l<7;l++)for(int m=l+1;m<8;m++)
   if(!fn((int)(bc[a[i]][1]+bc[a[j]][2]+bc[a[k]][3]+bc[a[l]][4]+bc[a[m]][5])))return false;
 return true;
}
bool add(U x){
 if(!fives(x,[](int q){return owner[q]<0;}))return false;
 int i=selected.size();selected.push_back(x);fives(x,[i](int q){owner[q]=i;return true;});return true;
}
void reset(const vector<U>&s){selected.clear();fill(owner.begin(),owner.end(),-1);for(U x:s){if(!add(x)){cerr<<"bad seed"<<endl;exit(2);}}}
void enumerate(int pos,int need,U x){
 if(!need){int r=popcount(x&tailmask);if(r>=2||golayset.count(x))pool.push_back(x);return;}
 for(int j=pos;j<=N-need;j++)enumerate(j+1,need-1,x|(1u<<j));
}
int score(U tails,int del){int n=0;for(U x:parent)if(!(x&(1u<<del))&&popcount(x&tails)>=2)n++;return n;}
void output(const string&p,const vector<U>&s){ofstream f(p);for(U x:s)f<<hex<<x<<"\n";}
int main(int argc,char**argv){
 string dir=argv[1];double budget=argc>2?stod(argv[2]):90;
 for(int n=0;n<=32;n++){bc[n][0]=1;for(int k=1;k<=8;k++)bc[n][k]=(n==0?0:bc[n-1][k-1]+bc[n-1][k]);}
 for(string fn:{"golay_octads.hex","a32.8.8.1713H"}){ifstream f(dir+"/"+fn);string line;while(getline(f,line))if(line.size()&&line[0]!='$'){U x=stoul(line,nullptr,16);if(fn[0]=='g'){golay.push_back(x);golayset.insert(x);}else parent.push_back(x);}}
 int bestscore=-1,bestdel=-1;U besttails=0;
 for(int del=0;del<32;del++)for(int start=0;start<6;start++){
   vector<int> ps;for(int j=0;j<32;j++)if(j!=del)ps.push_back(j);shuffle(ps.begin(),ps.end(),rng);U tails=0;for(int j=0;j<7;j++)tails|=1u<<ps[j];int cur=score(tails,del);
   for(int iter=0;iter<20;iter++){
     int bs=cur;U bt=tails;
     for(int a=0;a<32;a++)if(tails&(1u<<a))for(int b=0;b<32;b++)if(b!=del&&!(tails&(1u<<b))){U t=tails^(1u<<a)^(1u<<b);int s=score(t,del);if(s>bs){bs=s;bt=t;}}
     if(bs==cur)break;tails=bt;cur=bs;
   }
   if(cur>bestscore){bestscore=cur;bestdel=del;besttails=tails;}
 }
 cout<<"parent_best "<<bestscore<<" deleted "<<bestdel<<" tails "<<hex<<besttails<<dec<<endl;
 vector<int> perm(32,-1);int ci=0,ti=24;for(int j=0;j<32;j++)if(j!=bestdel)perm[j]=(besttails&(1u<<j))?ti++:ci++;
 vector<U> ps;for(U x:parent)if(!(x&(1u<<bestdel))&&popcount(x&besttails)>=2){U y=0;for(int j=0;j<32;j++)if(x&(1u<<j))y|=1u<<perm[j];ps.push_back(y);}
 owner.resize(bc[N][5],-1);enumerate(0,8,0);cout<<"pool "<<pool.size()<<endl;
 shuffle(pool.begin(),pool.end(),rng);reset(ps);for(U x:pool)add(x);vector<U>best=selected;cout<<"parent_fill "<<best.size()<<endl;output(dir+"/support31_best.hex",best);
 reset(golay);shuffle(pool.begin(),pool.end(),rng);for(U x:pool)add(x);cout<<"golay_fill "<<selected.size()<<endl;if(selected.size()>best.size()){best=selected;output(dir+"/support31_best.hex",best);}
 double start=now();int rounds=0;
 while(now()-start<budget){
   vector<U> base=best;double frac=.02+.13*(rng()%100)/100.;int mode=rng()%4;
   if(mode==0){U pair=0;while(popcount(pair)<2)pair|=1u<<(rng()%N);base.erase(remove_if(base.begin(),base.end(),[&](U x){return (x&pair)==pair;}),base.end());}
   else{shuffle(base.begin(),base.end(),rng);base.resize((size_t)(base.size()*(1-frac)));}
   reset(base);shuffle(pool.begin(),pool.end(),rng);for(U x:pool)add(x);rounds++;
   if(selected.size()>best.size()){best=selected;output(dir+"/support31_best.hex",best);cout<<"best "<<best.size()<<" round "<<rounds<<" seconds "<<now()-start<<endl;}
 }
 int bad=0;for(size_t i=0;i<best.size();i++)for(size_t j=0;j<i;j++)if(popcount(best[i]&best[j])>4)bad++;
 cout<<"final "<<best.size()<<" kissing_count "<<98304+1860+128*best.size()<<" rounds "<<rounds<<" seconds "<<now()-start<<" violations "<<bad<<endl;
 return bad?1:0;
}
