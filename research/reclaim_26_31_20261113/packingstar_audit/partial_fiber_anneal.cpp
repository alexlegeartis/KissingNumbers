// Individual-head shared-deletion search; new/new conflicts are evaluated
// directly against the active set instead of storing a quadratic graph.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace std;
template<class T> vector<T> readv(ifstream& f,size_t n){vector<T> v(n);f.read((char*)v.data(),n*sizeof(T));return v;}
int main(int argc,char**argv){
 if(argc!=6)return 2;
 ifstream f(argv[1],ios::binary);auto h=readv<int64_t>(f,3);
 int n=h[0],m=h[1],d=h[2];auto C=readv<int16_t>(f,(size_t)n*d);auto B=readv<int32_t>(f,(size_t)n*4);
 if(!f||d!=14)return 3;
 double seconds=stod(argv[2]);mt19937_64 rng(stoull(argv[3]));bool continuation=string(argv[5])=="continuation";
 vector<vector<int>> rev(m);for(int c=0;c<n;c++)for(int j=0;j<4;j++)rev[B[4*c+j]].push_back(c);
 vector<int> cover(m),selected(n),active,position(n,-1),touched,change(m),stamp(m),remove;
 vector<int> best,nonempty;int gain=0,victims=0,bestgain=0,bestnon=-1000000,maxselected=0,maxvictims=0,epoch=0;
 uint64_t proposals=0,accepted=0,restarts=0;
 auto start=chrono::steady_clock::now();auto elapsed=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};
 auto save=[&](){ofstream q(argv[4]);q<<"{\"gain\":"<<bestgain<<",\"best_nonempty_gain\":"<<bestnon<<",\"max_selected\":"<<maxselected<<",\"victims_at_max_selected\":"<<maxvictims<<",\"proposals\":"<<proposals<<",\"accepted\":"<<accepted<<",\"restarts\":"<<restarts<<",\"seconds\":"<<elapsed()<<",\"selected\":[";for(size_t i=0;i<best.size();i++){if(i)q<<',';q<<best[i];}q<<"],\"best_nonempty_selected\":[";for(size_t i=0;i<nonempty.size();i++){if(i)q<<',';q<<nonempty[i];}q<<"]}\n";};
 while(elapsed()<seconds){
  fill(cover.begin(),cover.end(),0);fill(selected.begin(),selected.end(),0);fill(position.begin(),position.end(),-1);active.clear();gain=0;victims=0;restarts++;
  double rs=elapsed(),duration=min(12.0,seconds-rs);
  while(elapsed()-rs<duration){
   proposals++;double phase=(elapsed()-rs)/max(0.01,duration);
   double temperature=(2.0*pow(2.0,double((restarts-1)%5)))*pow(0.003,phase);
   int c=rng()%n;bool dropping=false;
   if(!active.empty()&&rng()%100<12){c=active[rng()%active.size()];dropping=true;}
   else if(!active.empty()&&rng()%100<85){int a=active[rng()%active.size()];int v=B[4*a+rng()%4];c=rev[v][rng()%rev[v].size()];}
   if(!dropping&&selected[c])continue;
   remove.clear();if(dropping)remove.push_back(c);else{
    for(int a:active){int dot=0;for(int j=0;j<d;j++)dot+=int(C[d*c+j])*int(C[d*a+j]);if(dot>16)remove.push_back(a);}
   }
   if(remove.size()>16)continue;
   ++epoch;touched.clear();
   auto touch=[&](int v,int a){if(stamp[v]!=epoch){stamp[v]=epoch;change[v]=0;touched.push_back(v);}change[v]+=a;};
   for(int a:remove)for(int j=0;j<4;j++)touch(B[4*a+j],-1);
   if(!dropping)for(int j=0;j<4;j++)touch(B[4*c+j],1);
   int dv=0;for(int v:touched)dv+=(cover[v]+change[v]>0)-(cover[v]>0);
   int delta=(!dropping)-int(remove.size())-dv;
   double penalty=continuation?0.02+0.98*pow(phase,2.0):1.0;
   double merit=(!dropping)-int(remove.size())-penalty*dv;
   if(continuation)temperature=.25*pow(.02,phase);
   double u=(rng()+.5)/(double(UINT64_MAX)+1.0);
   if(merit<0&&u>=exp(merit/temperature))continue;
   accepted++;gain+=delta;victims+=dv;for(int v:touched)cover[v]+=change[v];
   for(int a:remove){int i=position[a],last=active.back();active[i]=last;position[last]=i;active.pop_back();position[a]=-1;selected[a]=0;}
   if(!dropping){position[c]=active.size();active.push_back(c);selected[c]=1;}
   if(int(active.size())>maxselected){maxselected=active.size();maxvictims=victims;}
   if(!active.empty()&&gain>bestnon){bestnon=gain;nonempty=active;}
   if(gain>bestgain){bestgain=gain;best=active;save();cout<<"gain "<<gain<<" selected "<<active.size()<<" victims "<<victims<<endl;}
  }
  save();cout<<"restart "<<restarts<<" best_gain "<<bestgain<<" best_nonempty "<<bestnon<<" max_selected "<<maxselected<<" proposals "<<proposals<<endl;
 }
 save();return 0;
}
