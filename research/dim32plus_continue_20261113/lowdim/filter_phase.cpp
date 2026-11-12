#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int maximum_dot=6; using O=array<int,5>;int8_t T[253][253][3];
bool conflict(O u,O v){
 for(int t=0;t<3;t++)if(T[u[0]][v[0]][t]+T[u[1]][v[1]][(v[3]-u[3]+t+3)%3]+T[u[2]][v[2]][(v[4]-u[4]+t+3)%3]>maximum_dot)return true;
 return false;
}
int main(int argc,char**argv){
 if(argc!=3){cerr<<"usage: filter maximum_blockers maximum_dot\n";return 2;} maximum_dot=stoi(argv[2]);
 ifstream tab("T.bin",ios::binary);tab.read((char*)T,sizeof(T));
 vector<O> seed; ifstream s("selected.bin",ios::binary);O u;while(s.read((char*)u.data(),sizeof(u)))seed.push_back(u);
 vector<int> rows[3][253];for(int i=0;i<(int)seed.size();i++)for(int k=0;k<3;k++)rows[k][seed[i][k]].push_back(i);
 int original=seed.size(),limit=stoi(argv[1]);ifstream pool("pool.bin",ios::binary);ofstream small("small.bin",ios::binary);long long scanned=0,kept=0,hist[50]={};
 while(pool.read((char*)u.data(),sizeof(u))){
  int count=0;
  for(int k=0;k<3;k++){
   for(int index:rows[k][u[k]]){
    O v=seed[index];if(k>0&&u[0]==v[0])continue;if(k>1&&u[1]==v[1])continue;
    if(u==v)continue;
    if(conflict(u,v))count++;
    if(count>limit)break;
   }
   if(count>limit)break;
  }
  hist[min(count,49)]++;scanned++;
  if(count==0){
   bool existing=false;for(int j:rows[0][u[0]])if(seed[j]==u){existing=true;break;}
   if(!existing){int i=seed.size();seed.push_back(u);for(int k=0;k<3;k++)rows[k][u[k]].push_back(i);cerr<<"ADDED "<<i-original+1<<" at "<<scanned<<"\n";}
  }
  if(count<=limit){small.write((char*)u.data(),sizeof(u));kept++;}
 }
 ofstream out("expanded.bin",ios::binary);for(auto v:seed)out.write((char*)v.data(),sizeof(v));
 cerr<<"scanned "<<scanned<<" kept "<<kept<<" added "<<seed.size()-original<<" histogram";for(int i=0;i<=limit+1;i++)cerr<<" "<<hist[i];cerr<<"\n";
}
