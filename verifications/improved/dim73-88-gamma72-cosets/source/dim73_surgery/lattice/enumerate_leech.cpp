// 24-dimensional fixed-parity list enumeration in the Leech lattice.
// Interval mode encloses an exact rational LDL; fast mode is exploratory.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace std;
constexpr int N=24;
long long gram[N][N];
int parity[N],q[N],v[N];
double dl[N],dh[N],ml[N][N],mh[N][N],mu[N][N],d[N];
double radius2,seconds;
uint64_t nodes=0,leaves=0,hist[256]={},levels[N]={};
bool stopped=false,fast=false;
ofstream records;
ofstream selected_vectors;
vector<uint32_t> wanted;
size_t wanted_position=0;
chrono::steady_clock::time_point started;
inline double down(double x){return nextafter(x,-numeric_limits<double>::infinity());}
inline double up(double x){return nextafter(x,numeric_limits<double>::infinity());}
void leaf(){
    long long norm=0;
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)norm+=(long long)v[i]*gram[i][j]*v[j];
    if(norm>radius2)return;
    if(norm<0 || norm%4 || norm/4>255){cerr<<"Unexpected norm "<<norm<<endl;exit(5);}
    uint32_t key=0;
    for(int i=0;i<N;++i)key|=uint32_t(((v[i]-q[i])/2)&1)<<i;
    uint32_t record=key|uint32_t(norm/4)<<24;
    records.write((char*)&record,4);
    if(wanted_position<wanted.size() && wanted[wanted_position]==leaves){
        uint32_t index=(uint32_t)leaves;
        selected_vectors.write((char*)&index,4);
        for(int i=0;i<N;++i){int32_t coord=v[i];selected_vectors.write((char*)&coord,4);}
        ++wanted_position;
    }
    ++leaves;++hist[norm/4];
}
inline bool tick(int k){
    ++levels[k];
    if(((++nodes)&((1ULL<<20)-1))==0){
        if(chrono::duration<double>(chrono::steady_clock::now()-started).count()>seconds){
            stopped=true;return false;
        }
    }
    return true;
}
void visit_fast(int k,double partial){
    if(k<0){leaf();return;}
    if(!tick(k))return;
    double c=0;for(int j=k+1;j<N;++j)c+=mu[j][k]*v[j];
    double bound=sqrt(max(0.0,(radius2+1e-7-partial)/d[k]));
    int low=(int)ceil(-c-bound-1e-8),high=(int)floor(-c+bound+1e-8);
    if((low&1)!=parity[k])++low;
    for(int z=low;z<=high;z+=2){
        double next=partial+d[k]*(z+c)*(z+c);
        if(next>radius2+1e-7)continue;
        v[k]=z;visit_fast(k-1,next);if(stopped)return;
    }
}
void visit_interval(int k,double lower){
    if(k<0){leaf();return;}
    if(!tick(k))return;
    double cl=0,ch=0;
    for(int j=k+1;j<N;++j){
        double a,b;
        if(v[j]>=0){a=down(ml[j][k]*v[j]);b=up(mh[j][k]*v[j]);}
        else{a=down(mh[j][k]*v[j]);b=up(ml[j][k]*v[j]);}
        cl=down(cl+a);ch=up(ch+b);
    }
    double rem=up(radius2-lower);if(rem<0)return;
    double bound=up(sqrt(up(rem/dl[k])));
    int low=(int)ceil(down(-ch-bound)),high=(int)floor(up(-cl+bound));
    if((low&1)!=parity[k])++low;
    for(int z=low;z<=high;z+=2){
        double zl=down(z+cl),zh=up(z+ch),sq=0;
        if(zl>0)sq=down(zl*zl);else if(zh<0)sq=down(zh*zh);
        double next=down(lower+down(dl[k]*max(0.0,sq)));
        if(next>radius2)continue;
        v[k]=z;visit_interval(k-1,max(0.0,next));if(stopped)return;
    }
}
int main(int argc,char**argv){
    if(argc!=7&&argc!=8){cerr<<"ldl q radius seconds output-prefix fast|interval [wanted-indices]\n";return 2;}
    int n;ifstream in(argv[1]),p(argv[2]);in>>n;if(n!=N)return 3;
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)in>>gram[i][j];
    for(int i=0;i<N;++i){in>>dl[i]>>dh[i];d[i]=(dl[i]+dh[i])/2;}
    for(int i=0;i<N;++i)for(int j=0;j<N;++j){in>>ml[i][j]>>mh[i][j];mu[i][j]=(ml[i][j]+mh[i][j])/2;}
    for(int i=0;i<N;++i){p>>q[i];parity[i]=q[i]&1;}
    if(!in||!p)return 4;
    radius2=stod(argv[3]);seconds=stod(argv[4]);fast=string(argv[6])=="fast";
    records.open(string(argv[5])+".bin",ios::binary);
    if(argc==8){
        ifstream f(argv[7],ios::binary|ios::ate);size_t size=f.tellg();f.seekg(0);
        wanted.resize(size/4);f.read((char*)wanted.data(),size);if(!f)return 6;
        selected_vectors.open(string(argv[5])+".selected.bin",ios::binary);
    }
    started=chrono::steady_clock::now();
    if(fast)visit_fast(N-1,0);else visit_interval(N-1,0);
    if(!stopped && wanted_position!=wanted.size()){cerr<<"missing selected vectors\n";return 7;}
    double elapsed=chrono::duration<double>(chrono::steady_clock::now()-started).count();
    ofstream out(string(argv[5])+".json");
    out<<"{\n  \"complete\": "<<(stopped?"false":"true")<<",\n  \"certified_intervals\": "<<(fast?"false":"true")
       <<",\n  \"radius_squared_leech\": "<<radius2/4<<",\n  \"nodes\": "<<nodes
       <<",\n  \"seconds\": "<<elapsed<<",\n  \"histogram\": {";
    bool first=true;for(int i=0;i<256;++i)if(hist[i]){out<<(first?"":", ")<<'"'<<i<<"\": "<<hist[i];first=false;}
    out<<"}\n}\n";
    cout<<"mode="<<(fast?"fast":"interval")<<" complete="<<!stopped<<" nodes="<<nodes
        <<" records="<<leaves<<" seconds="<<elapsed<<" histogram=";
    for(int i=0;i<256;++i)if(hist[i])cout<<i<<":"<<hist[i]<<" ";cout<<endl;
}
