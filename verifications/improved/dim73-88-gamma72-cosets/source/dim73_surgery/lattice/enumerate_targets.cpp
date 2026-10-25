// Enumerate short vectors of prescribed residue modulo 4 in the Leech lattice.
// Outward rounded intervals enclose the exact rational LDL decomposition.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
using namespace std;
constexpr int N=24;
int64_t gram[N][N];
int32_t residue[N],v[N],case_id;
double dl[N],dh[N],ml[N][N],mh[N][N];
int radius2;
uint64_t nodes=0,leaves=0;
ofstream output;
double down(double x){return nextafter(x,-numeric_limits<double>::infinity());}
double up(double x){return nextafter(x,numeric_limits<double>::infinity());}
void visit(int k,double lower){
    ++nodes;
    if(k<0){
        int64_t norm=0;
        for(int i=0;i<N;++i)for(int j=0;j<N;++j)norm+=int64_t(v[i])*gram[i][j]*v[j];
        if(norm>radius2)return;
        if(norm<radius2){cerr<<"SHORTER vector in case "<<case_id<<": "<<norm<<endl;exit(5);}
        output.write((char*)&case_id,4);output.write((char*)v,4*N);++leaves;return;
    }
    double cl=0,ch=0;
    for(int j=k+1;j<N;++j){
        double lo,hi;
        if(v[j]>=0){lo=down(ml[j][k]*v[j]);hi=up(mh[j][k]*v[j]);}
        else{lo=down(mh[j][k]*v[j]);hi=up(ml[j][k]*v[j]);}
        cl=down(cl+lo);ch=up(ch+hi);
    }
    double rem=up(radius2-lower);if(rem<0)return;
    double bound=up(sqrt(up(rem/dl[k])));
    int low=int(ceil(down(-ch-bound))),high=int(floor(up(-cl+bound)));
    low+=((residue[k]-low)%4+4)%4;
    for(int z=low;z<=high;z+=4){
        double zl=down(z+cl),zh=up(z+ch),sq=0;
        if(zl>0)sq=down(zl*zl);else if(zh<0)sq=down(zh*zh);
        double next=down(lower+down(dl[k]*max(0.0,sq)));
        if(next>radius2)continue;
        v[k]=z;visit(k-1,max(0.0,next));
    }
}
int main(int argc,char**argv){
    if(argc!=4){cerr<<"ldl cases output-prefix\n";return 2;}
    ifstream in(argv[1]),cases(argv[2]);int n,count;in>>n;if(n!=N)return 3;
    for(auto&row:gram)for(auto&x:row)in>>x;
    for(int i=0;i<N;++i)in>>dl[i]>>dh[i];
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)in>>ml[i][j]>>mh[i][j];
    cases>>count;if(!in||!cases)return 4;
    output.open(string(argv[3])+".vectors.bin",ios::binary);
    auto started=chrono::steady_clock::now();
    for(case_id=0;case_id<count;++case_id){
        cases>>radius2;for(auto&x:residue)cases>>x;if(!cases)return 4;
        visit(23,0);
    }
    double elapsed=chrono::duration<double>(chrono::steady_clock::now()-started).count();
    ofstream report(string(argv[3])+".enum.json");
    report<<"{\"complete\":true,\"cases\":"<<count<<",\"nodes\":"<<nodes<<",\"vectors\":"<<leaves<<",\"seconds\":"<<elapsed<<"}\n";
    cout<<"cases="<<count<<" nodes="<<nodes<<" vectors="<<leaves<<" seconds="<<elapsed<<endl;
}
