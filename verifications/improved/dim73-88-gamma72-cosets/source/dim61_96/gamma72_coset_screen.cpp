// Exact integer screening of all 2^24 third-block cosets from Leech lists.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
constexpr uint32_t MASK=(1U<<24)-1;
uint32_t amap[3][256],nmap[3][256],arows[24];
int pivots[12];
uint32_t apply(uint32_t table[3][256],uint32_t z){return table[0][z&255]^table[1][(z>>8)&255]^table[2][z>>16];}
vector<uint32_t> readbin(const string&path){
    ifstream f(path,ios::binary|ios::ate);size_t size=f.tellg();f.seekg(0);
    vector<uint32_t> v(size/4);f.read((char*)v.data(),size);if(!f)throw runtime_error("read "+path);return v;
}
void maps(const string&path,uint32_t&mask,int&base){
    ifstream f(path);for(auto&x:arows)f>>x;f>>mask>>base;if(!f)throw runtime_error("map read");
    uint32_t basis[24]={};
    for(int i=0;i<24;++i){uint32_t r=arows[i]^(1U<<i);for(int j=0;j<24;++j)if(r>>j&1){if(basis[j])r^=basis[j];else{basis[j]=r;break;}}}
    int count=0;for(int j=0;j<24;++j)if(basis[j])pivots[count++]=j;
    if(count!=12)throw runtime_error("N rank");
    // Row-echelon pivot coordinates are not coefficient coordinates, but their
    // restriction to this subspace is injective, which is all grouping needs.
    for(int part=0;part<3;++part)for(int z=0;z<256;++z){
        uint32_t a=0;for(int bit=0;bit<8;++bit)if(z>>bit&1)a^=arows[8*part+bit];
        amap[part][z]=a;uint32_t b=a^(uint32_t(z)<<(8*part)),r=0;
        for(int j=0;j<12;++j)r|=((b>>pivots[j])&1)<<j;
        nmap[part][z]=r;
    }
}
int main(int argc,char**argv){
    if(argc!=6){cerr<<"map list1 list2 list3 output-prefix\n";return 2;}
    uint32_t mask;int base;maps(argv[1],mask,base);
    vector<uint32_t> raw[3];for(int i=0;i<3;++i)raw[i]=readbin(argv[2+i]);
    vector<uint32_t> grouped[2][6][4096],third[6];
    for(int i=0;i<2;++i)for(uint32_t rec:raw[i]){
        int cost=rec>>24;uint32_t z=rec&MASK;
        if(cost>=2&&cost<=5)grouped[i][cost][apply(nmap,z)].push_back(z);
    }
    for(uint32_t rec:raw[2])if((rec>>24)>=2&&(rec>>24)<=5)third[rec>>24].push_back(rec&MASK);
    vector<uint8_t> minimum(1<<24,255);
    vector<uint32_t> count(1<<24,0);
    uint64_t pairs=0,triples=0;auto started=chrono::steady_clock::now();
    for(int sum=4;sum<=7;++sum){
        unordered_map<uint32_t,uint32_t> p;
        for(int c1=2;c1<=5;++c1){int c2=sum-c1;if(c2<2||c2>5)continue;
            for(int r=0;r<4096;++r)for(uint32_t z1:grouped[0][c1][r]){
                uint32_t az1=apply(amap,z1);
                for(uint32_t z2:grouped[1][c2][r]){++p[az1^z2];++pairs;}
            }
        }
        for(int c3=2;c3<=9-sum;++c3)for(auto item:p)for(uint32_t z3:third[c3]){
            uint32_t h=item.first^z3;int total=sum+c3;
            if(minimum[h]>total){minimum[h]=total;count[h]=item.second;}
            else if(minimum[h]==total)count[h]+=item.second;
            ++triples;
        }
        cout<<"pair norm="<<sum<<" distinct="<<p.size()<<" cumulative triples="<<triples<<endl;
    }
    map<int,uint64_t> hist;uint64_t candidates=0;ofstream cand(string(argv[5])+".candidates.bin",ios::binary);
    for(uint32_t h=0;h<=MASK;++h)if((__builtin_parity(h&mask)^base)==1){
        ++hist[minimum[h]];
        if(minimum[h]==255){++candidates;cand.write((char*)&h,4);}
    }
    ofstream vals(string(argv[5])+".minimum.bin",ios::binary);vals.write((char*)minimum.data(),minimum.size());
    ofstream counts(string(argv[5])+".counts.bin",ios::binary);counts.write((char*)count.data(),4*count.size());
    double elapsed=chrono::duration<double>(chrono::steady_clock::now()-started).count();
    ofstream out(string(argv[5])+".json");out<<"{\n  \"anisotropic_cosets\": 8388608,\n  \"candidates_above_9\": "<<candidates
        <<",\n  \"pair_products\": "<<pairs<<",\n  \"triple_products\": "<<triples<<",\n  \"seconds\": "<<elapsed
        <<",\n  \"minimum_histogram_unscaled\": {";
    bool first=true;for(auto item:hist){out<<(first?"":", ")<<'"'<<item.first<<"\": "<<item.second;first=false;}
    out<<"}\n}\n";cout<<"candidates="<<candidates<<" seconds="<<elapsed<<endl;
    return 0;
}
