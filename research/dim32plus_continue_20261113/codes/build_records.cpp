#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;
struct Rec { uint64_t a; array<uint32_t,8> old; };
vector<uint64_t> read(const char* path) {
  ifstream f(path,ios::binary|ios::ate); assert(f);
  auto bytes=f.tellg(); assert(bytes%8==0); f.seekg(0);
  vector<uint64_t> a(size_t(bytes)/8); f.read((char*)a.data(),bytes); assert(f); return a;
}
int main(int argc,char**argv) {
  assert(argc==4 && sizeof(Rec)==40);
  auto candidates=read(argv[1]), old=read(argv[2]);
  vector<Rec> records(candidates.size());
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<candidates.size();i++) {
    Rec r; r.a=candidates[i]; int count=0;
    for(size_t j=0;j<old.size();j++) {
      int k=__builtin_popcountll(r.a&old[j]);
      if(k>8) { assert(k==10 && count<8); r.old[count++]=j; }
    }
    assert(count==8); records[i]=r;
  }
  ofstream out(argv[3],ios::binary); out.write((char*)records.data(),records.size()*sizeof(Rec));
  cout<<"records "<<records.size()<<" old "<<old.size()<<" exact pair checks "<<records.size()*old.size()<<endl;
}
