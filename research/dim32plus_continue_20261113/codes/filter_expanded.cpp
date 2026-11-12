#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using U=uint64_t;
struct Rec {U a;std::array<uint32_t,8> old;};
template<class T>std::vector<T> read(const char*p){
 std::ifstream f(p,std::ios::binary|std::ios::ate);assert(f);size_t n=f.tellg();assert(n%sizeof(T)==0);f.seekg(0);
 std::vector<T>a(n/sizeof(T));assert(f.read((char*)a.data(),n));return a;
}
int main(int argc,char**argv){
 assert(argc==4&&sizeof(Rec)==40);auto pool=read<Rec>(argv[1]);auto extra=read<U>(argv[2]);std::vector<Rec>out;
 for(auto&r:pool){bool good=true;for(U b:extra)if(__builtin_popcountll(r.a&b)>8){good=false;break;}if(good)out.push_back(r);}
 std::ofstream f(argv[3],std::ios::binary);f.write((char*)out.data(),out.size()*sizeof(Rec));assert(f);
 std::cout<<"candidates "<<pool.size()<<" extra "<<extra.size()<<" retained "<<out.size()<<"\n";
}
