// Exact intersection screening only; no heuristic floating-point tests.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
struct Support { uint64_t lo,hi; };
int main(int argc,char**argv) {
 if(argc!=3) return 2;
 std::ofstream proof(argv[2],std::ios::binary); if(!proof) return 4;
 std::ifstream f(argv[1],std::ios::binary); std::vector<Support> pool;
 Support s; while(f.read(reinterpret_cast<char*>(&s),sizeof(s))) pool.push_back(s);
 if(pool.size()!=26022) return 3;
 uint64_t total=0,survive=0, checks=0;
 auto test=[&](Support x) {
  ++total;
  for(uint32_t i=0;i<pool.size();i++) { auto y=pool[i]; ++checks; if(__builtin_popcountll(x.lo&y.lo)+__builtin_popcountll(x.hi&y.hi)>12) return i; }
  ++survive; std::cout<<"survivor "<<x.lo<<" "<<x.hi<<"\n"; return UINT32_MAX;
 };
 auto add=[](Support& x,int p) { if(p<64) x.lo|=1ULL<<p; else x.hi|=1ULL<<(p-64); };
 std::array<int,6> a={0,1,2,3,4,5};
 do {
  Support x={0,0}; for(int r:a) for(int c=0;c<4;c++) add(x,4*r+c); uint32_t witness=test(x); proof.write(reinterpret_cast<char*>(&witness),sizeof(witness));
  int k=5; while(k>=0 && a[k]==18+k) --k; if(k<0) break;
  ++a[k]; for(int j=k+1;j<6;j++) a[j]=a[j-1]+1;
 } while(true);
 std::cout<<"full_rows tested="<<total<<" survivors="<<survive<<" checks="<<checks<<"\n";
 std::mt19937_64 rng(20261104); std::array<int,24> rows; for(int j=0;j<24;j++) rows[j]=j;
 total=survive=checks=0;
 for(int t=0;t<1000000;t++) {
  std::shuffle(rows.begin(),rows.end(),rng); Support x={0,0};
  // Six triple rows and six singleton rows: 24 cells on exactly 12 rows.
  for(int j=0;j<12;j++) {int c=rng()%4; for(int d=0;d<4;d++) if(j<6 ? d!=c : d==c) add(x,4*rows[j]+d);}
  test(x);
 }
 std::cout<<"triple_single tested="<<total<<" survivors="<<survive<<" checks="<<checks<<"\n";
}
