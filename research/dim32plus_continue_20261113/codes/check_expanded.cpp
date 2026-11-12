// Independent scalar pair proof for every retained 73-bit weight-sixteen support.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
using U=uint64_t;
void require(bool b){if(!b)throw std::runtime_error("expanded support certificate rejected");}
int wt(U x){return __builtin_popcountll(x);}
struct Word{U lo,hi;};
int dot(Word a,Word b){return wt(a.lo&b.lo)+wt(a.hi&b.hi);}
int main(int argc,char**argv){try{
 require(argc==2);std::ifstream f(argv[1]);size_t n,k;f>>n>>k;require(bool(f)&&(n==34028||n==40428));
 std::vector<int>head(64);std::set<int>positions;Word body{0,0};
 for(int&c:head){f>>c;require(c>=0&&c<73&&positions.insert(c).second);if(c<64)body.lo|=U(1)<<c;else body.hi|=U(1)<<(c-64);}
 std::vector<Word>parent(n),extra;std::set<std::pair<U,U>>seen;
 for(auto&w:parent){f>>w.lo>>w.hi;require(bool(f)&&w.hi<512&&wt(w.lo)+wt(w.hi)==16&&seen.emplace(w.lo,w.hi).second);if((w.lo&~body.lo)||(w.hi&~body.hi))extra.push_back(w);}
 require(extra.size()==n-30828);U pairs=0,added_pairs=0;
 for(size_t i=0;i<n;i++)for(size_t j=0;j<i;j++){require(dot(parent[i],parent[j])<=8);pairs++;}
 for(size_t i=0;i<k;i++){U s;f>>s;require(bool(f)&&wt(s)==16);Word w{0,0};for(int j=0;j<64;j++)if(s>>j&1){int c=head[j];if(c<64)w.lo|=U(1)<<c;else w.hi|=U(1)<<(c-64);}for(auto b:extra){require(dot(w,b)<=8);added_pairs++;}}
 std::cout<<"{\"status\":\"PASS\",\"parent_support_pairs\":"<<pairs<<",\"new_extra_pairs\":"<<added_pairs<<",\"extra_supports\":"<<extra.size()<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
